#!/usr/bin/env python3
"""Fig. 5c sparse adaptive short-time drift experiment."""

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from hsr_sim.adaptive import relocation_attempts_per_mcs, relocation_proposal  # noqa: E402
from hsr_sim.adaptive_theory import beta_c, d_beta, p_ad_star, p_c, p_static_star  # noqa: E402
from hsr_sim.engine import C, D, SimParams, db_update, graph_to_neighbors, precompute_neighbor_combos, steering_count  # noqa: E402


DRIFT_FIELDS = [
    "exp_id", "replicate", "seed", "N", "k", "n", "M", "r0", "r", "R",
    "omega", "c", "p_target", "p_actual", "x_target", "x_start", "x_end",
    "beta", "gamma", "relax_mcs", "tau_mcs", "delta_x", "drift", "attempts",
    "relocation_accepted", "x_after_relax", "z", "eta0", "d_beta_theory",
    "p_c_theory",
]

PC_FIELDS = ["exp_id", "gamma", "beta", "x_target", "p_c_theory", "p_c_empirical", "p_c_ci95_low", "p_c_ci95_high", "num_p", "num_replicates"]
GLOBAL_FIELDS = ["exp_id", "beta", "beta_c", "p_static_star", "theory_p_ad_star", "theory_x_bottleneck", "theory_branch", "empirical_p_ad_star", "empirical_ci95_low", "empirical_ci95_high", "normalized_theory", "normalized_empirical"]
GAMMA_FIELDS = ["exp_id", "gamma", "beta", "x_target", "theory_p_c", "empirical_p_c", "num_replicates"]


def crossing_from_means(means: list[tuple[float, float]]) -> float:
    """Linear interpolation at the first observed drift sign change."""
    for (p0, g0), (p1, g1) in zip(means[:-1], means[1:]):
        if g0 == 0:
            return float(p0)
        if g0 * g1 <= 0:
            return float(p0 + (0.0 - g0) * (p1 - p0) / (g1 - g0)) if g1 != g0 else float(p0)
    return float("nan")


def bootstrap_crossing(
    grouped: dict[float, list[float]],
    bootstrap: int,
    seed: int,
) -> tuple[float, float, float]:
    """Return crossing and percentile CI from replicate-level drift values."""
    means = [(p, float(np.mean(values))) for p, values in sorted(grouped.items())]
    estimate = crossing_from_means(means)
    if bootstrap <= 0 or not np.isfinite(estimate):
        return estimate, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    roots = []
    for _ in range(bootstrap):
        sampled = [(p, float(np.mean(rng.choice(values, size=len(values), replace=True)))) for p, values in sorted(grouped.items())]
        root = crossing_from_means(sampled)
        if np.isfinite(root):
            roots.append(root)
    if len(roots) < max(20, bootstrap // 20):
        return estimate, float("nan"), float("nan")
    low, high = np.percentile(np.asarray(roots), [2.5, 97.5])
    return estimate, float(low), float(high)


def bootstrap_global_envelope(
    by_x: dict[float, dict[float, list[float]]],
    bootstrap: int,
    seed: int,
) -> tuple[float, float, float]:
    """Bootstrap max_x p_c(x) from nested x -> p -> replicate drifts."""
    point_roots = []
    for grouped in by_x.values():
        root, _, _ = bootstrap_crossing(grouped, 0, seed)
        if np.isfinite(root):
            point_roots.append(root)
    estimate = max(point_roots) if point_roots else float("nan")
    if bootstrap <= 0 or not by_x:
        return estimate, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    envelopes = []
    for _ in range(bootstrap):
        roots = []
        for grouped in by_x.values():
            sampled = {
                p: list(rng.choice(values, size=len(values), replace=True))
                for p, values in grouped.items()
            }
            root = crossing_from_means([(p, float(np.mean(v))) for p, v in sorted(sampled.items())])
            if np.isfinite(root):
                roots.append(root)
        if roots:
            envelopes.append(max(roots))
    if len(envelopes) < max(20, bootstrap // 20):
        return estimate, float("nan"), float("nan")
    low, high = np.percentile(np.asarray(envelopes), [2.5, 97.5])
    return estimate, float(low), float(high)


def exact_state(N: int, p: float, x: float, rng: np.random.Generator):
    nS = steering_count(N, p)
    mask = np.zeros(N, dtype=bool)
    if nS:
        mask[rng.choice(N, size=nS, replace=False)] = True
    strategy = np.zeros(N, dtype=np.int8); strategy[mask] = C
    free = np.flatnonzero(~mask); nC = int(round(x * free.size))
    if nC:
        strategy[rng.choice(free, size=nC, replace=False)] = C
    return strategy, mask


def exposure(strategy, neighbors):
    return np.asarray([np.mean(strategy[nb] == D) for nb in neighbors], dtype=float)


def neutral_relocate(strategy, mask, neighbors, combos, params, beta, gamma, relax_mcs, rng):
    neutral = SimParams(**{**params.__dict__, "omega": 0.0})
    steps = int(round(relax_mcs * params.N))
    for step in range(steps):
        db_update(strategy, mask, neighbors, combos, neutral, rng)
        # Preserve the requested x while allowing neutral pair organization.
        free = np.flatnonzero(~mask); c_nodes = free[strategy[free] == C]; d_nodes = free[strategy[free] == D]
        target = int(round(params.x0 * free.size)); delta = target - c_nodes.size
        if delta > 0 and d_nodes.size: strategy[rng.choice(d_nodes, size=min(delta, d_nodes.size), replace=False)] = C
        elif delta < 0 and c_nodes.size: strategy[rng.choice(c_nodes, size=min(-delta, c_nodes.size), replace=False)] = D
        # One DB update is 1/N evolutionary time; its relocation count has
        # Poisson mean gamma*p/N?  Across N updates this is gamma*p*N per MCS.
        nrel = int(rng.poisson(gamma * float(mask.mean())))
        h = exposure(strategy, neighbors)
        for _ in range(nrel): relocation_proposal(strategy, mask, neighbors, beta, rng, exposure=h)


def run_one(task):
    cfg, rep = task; seed = int(cfg["seed0"] + 1_000_000 * rep + cfg["p_index"] * 10_000 + cfg["beta_index"] * 100)
    rng = np.random.default_rng(seed); N, k, n = cfg["N"], cfg["k"], cfg["n"]
    r0 = k * n / (k + n - 1); M = math.comb(k, n - 1); r = cfg["R"] * r0
    params = SimParams(N=N, k=k, n=n, r=r, omega=cfg["omega"], c=cfg["c"], M_mode="all", p_s=cfg["p"], x0=cfg["x"], T_mcs=1)
    G = nx.random_regular_graph(k, N, seed=int(rng.integers(0, 2**32 - 1))); neighbors = graph_to_neighbors(G, N); combos = precompute_neighbor_combos(neighbors, n)
    strategy, mask = exact_state(N, cfg["p"], cfg["x"], rng); neutral_relocate(strategy, mask, neighbors, combos, params, cfg["beta"], cfg["gamma"], cfg["relax_mcs"], rng)
    x0 = float(np.mean(strategy[~mask] == C)); h = exposure(strategy, neighbors); eta0 = float(mask.mean() / (1-mask.mean()) * np.mean(h[mask])) if mask.any() else 0.0
    z = float(sum((strategy[u] == C and strategy[v] == C) for u, v in G.edges()) * 2 / (k*N))
    steps = max(1, int(round(cfg["tau_mcs"] * N))); accepted = 0; attempts = 0
    for _ in range(steps):
        db_update(strategy, mask, neighbors, combos, params, rng)
        h = exposure(strategy, neighbors)
        # Distribute the continuous-time relocation process across DB events:
        # N DB updates per MCS and mean gamma*p proposals per DB update.
        nrel = int(rng.poisson(cfg["gamma"] * float(mask.mean())))
        for _ in range(nrel):
            result = relocation_proposal(strategy, mask, neighbors, cfg["beta"], rng, exposure=h); accepted += int(result.accepted); attempts += 1
    x1 = float(np.mean(strategy[~mask] == C)); p_actual = float(mask.mean())
    return {"exp_id": cfg["exp_id"], "replicate": rep, "seed": seed, "N": N, "k": k, "n": n, "M": M, "r0": r0, "r": r, "R": cfg["R"], "omega": cfg["omega"], "c": cfg["c"], "p_target": cfg["p"], "p_actual": p_actual, "x_target": cfg["x"], "x_start": x0, "x_end": x1, "beta": cfg["beta"], "gamma": cfg["gamma"], "relax_mcs": cfg["relax_mcs"], "tau_mcs": cfg["tau_mcs"], "delta_x": x1-x0, "drift": (x1-x0)/cfg["tau_mcs"], "attempts": attempts, "relocation_accepted": accepted, "x_after_relax": x0, "z": z, "eta0": eta0, "d_beta_theory": d_beta(cfg["x"], cfg["beta"], k), "p_c_theory": p_c(cfg["x"], cfg["beta"], k, cfg["omega"], M, cfg["c"], r, r0)}


def floats(text): return [float(x) for x in text.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--out-dir", default="results/fig5c_threshold"); ap.add_argument("--exp-id", default="Fig5cThresholdPilot"); ap.add_argument("--N", type=int, default=1000); ap.add_argument("--R", type=float, default=0.5); ap.add_argument("--omega", type=float, default=1e-3); ap.add_argument("--c", type=float, default=1.0); ap.add_argument("--p-grid", default="0.04,0.06,0.08,0.10,0.12"); ap.add_argument("--x-grid", default="0.8,0.9,0.95"); ap.add_argument("--beta-grid", default="0,2,5"); ap.add_argument("--gamma", type=float, default=10); ap.add_argument("--gamma-grid", default="", help="Optional comma-separated gamma values for convergence."); ap.add_argument("--relax-mcs", type=float, default=5); ap.add_argument("--tau-mcs", type=float, default=2); ap.add_argument("--replicates", type=int, default=10); ap.add_argument("--bootstrap", type=int, default=500); ap.add_argument("--seed0", type=int, default=20260904); ap.add_argument("--workers", type=int, default=1); ap.add_argument("--progress-every", type=int, default=10); ap.add_argument("--no-resume", action="store_true"); args=ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True); ps, xs, bs = floats(args.p_grid), floats(args.x_grid), floats(args.beta_grid); gs = floats(args.gamma_grid) if args.gamma_grid else [args.gamma]
    tasks = [{"exp_id": args.exp_id, "N": args.N, "k": 10, "n": 5, "R": args.R, "omega": args.omega, "c": args.c, "p": p, "x": x, "beta": b, "gamma": g, "relax_mcs": args.relax_mcs, "tau_mcs": args.tau_mcs, "seed0": args.seed0, "p_index": ip, "beta_index": ib, "gamma_index": ig} for ig,g in enumerate(gs) for ip,p in enumerate(ps) for ib,b in enumerate(bs) for x in xs]
    jobs = [(cfg, rep) for cfg in tasks for rep in range(args.replicates)]
    metadata = {"experiment": args.exp_id, "model": "sparse adaptive short-time generator drift", "parameters": {"N": args.N, "k": 10, "n": 5, "R": args.R, "omega": args.omega, "c": args.c, "p_grid": ps, "x_grid": xs, "beta_grid": bs, "gamma_grid": gs, "relax_mcs": args.relax_mcs, "tau_mcs": args.tau_mcs, "replicates": args.replicates}, "rate_definition": "gamma*p*N relocation proposals per MCS in expectation", "resume_safe": True}
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    checkpoint = out / "checkpoint.json"; completed = set()
    if checkpoint.exists() and not args.no_resume:
        completed = set(json.loads(checkpoint.read_text(encoding="utf-8")).get("completed", []))
    # A checkpoint may have been written before the raw CSV flush completed.
    # Trust only task keys that are actually present in the raw output.
    raw_existing = out / "fig5c_local_drift.csv"
    if completed and raw_existing.exists():
        raw_keys = set()
        with raw_existing.open("r", newline="") as handle:
            for raw in csv.DictReader(handle):
                raw_keys.add(f"{raw['gamma']}_{raw['p_target']}_{raw['beta']}_{raw['x_target']}_{raw['replicate']}")
        completed &= raw_keys
    jobs = [j for j in jobs if f"{j[0]['gamma']}_{j[0]['p']}_{j[0]['beta']}_{j[0]['x']}_{j[1]}" not in completed]
    total_jobs = len(jobs)
    print(f"Fig. 5c: {total_jobs} tasks pending; workers={args.workers}; bootstrap={args.bootstrap}", flush=True)
    if total_jobs == 0:
        print("Fig. 5c: no new simulation tasks; rebuilding aggregates from raw CSV", flush=True)
    rows = []
    if args.workers > 1:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            iterator = pool.imap_unordered(run_one, jobs, chunksize=1)
            for done, row in enumerate(iterator, 1):
                rows.append(row)
                completed.add(f"{row['gamma']}_{row['p_target']}_{row['beta']}_{row['x_target']}_{row['replicate']}")
                if done == 1 or done % max(1, args.progress_every) == 0 or done == total_jobs:
                    print(f"Fig. 5c progress: {done}/{total_jobs} ({100*done/max(1,total_jobs):.1f}%)", flush=True)
                    checkpoint.write_text(json.dumps({"completed": sorted(completed)}, indent=2), encoding="utf-8")
    else:
        for done, job in enumerate(jobs, 1):
            row = run_one(job); rows.append(row)
            completed.add(f"{row['gamma']}_{row['p_target']}_{row['beta']}_{row['x_target']}_{row['replicate']}")
            if done == 1 or done % max(1, args.progress_every) == 0 or done == total_jobs:
                print(f"Fig. 5c progress: {done}/{total_jobs} ({100*done/max(1,total_jobs):.1f}%)", flush=True)
                checkpoint.write_text(json.dumps({"completed": sorted(completed)}, indent=2), encoding="utf-8")
    drift_path = out/"fig5c_local_drift.csv"
    append_mode = drift_path.exists() and bool(completed)
    with drift_path.open("a" if append_mode else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=DRIFT_FIELDS)
        if not append_mode:
            w.writeheader()
        w.writerows(rows)
    completed.update(f"{r['gamma']}_{r['p_target']}_{r['beta']}_{r['replicate']}" for r in rows)
    checkpoint.write_text(json.dumps({"completed": sorted(completed)}, indent=2), encoding="utf-8")
    # On resume, aggregate the complete raw file, not only rows produced in
    # the current invocation.
    rows = pd.read_csv(drift_path).to_dict("records")
    # Linear drift-vs-p estimate for each (beta,x), retaining all replicate-level rows.
    pc_rows=[]
    for g in gs:
     for b in bs:
      for x in xs:
        group=[r for r in rows if r["gamma"]==g and r["beta"]==b and r["x_target"]==x]; grouped={}
        for p in ps:
          grouped[p] = [float(r["drift"]) for r in group if abs(float(r["p_target"])-p)<1e-12]
        root, ci_low, ci_high = bootstrap_crossing(grouped, args.bootstrap, args.seed0 + int(g*1000) + int(b*100) + int(x*10))
        theory=p_c(x,b,10,args.omega,math.comb(10,4),args.c,args.R*(10*5/(10+5-1)),10*5/(10+5-1))
        pc_rows.append({"exp_id":args.exp_id,"gamma":g,"beta":b,"x_target":x,"p_c_theory":theory,"p_c_empirical":root,"p_c_ci95_low":ci_low,"p_c_ci95_high":ci_high,"num_p":len(ps),"num_replicates":args.replicates})
    with (out/"fig5c_pc_by_x.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=PC_FIELDS); w.writeheader(); w.writerows(pc_rows)
    r0=50/14; M=math.comb(10,4); pstar=p_static_star(args.omega,M,args.c,args.R*r0,r0); glob=[]
    reference_gamma = args.gamma if args.gamma in gs else max(gs)
    for b in bs:
      th,xstar,branch=p_ad_star(b,10,args.omega,M,args.c,args.R*r0,r0); vals=[r["p_c_empirical"] for r in pc_rows if r["gamma"]==reference_gamma and r["beta"]==b and np.isfinite(r["p_c_empirical"])]
      by_x = {}
      for x in xs:
        grouped = {}
        for p in ps:
          grouped[p] = [float(r["drift"]) for r in rows if r["gamma"]==reference_gamma and r["beta"]==b and r["x_target"]==x and abs(float(r["p_target"])-p)<1e-12]
        by_x[x] = grouped
      emp, ci_low, ci_high = bootstrap_global_envelope(by_x, args.bootstrap, args.seed0 + int(b*1000) + 991)
      glob.append({"exp_id":args.exp_id,"beta":b,"beta_c":beta_c(10),"p_static_star":pstar,"theory_p_ad_star":th,"theory_x_bottleneck":xstar,"theory_branch":branch,"empirical_p_ad_star":emp,"empirical_ci95_low":ci_low,"empirical_ci95_high":ci_high,"normalized_theory":th/pstar if pstar else float("nan"),"normalized_empirical":emp/pstar if pstar and np.isfinite(emp) else float("nan")})
    with (out/"fig5c_global_threshold.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=GLOBAL_FIELDS); w.writeheader(); w.writerows(glob)
    gamma_rows=[]
    for g in gs:
      for b in bs:
        vals=[r["p_c_empirical"] for r in pc_rows if r["gamma"]==g and r["beta"]==b and np.isfinite(r["p_c_empirical"])]
        th,_,_=p_ad_star(b,10,args.omega,M,args.c,args.R*r0,r0)
        gamma_rows.append({"exp_id":args.exp_id,"gamma":g,"beta":b,"x_target":"all","theory_p_c":th,"empirical_p_c":max(vals) if vals else float("nan"),"num_replicates":args.replicates})
    with (out/"fig5c_gamma_convergence.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=GAMMA_FIELDS); w.writeheader(); w.writerows(gamma_rows)
    print(f"Wrote {len(rows)} drift rows to {out}")

if __name__ == "__main__": main()
