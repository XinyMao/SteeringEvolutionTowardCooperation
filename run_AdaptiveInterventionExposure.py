#!/usr/bin/env python3
"""Fig. 5b: relocation-only steering-defector exposure.

Strategies are frozen for each replicate.  Only the steering mask evolves,
with total proposal rate ``gamma * p * N`` per MCS.  The output is a tidy CSV
that keeps raw exposure, normalized exposure, theory values, and simple
stationarity diagnostics together.
"""

from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from hsr_sim.adaptive import local_defector_exposure, relocation_proposal  # noqa: E402
from hsr_sim.engine import C, graph_to_neighbors, make_graph, steering_count  # noqa: E402


FIELDS = [
    "exp_id", "replicate", "seed", "N", "k", "n", "M", "p_target",
    "p_actual", "x_target", "x_actual", "beta", "gamma", "burn_in_mcs",
    "sample_mcs", "samples", "eta_mean", "eta_sd", "eta_ci95_low",
    "eta_ci95_high", "eta_baseline", "eta_normalized", "theory_eta",
    "theory_eta_normalized", "z", "a", "attempts", "accepted",
    "acceptance_rate", "tail_first_half_mean", "tail_second_half_mean",
    "tail_drift",
]


def theory_values(p: float, x: float, beta: float, k: int, z: float) -> tuple[float, float, float, float]:
    a = p + (1.0 - p) * x
    eta0 = p * (a - z) / ((1.0 - p) * a) if p > 0 and a > 0 else 0.0
    # Weak-targeting pair prediction from Adaptive_SI Eq. (83).
    weak = eta0 + beta * p * x * z * (a - z) / (k * a**3) if a > 0 else eta0
    if eta0 > 0:
        normalized = weak / eta0
    else:
        normalized = float("nan")
    return eta0, weak, normalized, a


def run_one(spec: tuple[dict, int]) -> dict:
    cfg, replicate = spec
    seed = int(cfg["seed0"] + 1_000_000 * replicate)
    rng = np.random.default_rng(seed)
    N, k, n = int(cfg["N"]), int(cfg["k"]), int(cfg["n"])
    p_target, x_target = float(cfg["p"]), float(cfg["x"])
    graph_seed = int(rng.integers(0, 2**32 - 1))
    import networkx as nx
    G = nx.random_regular_graph(k, N, seed=graph_seed)
    neighbors = graph_to_neighbors(G, N)
    nS = steering_count(N, p_target)
    steering = np.zeros(N, dtype=bool)
    steering[rng.choice(N, size=nS, replace=False)] = True
    strategy = np.zeros(N, dtype=np.int8)
    # Steering positions are expressed cooperators by definition in the
    # adaptive theory (A = S union C).  Relocation changes only the steering
    # identity, so both the source and target remain cooperative after a
    # swap; keeping these nodes C is essential for the beta=0 normalization.
    strategy[steering] = C
    free = np.flatnonzero(~steering)
    nC = int(round(x_target * free.size))
    if nC:
        strategy[rng.choice(free, size=nC, replace=False)] = C
    x_actual = float(np.mean(strategy[~steering] == C))
    p_actual = float(np.mean(steering))
    exposure = local_defector_exposure(strategy, neighbors)
    # z is unchanged during relocation because both source and target express C.
    A = strategy == C
    z = float(sum(A[u] and A[v] for u, v in G.edges()) * 2.0 / (k * N))
    eta0, theory_eta, theory_norm, a = theory_values(p_actual, x_actual, beta=float(cfg["beta"]), k=k, z=z)
    burn_attempts = int(round(float(cfg["burn_in_mcs"]) * cfg["gamma"] * p_actual * N))
    sample_attempts = int(round(float(cfg["sample_mcs"]) * cfg["gamma"] * p_actual * N))
    accepted = 0
    attempts = 0
    for _ in range(burn_attempts):
        result = relocation_proposal(strategy, steering, neighbors, float(cfg["beta"]), rng, exposure=exposure)
        accepted += int(result.accepted); attempts += 1
    sample_values: list[float] = []
    for _ in range(max(sample_attempts, 1)):
        result = relocation_proposal(strategy, steering, neighbors, float(cfg["beta"]), rng, exposure=exposure)
        accepted += int(result.accepted); attempts += 1
        # eta = E_SD / [k N (1-p)] = p/(1-p) * mean(h | steering).
        sample_values.append(
            float(p_actual / (1.0 - p_actual) * np.mean(exposure[steering]))
            if nS and p_actual < 1.0 else 0.0
        )
    vals = np.asarray(sample_values, dtype=float)
    mean = float(vals.mean())
    sd = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
    ci = 1.96 * sd / math.sqrt(vals.size) if vals.size else float("nan")
    half = max(1, vals.size // 2)
    first, second = float(vals[:half].mean()), float(vals[-half:].mean())
    return {
        "exp_id": cfg["exp_id"], "replicate": replicate, "seed": seed,
        "N": N, "k": k, "n": n, "M": math.comb(k, n - 1),
        "p_target": p_target, "p_actual": p_actual, "x_target": x_target,
        "x_actual": x_actual, "beta": cfg["beta"], "gamma": cfg["gamma"],
        "burn_in_mcs": cfg["burn_in_mcs"], "sample_mcs": cfg["sample_mcs"],
        "samples": vals.size, "eta_mean": mean, "eta_sd": sd,
        "eta_ci95_low": mean - ci, "eta_ci95_high": mean + ci,
        "eta_baseline": eta0, "eta_normalized": mean / eta0 if eta0 > 0 else float("nan"),
        "theory_eta": theory_eta, "theory_eta_normalized": theory_norm,
        "z": z, "a": a, "attempts": attempts, "accepted": accepted,
        "acceptance_rate": accepted / attempts if attempts else float("nan"),
        "tail_first_half_mean": first, "tail_second_half_mean": second,
        "tail_drift": second - first,
    }


def parse_grid(text: str, cast=float) -> list:
    return [cast(v.strip()) for v in text.split(",") if v.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results/fig5b_exposure/fig5b_summary.csv")
    parser.add_argument("--exp-id", default="Fig5bExposurePilot")
    parser.add_argument("--N", type=int, default=1000)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--n", type=int, default=4)
    parser.add_argument("--p-grid", default="0.02,0.05,0.10")
    parser.add_argument("--x-grid", default="0.2,0.4,0.6,0.8")
    parser.add_argument("--beta-grid", default="0,0.5,1,2")
    parser.add_argument("--gamma", type=float, default=5.0)
    parser.add_argument("--burn-in-mcs", type=float, default=50.0)
    parser.add_argument("--sample-mcs", type=float, default=50.0)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--seed0", type=int, default=20260904)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=20)
    args = parser.parse_args()
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    configs = [{"exp_id": args.exp_id, "N": args.N, "k": args.k, "n": args.n,
                "p": p, "x": x, "beta": beta, "gamma": args.gamma,
                "burn_in_mcs": args.burn_in_mcs, "sample_mcs": args.sample_mcs,
                "seed0": args.seed0}
               for p in parse_grid(args.p_grid) for x in parse_grid(args.x_grid)
               for beta in parse_grid(args.beta_grid)]
    tasks = [(cfg, rep) for cfg in configs for rep in range(args.replicates)]
    ctx = mp.get_context("spawn")
    rows = []
    rows = []
    total = len(tasks)
    print(f"Fig. 5b: {total} tasks pending; workers={args.workers}", flush=True)
    if args.workers > 1:
        with ctx.Pool(args.workers) as pool:
            for done, row in enumerate(pool.imap_unordered(run_one, tasks, chunksize=1), 1):
                rows.append(row)
                if done == 1 or done % max(1, args.progress_every) == 0 or done == total:
                    print(f"Fig. 5b progress: {done}/{total} ({100*done/max(1,total):.1f}%)", flush=True)
    else:
        for done, task in enumerate(tasks, 1):
            rows.append(run_one(task))
            if done == 1 or done % max(1, args.progress_every) == 0 or done == total:
                print(f"Fig. 5b progress: {done}/{total} ({100*done/max(1,total):.1f}%)", flush=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS); writer.writeheader(); writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
