"""Formal sparse-steering drift experiments.

This runner targets the sparse theory

    dx/dt = (k - 2) / (k - 1) * (1 - x) * [p - omega * M * c * x * (1 - r/r0)]

and estimates conditional short-window drift after a local relaxation period.
It writes MATLAB-readable raw and aggregated CSV/MAT files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from hsr_sim.engine import (  # noqa: E402
    C,
    OutputConfig,
    SimParams,
    append_csv,
    graph_to_neighbors,
    load_checkpoint,
    make_graph,
    precompute_neighbor_combos,
    steering_count,
    update_checkpoint,
    write_metadata,
    db_update,
)


RAW_FIELDS = [
    "exp_id", "run_id", "replicate_id", "seed", "accepted",
    "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
    "p_s_target", "num_steering", "p_s_actual", "P_actual", "p_star",
    "p_over_p_star", "x0_target", "x_start", "x_after_relax", "x_end",
    "x_window", "T_relax_mcs", "tau_mcs", "delta_x", "drift",
    "G_hat", "G_theory_start", "topology", "topology_param",
]

AGG_FIELDS = [
    "exp_id", "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
    "p_s_target", "mean_p_s_actual", "mean_P_actual", "p_star",
    "p_over_p_star", "x0_target", "x_window", "T_relax_mcs", "tau_mcs",
    "topology", "topology_param", "num_runs", "num_accepted",
    "acceptance_rate", "mean_x_after_relax", "mean_drift", "std_drift",
    "sem_drift", "ci95_low", "ci95_high", "mean_G_hat", "sem_G_hat",
    "mean_G_theory_start",
]

THRESHOLD_FIELDS = [
    "exp_id", "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
    "topology", "topology_param", "x0_target", "T_relax_mcs", "tau_mcs",
    "p_star_theory", "p_star_hat", "P_star_theory", "P_star_hat",
    "slope", "intercept", "r_squared", "num_p_values", "num_accepted",
]

UNIFORM_THRESHOLD_FIELDS = [
    "exp_id", "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
    "topology", "topology_param", "T_relax_mcs", "tau_mcs",
    "p_star_theory", "p_star_hat_uniform", "P_star_theory",
    "P_star_hat_uniform", "worst_x0_target", "num_x0_values", "num_accepted",
]

POOLED_THRESHOLD_FIELDS = [
    "exp_id", "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
    "topology", "topology_param", "T_relax_mcs", "tau_mcs",
    "p_star_theory", "p_star_hat_x1", "P_star_theory", "P_star_hat_x1",
    "intercept", "coef_p", "coef_x", "r_squared", "num_cells",
    "num_accepted", "x_min", "x_max",
]


def parse_csv_floats(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def default_r0(k: int, n: int) -> float:
    return float(k * n / (k + n - 1))


def default_M(k: int, n: int) -> int:
    return int(math.comb(k, n - 1))


def _is_number(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def parse_cases(values: Sequence[str]) -> list[dict[str, float | int | str | None]]:
    """Parse theory cases.

    Supported forms:
        k,n
        k,n,topology
        k,n,r0,topology
        k,n,r0,M,topology

    r0 and M are intentionally allowed as an explicit pair because future
    theory variants may change both together.
    """
    cases = []
    for raw in values:
        parts = [p.strip() for p in raw.split(",")]
        if len(parts) not in (2, 3, 4, 5):
            raise ValueError("--case must be k,n[,topology] or k,n,r0[,M],topology")
        k, n = int(parts[0]), int(parts[1])
        r0: float | None = None
        M: int | None = None
        topology = "random_regular"
        if len(parts) == 3:
            if _is_number(parts[2]):
                r0 = float(parts[2])
            else:
                topology = parts[2]
        elif len(parts) == 4:
            r0 = float(parts[2])
            topology = parts[3]
        elif len(parts) == 5:
            r0 = float(parts[2])
            M = int(parts[3])
            topology = parts[4]
        cases.append({"k": k, "n": n, "r0": r0, "M": M, "topology": topology})
    return cases


def resolved_M(k: int, n: int, m_mode: str, m_value: int | None) -> int:
    total = math.comb(k, n - 1)
    if m_mode == "all":
        return total
    if m_value is None or not (1 <= m_value <= total):
        raise ValueError(f"--M must be in [1, {total}] when --M-mode sampled")
    return int(m_value)


def exact_state_at_x0(params: SimParams, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    N = params.N
    nS = steering_count(N, params.p_s)
    nodes = np.arange(N)
    steering_nodes = rng.choice(nodes, size=nS, replace=False) if nS > 0 else np.array([], dtype=np.int64)
    is_steering = np.zeros(N, dtype=bool)
    is_steering[steering_nodes] = True
    strategy = np.zeros(N, dtype=np.int8)
    strategy[is_steering] = C
    free_nodes = nodes[~is_steering]
    nC_free = int(math.floor(params.x0 * len(free_nodes)))
    if nC_free > 0:
        strategy[rng.choice(free_nodes, size=nC_free, replace=False)] = C
    return strategy, is_steering


def x_free(strategy: np.ndarray, is_steering: np.ndarray) -> float:
    free = ~is_steering
    return float(np.mean(strategy[free] == C)) if np.any(free) else float("nan")


def sparse_run_id(exp_id: str, params: SimParams, r0: float, epsilon: float, ratio: float, rep: int, relax: int, tau: int) -> str:
    def fmt(x: float) -> str:
        return f"{x:.6g}".replace(".", "p").replace("-", "m")
    return (
        f"{exp_id}_N{params.N}_k{params.k}_n{params.n}_M{params.resolved_M()}_"
        f"r0{fmt(r0)}_R{fmt(params.r / r0)}_eps{fmt(epsilon)}_"
        f"q{fmt(ratio)}_x{fmt(params.x0)}_{params.topology}_rel{relax}_tau{tau}_rep{rep:06d}"
    )


def build_params(args: argparse.Namespace) -> list[tuple[SimParams, float, float, float, float]]:
    cases = parse_cases(args.case)
    R_values = parse_csv_floats(args.R_grid)
    eps_values = parse_csv_floats(args.epsilon_grid)
    ratios = parse_csv_floats(args.p_ratio_grid)
    params = []
    for case in cases:
        k = int(case["k"])
        n = int(case["n"])
        r0 = default_r0(k, n) if case["r0"] is None else float(case["r0"])
        topology = str(case["topology"])
        case_M = case["M"]
        if case_M is not None:
            total_M = default_M(k, n)
            if int(case_M) > total_M:
                raise ValueError(
                    f"Explicit case M={case_M} exceeds C(k,n-1)={total_M}; "
                    "the current engine can only simulate all groups or a sampled subset."
                )
            M = int(case_M)
            M_mode = "all" if M == total_M else "sampled"
            M_value = None if M_mode == "all" else M
        else:
            M = resolved_M(k, n, args.M_mode, args.M)
            M_mode = args.M_mode
            M_value = args.M
        for R in R_values:
            p_star = None
            for epsilon in eps_values:
                omega = epsilon / (M * args.c)
                r = R * r0
                p_star = epsilon * (1.0 - R)
                N = args.N
                if args.auto_N:
                    if p_star <= 0:
                        raise ValueError("auto_N requires R < 1.")
                    N = max(args.N, int(math.ceil(args.min_steering_at_threshold / p_star)))
                    if topology in {"square_lattice", "moore_lattice"}:
                        L = int(math.ceil(math.sqrt(N)))
                        N = L * L
                if args.max_N > 0 and N > args.max_N:
                    raise ValueError(f"Required N={N} exceeds --max-N={args.max_N}.")
                for ratio in ratios:
                    p = p_star * ratio if p_star is not None else 0.0
                    if p > args.max_p:
                        continue
                    for x0 in parse_csv_floats(args.x_grid):
                        params.append((
                            SimParams(
                                N=N, k=k, n=n, r=r, omega=omega, c=args.c,
                                M_mode=M_mode, M=M_value, p_s=p, x0=x0,
                                topology=topology, topology_param=args.topology_param,
                                T_mcs=args.tau_mcs, record_interval_mcs=1,
                                seed=args.seed0,
                            ),
                            r0, R, epsilon, ratio,
                        ))
    return params


def run_one(task: tuple[SimParams, float, float, float, float, int, int, int, float, str, int]) -> dict:
    params, r0, R, epsilon, ratio, rep, T_relax, tau, x_window, exp_id, seed = task
    rng = np.random.default_rng(seed)
    params = SimParams(**{**asdict(params), "seed": seed})
    G = make_graph(params, rng)
    neighbors = graph_to_neighbors(G, params.N)
    combos = None if params.M_mode == "all" else precompute_neighbor_combos(neighbors, params.n)
    strategy, is_steering = exact_state_at_x0(params, rng)

    x_start = x_free(strategy, is_steering)
    for _ in range(T_relax * params.N):
        db_update(strategy, is_steering, neighbors, combos, params, rng)
    x_after_relax = x_free(strategy, is_steering)
    accepted = abs(x_after_relax - params.x0) <= x_window

    x_end = float("nan")
    dx = float("nan")
    drift = float("nan")
    G_hat = float("nan")
    if accepted:
        for _ in range(tau * params.N):
            db_update(strategy, is_steering, neighbors, combos, params, rng)
        x_end = x_free(strategy, is_steering)
        dx = x_end - x_after_relax
        drift = dx / float(tau)
        if x_after_relax < 1.0:
            G_hat = (params.k - 1.0) / (params.k - 2.0) * drift / (1.0 - x_after_relax)

    p_star = epsilon * (1.0 - R)
    p_actual = float(np.mean(is_steering))
    P_actual = p_actual / epsilon if epsilon > 0 else float("nan")
    G_theory = p_actual - epsilon * x_after_relax * (1.0 - R)
    run_id = sparse_run_id(exp_id, params, r0, epsilon, ratio, rep, T_relax, tau)
    return {
        "exp_id": exp_id,
        "run_id": run_id,
        "replicate_id": rep,
        "seed": seed,
        "accepted": bool(accepted),
        "N": params.N,
        "k": params.k,
        "n": params.n,
        "M": params.resolved_M(),
        "r0": r0,
        "R": R,
        "r": params.r,
        "omega": params.omega,
        "epsilon": epsilon,
        "c": params.c,
        "p_s_target": params.p_s,
        "num_steering": int(np.sum(is_steering)),
        "p_s_actual": p_actual,
        "P_actual": P_actual,
        "p_star": p_star,
        "p_over_p_star": ratio,
        "x0_target": params.x0,
        "x_start": x_start,
        "x_after_relax": x_after_relax,
        "x_end": x_end,
        "x_window": x_window,
        "T_relax_mcs": T_relax,
        "tau_mcs": tau,
        "delta_x": dx,
        "drift": drift,
        "G_hat": G_hat,
        "G_theory_start": G_theory,
        "topology": params.topology,
        "topology_param": params.topology_param,
    }


def append_many(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    rows = list(rows)
    if not rows:
        return
    new_file = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerows(rows)


def aggregate_outputs(out: Path, exp_id: str, save_mat: bool = True) -> None:
    import pandas as pd

    raw_path = out / "sparse_drift_raw.csv"
    if not raw_path.exists():
        return
    df = pd.read_csv(raw_path).drop_duplicates(subset=["run_id"], keep="last")
    keys = [
        "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
        "p_s_target", "p_over_p_star", "x0_target", "x_window",
        "T_relax_mcs", "tau_mcs", "topology", "topology_param",
    ]
    accepted = df[df["accepted"] == True].copy()  # noqa: E712
    grouped_all = df.groupby(keys, as_index=False).agg(num_runs=("run_id", "count"))
    if accepted.empty:
        return
    grouped_acc = accepted.groupby(keys, as_index=False).agg(
        num_accepted=("run_id", "count"),
        mean_p_s_actual=("p_s_actual", "mean"),
        mean_P_actual=("P_actual", "mean"),
        p_star=("p_star", "mean"),
        mean_x_after_relax=("x_after_relax", "mean"),
        mean_drift=("drift", "mean"),
        std_drift=("drift", "std"),
        mean_G_hat=("G_hat", "mean"),
        std_G_hat=("G_hat", "std"),
        mean_G_theory_start=("G_theory_start", "mean"),
    )
    agg = grouped_acc.merge(grouped_all, on=keys, how="left")
    agg["acceptance_rate"] = agg["num_accepted"] / agg["num_runs"]
    agg["sem_drift"] = agg["std_drift"] / np.sqrt(agg["num_accepted"])
    agg["ci95_low"] = agg["mean_drift"] - 1.96 * agg["sem_drift"]
    agg["ci95_high"] = agg["mean_drift"] + 1.96 * agg["sem_drift"]
    agg["sem_G_hat"] = agg["std_G_hat"] / np.sqrt(agg["num_accepted"])
    agg.insert(0, "exp_id", exp_id)
    agg.to_csv(out / "sparse_drift_aggregate.csv", index=False, columns=AGG_FIELDS)

    threshold_rows = []
    tkeys = [
        "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
        "topology", "topology_param", "x0_target", "T_relax_mcs", "tau_mcs",
    ]
    for tvals, g in agg.groupby(tkeys, dropna=False):
        if len(g) < 3:
            continue
        # Several target p values can collapse onto the same actual p because
        # the steering set size is integer-valued. Fit against distinct actual
        # p values, using accepted-sample counts as weights.
        g_fit = g.groupby("mean_p_s_actual", as_index=False).apply(
            lambda h: pd.Series({
                "mean_G_hat": np.average(h["mean_G_hat"], weights=h["num_accepted"]),
                "num_accepted": h["num_accepted"].sum(),
            }),
            include_groups=False,
        )
        x = g_fit["mean_p_s_actual"].to_numpy(dtype=float)
        y = g_fit["mean_G_hat"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if np.sum(ok) < 3:
            continue
        slope, intercept = np.polyfit(x[ok], y[ok], 1)
        yhat = slope * x[ok] + intercept
        ss_res = float(np.sum((y[ok] - yhat) ** 2))
        ss_tot = float(np.sum((y[ok] - np.mean(y[ok])) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        p_hat = -intercept / slope if slope != 0 else float("nan")
        row = dict(zip(tkeys, tvals))
        p_theory_x = (
            row["epsilon"] * row["x0_target"] * (1.0 - row["R"])
        )
        row.update({
            "exp_id": exp_id,
            "p_star_theory": p_theory_x,
            "p_star_hat": p_hat,
            "P_star_theory": (
                p_theory_x / row["epsilon"]
                if row["epsilon"] > 0
                else float("nan")
            ),
            "P_star_hat": p_hat / row["epsilon"] if row["epsilon"] > 0 else float("nan"),
            "slope": slope,
            "intercept": intercept,
            "r_squared": r_squared,
            "num_p_values": int(np.sum(ok)),
            "num_accepted": int(g_fit["num_accepted"].sum()),
        })
        threshold_rows.append(row)
    if threshold_rows:
        threshold_df = pd.DataFrame(threshold_rows)
        threshold_df.to_csv(out / "sparse_thresholds.csv", index=False, columns=THRESHOLD_FIELDS)
        uniform_rows = []
        ukeys = [
            "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
            "topology", "topology_param", "T_relax_mcs", "tau_mcs",
        ]
        for uvals, ug in threshold_df.groupby(ukeys, dropna=False):
            finite = ug[np.isfinite(ug["p_star_hat"])]
            if finite.empty:
                continue
            worst = finite.loc[finite["p_star_hat"].idxmax()]
            row = dict(zip(ukeys, uvals))
            p_theory_uniform = row["epsilon"] * (1.0 - row["R"])
            row.update({
                "exp_id": exp_id,
                "p_star_theory": p_theory_uniform,
                "p_star_hat_uniform": float(worst["p_star_hat"]),
                "P_star_theory": (
                    p_theory_uniform / row["epsilon"]
                    if row["epsilon"] > 0
                    else float("nan")
                ),
                "P_star_hat_uniform": float(worst["P_star_hat"]),
                "worst_x0_target": float(worst["x0_target"]),
                "num_x0_values": int(len(finite)),
                "num_accepted": int(finite["num_accepted"].sum()),
            })
            uniform_rows.append(row)
        if uniform_rows:
            pd.DataFrame(uniform_rows).to_csv(
                out / "sparse_uniform_thresholds.csv",
                index=False,
                columns=UNIFORM_THRESHOLD_FIELDS,
            )

    pooled_rows = []
    pkeys = [
        "N", "k", "n", "M", "r0", "R", "r", "omega", "epsilon", "c",
        "topology", "topology_param", "T_relax_mcs", "tau_mcs",
    ]
    for pvals, pg in agg.groupby(pkeys, dropna=False):
        x = pg["mean_x_after_relax"].to_numpy(dtype=float)
        ps = pg["mean_p_s_actual"].to_numpy(dtype=float)
        y = pg["mean_G_hat"].to_numpy(dtype=float)
        w = pg["num_accepted"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(ps) & np.isfinite(y) & np.isfinite(w) & (w > 0)
        if np.sum(ok) < 6:
            continue
        X = np.column_stack([np.ones(np.sum(ok)), ps[ok], x[ok]])
        sw = np.sqrt(w[ok])
        coef = np.linalg.lstsq(X * sw[:, None], y[ok] * sw, rcond=None)[0]
        intercept, coef_p, coef_x = coef
        yhat = X @ coef
        ybar = np.average(y[ok], weights=w[ok])
        ss_res = float(np.sum(w[ok] * (y[ok] - yhat) ** 2))
        ss_tot = float(np.sum(w[ok] * (y[ok] - ybar) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        row = dict(zip(pkeys, pvals))
        p_hat_x1 = -(intercept + coef_x) / coef_p if coef_p != 0 else float("nan")
        row.update({
            "exp_id": exp_id,
            "p_star_theory": float(pg["p_star"].mean()),
            "p_star_hat_x1": p_hat_x1,
            "P_star_theory": float(pg["p_star"].mean() / row["epsilon"]),
            "P_star_hat_x1": p_hat_x1 / row["epsilon"] if row["epsilon"] > 0 else float("nan"),
            "intercept": intercept,
            "coef_p": coef_p,
            "coef_x": coef_x,
            "r_squared": r_squared,
            "num_cells": int(np.sum(ok)),
            "num_accepted": int(np.sum(w[ok])),
            "x_min": float(np.min(x[ok])),
            "x_max": float(np.max(x[ok])),
        })
        pooled_rows.append(row)
    pooled_df = pd.DataFrame(pooled_rows)
    if pooled_rows:
        pooled_df.to_csv(out / "sparse_pooled_thresholds.csv", index=False, columns=POOLED_THRESHOLD_FIELDS)

    if not save_mat:
        return

    try:
        from scipy.io import savemat
        payload = {"raw": df.to_records(index=False), "aggregate": agg.to_records(index=False)}
        if threshold_rows:
            payload["thresholds"] = threshold_df.to_records(index=False)
            if uniform_rows:
                payload["uniform_thresholds"] = pd.DataFrame(uniform_rows).to_records(index=False)
        if pooled_rows:
            payload["pooled_thresholds"] = pooled_df.to_records(index=False)
        savemat(out / "sparse_drift.mat", payload)
    except Exception as exc:
        print(f"MAT export skipped: {exc}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Formal sparse-steering conditional drift experiment.")
    parser.add_argument("--out-dir", default="results/exp03_sparse_drift")
    parser.add_argument("--exp-id", default="E03sparse")
    parser.add_argument(
        "--case",
        action="append",
        default=None,
        help=(
            "Repeatable theory case. Forms: k,n; k,n,topology; "
            "k,n,r0,topology; k,n,r0,M,topology."
        ),
    )
    parser.add_argument("--N", type=int, default=400)
    parser.add_argument("--auto-N", action="store_true", help="Increase N so N*p_star reaches --min-steering-at-threshold.")
    parser.add_argument("--max-N", type=int, default=0, help="Abort auto-N above this value; 0 disables.")
    parser.add_argument("--min-steering-at-threshold", type=int, default=20)
    parser.add_argument("--c", type=float, default=1.0)
    parser.add_argument("--M-mode", choices=["all", "sampled"], default="all")
    parser.add_argument("--M", type=int, default=None)
    parser.add_argument("--R-grid", default="0.3,0.5,0.7,0.85,0.95")
    parser.add_argument("--epsilon-grid", default="0.0025,0.005,0.01")
    parser.add_argument("--p-ratio-grid", default="0,0.5,0.75,0.9,1.0,1.1,1.25,1.5")
    parser.add_argument("--x-grid", default="0.2,0.5,0.8,0.95")
    parser.add_argument("--max-p", type=float, default=0.5)
    parser.add_argument("--topology-param", type=float, default=1.0)
    parser.add_argument("--T-relax-mcs", type=int, default=20)
    parser.add_argument("--tau-mcs", type=int, default=3)
    parser.add_argument("--x-window", type=float, default=0.03)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed0", type=int, default=20260720)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument(
        "--mp-start-method",
        choices=["auto", "fork", "spawn", "forkserver"],
        default="auto",
        help="Multiprocessing start method. auto uses fork on Linux and spawn elsewhere.",
    )
    parser.add_argument("--progress-every", type=int, default=200)
    parser.add_argument("--no-mat", action="store_true", help="Skip MAT export after aggregation; CSV outputs are still written.")
    parser.add_argument("--no-checkpoint", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.case is None:
        args.case = ["10,5,random_regular"]
    if args.smoke:
        args.out_dir = "/private/tmp/hsr_sparse_drift_smoke"
        args.N = 60
        args.case = ["10,5,random_regular"]
        args.R_grid = "0.6"
        args.epsilon_grid = "0.01"
        args.p_ratio_grid = "0,1.0,1.5"
        args.x_grid = "0.5,0.8"
        args.T_relax_mcs = 1
        args.tau_mcs = 1
        args.x_window = 0.2
        args.replicates = 2
        args.workers = 1

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cfg = OutputConfig(out_dir=args.out_dir, exp_id=args.exp_id, save_timeseries=False, use_checkpoint=not args.no_checkpoint)
    params = build_params(args)
    if not params:
        raise RuntimeError("No parameter combinations to run.")
    write_metadata(out, cfg, params[0][0], note="Formal sparse-steering conditional drift experiment.")
    with (out / "sparse_design.json").open("w", encoding="utf-8") as f:
        json.dump({"created_at": time.strftime("%Y-%m-%d %H:%M:%S"), "args": vars(args)}, f, indent=2)

    completed = load_checkpoint(out / "checkpoint.json") if cfg.use_checkpoint else set()
    tasks = []
    pidx = 0
    for base, r0, R, epsilon, ratio in params:
        for rep in range(args.replicates):
            seed = int(args.seed0 + 1000000 * pidx + rep)
            run_id = sparse_run_id(args.exp_id, base, r0, epsilon, ratio, rep, args.T_relax_mcs, args.tau_mcs)
            if run_id not in completed:
                tasks.append((base, r0, R, epsilon, ratio, rep, args.T_relax_mcs, args.tau_mcs, args.x_window, args.exp_id, seed))
        pidx += 1

    total = len(tasks)
    print(f"Sparse drift: starting {total} runs")
    print(f"Output directory: {out}")
    done = 0
    raw_path = out / "sparse_drift_raw.csv"
    checkpoint_path = out / "checkpoint.json"

    if args.workers <= 1:
        for task in tasks:
            row = run_one(task)
            append_csv(raw_path, row, RAW_FIELDS)
            completed.add(row["run_id"])
            done += 1
            if args.progress_every > 0 and (done == total or done % args.progress_every == 0):
                print(f"Sparse drift: {done}/{total} runs complete")
                if cfg.use_checkpoint:
                    update_checkpoint(checkpoint_path, completed)
    else:
        available_methods = mp.get_all_start_methods()
        if args.mp_start_method == "auto":
            start_method = "fork" if sys.platform.startswith("linux") and "fork" in available_methods else "spawn"
        else:
            start_method = args.mp_start_method
        if start_method not in available_methods:
            raise ValueError(
                f"Multiprocessing start method {start_method!r} is unavailable; "
                f"choose from {available_methods}."
            )
        ctx = mp.get_context(start_method)
        print(f"Sparse drift: using {args.workers} workers ({start_method})")
        with ctx.Pool(processes=args.workers) as pool:
            batch = []
            for row in pool.imap_unordered(run_one, tasks, chunksize=args.chunk_size):
                batch.append(row)
                completed.add(row["run_id"])
                done += 1
                if len(batch) >= args.chunk_size:
                    append_many(raw_path, batch, RAW_FIELDS)
                    batch.clear()
                if args.progress_every > 0 and (done == total or done % args.progress_every == 0):
                    if batch:
                        append_many(raw_path, batch, RAW_FIELDS)
                        batch.clear()
                    print(f"Sparse drift: {done}/{total} runs complete")
                    if cfg.use_checkpoint:
                        update_checkpoint(checkpoint_path, completed)
            if batch:
                append_many(raw_path, batch, RAW_FIELDS)

    if cfg.use_checkpoint:
        update_checkpoint(checkpoint_path, completed)
    aggregate_outputs(out, args.exp_id, save_mat=not args.no_mat)
    print("Sparse drift completed.")


if __name__ == "__main__":
    main()
