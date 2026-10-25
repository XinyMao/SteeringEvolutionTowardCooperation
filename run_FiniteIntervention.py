"""Finite-steering Monte Carlo experiments.

This runner is designed for the manuscript subsection
"Finite steering creates state-dependent intervention opportunities".

Implemented experiment modes:
  drift     State-conditioned short-time drift, Fig. 3a.
  boundary  Theory-guided local r*(x,p) sign diagram, Fig. 3b.
  global    Long-time high-cooperation probability, Fig. 3c.
  exposure  Fixed-budget placement/exposure mechanism, Fig. 3d.

The runner writes raw trajectory rows plus aggregate CSV/MAT files.  Rows retain
network_id, steering_id, snapshot_id, and trajectory_id so uncertainty can be
computed with hierarchical resampling instead of treating all trajectories as
exchangeable.
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
from hsr_sim.finite import (  # noqa: E402
    classical_r0,
    constrained_neutral_relaxation,
    eta0_pair_closure,
    finite_drift_theory,
    finite_local_stats,
    initialize_exact_state,
    local_r_star,
    run_short_drift_from_snapshot,
    select_steering_nodes,
)


RAW_FIELDS = [
    "exp_id", "run_id", "mode", "accepted",
    "network_id", "steering_id", "snapshot_id", "trajectory_id",
    "network_seed", "steering_seed", "strategy_seed", "relax_seed", "trajectory_seed",
    "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
    "p_target", "num_steering", "p_actual",
    "x_target", "x_initial", "x_after_relax", "x_final",
    "delta_t_mcs", "relax_mcs", "delta_x", "drift_empirical",
    "drift_theoretical", "drift_theoretical_closure_eta",
    "r_star_local", "r_star_closure_eta", "r_star_global_closure_eta",
    "eta0", "eta0_closure", "q_S_given_D", "q_C_given_D", "q_C_given_S", "q_S_given_S",
    "E_DS", "E_SS", "density_DS", "density_SS",
    "placement_protocol", "placement_bias", "metropolis_steps",
    "topology", "topology_param",
]

AGG_FIELDS = [
    "exp_id", "mode", "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
    "p_target", "mean_p_actual", "x_target", "delta_t_mcs", "relax_mcs",
    "placement_protocol", "placement_bias", "topology", "topology_param",
    "num_networks", "num_steering_sets", "num_snapshots", "num_trajectories",
    "mean_x_after_relax", "mean_eta0", "mean_q_S_given_D",
    "mean_drift_empirical", "std_drift_empirical", "sem_naive",
    "ci95_low_boot", "ci95_high_boot",
    "mean_drift_theoretical", "mean_drift_theoretical_closure_eta",
    "mean_r_star_local", "mean_r_star_closure_eta", "sign_accuracy",
    "rmse", "nrmse", "mean_E_DS", "mean_E_SS",
]

BOUNDARY_FIELDS = [
    "exp_id", "N", "k", "n", "M", "M_mode", "p_target", "mean_p_actual",
    "x_target", "omega", "c", "topology", "topology_param",
    "r_star_theory_mean", "r_star_hat", "ci95_low_boot", "ci95_high_boot",
    "slope", "intercept", "r_squared", "num_r_values", "num_trajectories",
    "mae_vs_classical_r0", "abs_error_vs_theory",
]

GLOBAL_RAW_FIELDS = [
    "exp_id", "run_id", "mode", "accepted",
    "network_id", "steering_id", "snapshot_id", "trajectory_id",
    "network_seed", "steering_seed", "strategy_seed", "relax_seed", "trajectory_seed",
    "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
    "p_target", "num_steering", "p_actual",
    "x_target", "x_initial", "x_after_relax", "x_final", "mean_x_tail",
    "global_T_mcs", "global_record_interval_mcs", "global_tail_fraction",
    "high_coop_threshold", "high_coop",
    "r_star_global_closure_eta", "eta0", "eta0_closure",
    "q_S_given_D", "q_C_given_D", "q_C_given_S", "q_S_given_S",
    "E_DS", "E_SS", "density_DS", "density_SS",
    "placement_protocol", "placement_bias", "metropolis_steps",
    "topology", "topology_param",
    "x0_requested", "x0_actual", "initial_state_seed", "valid_theory_state",
    "x_domain_min", "x_domain_max", "high_coop_indicator",
    "high_coop_probability_group", "final_x_free",
    "hitting_time_to_high_cooperation", "reached_high_cooperation",
    "initial_eta0", "initial_q_S_given_D", "initial_q_C_given_D",
    "initial_q_C_given_S", "initial_q_S_given_S",
]

GLOBAL_AGG_FIELDS = [
    "exp_id", "mode", "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
    "p_target", "mean_p_actual", "global_T_mcs", "global_record_interval_mcs",
    "global_tail_fraction", "high_coop_threshold",
    "placement_protocol", "placement_bias", "topology", "topology_param",
    "num_networks", "num_steering_sets", "num_snapshots", "num_trajectories",
    "num_x0_values", "high_coop_probability", "mean_final_x_free",
    "mean_tail_x_free", "sem_tail_x_free", "worst_x0_target",
    "min_high_coop_probability_by_x", "r_star_global_closure_eta",
    "mean_eta0", "mean_eta0_closure", "mean_E_DS", "mean_E_SS",
]


def parse_csv_floats(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_csv_strings(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def parse_cases(values: Sequence[str] | None, default_k: int, default_n: int, default_topology: str) -> list[dict]:
    if not values:
        return [{"k": default_k, "n": default_n, "topology": default_topology}]
    cases = []
    for raw in values:
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        if len(parts) not in (2, 3):
            raise ValueError("--case must be k,n or k,n,topology")
        cases.append({
            "k": int(parts[0]),
            "n": int(parts[1]),
            "topology": parts[2] if len(parts) == 3 else default_topology,
        })
    return cases


def fmt(x: float) -> str:
    return f"{x:.6g}".replace(".", "p").replace("-", "m")


def run_id_for(row: dict) -> str:
    return (
        f"{row['exp_id']}_{row['mode']}_N{row['N']}_k{row['k']}_n{row['n']}_"
        f"p{fmt(row['p_target'])}_x{fmt(row['x_target'])}_r{fmt(row['r'])}_"
        f"{row['placement_protocol']}{fmt(row['placement_bias'])}_"
        f"g{row['network_id']:04d}_s{row['steering_id']:04d}_"
        f"z{row['snapshot_id']:04d}_tr{row['trajectory_id']:04d}"
    )


def task_id_for(row: dict) -> str:
    return (
        f"{row['exp_id']}_{row['mode']}_N{row['N']}_k{row['k']}_n{row['n']}_"
        f"p{fmt(row['p_target'])}_x{fmt(row['x_target'])}_r{fmt(row['r'])}_"
        f"{row['placement_protocol']}{fmt(row['placement_bias'])}_"
        f"g{row['network_id']:04d}_s{row['steering_id']:04d}_z{row['snapshot_id']:04d}"
    )


def load_completed_run_ids(path: Path) -> set[str]:
    """Recover completed trajectory IDs directly from an existing raw CSV.

    The raw file is authoritative.  In particular, a checkpoint created by a
    smaller pilot must not suppress trajectory IDs added by a later full run.
    This also makes resume safe if a job is interrupted after rows are flushed
    but before checkpoint.json is updated, or if a raw file is incomplete.
    """
    if not path.exists():
        return set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "run_id" not in (reader.fieldnames or []):
            raise RuntimeError(f"Existing raw CSV has no run_id column: {path}")
        return {str(row["run_id"]) for row in reader if row.get("run_id")}


def append_many(path: Path, rows: Iterable[dict], fields: Sequence[str]) -> None:
    rows = list(rows)
    if not rows:
        return
    new_file = not path.exists()
    if not new_file:
        with path.open("r", newline="", encoding="utf-8") as f:
            existing = next(csv.reader(f), [])
        if existing != list(fields):
            raise RuntimeError(
                f"Existing CSV schema at {path} does not match this runner. "
                "Use a fresh --out-dir or archive the old output before resuming."
            )
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerows(rows)


def run_global_from_snapshot(
    params: SimParams,
    neighbors: list[np.ndarray],
    combos: list[np.ndarray],
    strategy0: np.ndarray,
    is_steering: np.ndarray,
    T_mcs: float,
    record_interval_mcs: float,
    tail_fraction: float,
    high_threshold: float,
    seed: int,
) -> tuple[float, float, bool, float, bool]:
    rng = np.random.default_rng(seed)
    strategy = np.array(strategy0, copy=True)
    steps = int(round(T_mcs * params.N))
    record_every = max(1, int(round(record_interval_mcs * params.N)))
    xs: list[float] = []
    hitting_time = float("nan")
    for step in range(1, steps + 1):
        db_update(strategy, is_steering, neighbors, combos, params, rng)
        if step % record_every == 0 or step == steps:
            current_x = float(np.mean(strategy[~is_steering] == 1))
            xs.append(current_x)
            if not np.isfinite(hitting_time) and current_x >= high_threshold:
                hitting_time = float(step / params.N)
    if not xs:
        xs.append(float(np.mean(strategy[~is_steering] == 1)))
    tail_n = max(1, int(math.ceil(tail_fraction * len(xs))))
    mean_tail = float(np.mean(xs[-tail_n:]))
    final_x = float(xs[-1])
    reached = bool(np.isfinite(hitting_time))
    return final_x, mean_tail, bool(mean_tail >= high_threshold), hitting_time, reached


def finite_global_r_star_closure_eta(params: SimParams, p_actual: float, x_grid: Sequence[float]) -> float:
    values = []
    for x in x_grid:
        eta0 = eta0_pair_closure(x, p_actual, params.k)
        val = local_r_star(x, p_actual, eta0, params.k, params.n, params.resolved_M(), params.c, params.omega)
        if np.isfinite(val):
            values.append(val)
    return float(np.max(values)) if values else float("nan")


def make_bias_strategy(N: int, x_target: float, rng: np.random.Generator) -> np.ndarray:
    strategy = np.zeros(N, dtype=np.int8)
    n_c = int(round(x_target * N))
    if n_c > 0:
        strategy[rng.choice(np.arange(N), size=n_c, replace=False)] = 1
    return strategy


def run_snapshot_task(task: tuple[SimParams, dict, int, int, int, Sequence[int], float, float, str]) -> list[dict]:
    base, spec, network_id, steering_id, snapshot_id, trajectory_ids, relax_mcs, tau_mcs, exp_id = task
    network_seed = int(spec["network_seed"])
    steering_seed = int(spec["steering_seed"])
    strategy_seed = int(spec["strategy_seed"])
    relax_seed = int(spec["relax_seed"])
    mode = str(spec["mode"])
    protocol = str(spec["placement_protocol"])
    placement_bias = float(spec["placement_bias"])
    metropolis_steps = int(spec["metropolis_steps"])

    rng_g = np.random.default_rng(network_seed)
    G = make_graph(base, rng_g)
    neighbors = graph_to_neighbors(G, base.N)
    combos = precompute_neighbor_combos(neighbors, base.n)

    num_steer = steering_count(base.N, base.p_s)
    bias_strategy = None
    if mode == "exposure" or protocol != "random":
        bias_strategy = make_bias_strategy(base.N, base.x0, np.random.default_rng(strategy_seed))
    steering_nodes = select_steering_nodes(
        G,
        bias_strategy,
        num_steer,
        np.random.default_rng(steering_seed),
        protocol=protocol,
        bias=placement_bias,
        metropolis_steps=metropolis_steps,
    )
    strategy, is_steering, target_count = initialize_exact_state(
        base.N,
        steering_nodes,
        base.x0,
        np.random.default_rng(strategy_seed),
    )
    x_initial = float(np.mean(strategy[~is_steering] == 1))

    constrained_neutral_relaxation(
        base,
        neighbors,
        combos,
        strategy,
        is_steering,
        target_count,
        np.random.default_rng(relax_seed),
        relax_mcs=relax_mcs,
    )
    start_stats = finite_local_stats(G, neighbors, strategy, is_steering)
    x_relaxed = float(start_stats["x"])
    p_actual = float(start_stats["p_actual"])
    eta0 = float(start_stats["eta0"])
    eta0_closed = eta0_pair_closure(x_relaxed, p_actual, base.k)
    M = base.resolved_M()
    r0 = classical_r0(base.k, base.n)
    r_star = local_r_star(x_relaxed, p_actual, eta0, base.k, base.n, M, base.c, base.omega)
    r_star_closed = local_r_star(x_relaxed, p_actual, eta0_closed, base.k, base.n, M, base.c, base.omega)
    r_global_closed = finite_global_r_star_closure_eta(base, p_actual, parse_csv_floats(str(spec["global_x_grid"])))
    r_values = [float(r) for r in spec.get("r_values", [base.r])]

    if mode == "global":
        rows = []
        T_global = float(spec["global_T_mcs"])
        record_interval = float(spec["global_record_interval_mcs"])
        tail_fraction = float(spec["global_tail_fraction"])
        high_threshold = float(spec["high_coop_threshold"])
        for ridx, r in enumerate(r_values):
            params_r = SimParams(**{**asdict(base), "r": r})
            for tr in trajectory_ids:
                trajectory_seed = int(spec["trajectory_seed0"]) + 1_000_000 * ridx + tr
                x_final, mean_tail, high_coop, hitting_time, reached_high = run_global_from_snapshot(
                    params_r,
                    neighbors,
                    combos,
                    strategy,
                    is_steering,
                    T_mcs=T_global,
                    record_interval_mcs=record_interval,
                    tail_fraction=tail_fraction,
                    high_threshold=high_threshold,
                    seed=trajectory_seed,
                )
                row = {
                    "exp_id": exp_id,
                    "mode": mode,
                    "accepted": True,
                    "network_id": network_id,
                    "steering_id": steering_id,
                    "snapshot_id": snapshot_id,
                    "trajectory_id": tr,
                    "network_seed": network_seed,
                    "steering_seed": steering_seed,
                    "strategy_seed": strategy_seed,
                    "relax_seed": relax_seed,
                    "trajectory_seed": trajectory_seed,
                    "N": base.N,
                    "k": base.k,
                    "n": base.n,
                    "M": M,
                    "M_mode": base.M_mode,
                    "r0": r0,
                    "r": r,
                    "omega": base.omega,
                    "c": base.c,
                    "p_target": base.p_s,
                    "num_steering": int(start_stats["num_steering"]),
                    "p_actual": p_actual,
                    "x_target": float(spec.get("x0_requested", base.x0)),
                    "x_initial": x_initial,
                    "x_after_relax": x_relaxed,
                    "x_final": x_final,
                    "mean_x_tail": mean_tail,
                    "global_T_mcs": T_global,
                    "global_record_interval_mcs": record_interval,
                    "global_tail_fraction": tail_fraction,
                    "high_coop_threshold": high_threshold,
                    "high_coop": high_coop,
                    "r_star_global_closure_eta": r_global_closed,
                    "eta0": eta0,
                    "eta0_closure": eta0_closed,
                    "q_S_given_D": start_stats["q_S_given_D"],
                    "q_C_given_D": start_stats["q_C_given_D"],
                    "q_C_given_S": start_stats["q_C_given_S"],
                    "q_S_given_S": start_stats["q_S_given_S"],
                    "E_DS": start_stats["E_DS"],
                    "E_SS": start_stats["E_SS"],
                    "density_DS": start_stats["density_DS"],
                    "density_SS": start_stats["density_SS"],
                    "placement_protocol": protocol,
                    "placement_bias": placement_bias,
                    "metropolis_steps": metropolis_steps,
                    "topology": base.topology,
                    "topology_param": base.topology_param,
                    "x0_requested": float(spec.get("x0_requested", base.x0)),
                    "x0_actual": float(base.x0),
                    "initial_state_seed": strategy_seed,
                    "valid_theory_state": bool(spec.get("valid_theory_state", True)),
                    "x_domain_min": float(spec.get("x_domain_min", float("nan"))),
                    "x_domain_max": float(spec.get("x_domain_max", float("nan"))),
                    "high_coop_indicator": high_coop,
                    "high_coop_probability_group": float("nan"),
                    "final_x_free": x_final,
                    "hitting_time_to_high_cooperation": hitting_time,
                    "reached_high_cooperation": reached_high,
                    "initial_eta0": eta0,
                    "initial_q_S_given_D": start_stats["q_S_given_D"],
                    "initial_q_C_given_D": start_stats["q_C_given_D"],
                    "initial_q_C_given_S": start_stats["q_C_given_S"],
                    "initial_q_S_given_S": start_stats["q_S_given_S"],
                    "_task_id": spec.get("task_id", ""),
                }
                row["run_id"] = run_id_for(row)
                rows.append(row)
        return rows

    rows = []
    for ridx, r in enumerate(r_values):
        params_r = SimParams(**{**asdict(base), "r": r})
        theory = finite_drift_theory(x_relaxed, p_actual, eta0, base.k, base.n, M, base.c, base.omega, r)
        theory_closed = finite_drift_theory(x_relaxed, p_actual, eta0_closed, base.k, base.n, M, base.c, base.omega, r)
        for tr in trajectory_ids:
            trajectory_seed = int(spec["trajectory_seed0"]) + 1_000_000 * ridx + tr
            _, x_final, drift = run_short_drift_from_snapshot(
                params_r, neighbors, combos, strategy, is_steering, tau_mcs=tau_mcs, seed=trajectory_seed
            )
            row = {
                "exp_id": exp_id,
                "mode": mode,
                "accepted": True,
                "network_id": network_id,
                "steering_id": steering_id,
                "snapshot_id": snapshot_id,
                "trajectory_id": tr,
                "network_seed": network_seed,
                "steering_seed": steering_seed,
                "strategy_seed": strategy_seed,
                "relax_seed": relax_seed,
                "trajectory_seed": trajectory_seed,
                "N": base.N,
                "k": base.k,
                "n": base.n,
                "M": M,
                "M_mode": base.M_mode,
                "r0": r0,
                "r": r,
                "omega": base.omega,
                "c": base.c,
                "p_target": base.p_s,
                "num_steering": int(start_stats["num_steering"]),
                "p_actual": p_actual,
                "x_target": base.x0,
                "x_initial": x_initial,
                "x_after_relax": x_relaxed,
                "x_final": x_final,
                "delta_t_mcs": tau_mcs,
                "relax_mcs": relax_mcs,
                "delta_x": x_final - x_relaxed,
                "drift_empirical": drift,
                "drift_theoretical": theory,
                "drift_theoretical_closure_eta": theory_closed,
                "r_star_local": r_star,
                "r_star_closure_eta": r_star_closed,
                "r_star_global_closure_eta": r_global_closed,
                "eta0": eta0,
                "eta0_closure": eta0_closed,
                "q_S_given_D": start_stats["q_S_given_D"],
                "q_C_given_D": start_stats["q_C_given_D"],
                "q_C_given_S": start_stats["q_C_given_S"],
                "q_S_given_S": start_stats["q_S_given_S"],
                "E_DS": start_stats["E_DS"],
                "E_SS": start_stats["E_SS"],
                "density_DS": start_stats["density_DS"],
                "density_SS": start_stats["density_SS"],
                "placement_protocol": protocol,
                "placement_bias": placement_bias,
                "metropolis_steps": metropolis_steps,
                "topology": base.topology,
                "topology_param": base.topology_param,
                "_task_id": spec.get("task_id", ""),
            }
            row["run_id"] = run_id_for(row)
            rows.append(row)
    return rows


def hierarchical_bootstrap_mean(
    df,
    value_col: str,
    n_boot: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    networks = df["network_id"].drop_duplicates().to_numpy()
    if networks.size == 0:
        return float("nan"), float("nan")
    vals = []
    for _ in range(n_boot):
        pieces = []
        sampled_g = rng.choice(networks, size=networks.size, replace=True)
        for g in sampled_g:
            dg = df[df["network_id"] == g]
            steering = dg["steering_id"].drop_duplicates().to_numpy()
            sampled_s = rng.choice(steering, size=steering.size, replace=True)
            for s in sampled_s:
                ds = dg[dg["steering_id"] == s]
                snapshots = ds["snapshot_id"].drop_duplicates().to_numpy()
                sampled_z = rng.choice(snapshots, size=snapshots.size, replace=True)
                for z in sampled_z:
                    dz = ds[ds["snapshot_id"] == z]
                    idx = rng.choice(dz.index.to_numpy(), size=len(dz), replace=True)
                    pieces.append(dz.loc[idx, value_col].to_numpy(dtype=float))
        if pieces:
            vals.append(float(np.nanmean(np.concatenate(pieces))))
    if not vals:
        return float("nan"), float("nan")
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def aggregate_outputs(out: Path, exp_id: str, n_boot: int, seed: int, save_mat: bool = True) -> None:
    import pandas as pd

    raw_path = out / "finite_drift_raw.csv"
    if not raw_path.exists():
        return
    df = pd.read_csv(raw_path).drop_duplicates(subset=["run_id"], keep="last")
    df["drift_theoretical"] = [
        finite_drift_theory(row.x_after_relax, row.p_actual, row.eta0, int(row.k), int(row.n), int(row.M), row.c, row.omega, row.r)
        for row in df.itertuples(index=False)
    ]
    if "eta0_closure" in df.columns:
        df["drift_theoretical_closure_eta"] = [
            finite_drift_theory(
                row.x_after_relax,
                row.p_actual,
                row.eta0_closure,
                int(row.k),
                int(row.n),
                int(row.M),
                row.c,
                row.omega,
                row.r,
            )
            for row in df.itertuples(index=False)
        ]
    keys = [
        "mode", "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
        "p_target", "x_target", "delta_t_mcs", "relax_mcs",
        "placement_protocol", "placement_bias", "topology", "topology_param",
    ]
    rows = []
    for vals, g in df.groupby(keys, dropna=False):
        drift = g["drift_empirical"].to_numpy(dtype=float)
        theory = g["drift_theoretical"].to_numpy(dtype=float)
        diff = drift - theory
        rmse = float(np.sqrt(np.nanmean(diff * diff)))
        denom = float(np.nanmax(theory) - np.nanmin(theory)) if len(theory) else float("nan")
        nrmse = rmse / denom if np.isfinite(denom) and denom > 0 else float("nan")
        ci_low, ci_high = hierarchical_bootstrap_mean(g, "drift_empirical", n_boot, seed)
        row = dict(zip(keys, vals))
        row.update({
            "exp_id": exp_id,
            "mean_p_actual": float(g["p_actual"].mean()),
            "num_networks": int(g["network_id"].nunique()),
            "num_steering_sets": int(g[["network_id", "steering_id"]].drop_duplicates().shape[0]),
            "num_snapshots": int(g[["network_id", "steering_id", "snapshot_id"]].drop_duplicates().shape[0]),
            "num_trajectories": int(len(g)),
            "mean_x_after_relax": float(g["x_after_relax"].mean()),
            "mean_eta0": float(g["eta0"].mean()),
            "mean_q_S_given_D": float(g["q_S_given_D"].mean()),
            "mean_drift_empirical": float(np.nanmean(drift)),
            "std_drift_empirical": float(np.nanstd(drift, ddof=1)) if len(g) > 1 else 0.0,
            "sem_naive": float(np.nanstd(drift, ddof=1) / np.sqrt(len(g))) if len(g) > 1 else 0.0,
            "ci95_low_boot": ci_low,
            "ci95_high_boot": ci_high,
            "mean_drift_theoretical": float(np.nanmean(theory)),
            "mean_drift_theoretical_closure_eta": float(g["drift_theoretical_closure_eta"].mean())
            if "drift_theoretical_closure_eta" in g.columns else float("nan"),
            "mean_r_star_local": float(g["r_star_local"].mean()),
            "mean_r_star_closure_eta": float(g["r_star_closure_eta"].mean())
            if "r_star_closure_eta" in g.columns else float("nan"),
            "sign_accuracy": float(np.nanmean(np.sign(drift) == np.sign(theory))),
            "rmse": rmse,
            "nrmse": nrmse,
            "mean_E_DS": float(g["E_DS"].mean()),
            "mean_E_SS": float(g["E_SS"].mean()),
        })
        rows.append(row)
    agg = pd.DataFrame(rows)
    if not agg.empty:
        agg.to_csv(out / "finite_drift_aggregate.csv", index=False, columns=AGG_FIELDS)

    boundary_rows = []
    bdf = agg[agg["mode"] == "boundary"].copy()
    bkeys = ["N", "k", "n", "M", "M_mode", "p_target", "x_target", "omega", "c", "topology", "topology_param"]
    for vals, g in bdf.groupby(bkeys, dropna=False):
        x = g["r"].to_numpy(dtype=float)
        y = g["mean_drift_empirical"].to_numpy(dtype=float)
        w = g["num_trajectories"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y) & (w > 0)
        if np.sum(ok) < 3:
            continue
        X = np.column_stack([np.ones(np.sum(ok)), x[ok]])
        sw = np.sqrt(w[ok])
        coef = np.linalg.lstsq(X * sw[:, None], y[ok] * sw, rcond=None)[0]
        intercept, slope = float(coef[0]), float(coef[1])
        yhat = X @ coef
        ybar = float(np.average(y[ok], weights=w[ok]))
        ss_res = float(np.sum(w[ok] * (y[ok] - yhat) ** 2))
        ss_tot = float(np.sum(w[ok] * (y[ok] - ybar) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        r_hat = -intercept / slope if slope != 0 else float("nan")
        row = dict(zip(bkeys, vals))
        theory_col = "mean_r_star_closure_eta" if "mean_r_star_closure_eta" in g.columns else "mean_r_star_local"
        theory_mean = float(g[theory_col].mean())
        row.update({
            "exp_id": exp_id,
            "mean_p_actual": float(g["mean_p_actual"].mean()),
            "r_star_theory_mean": theory_mean,
            "r_star_hat": r_hat,
            "ci95_low_boot": float("nan"),
            "ci95_high_boot": float("nan"),
            "slope": slope,
            "intercept": intercept,
            "r_squared": r2,
            "num_r_values": int(np.sum(ok)),
            "num_trajectories": int(g["num_trajectories"].sum()),
            "mae_vs_classical_r0": abs(r_hat - float(g["r0"].mean())) if np.isfinite(r_hat) else float("nan"),
            "abs_error_vs_theory": abs(r_hat - theory_mean) if np.isfinite(r_hat) else float("nan"),
        })
        boundary_rows.append(row)
    if boundary_rows:
        pd.DataFrame(boundary_rows).to_csv(out / "finite_boundary_estimates.csv", index=False, columns=BOUNDARY_FIELDS)

    if not save_mat:
        return
    try:
        from scipy.io import savemat

        payload = {"raw": df.to_records(index=False)}
        if not agg.empty:
            payload["aggregate"] = agg.to_records(index=False)
        if boundary_rows:
            payload["boundary"] = pd.DataFrame(boundary_rows).to_records(index=False)
        savemat(out / "finite_steering.mat", payload, long_field_names=True)
    except Exception as exc:
        print(f"MAT export skipped: {exc}")


def aggregate_global_outputs(out: Path, exp_id: str) -> None:
    import pandas as pd

    raw_path = out / "finite_global_raw.csv"
    if not raw_path.exists():
        return
    df = pd.read_csv(raw_path).drop_duplicates(subset=["run_id"], keep="last")
    keys = [
        "mode", "N", "k", "n", "M", "M_mode", "r0", "r", "omega", "c",
        "p_target", "global_T_mcs", "global_record_interval_mcs",
        "global_tail_fraction", "high_coop_threshold",
        "placement_protocol", "placement_bias", "topology", "topology_param",
    ]
    rows = []
    for vals, g in df.groupby(keys, dropna=False):
        by_x = g.groupby("x_target")["high_coop"].mean()
        row = dict(zip(keys, vals))
        row.update({
            "exp_id": exp_id,
            "mean_p_actual": float(g["p_actual"].mean()),
            "num_networks": int(g["network_id"].nunique()),
            "num_steering_sets": int(g[["network_id", "steering_id"]].drop_duplicates().shape[0]),
            "num_snapshots": int(g[["network_id", "steering_id", "snapshot_id"]].drop_duplicates().shape[0]),
            "num_trajectories": int(len(g)),
            "num_x0_values": int(g["x_target"].nunique()),
            "high_coop_probability": float(g["high_coop"].mean()),
            "mean_final_x_free": float(g["x_final"].mean()),
            "mean_tail_x_free": float(g["mean_x_tail"].mean()),
            "sem_tail_x_free": float(g["mean_x_tail"].std(ddof=1) / math.sqrt(len(g))) if len(g) > 1 else 0.0,
            "worst_x0_target": float(by_x.idxmin()) if len(by_x) else float("nan"),
            "min_high_coop_probability_by_x": float(by_x.min()) if len(by_x) else float("nan"),
            "r_star_global_closure_eta": float(g["r_star_global_closure_eta"].mean()),
            "mean_eta0": float(g["eta0"].mean()),
            "mean_eta0_closure": float(g["eta0_closure"].mean()),
            "mean_E_DS": float(g["E_DS"].mean()),
            "mean_E_SS": float(g["E_SS"].mean()),
        })
        rows.append(row)
    if rows:
        pd.DataFrame(rows).to_csv(out / "finite_global_aggregate.csv", index=False, columns=GLOBAL_AGG_FIELDS)


def build_specs(args: argparse.Namespace) -> list[tuple[SimParams, dict]]:
    p_values = parse_csv_floats(args.p_grid)
    x_values = parse_csv_floats(args.x_grid)
    protocols = parse_csv_strings(args.placement_protocols)
    biases = parse_csv_floats(args.placement_bias_grid)
    cases = parse_cases(args.case, args.k, args.n, args.topology)
    specs: list[tuple[SimParams, dict]] = []
    args._boundary_exclusions = []

    for case in cases:
        k = int(case["k"])
        n = int(case["n"])
        topology = str(case["topology"])
        base_kwargs = dict(
            N=args.N,
            k=k,
            n=n,
            c=args.c,
            omega=args.omega,
            M_mode=args.M_mode,
            M=args.M,
            topology=topology,
            topology_param=args.topology_param,
            T_mcs=1,
            record_interval_mcs=1,
            seed=args.seed0,
        )
        parameter_grid = None
        if args.parameter_grid_file:
            import pandas as pd

            parameter_grid = pd.read_csv(args.parameter_grid_file)
            if "requires_simulation" in parameter_grid:
                flag = parameter_grid["requires_simulation"]
                if flag.dtype == object:
                    flag = flag.astype(str).str.lower().isin(["true", "1", "yes"])
                parameter_grid = parameter_grid[flag.astype(bool)]
        if args.r_grid:
            r_values = parse_csv_floats(args.r_grid)
        else:
            r0 = classical_r0(k, n)
            r_values = [r0 - 0.35, r0 - 0.15, r0, r0 + 0.15]

        if args.mode == "boundary":
            r_offsets = parse_csv_floats(args.r_offset_grid)
            boundary_p_values = (
                parse_csv_floats(args.boundary_p_grid)
                if args.boundary_p_grid else [args.boundary_p]
            )
            for boundary_p in boundary_p_values:
                for x in x_values:
                    eta0 = eta0_pair_closure(x, boundary_p, k)
                    dummy = SimParams(**base_kwargs, p_s=boundary_p, x0=x, r=classical_r0(k, n))
                    center = local_r_star(
                        x, boundary_p, eta0, k, n, dummy.resolved_M(), args.c, args.omega
                    )
                    if not np.isfinite(center):
                        args._boundary_exclusions.append({
                            "p": boundary_p, "x": x, "r_theory": center,
                            "r": float("nan"), "reason": "invalid_finite_theory",
                        })
                        continue
                    if args.boundary_r_min is not None and center < args.boundary_r_min:
                        args._boundary_exclusions.append({
                            "p": boundary_p, "x": x, "r_theory": center,
                            "r": float("nan"),
                            "reason": "theoretical_crossing_below_physical_domain",
                        })
                        continue
                    candidate_r = [center + off for off in r_offsets]
                    r_values = []
                    for r_value in candidate_r:
                        if args.boundary_r_min is not None and r_value < args.boundary_r_min:
                            args._boundary_exclusions.append({
                                "p": boundary_p, "x": x, "r_theory": center,
                                "r": r_value, "reason": "r_below_physical_lower_bound",
                            })
                        else:
                            r_values.append(r_value)
                    r_values = sorted(set(round(value, 12) for value in r_values))
                    if len(r_values) < 3:
                        args._boundary_exclusions.append({
                            "p": boundary_p, "x": x, "r_theory": center,
                            "r": float("nan"), "reason": "fewer_than_three_valid_r_values",
                        })
                        continue
                    specs.append((SimParams(**base_kwargs, p_s=boundary_p, x0=x, r=center), {
                        "mode": args.mode,
                        "placement_protocol": "random",
                        "placement_bias": 0.0,
                        "metropolis_steps": args.metropolis_steps,
                        "global_x_grid": args.global_x_grid,
                        "r_values": r_values,
                    }))
            continue

        if args.mode == "global":
            for p in p_values:
                for x in x_values:
                    use_x = x
                    domain_info = {"x0_requested": x, "valid_theory_state": True,
                                   "x_domain_min": float("nan"), "x_domain_max": float("nan")}
                    if args.theory_global_file:
                        import pandas as pd

                        theory = pd.read_csv(args.theory_global_file)
                        match = theory[np.isclose(theory["p"], p)]
                        if not match.empty:
                            lo = float(match.iloc[0]["x_domain_min"])
                            hi = float(match.iloc[0]["x_domain_max"])
                            domain_info.update({"x_domain_min": lo, "x_domain_max": hi})
                            if x < lo or x > hi:
                                if args.adjust_outside_domain:
                                    use_x = min(max(x, lo + args.domain_margin), hi - args.domain_margin)
                                else:
                                    domain_info["valid_theory_state"] = False
                    if parameter_grid is not None:
                        p_rows = parameter_grid[np.isclose(parameter_grid["p"], p)]
                        use_r_values = sorted(p_rows["r"].astype(float).unique())
                    else:
                        use_r_values = r_values
                    for r in use_r_values:
                        specs.append((SimParams(**base_kwargs, p_s=p, x0=use_x, r=r), {
                            "mode": args.mode,
                            "placement_protocol": "random",
                            "placement_bias": 0.0,
                            "metropolis_steps": args.metropolis_steps,
                            "global_x_grid": args.global_x_grid,
                            "global_T_mcs": args.global_T_mcs,
                            "global_record_interval_mcs": args.global_record_interval_mcs,
                            "global_tail_fraction": args.global_tail_fraction,
                            "high_coop_threshold": args.high_coop_threshold,
                            **domain_info,
                        }))
            continue

        if args.mode == "exposure":
            protocols = protocols or ["random", "d_interface", "c_core", "clustered", "metropolis_exposure"]
            for protocol in protocols:
                use_biases = [0.0] if protocol in {"random", "clustered"} else biases
                for bias in use_biases:
                    for p in p_values:
                        for x in x_values:
                            for r in r_values:
                                specs.append((SimParams(**base_kwargs, p_s=p, x0=x, r=r), {
                                    "mode": args.mode,
                                    "placement_protocol": protocol,
                                    "placement_bias": bias,
                                    "metropolis_steps": args.metropolis_steps,
                                    "global_x_grid": args.global_x_grid,
                                }))
            continue

        for p in p_values:
            for x in x_values:
                for r in r_values:
                    specs.append((SimParams(**base_kwargs, p_s=p, x0=x, r=r), {
                        "mode": args.mode,
                        "placement_protocol": "random",
                        "placement_bias": 0.0,
                        "metropolis_steps": args.metropolis_steps,
                        "global_x_grid": args.global_x_grid,
                    }))
    return specs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite-steering state-conditioned Monte Carlo experiments.")
    parser.add_argument("--mode", choices=["drift", "boundary", "global", "exposure"], default="drift")
    parser.add_argument("--out-dir", default="results/exp04_finite_steering")
    parser.add_argument("--exp-id", default="E04finite")
    parser.add_argument("--N", type=int, default=2000)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--n", type=int, default=4)
    parser.add_argument(
        "--case",
        action="append",
        default=None,
        help="Repeatable case: k,n or k,n,topology. If omitted, --k/--n/--topology are used.",
    )
    parser.add_argument("--c", type=float, default=1.0)
    parser.add_argument("--omega", type=float, default=1e-3)
    parser.add_argument("--M-mode", choices=["all", "sampled"], default="all")
    parser.add_argument("--M", type=int, default=None)
    parser.add_argument("--topology", default="random_regular")
    parser.add_argument("--topology-param", type=float, default=1.0)
    parser.add_argument("--p-grid", default="0.02,0.05,0.10,0.15")
    parser.add_argument("--x-grid", default="0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,0.95")
    parser.add_argument("--r-grid", default="")
    parser.add_argument("--parameter-grid-file", default="")
    parser.add_argument("--theory-global-file", default="")
    parser.add_argument("--adjust-outside-domain", action="store_true")
    parser.add_argument("--domain-margin", type=float, default=1e-4)
    parser.add_argument("--boundary-p", type=float, default=0.02)
    parser.add_argument(
        "--boundary-p-grid", default="",
        help="Optional comma-separated boundary p values; preserves --boundary-p compatibility.",
    )
    parser.add_argument(
        "--boundary-r-min", type=float, default=None,
        help="Optional physical lower r bound. Crossings below it are audited and skipped.",
    )
    parser.add_argument("--r-offset-grid", default="-0.30,-0.20,-0.10,-0.05,0,0.05,0.10,0.20,0.30")
    parser.add_argument("--global-x-grid", default="0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,0.95")
    parser.add_argument("--placement-protocols", default="random,d_interface,c_core,clustered,metropolis_exposure")
    parser.add_argument("--placement-bias-grid", default="-4,-2,-1,0,1,2,4")
    parser.add_argument("--metropolis-steps", type=int, default=2000)
    parser.add_argument("--global-T-mcs", type=float, default=200.0)
    parser.add_argument("--global-record-interval-mcs", type=float, default=5.0)
    parser.add_argument("--global-tail-fraction", type=float, default=0.25)
    parser.add_argument("--high-coop-threshold", type=float, default=0.9)
    parser.add_argument("--relax-mcs", type=float, default=30.0)
    parser.add_argument("--tau-mcs", type=float, default=1.0)
    parser.add_argument("--networks", type=int, default=20)
    parser.add_argument("--steering-sets", type=int, default=10)
    parser.add_argument("--snapshots", type=int, default=1)
    parser.add_argument("--trajectories", type=int, default=50)
    parser.add_argument("--seed0", type=int, default=20260722)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--skip-aggregate", action="store_true", help="Write raw simulation output only; skip aggregate CSV/MAT export.")
    parser.add_argument("--no-mat", action="store_true", help="Write aggregate CSV files but skip MATLAB .mat export.")
    parser.add_argument("--no-checkpoint", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.smoke:
        if "--out-dir" not in sys.argv:
            args.out_dir = "/private/tmp/hsr_finite_smoke"
        args.N = 80
        args.k = 6
        args.n = 4
        args.case = ["4,3,random_regular", "6,4,random_regular"]
        args.M_mode = "sampled"
        args.M = 1
        args.p_grid = "0.05,0.10"
        args.x_grid = "0.3,0.7"
        r0 = classical_r0(args.k, args.n)
        args.r_grid = f"{r0 - 0.15},{r0 + 0.15}"
        args.relax_mcs = 1.0
        args.tau_mcs = 0.5
        args.networks = 2
        args.steering_sets = 2
        args.snapshots = 1
        args.trajectories = 2
        args.workers = 1
        args.bootstrap = 20
        args.progress_every = 2

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cfg = OutputConfig(out_dir=args.out_dir, exp_id=args.exp_id, save_timeseries=False, use_checkpoint=not args.no_checkpoint)
    spec_params = build_specs(args)
    if not spec_params:
        raise RuntimeError("No finite-steering parameter combinations to run.")

    write_metadata(out, cfg, spec_params[0][0], note="Finite-steering state-conditioned experiments.")
    with (out / "finite_design.json").open("w", encoding="utf-8") as f:
        json.dump({"created_at": time.strftime("%Y-%m-%d %H:%M:%S"), "args": vars(args)}, f, indent=2)
    if args.mode == "boundary":
        exclusions = list(getattr(args, "_boundary_exclusions", []))
        exclusion_path = out / "finite_boundary_excluded_cells.csv"
        with exclusion_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["p", "x", "r_theory", "r", "reason"]
            )
            writer.writeheader()
            writer.writerows(exclusions)
        print(f"Boundary exclusions: {len(exclusions)} -> {exclusion_path}")

    raw_path = out / ("finite_global_raw.csv" if args.mode == "global" else "finite_drift_raw.csv")
    completed = load_checkpoint(out / "checkpoint.json") if cfg.use_checkpoint else set()
    completed_run_ids = load_completed_run_ids(raw_path)
    tasks = []
    spec_idx = 0
    for base, spec in spec_params:
        for g in range(args.networks):
            network_seed = args.seed0 + 10_000_000 * spec_idx + 100_000 * g
            for s in range(args.steering_sets):
                steering_seed = network_seed + 1000 * s + 11
                for z in range(args.snapshots):
                    strategy_seed = network_seed + 1000 * s + 100 * z + 23
                    relax_seed = network_seed + 1000 * s + 100 * z + 37
                    trajectory_seed0 = network_seed + 1000 * s + 100 * z + 10_000
                    full_spec = {
                        **spec,
                        "network_seed": network_seed,
                        "steering_seed": steering_seed,
                        "strategy_seed": strategy_seed,
                        "relax_seed": relax_seed,
                        "trajectory_seed0": trajectory_seed0,
                    }
                    probe = {
                        "exp_id": args.exp_id,
                        "mode": args.mode,
                        "N": base.N,
                        "k": base.k,
                        "n": base.n,
                        "p_target": base.p_s,
                        "x_target": float(spec.get("x0_requested", base.x0)),
                        "r": base.r,
                        "placement_protocol": spec["placement_protocol"],
                        "placement_bias": spec["placement_bias"],
                        "network_id": g,
                        "steering_id": s,
                        "snapshot_id": z,
                        "trajectory_id": args.trajectories - 1,
                    }
                    snapshot_task_id = task_id_for(probe)
                    full_spec["task_id"] = snapshot_task_id
                    expected_run_ids: dict[int, str] = {}
                    for trajectory_id in range(args.trajectories):
                        trajectory_probe = dict(probe)
                        trajectory_probe["trajectory_id"] = trajectory_id
                        expected_run_ids[trajectory_id] = run_id_for(trajectory_probe)
                    missing_trajectory_ids = [
                        trajectory_id for trajectory_id, run_id in expected_run_ids.items()
                        if run_id not in completed_run_ids
                    ]
                    raw_task_complete = not missing_trajectory_ids
                    # A checkpoint records completion at snapshot-task level,
                    # whereas the raw CSV records individual trajectory IDs.
                    # The latter must win when a pilot (for example 20
                    # trajectories) is extended to a full design (30).
                    if not raw_task_complete:
                        tasks.append((base, full_spec, g, s, z, missing_trajectory_ids,
                                      args.relax_mcs, args.tau_mcs, args.exp_id))
        spec_idx += 1

    raw_fields = GLOBAL_RAW_FIELDS if args.mode == "global" else RAW_FIELDS
    checkpoint_path = out / "checkpoint.json"
    total = len(tasks)
    print(f"Finite steering ({args.mode}): starting {total} snapshot tasks")
    rows_per_task = args.trajectories
    if spec_params and "r_values" in spec_params[0][1]:
        rows_per_task *= len(spec_params[0][1]["r_values"])
    print(f"Each snapshot task emits {rows_per_task} trajectory rows")
    print(f"Output directory: {out}")

    done = 0
    if args.workers <= 1:
        for task in tasks:
            rows = run_snapshot_task(task)
            append_many(raw_path, rows, raw_fields)
            completed.add(str(rows[0].get("_task_id") or task_id_for(rows[0])))
            done += 1
            if args.progress_every > 0 and (done == total or done % args.progress_every == 0):
                print(f"Finite steering: {done}/{total} snapshot tasks complete")
                if cfg.use_checkpoint:
                    update_checkpoint(checkpoint_path, completed)
    else:
        ctx = mp.get_context("spawn")
        print(f"Finite steering: using {args.workers} workers")
        with ctx.Pool(processes=args.workers) as pool:
            batch: list[dict] = []
            for rows in pool.imap_unordered(run_snapshot_task, tasks, chunksize=args.chunk_size):
                batch.extend(rows)
                completed.add(str(rows[0].get("_task_id") or task_id_for(rows[0])))
                done += 1
                if len(batch) >= args.chunk_size * args.trajectories:
                    append_many(raw_path, batch, raw_fields)
                    batch.clear()
                if args.progress_every > 0 and (done == total or done % args.progress_every == 0):
                    if batch:
                        append_many(raw_path, batch, raw_fields)
                        batch.clear()
                    print(f"Finite steering: {done}/{total} snapshot tasks complete")
                    if cfg.use_checkpoint:
                        update_checkpoint(checkpoint_path, completed)
            if batch:
                append_many(raw_path, batch, raw_fields)

    if cfg.use_checkpoint:
        update_checkpoint(checkpoint_path, completed)
    if not args.skip_aggregate:
        if args.mode == "global":
            aggregate_global_outputs(out, args.exp_id)
        else:
            aggregate_outputs(out, args.exp_id, args.bootstrap, args.seed0 + 999, save_mat=not args.no_mat)
    print("Finite steering completed.")


if __name__ == "__main__":
    main()
