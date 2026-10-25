#!/usr/bin/env python3
"""Fig. 5d finite-steering coevolution and optional S-C-D snapshots.

The steering convention is explicit: steering positions are always expressed
cooperators (S), while x is cooperation among non-steering positions.  DB
updates reset a steering death site to C; relocation swaps only the steering
mask between a steering site and a cooperative non-steering site.
"""
from __future__ import annotations

import argparse, csv, json, math, multiprocessing as mp
from dataclasses import asdict
from pathlib import Path
import sys
import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hsr_sim.engine import (C, D, SimParams, db_update, graph_to_neighbors,
    make_graph, precompute_neighbor_combos, steering_count)
from hsr_sim.adaptive import local_defector_exposure, relocation_proposal

TRAJ_FIELDS = ["exp_id","replicate","time_mcs","N","k","n","p","x0","x","z","a",
    "d","dcrit","z_ad","z_rand","z_n_star","delta_ad","beta","gamma","R",
    "topology","seed"]
SUMMARY_FIELDS = ["exp_id","replicate","topology","seed","N","k","n","p","x0",
    "beta","gamma","R","min_delta_ad","fraction_delta_positive",
    "crossed_z_ad","x_final","z_final"]

def theory_row(x: float, p: float, k: int):
    a = p + (1.0-p)*x
    d = (a-a*a) / a if a > 0 else float("nan")
    dcrit = 2*a*(1-a)/(a*(1-p)+p) if a > 0 else float("nan")
    zad = a - 2*a*a*(1-a)/(a*(1-p)+p) if a > 0 else float("nan")
    zrand = a*a
    den = (k-1)*(a*(1+p)-p)
    zn = a*a*(1+(k-2)*a)/den if den > 0 else float("nan")
    return a, d, dcrit, zad, zrand, zn

def exact_state(N, p, x, rng):
    ns = steering_count(N, p)
    mask = np.zeros(N, dtype=bool)
    if ns: mask[rng.choice(N, size=ns, replace=False)] = True
    strategy = np.zeros(N, dtype=np.int8)
    # Critical convention: A = S union C, so all steering sites start as C.
    strategy[mask] = C
    free = np.flatnonzero(~mask)
    nc = int(round(x * free.size))
    if nc: strategy[rng.choice(free, size=nc, replace=False)] = C
    return strategy, mask

def pair_z(G, strategy, mask):
    # A means expressed cooperative: steering sites count as cooperative.
    A = (strategy == C) | mask
    return float(sum(A[u] and A[v] for u,v in G.edges()) * 2.0 / (G.number_of_nodes() * max(dict(G.degree()).values())))

def one_run(task):
    cfg, rep = task
    seed = int(cfg["seed0"] + 1_000_000 * rep)
    rng = np.random.default_rng(seed)
    N,k,n = cfg["N"],cfg["k"],cfg["n"]
    params = SimParams(N=N,k=k,n=n,r=cfg["R"]*k*n/(k+n-1),omega=cfg["omega"],
                       c=cfg["c"],M_mode="all",p_s=cfg["p"],x0=cfg["x0"],
                       topology=cfg["topology"],topology_param=cfg["topology_param"])
    G = make_graph(params, rng); neighbors = graph_to_neighbors(G,N)
    combos = precompute_neighbor_combos(neighbors,n)
    strategy, mask = exact_state(N,cfg["p"],cfg["x0"],rng)
    assert int(mask.sum()) == steering_count(N,cfg["p"])
    assert np.all(strategy[mask] == C)
    out_rows=[]; snap_dir = Path(cfg["snapshot_dir"]); snap_dir.mkdir(parents=True,exist_ok=True)
    total = int(round(cfg["T_mcs"]*N)); rec = max(1,int(round(cfg["record_mcs"]*N)))
    snap_every = max(1,int(round(cfg["snapshot_mcs"]*N))) if cfg["snapshot_mcs"] > 0 else 0
    min_delta=float("inf"); positive=0; count=0; crossed=False
    def record(step):
        nonlocal min_delta,positive,count,crossed
        p=float(mask.mean()); x=float(np.mean(strategy[~mask]==C)); z=pair_z(G,strategy,mask)
        a,d,dc,zad,zr,zn=theory_row(x,p,k); delta=z-zad
        min_delta=min(min_delta,delta); positive += int(delta>0); count += 1; crossed |= bool(delta<0)
        out_rows.append({"exp_id":cfg["exp_id"],"replicate":rep,"time_mcs":step/N,"N":N,"k":k,"n":n,
            "p":p,"x0":cfg["x0"],"x":x,"z":z,"a":a,"d":d,"dcrit":dc,"z_ad":zad,"z_rand":zr,
            "z_n_star":zn,"delta_ad":delta,"beta":cfg["beta"],"gamma":cfg["gamma"],
            "R":cfg["R"],"topology":cfg["topology"],"seed":seed})
        if snap_every and (step % snap_every == 0 or step == 0 or step == total):
            tag = f"{cfg['topology']}_p{cfg['p']:.5f}_x{cfg['x0']:.5f}_b{cfg['beta']:.3f}_rep{rep:04d}_t{step:07d}"
            np.savez_compressed(snap_dir / f"{tag}.npz",
                time_mcs=step/N, strategy=strategy.copy(), is_steering=mask.copy(),
                edges=np.asarray(list(G.edges()),dtype=np.int64), N=N,k=k,n=n,p=p,x=x,z=z,
                beta=cfg["beta"],gamma=cfg["gamma"],R=cfg["R"],seed=seed)
    record(0)
    for step in range(1,total+1):
        db_update(strategy,mask,neighbors,combos,params,rng)
        h=local_defector_exposure(strategy,neighbors)
        nrel=int(rng.poisson(cfg["gamma"]*float(mask.mean())))
        for _ in range(nrel): relocation_proposal(strategy,mask,neighbors,cfg["beta"],rng,exposure=h)
        if step % rec == 0 or step == total: record(step)
    summary={"exp_id":cfg["exp_id"],"replicate":rep,"topology":cfg["topology"],"seed":seed,
        "N":N,"k":k,"n":n,"p":cfg["p"],"x0":cfg["x0"],"beta":cfg["beta"],"gamma":cfg["gamma"],"R":cfg["R"],
        "min_delta_ad":min_delta,"fraction_delta_positive":positive/max(count,1),"crossed_z_ad":crossed,
        "x_final":out_rows[-1]["x"],"z_final":out_rows[-1]["z"]}
    return out_rows,summary

def floats(s): return [float(v.strip()) for v in s.split(',') if v.strip()]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir',default='results/fig5d_phase_space'); ap.add_argument('--exp-id',default='Fig5dPilot')
    ap.add_argument('--N',type=int,default=1000); ap.add_argument('--R',type=float,default=.5)
    ap.add_argument('--omega',type=float,default=1e-3); ap.add_argument('--c',type=float,default=1.)
    ap.add_argument('--p-grid',default='0.03'); ap.add_argument('--x-grid',default='0.2,0.6,0.9')
    ap.add_argument('--beta-grid',default='0,0.5'); ap.add_argument('--gamma',type=float,default=20.)
    ap.add_argument('--T-mcs',type=int,default=100); ap.add_argument('--record-mcs',type=float,default=1.)
    ap.add_argument('--snapshot-mcs',type=float,default=25.); ap.add_argument('--replicates',type=int,default=4)
    ap.add_argument('--topology',choices=['random_regular','square_lattice','moore_lattice'],default='random_regular')
    ap.add_argument('--topology-param',type=float,default=1.); ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--seed0',type=int,default=20260908); ap.add_argument('--progress-every',type=int,default=10)
    a=ap.parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True); (out/'snapshots').mkdir(exist_ok=True)
    if a.topology=='square_lattice' and int(math.isqrt(a.N))**2 != a.N: ap.error('square_lattice requires perfect-square N')
    if a.topology=='moore_lattice' and int(math.isqrt(a.N))**2 != a.N: ap.error('moore_lattice requires perfect-square N')
    if a.topology=='square_lattice' and a.N and a.N<4: ap.error('N too small')
    jobs=[{'exp_id':a.exp_id,'N':a.N,'k':6 if a.topology=='random_regular' else (4 if a.topology=='square_lattice' else 8),'n':4,
        'R':a.R,'omega':a.omega,'c':a.c,'p':p,'x0':x,'beta':b,'gamma':a.gamma,'T_mcs':a.T_mcs,
        'record_mcs':a.record_mcs,'snapshot_mcs':a.snapshot_mcs,'topology':a.topology,'topology_param':a.topology_param,
        'seed0':a.seed0,'exp_id':a.exp_id,'snapshot_dir':str(out/'snapshots')}
        for p in floats(a.p_grid) for x in floats(a.x_grid) for b in floats(a.beta_grid)]
    tasks=[(j,r) for j in jobs for r in range(a.replicates)]; print(f'Fig. 5d: {len(tasks)} trajectories; workers={a.workers}',flush=True)
    results=[]
    if a.workers>1:
        with mp.get_context('spawn').Pool(a.workers) as pool:
            it=pool.imap(one_run,tasks,chunksize=1)
            for i,v in enumerate(it,1):
                results.append(v)
                if i == 1 or i % max(1,a.progress_every) == 0 or i == len(tasks):
                    print(f'Fig. 5d progress: {i}/{len(tasks)}',flush=True)
    else:
        for i,t in enumerate(tasks,1): results.append(one_run(t)); print(f'Fig. 5d progress: {i}/{len(tasks)}',flush=True)
    rows=[r for rr,_ in results for r in rr]; sums=[s for _,s in results]
    with (out/'fig5d_trajectories.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=TRAJ_FIELDS); w.writeheader(); w.writerows(rows)
    with (out/'fig5d_summary.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=SUMMARY_FIELDS); w.writeheader(); w.writerows(sums)
    (out/'metadata.json').write_text(json.dumps({'experiment':a.exp_id,'model':'finite DB + adaptive relocation','parameters':vars(a),'steering_convention':'strategy[steering]=C; x=cooperation among non-steering'},indent=2))
    print(f'Wrote {len(rows)} trajectory rows and {len(sums)} summaries to {out}',flush=True)
if __name__=='__main__': main()
