#!/usr/bin/env python3
"""Referee round 2, major concern 10: correlation-aware inference for the rate ladder.

The published ladder (2026-08-26_rate_scan) re-extracted ONE Sim-1 dataset at
twelve clock compressions k and fitted ln|Delta| against ln mu with per-rung
errors taken as independent. The rungs share every event, so their errors are
strongly correlated. This script rebuilds the ladder SHARD BY SHARD, keeping the
per-shard delay histograms of every rung and both matching rules, so that a
bootstrap over shards (the independent simulation jobs, each with its own seed)
applies the SAME resample to all rungs and carries the full cross-rung
covariance into the exponent.

Per shard: singles are loaded once; for each k the decay epochs are compressed
exactly as rescale_rate.py does (t_new = t0/k + (t - t0) per eventID group);
clinical_per_shard (the published extractor, imported unchanged) builds triples
under both rules with the SAME timing-smear draws for every k and rule
(rng seeded per shard, re-created per call), so rungs differ only in rate.

Output: ladder_joint_hists.npz  hist[shard, k, rule, voi+1, bin]  (int32)
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import h5py

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX  # noqa: E402  (published extractor, unchanged)

SRC = Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5"
CFG = json.loads((RES / "2026-08-17_derivation_test/sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
KS = [0.15, 0.25, 0.35, 0.50, 0.70, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00]
RULES = ["closest", "all"]
SIGMA_T_NS = 0.0909
NB = len(SX.TAU_BINS) - 1
OUT = RES / "2026-10-08_c0_ladder/ladder_joint"


def rescale_times(t, eid, k):
    order = np.argsort(eid, kind="stable")
    es, ts = eid[order], t[order]
    newgrp = np.empty(len(es), dtype=bool); newgrp[0] = True; newgrp[1:] = es[1:] != es[:-1]
    gid = np.cumsum(newgrp) - 1
    base = np.minimum.reduceat(ts, np.flatnonzero(newgrp))[gid]
    out = np.empty_like(t)
    out[order] = base / k + (ts - base)
    return out


def one_shard(args):
    sh_id, a, b = args
    with h5py.File(SRC, "r") as f:
        g = f["singles"]
        t = g["time_ns"][a:b].astype(np.float64); eid = g["eventID"][a:b]
        in5 = g["in_511"][a:b]; inp = g["in_prompt"][a:b]
        gx = g["gx"][a:b]; gy = g["gy"][a:b]; gz = g["gz"][a:b]
    H = np.zeros((len(KS), 2, len(VOIS) + 1, 2, NB), dtype=np.int32)   # last-but-one: 0 signal, 1 randoms
    for ik, k in enumerate(KS):
        tk = t if k == 1.0 else rescale_times(t, eid, k)
        for ir, rule in enumerate(RULES):
            rng = np.random.default_rng([90900, int(sh_id)])   # identical smear for every k and rule
            out = SX.clinical_per_shard(tk, in5, inp, gx, gy, gz, SIGMA_T_NS, rng, VOIS, OUTLINE, matching=rule)
            if out is None:
                continue
            by_sig, by_rand = out
            for vi in range(-1, len(VOIS)):
                if by_sig[vi]:
                    H[ik, ir, vi + 1, 0] = np.histogram(np.concatenate(by_sig[vi]), bins=SX.TAU_BINS)[0]
                if by_rand[vi]:
                    H[ik, ir, vi + 1, 1] = np.histogram(np.concatenate(by_rand[vi]), bins=SX.TAU_BINS)[0]
    return sh_id, H


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    with h5py.File(SRC, "r") as f:
        sh = f["singles"]["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh)
    print(f"{len(bounds)} shards, {len(KS)} rungs, {len(RULES)} rules", flush=True)
    t0 = time.time()
    res = {}
    with Pool(24) as pool:
        for i, (sid, H) in enumerate(pool.imap_unordered(one_shard, bounds)):
            res[sid] = H
            if i % 8 == 0:
                print(f"  {i+1}/{len(bounds)} shards  {time.time()-t0:.0f} s", flush=True)
    ids = sorted(res)
    np.savez_compressed(OUT / "ladder_joint_hists.npz", hist=np.stack([res[s] for s in ids]),
                        shard_ids=np.array(ids), ks=np.array(KS), bins=SX.TAU_BINS,
                        voi_names=np.array(["bg"] + [v["name"] for v in VOIS]),
                        tau_cfg=np.array([np.nan] + [v["tau_cfg_ns"] for v in VOIS]))
    print(f"done {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
