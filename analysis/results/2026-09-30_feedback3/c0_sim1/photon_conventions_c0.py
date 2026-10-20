#!/usr/bin/env python3
"""Referee round 2, MC13: photon-1, photon-2 and symmetric delay conventions inside
the PRODUCTION clinical extractor (sim1_extract.clinical_per_shard, unchanged apart
from the delay function), full Sim-1 (128 shards), both matching rules.

Conventions (p_ann = TOF point, label-symmetric, and dp = prompt path from p_ann):
  photon-1   : (t1 - |p1-p_ann|/c) - (tp - dp/c)            [production]
  photon-2   : (t2 - |p2-p_ann|/c) - (tp - dp/c)
  symmetric  : ((t1+t2)/2 - |p1-p2|/(2c)) - (tp - dp/c)      [exchange-invariant]
For a TOF point between the two hits all three are algebraically identical. They
differ only when the TOF point falls outside the chord.

Output: photon_conventions.npz with hist[shard, conv, rule, voi+1, sig/rand, bin]
        and per-shard counts of accepted triples whose conventions differ.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np, h5py

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX  # noqa: E402

SRC = Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5"  # feedback round 3: corrected-yield arm C0
CFG = json.loads((RES / "2026-08-17_derivation_test/sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
CONVS = ["photon1", "photon2", "symmetric"]
RULES = ["closest", "all"]
C = SX.C_MM_NS
NB = len(SX.TAU_BINS) - 1
OUT = RES / "2026-09-30_feedback3/c0_sim1/photon_conventions_c0"
_ORIG = SX.per_event_tau
STATE = {"conv": "photon1", "n": 0, "ndiff": 0}


def tau_conv(p1, p2, pp, t1, t2, tp, p_ann):
    dp = np.linalg.norm(pp - p_ann, axis=-1)
    base = tp - dp / C
    a1 = (t1 - np.linalg.norm(p1 - p_ann, axis=-1) / C) - base
    a2 = (t2 - np.linalg.norm(p2 - p_ann, axis=-1) / C) - base
    sy = (0.5 * (t1 + t2) - np.linalg.norm(p1 - p2, axis=-1) / (2 * C)) - base
    STATE["n"] += len(a1)
    STATE["ndiff"] += int((np.abs(a1 - sy) > 1e-6).sum())
    return {"photon1": a1, "photon2": a2, "symmetric": sy}[STATE["conv"]]


SX.per_event_tau = tau_conv


def one_shard(args):
    sh_id, a, b = args
    with h5py.File(SRC, "r") as f:
        g = f["singles"]
        t = g["time_ns"][a:b].astype(np.float64)
        in5 = g["in_511"][a:b]; inp = g["in_prompt"][a:b]
        gx = g["gx"][a:b]; gy = g["gy"][a:b]; gz = g["gz"][a:b]
    H = np.zeros((len(CONVS), len(RULES), len(VOIS) + 1, 2, NB), dtype=np.int32)
    cnt = np.zeros((len(CONVS), len(RULES), 2), dtype=np.int64)
    for ic, conv in enumerate(CONVS):
        for ir, rule in enumerate(RULES):
            STATE.update(conv=conv, n=0, ndiff=0)
            rng = np.random.default_rng([90900, int(sh_id)])
            out = SX.clinical_per_shard(t, in5, inp, gx, gy, gz, 0.0909, rng, VOIS, OUTLINE, matching=rule)
            cnt[ic, ir] = (STATE["n"], STATE["ndiff"])
            if out is None:
                continue
            by_sig, by_rand = out
            for vi in range(-1, len(VOIS)):
                if by_sig[vi]:
                    H[ic, ir, vi + 1, 0] = np.histogram(np.concatenate(by_sig[vi]), bins=SX.TAU_BINS)[0]
                if by_rand[vi]:
                    H[ic, ir, vi + 1, 1] = np.histogram(np.concatenate(by_rand[vi]), bins=SX.TAU_BINS)[0]
    return sh_id, H, cnt


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    with h5py.File(SRC, "r") as f:
        sh = f["singles"]["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh)
    t0 = time.time()
    res = {}
    with Pool(24) as pool:
        for i, (sid, H, cnt) in enumerate(pool.imap_unordered(one_shard, bounds)):
            res[sid] = (H, cnt)
            if i % 16 == 0:
                print(f"  {i+1}/{len(bounds)}  {time.time()-t0:.0f} s", flush=True)
    ids = sorted(res)
    np.savez_compressed(OUT / "photon_conventions.npz",
                        hist=np.stack([res[s][0] for s in ids]), counts=np.stack([res[s][1] for s in ids]),
                        shard_ids=np.array(ids), convs=np.array(CONVS), rules=np.array(RULES),
                        bins=SX.TAU_BINS, voi_names=np.array(["bg"] + [v["name"] for v in VOIS]),
                        tau_cfg=np.array([np.nan] + [v["tau_cfg_ns"] for v in VOIS]))
    print(f"done {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
