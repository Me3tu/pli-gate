#!/usr/bin/env python3
"""Review round 13 (2026-10-14), referee M1: the timing-resolution sweep of the closest-prompt bias, repeated on the
corrected-yield zero-delay eight-tube stream C0 with the production clinical extractor (sim1_extract.clinical_per_shard,
imported unchanged, photon-1 delay convention), both matching rules.

The published sweep (Fig. 6a) used the discovery stream and the extractor of that stage. Here the single-photon timing
sigma takes the values of the published sweep, 90.9, 68.5, 30.0 and 15.0 ps (coincidence resolution 303, 228, 100 and
50 ps FWHM), plus 0 (noiseless timing). Smear draws per shard rng [round(sigma_ps*1000), shard]. The 90.9 ps arm uses
[90900, shard] and therefore reproduces photon_conventions_c0 (photon-1) bin for bin, which is checked.
Output crt_sweep_c0_hists.npz: hist[shard, sigma, rule, voi+1, sig/rnd, bin]."""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np, h5py

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX  # noqa: E402

SRC = str(Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5")
CFG = json.loads((RES / "2026-08-17_derivation_test/sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
SIG_PS = [90.9, 68.5, 30.0, 15.0, 0.0]
CRT_PS = [303, 228, 100, 50, 0]
RULES = ["closest", "all"]
NB = len(SX.TAU_BINS) - 1
OUT = RES / "2026-10-14_known_randoms"


def one_shard(args):
    sh_id, a, b, src = args
    with h5py.File(src, "r") as f:
        g = f["singles"]
        t = g["time_ns"][a:b].astype(np.float64)
        in5 = g["in_511"][a:b]; inp = g["in_prompt"][a:b]
        gx = g["gx"][a:b]; gy = g["gy"][a:b]; gz = g["gz"][a:b]
    H = np.zeros((len(SIG_PS), len(RULES), len(VOIS) + 1, 2, NB), dtype=np.int32)
    for js, sp in enumerate(SIG_PS):
        for ir, rule in enumerate(RULES):
            rng = np.random.default_rng([int(round(sp * 1000)), int(sh_id)])
            out = SX.clinical_per_shard(t, in5, inp, gx, gy, gz, sp / 1000.0, rng, VOIS, OUTLINE, matching=rule)
            if out is None: continue
            by_sig, by_rand = out
            for vi in range(-1, len(VOIS)):
                if by_sig[vi]: H[js, ir, vi + 1, 0] = np.histogram(np.concatenate(by_sig[vi]), bins=SX.TAU_BINS)[0]
                if by_rand[vi]: H[js, ir, vi + 1, 1] = np.histogram(np.concatenate(by_rand[vi]), bins=SX.TAU_BINS)[0]
    return sh_id, H


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    with h5py.File(SRC, "r") as f: sh = f["singles"]["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh); t0 = time.time(); res = {}
    with Pool(24) as pool:
        for i, (sid, H) in enumerate(pool.imap_unordered(one_shard, [(s, a, b, SRC) for s, a, b in bounds])):
            res[sid] = H
            if i % 16 == 0: print(f"  {i+1}/{len(bounds)}  {time.time()-t0:.0f} s", flush=True)
    ids = sorted(res); HH = np.stack([res[s] for s in ids])
    np.savez_compressed(OUT / "crt_sweep_c0_hists.npz", hist=HH, shard_ids=np.array(ids), sigma_ps=np.array(SIG_PS),
                        crt_ps=np.array(CRT_PS), rules=np.array(RULES), bins=SX.TAU_BINS,
                        voi_names=np.array(["bg"] + [v["name"] for v in VOIS]), tau_cfg=np.array([np.nan] + [v["tau_cfg_ns"] for v in VOIS]))
    P = np.load(RES / "2026-09-30_feedback3/c0_sim1/photon_conventions_c0/photon_conventions.npz")
    same = bool(np.array_equal(P["hist"][:, 0], HH[:, 0]))
    print(f"90.9 ps arm identical to photon_conventions_c0 photon-1: {same}")
    json.dump(dict(reproduces_photon_conventions_303ps=same, wall_s=time.time() - t0), open(OUT / "crt_sweep_c0_validation.json", "w"), indent=1)


if __name__ == "__main__":
    main()
