#!/usr/bin/env python3
"""Follow-up to factorization_test.py (audit of 2026-10-07): same builder and selection, with the
non-true candidates split by the origin of their prompt.  Classes per candidate triple:
  0 true    : all three singles share one simulated eventID
  1 same    : not true, and the prompt single's decay position (sx, sy, sz) lies in the same tube
              region as the pair's TOF annihilation point
  2 other   : not true, prompt decay elsewhere
Histograms H[shard, arm(0 signal, 1 delayed +50 ns), class, sel(0 |s|<=9, 1 |s|<=5), voi, bin]."""
import json, os, sys, time
from pathlib import Path
import h5py
import numpy as np

S = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-17_derivation_test"
sys.path.insert(0, str(S))
import sim1_extract as X  # noqa: E402

OUT = Path(__file__).resolve().parent
SINGLES = Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5"
CFG = json.loads((S / "sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
SIG_T_NS = 0.0909
W_WIDE = 9.0
EDGES = np.arange(-6.0, 9.0001, 0.05)
NB = len(EDGES) - 1
NV = len(VOIS)


def shard_hists(t_true, in5, inp, gx, gy, gz, eid, sx, sy, sz, rng):
    H = np.zeros((2, 3, 2, NV, NB), np.int64)
    n = len(t_true)
    t = t_true + rng.normal(0.0, SIG_T_NS, size=n)
    i5 = np.flatnonzero(in5.astype(bool)); ip = np.flatnonzero(inp.astype(bool))
    if len(i5) < 2 or len(ip) < 1:
        return H
    o5 = np.argsort(t[i5], kind="stable"); op = np.argsort(t[ip], kind="stable")
    i5 = i5[o5]; ip = ip[op]; t5 = t[i5]; tp = t[ip]
    ok = np.diff(t5) <= X.W_PAIR_NS
    pa = i5[:-1][ok]; pb = i5[1:][ok]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    pann = X.tof_p_ann(p1, p2, t[pa], t[pb])
    geom = (X.lor_perp(pann, p1, p2) <= OUTLINE.get("lor_tol_mm", 10.0)) & X.in_phantom_outline(
        pann[:, 0], pann[:, 1], pann[:, 2], OUTLINE)
    pa, pb, pann, p1, p2 = pa[geom], pb[geom], pann[geom], p1[geom], p2[geom]
    tbar = 0.5 * (t[pa] + t[pb])
    masks = X.voi_masks(pann, VOIS)
    for arm, shift in ((0, 0.0), (1, X.DELAY_NS)):
        sel, pick = X.all_prompt_indices(tbar + shift, tp, W_WIDE)
        if len(sel) == 0:
            continue
        q = ip[pick]
        pp = np.stack([gx[q], gy[q], gz[q]], -1)
        tau = X.per_event_tau(p1[sel], p2[sel], pp, t[pa[sel]], t[pb[sel]], t[q], pann[sel]) + shift
        s = t[q] - (tbar[sel] + shift)
        inside = np.abs(s) <= X.W_PROMPT_NS
        true = (eid[pa[sel]] == eid[pb[sel]]) & (eid[pa[sel]] == eid[q])
        src_q = np.stack([sx[q], sy[q], sz[q]], -1).astype(float)
        qmasks = X.voi_masks(src_q, VOIS)
        b = np.searchsorted(EDGES, tau, side="right") - 1
        good = (b >= 0) & (b < NB)
        for v in range(NV):
            mv = masks[v][sel] & good
            same = qmasks[v]
            cls = np.where(true, 0, np.where(same, 1, 2))
            for c in range(3):
                m = mv & (cls == c)
                H[arm, c, 1, v] += np.bincount(b[m & inside], minlength=NB)
                H[arm, c, 0, v] += np.bincount(b[m], minlength=NB)
    return H


def main():
    if not Path(SINGLES).exists(): raise SystemExit(f"{SINGLES}: raw data not included in the release, see analysis/README.md")
    t0 = time.time()
    f = h5py.File(SINGLES, "r"); g = f["singles"]
    bounds = X.find_shard_boundaries(g["shard_id"][:])
    rng = np.random.default_rng(90900)
    HS = []
    for k, (sh, s, e) in enumerate(bounds):
        HS.append(shard_hists(g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e], g["gx"][s:e], g["gy"][s:e],
                              g["gz"][s:e], g["eventID"][s:e], g["sx"][s:e], g["sy"][s:e], g["sz"][s:e], rng))
        if k % 16 == 0:
            print(f"shard {k + 1}/{len(bounds)}  {time.time() - t0:.0f}s", flush=True)
    f.close()
    np.savez_compressed(OUT / "factorization_hists2.npz", H=np.array(HS), edges=EDGES, voi=[v["name"] for v in VOIS],
                        tau_cfg=[v.get("tau_cfg_ns") for v in VOIS])
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
