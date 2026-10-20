#!/usr/bin/env python3
"""Review round 18, follow-up to part B of r18_const.py. The fitted lower accidental edge of the +-5 ns search lies
INSIDE -5 ns (d = -38 to -56 ps), whereas the randoms-arm edge at +5 ns lies beyond it (+44 to +64 ps, Note S17): the
accidental window is translated by about +50 ps, not widened, so a mirrored edge cannot describe it. Here the lower
edge (width s, position -5 - d) is translated to the upper end, A_tr(t) = Phi((5 - d - t)/s) Phi((t + 5 + d)/s), the
constant is the sideband plateau under A_tr, and the readout is the a priori acceptance of the true triples.
Same histograms, edges (from r18_const.json), bootstrap resamples (B = 1000, seed 20261019) as r18_const.py part A.
Writes r18_translate.json."""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = Path(__file__).resolve().parent
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
sd = lambda a: 0.5 * (np.nanpercentile(a, 84, axis=0) - np.nanpercentile(a, 16, axis=0))
tg = np.linspace(0.3, 12.0, 6000); tt = np.linspace(1.0, 5.0, 8001)
SIG_E = 0.35897
A359 = lambda t: Phi((5 - t) / SIG_E) * Phi((t + 5) / SIG_E)
E = np.exp(-tt[None, :] / tg[:, None]) * A359(tt)[None, :]; CUR = (E * tt).sum(1) / E.sum(1)
Z = np.load(RES / "2026-10-19_review16/rules16_hists.npz")
H = Z["hist"].astype(np.float64)
KS, ARMS, CLS = [float(k) for k in Z["ks"]], [str(a) for a in Z["arms"]], [str(c) for c in Z["classes"]]
NAMES, TCFG, BINS = [str(n) for n in Z["voi_names"]], Z["tau_cfg"].astype(float), Z["bins"]
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
bc = 0.5 * (BINS[:-1] + BINS[1:]); fw = (bc > 1.0) & (bc < 5.0); x = bc[fw]; sbm = (bc > -4.5) & (bc < -2.0)
iS, iG, iPc = (CLS.index(c) for c in ("sig", "G", "PA_corr"))
EF = json.loads((OUT / "r18_const.json").read_text())["mirror_edge_fit"]


def tau(spec):
    n = spec.sum(-1); m = (spec * x).sum(-1) / np.where(n > 0, n, np.nan)
    return np.interp(m, CUR, tg, left=np.nan, right=np.nan)


def dR_tr(Hs):
    res = np.zeros(Hs.shape[:-5] + (len(KS), len(ARMS), len(NAMES)))
    for ik, k in enumerate(KS):
        for ia, arm in enumerate(ARMS):
            e = EF[f"k={k}|{arm}"]; s, d = e["s_ps"] / 1000, e["d_ps"] / 1000
            A = Phi((5 - d - bc) / s) * Phi((bc + 5 + d) / s)
            hh = Hs[..., ik, ia, :, :, :]
            c = hh[..., iS, :][..., sbm].sum(-1) / A[sbm].sum()
            h = hh[..., fw]
            res[..., ik, ia, :] = 100 * (tau(h[..., iS, :] - c[..., None] * A[fw]) - tau(h[..., iG, :] + h[..., iPc, :])) / TCFG
    return res


S = H.shape[0]; flat = H.reshape(S, -1)
P0 = dR_tr(H.sum(0)); rng = np.random.default_rng(20261019); PB = []
for _ in range(20):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    PB.append(dR_tr((w @ flat).reshape((50,) + H.shape[1:])))
PB = np.concatenate(PB)
out = {"doc": __doc__, "ks": KS, "arms": {}}
for ia, arm in enumerate(ARMS):
    y0 = np.nanmean(P0[:, ia, TI], -1); yb = np.nanmean(PB[:, :, ia, TI], -1)
    out["arms"][arm] = [[float(y0[i]), float(sd(yb[:, i]))] for i in range(len(KS))]
    print(f"{arm:11s} translated-edge dR " + " ".join(f"{v:+.2f}+-{s:.2f}" for v, s in out["arms"][arm]))
json.dump(out, open(OUT / "r18_translate.json", "w"), indent=1)
