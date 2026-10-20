#!/usr/bin/env python3
"""Review round 17, verification follow-up: delay profile of the labeled accidental content (AP + PA_unc) in three
bands of the fit window, per unit a priori acceptance and relative to the sideband (-4.5 to -2 ns), for the
adjacent-pair and all-pairs builders. Separates the 1-3 ns deficit (pair breaking) from the 4.5-5 ns edge excess.
Pooled over the eight tubes, shard bootstrap (500 resamples). Input rules16_hists.npz, output rules16_flat_bands.json."""
import json
from pathlib import Path
import numpy as np
from scipy.special import erfc
D = Path(__file__).resolve().parent
Z = np.load(D / "rules16_hists.npz"); H = Z["hist"].astype(float)
KS, ARMS, CLS, NAMES, BINS = list(Z["ks"]), list(Z["arms"]), list(Z["classes"]), list(Z["voi_names"]), Z["bins"]
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
bc = 0.5 * (BINS[:-1] + BINS[1:]); Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
A = Phi((5 - bc) / 0.35897) * Phi((bc + 5) / 0.35897)
acc = H[..., CLS.index("AP"), :] + H[..., CLS.index("PA_unc"), :]
sb = (bc > -4.5) & (bc < -2); S = H.shape[0]; rng = np.random.default_rng(7)
bands = {"1-3": (1, 3), "3-4.5": (3, 4.5), "4.5-5": (4.5, 5)}
def prof(a):
    a = a[..., TI, :].sum(-2); ref = a[..., sb].sum(-1) / A[sb].sum()
    return np.stack([a[..., (bc > lo) & (bc < hi)].sum(-1) / A[(bc > lo) & (bc < hi)].sum() / ref - 1
                     for lo, hi in bands.values()], -1)
P0 = prof(acc.sum(0)); PB = []
for _ in range(10):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    PB.append(prof(np.tensordot(w, acc, axes=(1, 0))))
sd = np.concatenate(PB).std(0, ddof=1); out = {}
for arm in ["all", "ap_all", "closest", "ap_closest"]:
    ia = ARMS.index(arm); out[arm] = {}
    for ik, k in enumerate(KS):
        out[arm][str(k)] = {b: [100 * P0[ik, ia, j], 100 * sd[ik, ia, j]] for j, b in enumerate(bands)}
        print(arm, k, " ".join(f"{b}: {100*P0[ik,ia,j]:+.2f}+-{100*sd[ik,ia,j]:.2f}" for j, b in enumerate(bands)))
json.dump(out, open(D / "rules16_flat_bands.json", "w"), indent=1)
