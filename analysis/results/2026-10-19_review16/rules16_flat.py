#!/usr/bin/env python3
"""Review round 17 diagnostic: flatness of the true uncorrelated accidental content (AP + PA_unc, event labels)
per prompt rule, from rules16_hists.npz. F = density in the fit window (1-5 ns) / density on the sideband
(-4.5 to -2 ns), each divided by the a priori acceptance. F - 1 is the departure from flatness that a sideband
constant cannot see. Also prints the profile of the all-prompts arm in 1 ns steps. Pooled over the eight tubes."""
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
iAP, iPu, iG = CLS.index("AP"), CLS.index("PA_unc"), CLS.index("G")
fit = (bc > 1) & (bc < 5); sb = (bc > -4.5) & (bc < -2)
S = H.shape[0]; rng = np.random.default_rng(20261019)
acc = H[..., iAP, :] + H[..., iPu, :]                 # [shard, k, arm, voi, bin]
def F(a):
    a = a[..., TI, :].sum(-2)
    return (a[..., fit].sum(-1) / A[fit].sum()) / (a[..., sb].sum(-1) / A[sb].sum()) - 1
F0 = F(acc.sum(0)); FB = []
for _ in range(10):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    FB.append(F(np.tensordot(w, acc, axes=(1, 0))))
FB = np.concatenate(FB); sd = FB.std(0, ddof=1)
out = {}
for ia, arm in enumerate(ARMS):
    print(f"{arm:11s} " + "  ".join(f"k={k}: {100*F0[ik,ia]:+6.2f}+-{100*sd[ik,ia]:.2f}%" for ik, k in enumerate(KS)))
    out[arm] = {str(k): [100 * F0[ik, ia], 100 * sd[ik, ia]] for ik, k in enumerate(KS)}
ia = ARMS.index("all")
for ik, k in enumerate(KS):
    a = acc.sum(0)[ik, ia][TI].sum(0); ref = a[sb].sum() / A[sb].sum()
    prof = [(lo, a[(bc > lo) & (bc < lo + 1)].sum() / A[(bc > lo) & (bc < lo + 1)].sum() / ref - 1) for lo in range(-5, 5)]
    print(f"all k={k} profile (density/sideband - 1, %): " + " ".join(f"[{lo},{lo+1}] {100*v:+.1f}" for lo, v in prof))
json.dump(out, open(D / "rules16_flat.json", "w"), indent=1)
