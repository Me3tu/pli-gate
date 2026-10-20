#!/usr/bin/env python3
"""Shard independence of the eight-tube stream (Methods), written for release v1.1.0.

The published shard_independence.json had no script on disk (inventory flag F2). This script reproduces it.
Data: genuine (event-identity) triples of the corrected-yield zero-delay stream C0, all-prompts rule, rung k = 1.00
(the production rate), 128 shards, from 2026-10-14_known_randoms/known_randoms_hists.npz (class G).
Readout per tube: windowed mean on (1, 5) ns inverted with the sharp-edge moment curve (gridded as in
acceptance_tm_c0.py). Uncertainties of that readout:
  analytic   counting: sqrt(sample variance on (1, 5) ns / n) divided by the slope d(mean)/d(tau) of the moment curve
  boot       shard bootstrap, B = 1000, seed 1, standard deviation (ddof 1)
  ratio      boot / analytic, and its mean over the eight tubes
Half split: tubes read separately from the even and the odd shards (shard index 0, 2, ... against 1, 3, ...),
z = (tau_even - tau_odd) / (2 sigma_analytic of the full sample), chi2 = sum z^2 over the eight tubes, 8 dof.
Writes shard_independence.json next to this script (p_value is added to the published fields)."""
import json, os
from pathlib import Path
import numpy as np
from scipy.stats import chi2

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
Z = np.load(RES / "2026-10-14_known_randoms/known_randoms_hists.npz")
ik, ir, ic = list(Z["ks"]).index(1.0), list(Z["rules"]).index("all"), list(Z["classes"]).index("G")
NAMES = list(Z["voi_names"]); TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
G = Z["hist"][:, ik, ir][:, TI, ic].astype(float)              # [shard, tube, bin]
B_ = Z["bins"]; bc = 0.5 * (B_[:-1] + B_[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
G = G[..., sel]; S = G.shape[0]
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)
E = np.exp(-t[None, :] / tg[:, None]); CV = (E * t).sum(1) / E.sum(1); DCV = np.gradient(CV, tg)


def readout(h):
    n = h.sum(-1); m = (h * x).sum(-1) / n
    tau = np.interp(m, CV, tg)
    sd = np.sqrt((h * (x - m[..., None]) ** 2).sum(-1) / n) / np.sqrt(n) / np.interp(tau, tg, DCV)
    return tau, sd, n


tau, sd, n = readout(G.sum(0))
rng = np.random.default_rng(1); B = 1000
W = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(B)]).astype(float)
TB = readout((W @ G.reshape(S, -1)).reshape((B,) + G.shape[1:]))[0]
boot = TB.std(0, ddof=1)
even = np.arange(S) % 2 == 0
z = (readout(G[even].sum(0))[0] - readout(G[~even].sum(0))[0]) / (2 * sd)
c2 = float((z ** 2).sum())
out = dict(per_tube=[dict(tube=NAMES[i], n=int(n[j]), boot_ps=1000 * float(boot[j]), analytic_ps=1000 * float(sd[j]),
                          ratio=float(boot[j] / sd[j])) for j, i in enumerate(TI)],
           mean_ratio=float((boot / sd).mean()), half_split_z=[float(v) for v in z], chi2=c2,
           p_value=float(chi2.sf(c2, len(TI))))
json.dump(out, open(Path(__file__).resolve().parent / "shard_independence.json", "w"), indent=1)
print(f"boot/analytic {(boot / sd).min():.2f} to {(boot / sd).max():.2f}, mean {(boot / sd).mean():.2f}; "
      f"odd/even chi2 {c2:.1f} for {len(TI)} dof, p = {out['p_value']:.2f}")
