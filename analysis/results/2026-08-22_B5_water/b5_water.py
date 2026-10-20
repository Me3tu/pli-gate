#!/usr/bin/env python3
"""B5 - the water cylinder under both timing models, from data already on disk.

Comparison axis: recovered o-Ps lifetime in the VOI_15x25 water region, against
the configured 1.84 ns and Steinberger et al. 2024's published 1.815 ns. The two
timing models are the production one (direct annihilation at exactly 0 ps) and
Mercolli's fitted tau_2 = 388 ps. Both matching rules are shown for each. Nothing
is re-simulated: arm A is the archived O4 water run, arm B the PLI_DIRECT_TAU_PS
=388 run of 2026-08-19.
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
V2 = RES / "2026-08-19_V2_earlytime"
V4 = RES / "2026-08-19_V4_stats"
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa: E402

BC = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:])
SEL = (BC > X.TAU_MIN_NS) & (BC < X.W_PROMPT_NS)
W = X.W_PROMPT_NS - X.TAU_MIN_NS
VOI = "VOI_15x25"
CFG, PUB = 1.84, 1.815


def tau_hat(hs, hr):
    h = hs[SEL] - hr[SEL]; n = h.sum()
    if n <= 0: return float("nan")
    return X.truncation_correct(float((BC[SEL] * h).sum() / n) - X.TAU_MIN_NS, W)


def arm(npz, label):
    z = np.load(npz, allow_pickle=True)
    names = [str(s) for s in z["voi_names"]]
    i = names.index(VOI)
    out = {}
    for rule, ks, kr in (("closest", "closest_sig", "closest_rnd"),
                         ("all", "all_sig", "all_rnd")):
        S = z[ks][:, i, :]; R = z[kr][:, i, :]
        hs, hr = S.sum(0), R.sum(0)
        ns = S.shape[0]
        rng = np.random.default_rng(20260823)
        bs = np.array([tau_hat(S[k].sum(0), R[k].sum(0))
                       for k in (rng.integers(0, ns, ns) for _ in range(800))])
        out[rule] = dict(tau=tau_hat(hs, hr), sd=float(np.nanstd(bs, ddof=1)),
                         n_sig=int(hs[SEL].sum()), n_rnd=int(hr[SEL].sum()),
                         spectrum_sig=hs.tolist(), spectrum_rnd=hr.tolist())
        print(f"  {label:10s} {rule:8s} tau = {out[rule]['tau']:.4f} +/- {out[rule]['sd']:.4f} ns"
              f"   signal {out[rule]['n_sig']:,}  randoms {out[rule]['n_rnd']:,}")
    return out


print(f"water {VOI}; configured {CFG} ns, Steinberger et al. 2024 published {PUB} ns\n")
res = {"bin_centres_ns": BC.tolist(), "fit_window_ns": [X.TAU_MIN_NS, X.W_PROMPT_NS],
       "voi": VOI, "tau_configured_ns": CFG, "tau_published_ns": PUB,
       "arms": {}}
res["arms"]["0 ps"] = arm(V4 / "hists_water.npz", "0 ps")
res["arms"]["388 ps"] = arm(V2 / "hists_water388.npz", "388 ps")
a0 = res["arms"]["0 ps"]; a3 = res["arms"]["388 ps"]
print(f"\n  effect of the realistic direct component, closest-prompt: "
      f"{1000*(a3['closest']['tau'] - a0['closest']['tau']):+.0f} ps")
print(f"  effect of the realistic direct component, all-prompts   : "
      f"{1000*(a3['all']['tau'] - a0['all']['tau']):+.0f} ps")
Path("b5_water.json").write_text(json.dumps(res))
print("\nwrote b5_water.json")
