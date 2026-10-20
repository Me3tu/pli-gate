#!/usr/bin/env python3
"""Binomial maximum-likelihood edge fits (replacing the least-squares ratio fit of analyze.py, whose
data-based errors biased the width), prompt-origin split, and the all-prompts over-subtraction check.
Input factorization_hists2.npz: H[shard, arm(0 signal,1 delayed), class(0 true,1 same-tube prompt,
2 other prompt), sel(0 |s|<=9, 1 |s|<=5), voi, bin]."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import erfc
from scipy.stats import norm

OUT = Path(__file__).resolve().parent
Z = np.load(OUT / "factorization_hists2.npz", allow_pickle=True)
H = Z["H"].astype(float); E = Z["edges"]; c = 0.5 * (E[1:] + E[:-1]); names = list(Z["voi"]); tcfg = [float(x) for x in Z["tau_cfg"]]
NS, NV = H.shape[0], H.shape[4]
# consistency with the first run: true-triple histograms must be identical (same builder, seed, shards)
H1 = np.load(OUT / "factorization_hists.npz")["H"]
assert np.array_equal(H1[:, 0, 1], H[:, 0, 0].astype(H1.dtype)), "true-triple histograms differ from run 1"
FIT = (c > 3.0) & (c < 7.0)


def ml_edge(wide, inside):
    m = FIT & (wide > 0); n = wide[m]; k = inside[m]; t = c[m]
    def nll(p):
        s, d = p
        if s <= 0.05: return 1e30
        A = np.clip(norm.cdf((5 + d - t) / s), 1e-12, 1 - 1e-12)
        return -float((k * np.log(A) + (n - k) * np.log(1 - A)).sum())
    r = minimize(nll, [0.35, 0.03], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-9, maxiter=4000))
    s, d = r.x
    z = (5 + d - t) / s; A = np.clip(norm.cdf(z), 1e-12, 1 - 1e-12); phi = norm.pdf(z)
    dA_ds = -phi * z / s; dA_dd = phi / s
    w = n / (A * (1 - A))
    I = np.array([[(w * dA_ds ** 2).sum(), (w * dA_ds * dA_dd).sum()], [(w * dA_ds * dA_dd).sum(), (w * dA_dd ** 2).sum()]])
    cov = np.linalg.inv(I)
    chi2 = float((((k - n * A) ** 2) / (n * A * (1 - A) + 1e-12)).sum())
    return dict(sigma_ps=1000 * s, sigma_err_ps=1000 * np.sqrt(cov[0, 0]), shift_ps=1000 * d, shift_err_ps=1000 * np.sqrt(cov[1, 1]),
                pearson_chi2=chi2, nbins=int(m.sum()), n=float(n.sum()))


def sumv(h, v):  # h [shard, ..., voi, bin] -> summed over shards and the chosen voi (or all)
    x = h.sum(0)
    return x.sum(-2) if v is None else x[..., v, :]


res = dict(_doc=__doc__, pooled={}, tubes={})
for v in [None] + list(range(NV)):
    key = "pooled" if v is None else names[v]
    r = {}
    tw, ti = sumv(H[:, 0, 0, 0], v), sumv(H[:, 0, 0, 1], v)
    r["true"] = ml_edge(tw, ti)
    if v is None:
        for arm, an in ((0, "signal_accidental"), (1, "delayed")):
            for cl, cn in ((1, "same_tube_prompt"), (2, "other_prompt"), (None, "all")):
                hw = H[:, arm, 1:, 0].sum(1) if cl is None else H[:, arm, cl, 0]
                hi = H[:, arm, 1:, 1].sum(1) if cl is None else H[:, arm, cl, 1]
                if arm == 1 and cl is None:
                    hw = H[:, 1, :, 0].sum(1); hi = H[:, 1, :, 1].sum(1)   # delayed arm: every candidate is accidental
                r[f"{an}_{cn}"] = ml_edge(sumv(hw, v), sumv(hi, v))
        # share of accidental prompts originating in the same tube as the pair (signal arm, inside window, 1-5 ns)
        m15 = (c > 1) & (c < 5)
        same = sumv(H[:, 0, 1, 1], v)[m15].sum(); oth = sumv(H[:, 0, 2, 1], v)[m15].sum()
        r["accidental_same_tube_share_1_5"] = float(same / (same + oth))
        # all-prompts over-subtraction: delayed-arm count minus accidental triples of the signal window, 1-5 ns
        dl = H[:, 1, :, 1].sum(1)[..., m15].sum(-1).sum(-1)      # per shard
        ac = H[:, 0, 1:, 1].sum(1)[..., m15].sum(-1).sum(-1)
        tr = H[:, 0, 0, 1][..., m15].sum(-1).sum(-1)
        ex = dl - ac
        rng = np.random.default_rng(11); bs = [ex[i].sum() for i in rng.integers(0, NS, size=(1000, NS))]
        r["oversubtraction_1_5"] = dict(excess=float(ex.sum()), excess_sd=float(np.std(bs, ddof=1)), delayed=float(dl.sum()),
                                        true=float(tr.sum()), excess_over_true=float(ex.sum() / tr.sum()))
    res["pooled" if v is None else "tubes"] = r if v is None else res["tubes"]
    if v is not None:
        res["tubes"][key] = r
sw = np.array([res["tubes"][n]["true"]["sigma_ps"] for n in names]); se = np.array([res["tubes"][n]["true"]["sigma_err_ps"] for n in names])
w = 1 / se ** 2; mw = float((w * sw).sum() / w.sum())
res["true_width_homogeneity"] = dict(ivw_ps=mw, chi2=float((w * (sw - mw) ** 2).sum()), dof=NV - 1,
                                    r_with_tau_cfg=float(np.corrcoef(sw, tcfg)[0, 1]), min_ps=float(sw.min()), max_ps=float(sw.max()))
json.dump(res, open(OUT / "factorization_result2.json", "w"), indent=1)
P = res["pooled"]
for k, x in P.items():
    if isinstance(x, dict) and "sigma_ps" in x:
        print(f"{k:32s} sigma {x['sigma_ps']:.1f} +- {x['sigma_err_ps']:.1f}  shift {x['shift_ps']:+.1f} +- {x['shift_err_ps']:.1f}  chi2 {x['pearson_chi2']:.1f}/{x['nbins']}  n {x['n']:.0f}")
print("same-tube share of accidental prompts", round(P["accidental_same_tube_share_1_5"], 3))
print("oversubtraction", P["oversubtraction_1_5"])
print("width homogeneity", res["true_width_homogeneity"])
for n in names:
    t = res["tubes"][n]["true"]
    print(f"{n:30s} {t['sigma_ps']:.0f}+-{t['sigma_err_ps']:.0f}  {t['shift_ps']:+.0f}+-{t['shift_err_ps']:.0f}  chi2 {t['pearson_chi2']:.0f}/{t['nbins']}")
