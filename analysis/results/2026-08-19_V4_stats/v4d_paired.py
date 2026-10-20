#!/usr/bin/env python3
"""V4(d) — paired-arm covariance and its effect on sigma(Delta tau).

The two matching rules run on the SAME shards, so tau_closest and tau_all are
positively correlated and the production uncertainty on their difference,
formed as hypot(sigma_c, sigma_a), IGNORES that correlation and is therefore
conservative. This measures the correlation by resampling SHARDS coherently:
one bootstrap draw selects a set of shards and both arms are recomputed from
exactly those shards.

    var(Delta) = var(c) + var(a) - 2 cov(c, a)

REGISTERED EXPECTATION (the authors', restated): positive correlation, so sigma of
the difference SHRINKS relative to the independent-error assumption. If it
widens instead, stop and report.
"""
import json, os, sys
from pathlib import Path
import numpy as np
S = Path(__file__).parent
sys.path.insert(0, str(Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa

def tau_from_hists(h_sig, h_rnd):
    """subtract_and_mle, reduced to its histogram arithmetic."""
    bc = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:])
    sel = (bc > X.TAU_MIN_NS) & (bc < X.W_PROMPT_NS)
    hs = h_sig[sel] - h_rnd[sel]
    n = hs.sum()
    if n <= 0:
        return np.nan
    obs = float((bc[sel] * hs).sum() / n) - X.TAU_MIN_NS
    return X.truncation_correct(obs, X.W_PROMPT_NS - X.TAU_MIN_NS)

def run(npz, label, B=4000):
    d = np.load(npz, allow_pickle=True)
    names = list(d["voi_names"]); cfg = d["tau_cfg"]
    cs, cr, as_, ar = d["closest_sig"], d["closest_rnd"], d["all_sig"], d["all_rnd"]
    ns = cs.shape[0]
    rng = np.random.default_rng(20260819)
    print(f"\n=== {label}: {ns} shards, {len(names)} VOIs, {B} paired draws ===")
    print(f"{'VOI':>18} {'rho(c,a)':>9} {'sd_c':>7} {'sd_a':>7} {'sd_D paired':>12} "
          f"{'sd_D indep':>11} {'ratio':>6}")
    rows = []
    for ki, nm in enumerate(names):
        if nm == "bg" or not np.isfinite(cfg[ki]):
            continue
        tc = np.empty(B); ta = np.empty(B)
        for b in range(B):
            i = rng.integers(0, ns, ns)
            tc[b] = tau_from_hists(cs[i, ki].sum(0), cr[i, ki].sum(0))
            ta[b] = tau_from_hists(as_[i, ki].sum(0), ar[i, ki].sum(0))
        ok = np.isfinite(tc) & np.isfinite(ta)
        tc, ta = tc[ok], ta[ok]
        sc, sa = tc.std(ddof=1), ta.std(ddof=1)
        rho = float(np.corrcoef(tc, ta)[0, 1])
        sd_paired = (tc - ta).std(ddof=1)
        sd_indep = float(np.hypot(sc, sa))
        rows.append(dict(voi=nm, rho=rho, sd_closest=float(sc), sd_all=float(sa),
                         sd_diff_paired=float(sd_paired), sd_diff_independent=sd_indep,
                         ratio=float(sd_paired / sd_indep), n_draws=int(ok.sum())))
        print(f"{nm:>18} {rho:+9.3f} {sc:7.4f} {sa:7.4f} {sd_paired:12.4f} "
              f"{sd_indep:11.4f} {sd_paired/sd_indep:6.3f}")
    return rows

out = {}
for npz, lab in ((S / "hists_water.npz", "Steinberger water"),
                 (S / "hists_disks.npz", "Steinberger disks")):
    if npz.exists():
        out[lab] = run(npz, lab)
    else:
        print(f"SKIP {lab}: {npz.name} not present")
if out:
    allr = [r for v in out.values() for r in v]
    rr = np.array([r["ratio"] for r in allr]); rho = np.array([r["rho"] for r in allr])
    print(f"\nacross {len(allr)} regions: rho {rho.min():+.3f} to {rho.max():+.3f} "
          f"(median {np.median(rho):+.3f})")
    print(f"paired/independent sigma ratio: {rr.min():.3f} to {rr.max():.3f} "
          f"(median {np.median(rr):.3f})")
    verdict = ("CONSERVATIVE as expected: positive correlation shrinks sigma(Delta)"
               if np.median(rr) < 1 else
               "WIDENS -- contrary to expectation; the independent-error assumption is NOT conservative")
    print(f"VERDICT: {verdict}")
    json.dump(dict(per_phantom=out, median_rho=float(np.median(rho)),
                   median_ratio=float(np.median(rr)), verdict=verdict),
              open(S / "v4d_paired.json", "w"), indent=1)
    print(f"wrote {S/'v4d_paired.json'}")
