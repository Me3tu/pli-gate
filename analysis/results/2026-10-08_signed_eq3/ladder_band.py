#!/usr/bin/env python3
"""Review round 6 (2026-10-08): pointwise bootstrap band of the fitted rate law for both
rate replays, carrying the prefactor and the exponent together (the band drawn so far was
the exponent interval pivoted at the data centroid).

Identical in method and random stream to 2026-09-30_referee2/ladder_boot.py (discovery,
seed 20260930) and 2026-10-08_c0_ladder/ladder_boot_c0.py (corrected yield, six common
tubes, seed 20261008): the same shard resamples are drawn in the same order, so the
exponent intervals are reproduced exactly; here the intercept of every replicate is kept
as well, and the 2.5/97.5 percentiles of ln(prefactor) + exponent * ln(mu) are taken on a
grid of mu. Also reported: the same fits on the SIGNED per-tube rule difference
b_all - b_closest, and how often a tube's all-prompts bias is positive (where the
absolute-value statistic folds). Post hoc (review), not registered.
"""
import json, os
from pathlib import Path
import numpy as np

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = RES / "2026-10-08_signed_eq3"
COMMON = json.loads((RES / "2026-08-26_rate_scan/ladder12.json").read_text())["common_tubes"]
STREAMS = {
    "discovery": dict(npz=RES / "2026-09-30_referee2/ladder_joint/ladder_joint_hists.npz",
                      mu=RES / "2026-08-26_rate_scan/mu_ladder.json", seed=20260930,
                      ref=RES / "2026-09-30_referee2/ladder_joint/ladder_boot.json"),
    "corrected": dict(npz=RES / "2026-10-08_c0_ladder/ladder_joint/ladder_joint_hists.npz",
                      mu=RES / "2026-10-08_c0_ladder/mu_ladder_c0.json", seed=20261008,
                      ref=RES / "2026-10-08_c0_ladder/ladder_boot_c0.json"),
}
TAU_MIN, W_PROMPT = 1.0, 5.0
W = W_PROMPT - TAU_MIN
lam_g = np.geomspace(0.02, 50.0, 400001); u = W / lam_g
obs_g = lam_g - W * np.exp(-u) / (-np.expm1(-u))
GRID = np.geomspace(0.003, 0.13, 200)
out = dict(note=__doc__, mu_grid=GRID.tolist(), streams={})
for name, cfg in STREAMS.items():
    Z = np.load(cfg["npz"]); H = Z["hist"].astype(np.float64)
    KS, BINS, NAMES, TCFG = Z["ks"], Z["bins"], list(Z["voi_names"]), Z["tau_cfg"]
    MU = json.loads(cfg["mu"].read_text()); mu = np.array([MU[f"{k:.2f}"]["mu"] for k in KS])
    bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > TAU_MIN) & (bc < W_PROMPT)
    ci = [NAMES.index(t) for t in COMMON]
    def tau_from(Hs):
        hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]; ns = hs.sum(-1)
        t = (hs * bc[sel]).sum(-1) / np.where(ns > 0, ns, np.nan) - TAU_MIN
        return np.where(t > 0, np.interp(t, obs_g, lam_g, left=np.nan, right=np.nan), np.nan)
    def stats(Hs):
        tau = tau_from(Hs); b = 100.0 * (tau[..., ci] - TCFG[ci]) / TCFG[ci]
        folded = (np.abs(b[:, 0]) - np.abs(b[:, 1])).mean(-1)
        signed = (b[:, 1] - b[:, 0]).mean(-1)
        return folded, signed, b
    X = np.c_[np.log(mu), np.ones_like(mu)]
    fit = lambda y: np.linalg.lstsq(X, np.log(y), rcond=None)[0]
    f0, s0, b0 = stats(H.sum(0))
    p0f, p0s = fit(f0), fit(s0)
    rng = np.random.default_rng(cfg["seed"]); S = H.shape[0]; B = 4000; flat = H.reshape(S, -1)
    reps_f, reps_s, npos = [], [], 0
    for c in range(0, B, 250):
        nb = min(250, B - c)
        w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(nb)]).astype(np.float64)
        Hb = (w @ flat).reshape((nb,) + H.shape[1:])
        for i in range(nb):
            yf, ys, bb = stats(Hb[i])
            if np.all(np.isfinite(yf)) and np.all(yf > 0):
                reps_f.append(fit(yf))
            if np.all(np.isfinite(ys)) and np.all(ys > 0):
                reps_s.append(fit(ys))
            npos += int((bb[:, 1] > 0).sum())
    reps_f, reps_s = np.array(reps_f), np.array(reps_s)
    ref = json.loads(cfg["ref"].read_text())
    ref_ci = ref["boot_ci95"] if "boot_ci95" in ref else ref["sets"]["six_common"]["boot_ci95"]
    def band(reps, p0):
        lny = reps[:, 1][:, None] + reps[:, 0][:, None] * np.log(GRID)[None, :]
        lo, hi = np.percentile(lny, [2.5, 97.5], axis=0)
        cen = p0[1] + p0[0] * np.log(GRID)
        return dict(center=np.exp(cen).tolist(), lo=np.exp(lo).tolist(), hi=np.exp(hi).tolist(),
                    rel_halfwidth_at_data=[float(x) for x in
                        0.5 * (np.exp(np.interp(np.log(mu), np.log(GRID), hi)) - np.exp(np.interp(np.log(mu), np.log(GRID), lo)))
                        / np.exp(p0[1] + p0[0] * np.log(mu))])
    ci_f = np.percentile(reps_f[:, 0], [2.5, 97.5]); ci_s = np.percentile(reps_s[:, 0], [2.5, 97.5])
    cov = np.cov(reps_f.T)
    out["streams"][name] = dict(
        mu=mu.tolist(), folded_mean=f0.tolist(), signed_mean=s0.tolist(),
        folded=dict(exponent=float(p0f[0]), prefactor=float(np.exp(p0f[1])), ci95=ci_f.tolist(),
                    reference_ci95=ref_ci, B=int(len(reps_f)),
                    corr_exponent_lnprefactor=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
                    band=band(reps_f, p0f)),
        signed=dict(exponent=float(p0s[0]), prefactor=float(np.exp(p0s[1])), ci95=ci_s.tolist(), B=int(len(reps_s)),
                    band=band(reps_s, p0s)),
        n_positive_all_prompts_point=int((b0[:, 1] > 0).sum()), n_tube_rungs=int(b0[:, 1].size),
        frac_positive_all_prompts_boot=float(npos / (B * b0[:, 1].size)),
        max_abs_folded_minus_signed=float(np.abs(f0 - s0).max()))
    o = out["streams"][name]
    print(f"[{name}] folded exponent {p0f[0]:.4f} CI [{ci_f[0]:.4f}, {ci_f[1]:.4f}] (reference {ref_ci[0]:.4f}, {ref_ci[1]:.4f}) "
          f"prefactor {np.exp(p0f[1]):.1f} corr(beta, ln A) {o['folded']['corr_exponent_lnprefactor']:+.3f}")
    print(f"   signed exponent {p0s[0]:.4f} CI [{ci_s[0]:.4f}, {ci_s[1]:.4f}] prefactor {np.exp(p0s[1]):.1f}; "
          f"max |folded - signed| {o['max_abs_folded_minus_signed']:.4f} pp; positive all-prompts biases at the point "
          f"{o['n_positive_all_prompts_point']}/{o['n_tube_rungs']}, in the bootstrap {100*o['frac_positive_all_prompts_boot']:.2f}%")
    print("   band rel half-width at the rungs:", np.round(o["folded"]["band"]["rel_halfwidth_at_data"], 4))
json.dump(out, open(OUT / "ladder_band.json", "w"), indent=1)
print("wrote", OUT / "ladder_band.json")
