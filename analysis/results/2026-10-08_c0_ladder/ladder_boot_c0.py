#!/usr/bin/env python3
"""C0 rate ladder: joint shard bootstrap of the exponent, identical in method to
2026-09-30_referee2/ladder_boot.py (truncated-mean readout, mean over the tube set of
|b_closest| - |b_all|, OLS of ln mean|DiD| on ln mu, B = 4000 shard resamples applied to all
rungs, GLS with the bootstrap covariance). mu from mu_ladder_c0.json (all shards).
Tube sets: the six tubes common to all discovery-stream rungs, and all eight."""
import json, os
from pathlib import Path
import numpy as np

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
D = RES / "2026-10-08_c0_ladder"
Z = np.load(D / "ladder_joint/ladder_joint_hists.npz")
H = Z["hist"].astype(np.float64)
KS, BINS, NAMES, TCFG = Z["ks"], Z["bins"], list(Z["voi_names"]), Z["tau_cfg"]
MU = json.loads((D / "mu_ladder_c0.json").read_text())
mu = np.array([MU[f"{k:.2f}"]["mu"] for k in KS])
DISC = json.loads((RES / "2026-08-26_rate_scan/ladder12.json").read_text())
BOOT_DISC = json.loads((RES / "2026-09-30_referee2/ladder_joint/ladder_boot.json").read_text())
TAU_MIN, W_PROMPT = 1.0, 5.0
W = W_PROMPT - TAU_MIN
bc = 0.5 * (BINS[:-1] + BINS[1:])
sel = (bc > TAU_MIN) & (bc < W_PROMPT)
lam_g = np.geomspace(0.02, 50.0, 400001); u = W / lam_g
obs_g = lam_g - W * np.exp(-u) / (-np.expm1(-u))
def trunc_inv(t):
    return np.where(t > 0, np.interp(t, obs_g, lam_g, left=np.nan, right=np.nan), np.nan)
def tau_from(Hs):
    hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]; ns = hs.sum(-1)
    return trunc_inv((hs * bc[sel]).sum(-1) / np.where(ns > 0, ns, np.nan) - TAU_MIN)
X = np.c_[np.log(mu), np.ones_like(mu)]
def slope(y):
    return np.linalg.lstsq(X, np.log(y), rcond=None)[0]

out = dict(ks=KS.tolist(), mu=mu.tolist(), mu_over_k=(mu / KS).tolist(), seed=20261008, sets={})
for label, tubes in (("six_common", DISC["common_tubes"]), ("all_eight", [n for n in NAMES if n != "bg"])):
    ci = [NAMES.index(t) for t in tubes]
    def did(Hs):
        tau = tau_from(Hs)
        b = 100.0 * (tau[..., ci] - TCFG[ci]) / TCFG[ci]
        per = np.abs(b[:, 0]) - np.abs(b[:, 1])
        return per.mean(-1), per, b
    d0, per0, b0 = did(H.sum(0))
    beta0, icpt0 = slope(d0) if np.all(d0 > 0) else (np.nan, np.nan)
    rng = np.random.default_rng(20261008); S = H.shape[0]; B = 4000
    flat = H.reshape(S, -1); betas, lny = [], []
    for c in range(0, B, 250):
        nb = min(250, B - c)
        w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(nb)]).astype(np.float64)
        Hb = (w @ flat).reshape((nb,) + H.shape[1:])
        for i in range(nb):
            y, _, _ = did(Hb[i])
            if np.all(np.isfinite(y)) and np.all(y > 0):
                betas.append(slope(y)[0]); lny.append(np.log(y))
    betas, lny = np.array(betas), np.array(lny)
    lo, hi = np.percentile(betas, [2.5, 97.5])
    C = np.cov(lny.T); Ci = np.linalg.inv(C); V = np.linalg.inv(X.T @ Ci @ X); bg = V @ X.T @ Ci @ np.log(d0)
    corr = C / np.sqrt(np.outer(np.diag(C), np.diag(C))); iu = np.triu_indices(len(KS), 1)
    pref = d0 / mu
    out["sets"][label] = dict(tubes=tubes, mean_did=d0.tolist(), mean_did_sd=np.std(np.exp(lny), axis=0, ddof=1).tolist(),
        per_tube_did=per0.tolist(), exponent_ols=float(beta0), prefactor_ols=float(np.exp(icpt0)),
        boot_B=int(len(betas)), boot_sd=float(betas.std(ddof=1)), boot_ci95=[float(lo), float(hi)],
        gls_exponent=float(bg[0]), gls_se=float(np.sqrt(V[0, 0])), corr_median=float(np.median(corr[iu])),
        did_over_mu=pref.tolist(), did_over_mu_rel_spread=float(pref.std() / pref.mean()))
    print(f"[{label}] n={len(tubes)}  OLS exponent {beta0:.4f}  prefactor {np.exp(icpt0):.1f}  "
          f"boot 95% [{lo:.4f}, {hi:.4f}] (B={len(betas)})  GLS {bg[0]:.4f}+/-{np.sqrt(V[0,0]):.4f}  corr med {np.median(corr[iu]):.3f}")
    for k, m, y in zip(KS, mu, d0):
        print(f"   k={k:4.2f} mu={m:.5f} mean|DiD|={y:7.3f}  /mu={y/m:7.1f}")
# discovery-stream comparison at matching rungs
dd = np.array(BOOT_DISC["mean_did"]); out["discovery_mean_did"] = dd.tolist()
out["ratio_c0_over_discovery_six"] = (np.array(out["sets"]["six_common"]["mean_did"]) / dd).tolist()
print("C0/discovery mean|DiD| per rung (six tubes):", np.round(out["ratio_c0_over_discovery_six"], 3))
json.dump(out, open(D / "ladder_boot_c0.json", "w"), indent=1)
print("wrote", D / "ladder_boot_c0.json")
