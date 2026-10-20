#!/usr/bin/env python3
"""Does the prompt-window acceptance explain the all-prompts residual of the windowed-moment
estimator?  The production truncation correction assumes a sharp upper edge at T = 5 ns in
the corrected delay. The window is applied to the uncorrected pair-prompt time, so in the
corrected delay the upper edge is softened by the geometric path spread (sigma_e).

For each Sim-1 tube (full 128 shards, 303 ps, both rules) the windowed mean is inverted
 (i) with the production sharp-edge moment (published),
 (ii) with the acceptance-weighted moment m_A(tau) = int t e^{-t/tau} A(t) / int e^{-t/tau} A(t)
     on (1, 5) ns, A(t) = Phi((5-t)/sigma_e) Phi((t+5)/sigma_e), for sigma_e = 359 ps
     (a-priori geometric+timing value used by the component fit) and 290.6 ps (measured on the
     flat randoms arm). No parameter is fitted to the tube data.
Uncertainties: shard bootstrap (B = 1000), paired across variants.
"""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc

D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2]) / "2026-09-30_feedback3/c0_sim1"  # feedback round 3: C0
Z = np.load(D / "photon_conventions_c0/photon_conventions.npz")
H = Z["hist"][:, 0].astype(float)          # photon-1 (identical to all conventions) [shard, rule, voi+1, 2, bin]
BINS, NAMES, TCFG = Z["bins"], list(Z["voi_names"]), Z["tau_cfg"]
bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
tg = np.linspace(0.3, 12.0, 6000)
t = np.linspace(1.0, 5.0, 8001)

def moment_curve(sig_e, d=0.0):
    A = np.ones_like(t) if sig_e is None else Phi((5 + d - t) / sig_e) * Phi((t + 5 - d) / sig_e)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)          # conditional mean on (1,5)

# acceptance measured on this dataset's own randoms arm (upper edge, all-prompts rule, all regions)
from scipy.optimize import curve_fit
_r = H.sum(0)[1, :, 1, :].sum(0); _m = (bc > 2.5) & (bc < 7.5)
_p, _c = curve_fit(lambda z, a, s, d: a * Phi((5 + d - z) / s), bc[_m], _r[_m],
                   p0=[_r[(bc > 1) & (bc < 4)].mean(), 0.3, 0.0], sigma=np.sqrt(np.clip(_r[_m], 1, None)))
print(f"randoms-arm upper edge: sigma_e {1000*_p[1]:.1f} +- {1000*np.sqrt(_c[1,1]):.1f} ps, shift {1000*_p[2]:+.1f} ps")
CURVES = {"sharp": moment_curve(None), "acc_measured": moment_curve(_p[1], _p[2]),
          "acc359": moment_curve(0.35897), "acc291": moment_curve(0.2906)}
for k, m in CURVES.items():
    assert np.all(np.diff(m) > 0), k

def taus(Hs):
    hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]
    mean = (hs * x).sum(-1) / hs.sum(-1)
    return {k: np.interp(mean, m, tg, left=np.nan, right=np.nan) for k, m in CURVES.items()}

T0 = taus(H.sum(0))
rng = np.random.default_rng(4); S = H.shape[0]; B = 1000
flat = H.reshape(S, -1)
boots = {k: [] for k in CURVES}
for c0 in range(0, B, 250):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(250)]).astype(float)
    tb = taus((w @ flat).reshape((250,) + H.shape[1:]))
    for k in CURVES: boots[k].append(tb[k])
boots = {k: np.concatenate(v) for k, v in boots.items()}
out = {}
for ir, rule in enumerate(["closest", "all"]):
    print(f"\nrule = {rule}")
    print(f"  {'tube':28s} {'cfg':>6} | {'sharp':>14} | {'acc measured':>14} | {'acc 359 ps':>14} | {'acc 291 ps':>14}")
    rows = []
    for iv, nm in enumerate(NAMES):
        if nm == "bg": continue
        r = {"tube": nm, "tau_cfg": float(TCFG[iv])}
        s = []
        for k in CURVES:
            b = 100 * (T0[k][ir, iv] - TCFG[iv]) / TCFG[iv]
            e = 100 * np.nanstd(boots[k][:, ir, iv], ddof=1) / TCFG[iv]
            r[k] = [float(b), float(e)]; s.append(f"{b:+6.2f}±{e:4.2f}")
        rows.append(r)
        print(f"  {nm:28s} {TCFG[iv]:6.3f} | " + " | ".join(f"{v:>14}" for v in s))
    for k in CURVES:
        v = np.array([r[k][0] for r in rows])
        bm = np.array([np.mean(np.abs(100 * (boots[k][:, ir, iv] - TCFG[iv]) / TCFG[iv]))
                       for iv, nm in enumerate(NAMES) if nm != "bg"])
        mabs = np.mean(np.abs(100 * (boots[k][:, ir, 1:] - TCFG[1:]) / TCFG[1:]), axis=1)
        print(f"  mean |bias| {k:7s}: {np.mean(np.abs(v)):.2f}%  (boot sd {np.std(mabs, ddof=1):.2f}),  mean signed {np.mean(v):+.2f}%")
        out[f"{rule}_{k}_mean_abs"] = float(np.mean(np.abs(v))); out[f"{rule}_{k}_mean_abs_sd"] = float(np.std(mabs, ddof=1))
        out[f"{rule}_{k}_mean_signed"] = float(np.mean(v))
    out[rule] = rows
out["sigma_e_measured_ps"]=float(1000*_p[1]); out["sigma_e_measured_err_ps"]=float(1000*np.sqrt(_c[1,1])); out["edge_shift_ps"]=float(1000*_p[2])
json.dump(out, open(D / "acceptance_tm_c0.json", "w"), indent=1)
print("wrote acceptance_tm.json")
