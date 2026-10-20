#!/usr/bin/env python3
"""Analysis of rules16_hists.npz (see rules16.py), with the readouts of 2026-10-14_known_randoms/known_randoms_analysis.py.

Per rung (k = 0.5, 1.0, 2.0; mu from 2026-10-08_c0_ladder/mu_ladder_c0.json), arm and tube, fit window 1 < delay < 5 ns:
  bsub_sharp / bsub_acc359   residual of the randoms-subtracted windowed-moment readout against the configured lifetime
                             (sharp-edge production inversion, and inversion with the a priori acceptance 358.97 ps)
  dR_*                       randoms-correction error tau(sig - rnd) - tau(G + PA_corr), percent of configured
  bG_*                       genuine-triple residual
  phi_known                  (randoms arm - AP - PA_unc) / G, percent: departure of the delayed-window estimate from the
                             uncorrelated accidental content of the signal window
  rule_shift_*               signed difference b(all) - b(arm) on identical events (adjacent builder; for the all-pairs
                             arms, against ap_all), percentage points, the quantity of Eq. (4)
Pooled values: unweighted means over the eight tubes; mean absolute residual over the eight tubes.
Uncertainties: shard bootstrap B = 1000, seed 20261019, the same resamples for all arms and rungs (paired).
Writes rules16_summary.json."""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
D = RES / "2026-10-19_review16"
Z = np.load(D / "rules16_hists.npz")
H = Z["hist"].astype(np.float64)                      # [shard, k, arm, voi+1, cls, bin]
KS, ARMS, CLS = list(Z["ks"]), list(Z["arms"]), list(Z["classes"])
NAMES, TCFG, BINS = list(Z["voi_names"]), Z["tau_cfg"], Z["bins"]
MU = json.loads((RES / "2026-10-08_c0_ladder/mu_ladder_c0.json").read_text())
mu = np.array([MU[f"{k:.2f}"]["mu"] for k in KS])
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def moment_curve(sig_e):
    A = np.ones_like(t) if sig_e is None else Phi((5 - t) / sig_e) * Phi((t + 5) / sig_e)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


CURVES = {"sharp": moment_curve(None), "acc359": moment_curve(0.35897)}
iS, iR, iG, iAP, iPc, iPu = (CLS.index(c) for c in ("sig", "rnd", "G", "AP", "PA_corr", "PA_unc"))
iAll, iApAll = ARMS.index("all"), ARMS.index("ap_all")


def readouts(Hs):
    h = Hs[..., sel]

    def tau(spec, curve):
        n = spec.sum(-1); m = (spec * x).sum(-1) / np.where(n > 0, n, np.nan)
        return np.interp(m, CURVES[curve], tg, left=np.nan, right=np.nan)
    sub = h[..., iS, :] - h[..., iR, :]; gp = h[..., iG, :] + h[..., iPc, :]; g = h[..., iG, :]
    out = {}
    for cv in CURVES:
        ts, tgp, tgg = tau(sub, cv), tau(gp, cv), tau(g, cv)
        out[f"bsub_{cv}"] = 100 * (ts - TCFG) / TCFG
        out[f"dR_{cv}"] = 100 * (ts - tgp) / TCFG
        out[f"bG_{cv}"] = 100 * (tgg - TCFG) / TCFG
    nG = h[..., iG, :].sum(-1)
    out["phi_known"] = 100 * (h[..., iR, :].sum(-1) - h[..., iAP, :].sum(-1) - h[..., iPu, :].sum(-1)) / nG
    out["frac_PAc"] = 100 * h[..., iPc, :].sum(-1) / nG
    for cv in CURVES:
        b = out[f"bsub_{cv}"]                         # [..., arm, voi]
        ref = np.stack([b[..., iAll if not ARMS[a].startswith("ap_") else iApAll, :] for a in range(len(ARMS))], axis=-2)
        out[f"rule_shift_{cv}"] = ref - b
    return out


P0 = readouts(H.sum(0))                               # [k, arm, voi+1]
rng = np.random.default_rng(20261019); S = H.shape[0]; B = 1000; flat = H.reshape(S, -1)
keys = list(P0.keys()); boots = {k: [] for k in keys}
for c in range(0, B, 50):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    pb = readouts((w @ flat).reshape((50,) + H.shape[1:]))
    for k in keys: boots[k].append(pb[k])
boots = {k: np.concatenate(v) for k, v in boots.items()}
sd = lambda a: 0.5 * (np.nanpercentile(a, 84, axis=0) - np.nanpercentile(a, 16, axis=0))

out = dict(doc=__doc__, ks=KS, mu=mu.tolist(), arms=ARMS, tubes=[NAMES[i] for i in TI], B=B, per=[])
for ia, arm in enumerate(ARMS):
    for ik, k in enumerate(KS):
        e = dict(arm=arm, k=k, mu=float(mu[ik]))
        for key in keys:
            v0 = float(np.nanmean(P0[key][ik, ia, TI])); vb = np.nanmean(boots[key][:, ik, ia, TI], axis=-1)
            e[key] = [v0, float(sd(vb))]
            if key.startswith("bsub"):
                a0 = float(np.nanmean(np.abs(P0[key][ik, ia, TI]))); ab = np.nanmean(np.abs(boots[key][:, ik, ia, TI]), axis=-1)
                e["mae_" + key] = [a0, float(sd(ab))]
                e[key + "_range"] = [float(np.nanmin(P0[key][ik, ia, TI])), float(np.nanmax(P0[key][ik, ia, TI]))]
            if key.startswith("rule_shift"):
                e[key + "_range"] = [float(np.nanmin(P0[key][ik, ia, TI])), float(np.nanmax(P0[key][ik, ia, TI]))]
        e["per_tube_bsub_sharp"] = [[float(P0["bsub_sharp"][ik, ia, i]), float(sd(boots["bsub_sharp"][:, ik, ia, i]))] for i in TI]
        out["per"].append(e)
# rate slopes of the pooled randoms-correction error and of phi_known, per 0.01 in mu
out["slopes"] = {}
for ia, arm in enumerate(ARMS):
    for key in ("dR_sharp", "dR_acc359", "phi_known", "rule_shift_sharp"):
        y0 = np.nanmean(P0[key][:, ia, TI], axis=-1); yb = np.nanmean(boots[key][:, :, ia, TI], axis=-1)
        X = np.c_[mu / 0.01, np.ones_like(mu)]
        s0 = np.linalg.lstsq(X, y0, rcond=None)[0]; sb = np.array([np.linalg.lstsq(X, y, rcond=None)[0] for y in yb])
        out["slopes"].setdefault(arm, {})[key] = [float(s0[0]), float(sd(sb[:, 0]))]
json.dump(out, open(D / "rules16_summary.json", "w"), indent=1)
for e in out["per"]:
    print(f"{e['arm']:11s} k {e['k']:.1f} mu {e['mu']:.4f} | mae sharp {e['mae_bsub_sharp'][0]:5.2f} acc {e['mae_bsub_acc359'][0]:5.2f} "
          f"| signed acc {e['bsub_acc359'][0]:+5.2f}+-{e['bsub_acc359'][1]:.2f} | shift sharp {e['rule_shift_sharp'][0]:+6.2f}+-{e['rule_shift_sharp'][1]:.2f} "
          f"| dR acc {e['dR_acc359'][0]:+5.2f}+-{e['dR_acc359'][1]:.2f} | phi {e['phi_known'][0]:+5.2f}+-{e['phi_known'][1]:.2f} "
          f"| bG acc {e['bG_acc359'][0]:+5.2f} | PAc {e['frac_PAc'][0]:.3f}")
print(json.dumps(out["slopes"]))
