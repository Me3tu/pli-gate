#!/usr/bin/env python3
"""Review round 14: the analysis of known_randoms_analysis.py applied to the realistic-timing arm C2
(known_randoms_c2_hists.npz, known_randoms_c2.py). The accidental-prompt rate per rung is taken from the C0 ladder: the
two arms share every event and differ only in the direct delay of the annihilation photons, which leaves the prompt
singles, and so mu, unchanged (the production arms agree to 0.01%). The windowed-moment readout is biased under the
388 ps model by the direct-component leakage (Sec. 3.8), for every class alike, so the randoms-correction error dR, a
difference on the same readout, remains a test of the subtraction. Writes known_randoms_c2_summary.json.

Original docstring: Analysis of known_randoms_hists.npz (see known_randoms_c0.py).

Per rung (accidental-prompt rate mu from 2026-10-08_c0_ladder/mu_ladder_c0.json), rule, and tube, in the fit window
1 < delay < 5 ns (bin centers, as subtract_and_mle):
  counts: signal arm, randoms arm, and the labeled classes G, AP, PA_corr, PA_unc of the signal arm;
  delayed-window test: (randoms arm - AP - PA_unc) / G, the departure of the delayed-window estimate from the
      uncorrelated accidental content it is meant to estimate, and (randoms arm - all non-genuine) / G, the published
      over-subtraction phi (Note S17);
  lifetimes from the windowed mean, inverted with the sharp-edge moment (production) and with the a priori
      acceptance (sigma_e = 358.97 ps), for: signal minus randoms arm (the clinical readout), G + PA_corr (what an exact
      removal of the uncorrelated accidentals leaves), and G (genuine triples only);
  randoms-correction error  dR = tau(sig - rnd) - tau(G + PA_corr)   [percent of the configured lifetime]
  accidental-pair effect    dP = tau(G + PA_corr) - tau(G)
  genuine-triple residual   tau(G) - tau_cfg, with the a priori acceptance.
Uncertainties: shard bootstrap (B = 1000, seed 20261014), the same resample for every rung and rule; pooled values are
unweighted means over the eight tubes (and the six tubes common to all rungs of the discovery replay).
Writes known_randoms_summary.json."""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
D = RES / "2026-10-14_known_randoms"
Z = np.load(D / "known_randoms_c2_hists.npz")
H = Z["hist"].astype(np.float64)                     # [shard, k, rule, voi+1, cls, bin]
KS, RULES, CLS = list(Z["ks"]), list(Z["rules"]), list(Z["classes"])
NAMES, TCFG, BINS = list(Z["voi_names"]), Z["tau_cfg"], Z["bins"]
MU = json.loads((RES / "2026-10-08_c0_ladder/mu_ladder_c0.json").read_text()); mu = np.array([MU[f"{k:.2f}"]["mu"] for k in KS])
COMMON = json.loads((RES / "2026-08-26_rate_scan/ladder12.json").read_text())["common_tubes"]
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]; CI = [NAMES.index(t) for t in COMMON]
bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def moment_curve(sig_e):
    A = np.ones_like(t) if sig_e is None else Phi((5 - t) / sig_e) * Phi((t + 5) / sig_e)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


CURVES = {"sharp": moment_curve(None), "acc359": moment_curve(0.35897)}
iS, iR, iG, iAP, iPc, iPu = (CLS.index(c) for c in ("sig", "rnd", "G", "AP", "PA_corr", "PA_unc"))


def readouts(Hs):
    """Hs[..., voi+1, cls, bin] -> dict of arrays [..., voi+1]."""
    h = Hs[..., sel]
    def tau(spec, curve):
        n = spec.sum(-1); m = (spec * x).sum(-1) / np.where(n > 0, n, np.nan)
        return np.interp(m, CURVES[curve], tg, left=np.nan, right=np.nan)
    sub = h[..., iS, :] - h[..., iR, :]; gp = h[..., iG, :] + h[..., iPc, :]; g = h[..., iG, :]
    out = {}
    for cv in CURVES:
        ts, tgp, tgg = tau(sub, cv), tau(gp, cv), tau(g, cv)
        out[f"dR_{cv}"] = 100 * (ts - tgp) / TCFG; out[f"dP_{cv}"] = 100 * (tgp - tgg) / TCFG
        out[f"bG_{cv}"] = 100 * (tgg - TCFG) / TCFG; out[f"bsub_{cv}"] = 100 * (ts - TCFG) / TCFG
    nG = h[..., iG, :].sum(-1)
    out["phi_known"] = 100 * (h[..., iR, :].sum(-1) - h[..., iAP, :].sum(-1) - h[..., iPu, :].sum(-1)) / nG
    out["phi_pub"] = 100 * (h[..., iR, :].sum(-1) - (h[..., iS, :].sum(-1) - nG)) / nG
    out["frac_PAc"] = 100 * h[..., iPc, :].sum(-1) / nG
    out["frac_unc_of_rnd"] = 100 * (h[..., iAP, :].sum(-1) + h[..., iPu, :].sum(-1)) / h[..., iR, :].sum(-1)
    return out


P0 = readouts(H.sum(0))                               # [k, rule, voi+1]
rng = np.random.default_rng(20261014); S = H.shape[0]; B = 1000; flat = H.reshape(S, -1)
keys = list(P0.keys()); boots = {k: [] for k in keys}
for c in range(0, B, 50):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    Hb = (w @ flat).reshape((50,) + H.shape[1:])
    pb = readouts(Hb)
    for k in keys: boots[k].append(pb[k])
boots = {k: np.concatenate(v) for k, v in boots.items()}   # [B, k, rule, voi+1]
sd = lambda a: 0.5 * (np.nanpercentile(a, 84, axis=0) - np.nanpercentile(a, 16, axis=0))

out = dict(doc=__doc__, ks=KS, mu=mu.tolist(), rules=RULES, tubes=[NAMES[i] for i in TI], common=COMMON, B=B, per_rung={})
for ir, rule in enumerate(RULES):
    rr = []
    for ik in range(len(KS)):
        e = dict(k=KS[ik], mu=float(mu[ik]))
        for key in keys:
            for lab, idx in (("eight", TI), ("six", CI)):
                v0 = float(np.nanmean(P0[key][ik, ir, idx])); vb = np.nanmean(boots[key][:, ik, ir, idx], axis=-1)
                e[f"{key}_{lab}"] = [v0, float(sd(vb))]
            e[f"{key}_per_tube"] = [[float(P0[key][ik, ir, i]), float(sd(boots[key][:, ik, ir, i]))] for i in TI]
        rr.append(e)
    out["per_rung"][rule] = rr
# slope of the pooled randoms-correction error against mu (per 0.01 accidental prompts per window), bootstrap
for ir, rule in enumerate(RULES):
    for key in ("dR_sharp", "dR_acc359", "phi_known", "phi_pub", "dP_acc359"):
        y0 = np.nanmean(P0[key][:, ir, TI], axis=-1); yb = np.nanmean(boots[key][:, :, ir, TI], axis=-1)
        X = np.c_[mu / 0.01, np.ones_like(mu)]
        s0 = np.linalg.lstsq(X, y0, rcond=None)[0]; sb = np.array([np.linalg.lstsq(X, y, rcond=None)[0] for y in yb])
        out.setdefault("slopes", {}).setdefault(rule, {})[key] = dict(slope_per_0p01=float(s0[0]), slope_sd=float(sd(sb[:, 0])),
                                                                    intercept=float(s0[1]), intercept_sd=float(sd(sb[:, 1])))
json.dump(out, open(D / "known_randoms_c2_summary.json", "w"), indent=1)
for rule in RULES:
    print(f"== {rule}")
    for e in out["per_rung"][rule]:
        print(f"  k {e['k']:.2f} mu {e['mu']:.4f}  dR sharp {e['dR_sharp_eight'][0]:+6.2f}+-{e['dR_sharp_eight'][1]:.2f}  "
              f"dR acc {e['dR_acc359_eight'][0]:+6.2f}+-{e['dR_acc359_eight'][1]:.2f}  phi_known {e['phi_known_eight'][0]:+6.2f}+-{e['phi_known_eight'][1]:.2f}  "
              f"phi_pub {e['phi_pub_eight'][0]:+6.2f}  PAc {e['frac_PAc_eight'][0]:.2f}%  dP acc {e['dP_acc359_eight'][0]:+5.2f}  "
              f"bG acc {e['bG_acc359_eight'][0]:+5.2f}+-{e['bG_acc359_eight'][1]:.2f}  bsub acc {e['bsub_acc359_eight'][0]:+6.2f}")
    print("  slopes:", json.dumps(out["slopes"][rule]))
