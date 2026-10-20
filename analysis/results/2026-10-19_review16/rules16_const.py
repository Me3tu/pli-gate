#!/usr/bin/env python3
"""Review round 16b (2026-10-19, second report, Sec. 3): the six prompt rules of rules16 with the randoms estimated
as a constant fitted on negative delays (the Quadra procedure) instead of the delayed window.

Input: rules16_hists.npz (same events, same arms; signal histograms only, class 'sig').
Per rung (k = 0.5, 1.0, 2.0), arm and tube:
  c          mean signal counts per 50 ps bin on the sideband SB of path-corrected delay (primary [-4.5, -2.0] ns,
             check [-5.0, -2.0] ns), i.e. the constant background of the published analyses
  sharp      readout of sig - c on the fit window 1 < t < 5 ns with the sharp-edge inversion
  acc359     readout of sig - c*A(t) with the a priori prompt-window acceptance A(t), sigma_e = 358.97 ps, inverted
             with the same acceptance (the constant follows the acceptance like the accidental content)
  bsub_*     residual against the configured lifetime, % of tau_cfg
  dR_*       randoms-correction error tau(sig - background) - tau(G + PA_corr), % of tau_cfg
  phi_const  (estimated background in the fit window - AP - PA_unc) / G, %: departure of the constant from the
             uncorrelated accidental content of the fit window
  rule_shift b(all) - b(arm) on identical events (all-pairs arms against ap_all), percentage points
Pooled values: unweighted means over the eight tubes. Uncertainties: shard bootstrap B = 1000, seed 20261019 (the
same resamples as rules16_analysis.py). Writes rules16_const_summary.json and prints a table.
Revision b: the constant is the sideband plateau sum(sig)/sum(A) instead of the plain sideband mean."""
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
bc = 0.5 * (BINS[:-1] + BINS[1:])
sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
SIG_E = 0.35897
A_all = Phi((5 - bc) / SIG_E) * Phi((bc + 5) / SIG_E)
A = A_all[sel]
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def moment_curve(sig_e):
    Aa = np.ones_like(t) if sig_e is None else Phi((5 - t) / sig_e) * Phi((t + 5) / sig_e)
    E = np.exp(-t[None, :] / tg[:, None]) * Aa[None, :]
    return (E * t).sum(1) / E.sum(1)


CURVES = {"sharp": moment_curve(None), "acc359": moment_curve(SIG_E)}
iS, iG, iAP, iPc, iPu = (CLS.index(c) for c in ("sig", "G", "AP", "PA_corr", "PA_unc"))
iAll, iApAll = ARMS.index("all"), ARMS.index("ap_all")
SBS = {"sb45": (-4.5, -2.0), "sb50": (-5.0, -2.0)}


def tau(spec, curve):
    n = spec.sum(-1); m = (spec * x).sum(-1) / np.where(n > 0, n, np.nan)
    return np.interp(m, CURVES[curve], tg, left=np.nan, right=np.nan)


def readouts(Hs, sb):
    lo, hi = SBS[sb]
    sbm = (bc > lo) & (bc < hi)
    sig_all = Hs[..., iS, :]
    # plateau of the constant: sideband counts divided by the summed acceptance there (Poisson ML for c*A(t)), so
    # the softened window edge at -5 ns does not bias it. The sharp readout subtracts the same plateau flat.
    c = sig_all[..., sbm].sum(-1) / A_all[sbm].sum()
    h = Hs[..., sel]
    sig = h[..., iS, :]
    gp = h[..., iG, :] + h[..., iPc, :]; g = h[..., iG, :]
    out = {}
    for cv in CURVES:
        shape = np.ones_like(x) if cv == "sharp" else A
        sub = sig - c[..., None] * shape
        ts, tgp, tgg = tau(sub, cv), tau(gp, cv), tau(g, cv)
        out[f"bsub_{cv}"] = 100 * (ts - TCFG) / TCFG
        out[f"dR_{cv}"] = 100 * (ts - tgp) / TCFG
        out[f"bG_{cv}"] = 100 * (tgg - TCFG) / TCFG
    nG = h[..., iG, :].sum(-1)
    est = c * A.sum()
    out["phi_const"] = 100 * (est - h[..., iAP, :].sum(-1) - h[..., iPu, :].sum(-1)) / nG
    for cv in CURVES:
        b = out[f"bsub_{cv}"]
        ref = np.stack([b[..., iAll if not ARMS[a].startswith("ap_") else iApAll, :] for a in range(len(ARMS))], axis=-2)
        out[f"rule_shift_{cv}"] = ref - b
    return out


sd = lambda a: 0.5 * (np.nanpercentile(a, 84, axis=0) - np.nanpercentile(a, 16, axis=0))
result = dict(doc=__doc__, ks=KS, mu=mu.tolist(), arms=ARMS, sidebands={k: list(v) for k, v in SBS.items()}, B=1000)
S = H.shape[0]; flat = H.reshape(S, -1)
for sb in SBS:
    P0 = readouts(H.sum(0), sb)
    rng = np.random.default_rng(20261019); B = 1000
    keys = list(P0.keys()); boots = {k: [] for k in keys}
    for c0 in range(0, B, 50):
        w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
        pb = readouts((w @ flat).reshape((50,) + H.shape[1:]), sb)
        for k in keys: boots[k].append(pb[k])
    boots = {k: np.concatenate(v) for k, v in boots.items()}
    per = []
    for ia, arm in enumerate(ARMS):
        for ik, k in enumerate(KS):
            e = dict(arm=arm, k=k, mu=float(mu[ik]))
            for key in keys:
                v0 = float(np.nanmean(P0[key][ik, ia, TI])); vb = np.nanmean(boots[key][:, ik, ia, TI], axis=-1)
                e[key] = [v0, float(sd(vb))]
                if key.startswith("bsub"):
                    a0 = float(np.nanmean(np.abs(P0[key][ik, ia, TI]))); ab = np.nanmean(np.abs(boots[key][:, ik, ia, TI]), axis=-1)
                    e["mae_" + key] = [a0, float(sd(ab))]
                if key.startswith("rule_shift") or key.startswith("bsub"):
                    e[key + "_range"] = [float(np.nanmin(P0[key][ik, ia, TI])), float(np.nanmax(P0[key][ik, ia, TI]))]
            per.append(e)
    slopes = {}
    for ia, arm in enumerate(ARMS):
        for key in ("dR_acc359", "phi_const", "rule_shift_acc359", "rule_shift_sharp"):
            y0 = np.nanmean(P0[key][:, ia, TI], axis=-1); yb = np.nanmean(boots[key][:, :, ia, TI], axis=-1)
            X = np.c_[mu / 0.01, np.ones_like(mu)]
            s0 = np.linalg.lstsq(X, y0, rcond=None)[0]; sb_ = np.array([np.linalg.lstsq(X, y, rcond=None)[0] for y in yb])
            slopes.setdefault(arm, {})[key] = [float(s0[0]), float(sd(sb_[:, 0]))]
    result[sb] = dict(per=per, slopes=slopes)
    print(f"== sideband {SBS[sb]} ns")
    for e in per:
        if e["k"] != 1.0: continue
        print(f"{e['arm']:11s} mae acc {e['mae_bsub_acc359'][0]:5.2f} | signed acc {e['bsub_acc359'][0]:+6.2f}+-{e['bsub_acc359'][1]:.2f}"
              f" | shift acc {e['rule_shift_acc359'][0]:+6.2f}+-{e['rule_shift_acc359'][1]:.2f} sharp {e['rule_shift_sharp'][0]:+6.2f}+-{e['rule_shift_sharp'][1]:.2f}"
              f" | dR acc {e['dR_acc359'][0]:+6.2f}+-{e['dR_acc359'][1]:.2f} | phi {e['phi_const'][0]:+6.2f}+-{e['phi_const'][1]:.2f}")
    print(json.dumps(slopes))
json.dump(result, open(D / "rules16_const_summary.json", "w"), indent=1)
