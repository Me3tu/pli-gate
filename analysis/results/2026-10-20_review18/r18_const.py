#!/usr/bin/env python3
"""Review round 18 (third report on v2.17, 2026-10-20): the sideband constant, three follow-ups. No new simulation.

A  Rate dependence of the constant's randoms-correction error. From rules16_hists.npz (eight tubes, rate replay
   k = 0.5, 1, 2), the error dR = tau(sig - c*A) - tau(G + PA_corr) of the acceptance-aware readout (as in
   rules16_const.py rev. b, sideband -4.5 to -2.0 ns) at all three rungs, for every arm, with the joint shard
   bootstrap (B = 1000, seed 20261019, one resample shared by all rungs). Linearity: D21 = dR(2) - 2 dR(1),
   D10 = dR(1) - 2 dR(0.5), and the coefficient b of dR = a mu + b mu^2 (fit through the origin).
B  Mirrored accidental edge (optional item of the report). The lower edge of the +-5 ns prompt search lies in the
   signal histogram at -5 ns, where only accidentals contribute. Fit c*Phi((t + 5 + d)/s) to the signal histogram
   between -5.0 and -3.0 ns (all tubes pooled, per rung and arm, Poisson ML), mirror the fitted edge to +5 ns,
   A_acc(t) = Phi((5 + d - t)/s) Phi((t + 5 + d)/s), subtract c*A_acc with c from the sideband plateau under A_acc,
   and read out with the a priori acceptance of the true triples. Reports dR and the residual b against tau_cfg.
C  Published windows with the constant. From 2026-10-09_pair_window/c0/hists_pair_window.npz (Note S18; arm
   all4p2_t20 = all pairs within 4.2 ns, prompts within +-20 ns, and arm adj2 = production), the sharp readout
   (truncation moment on 1 to 5 ns) after subtracting either the delayed window or a flat constant from the
   sideband -4.5 to -2.0 ns (check -5.0 to -2.0), per tube, both rules. With prompts within +-20 ns the sideband and
   the fit window lie inside the flat part of the accidental content, so the constant is subtracted flat. Paired
   shard bootstrap (B = 400, seed 18) of b_const - b_delayed. For adj2 the constant is subtracted with the a priori
   acceptance and read out with it, as in Note S12.
D  Departure from flatness by delay band (1-3, 3-4.5, 4.5-5 ns) for every arm at k = 1 and 2 (as rules16_flat_bands).
Writes r18_const.json."""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc
from scipy.optimize import minimize

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = Path(__file__).resolve().parent
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
sd = lambda a, ax=0: 0.5 * (np.nanpercentile(a, 84, axis=ax) - np.nanpercentile(a, 16, axis=ax))
tg = np.linspace(0.3, 12.0, 6000); tt = np.linspace(1.0, 5.0, 8001)
SIG_E = 0.35897


def curve_from_A(Afun):
    Aa = np.ones_like(tt) if Afun is None else Afun(tt)
    E = np.exp(-tt[None, :] / tg[:, None]) * Aa[None, :]
    return (E * tt).sum(1) / E.sum(1)


A359 = lambda t: Phi((5 - t) / SIG_E) * Phi((t + 5) / SIG_E)
CURVES = {"sharp": curve_from_A(None), "acc359": curve_from_A(A359)}
result = {"doc": __doc__}

# ---------------------------------------------------------------- A, B, D: rules16 histograms
Z = np.load(RES / "2026-10-19_review16/rules16_hists.npz")
H = Z["hist"].astype(np.float64)                      # [shard, k, arm, voi+1, cls, bin]
KS, ARMS, CLS = [float(k) for k in Z["ks"]], [str(a) for a in Z["arms"]], [str(c) for c in Z["classes"]]
NAMES, TCFG, BINS = [str(n) for n in Z["voi_names"]], Z["tau_cfg"].astype(float), Z["bins"]
MU = json.loads((RES / "2026-10-08_c0_ladder/mu_ladder_c0.json").read_text())
mu = np.array([MU[f"{k:.2f}"]["mu"] for k in KS])
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
bc = 0.5 * (BINS[:-1] + BINS[1:]); fw = (bc > 1.0) & (bc < 5.0); x = bc[fw]
iS, iG, iAP, iPc, iPu = (CLS.index(c) for c in ("sig", "G", "AP", "PA_corr", "PA_unc"))
sbm = (bc > -4.5) & (bc < -2.0)
A_all = A359(bc)


def tau(spec, curve):
    n = spec.sum(-1); m = (spec * x).sum(-1) / np.where(n > 0, n, np.nan)
    return np.interp(m, curve, tg, left=np.nan, right=np.nan)


def dR_const(Hs):
    sig = Hs[..., iS, :]
    c = sig[..., sbm].sum(-1) / A_all[sbm].sum()
    h = Hs[..., fw]
    sub = h[..., iS, :] - c[..., None] * A_all[fw]
    ref = h[..., iG, :] + h[..., iPc, :]
    return 100 * (tau(sub, CURVES["acc359"]) - tau(ref, CURVES["acc359"])) / TCFG


# edge fit for B on the full-sample histograms (per rung and arm, tubes pooled)
lo = (bc > -5.0) & (bc < -3.0)
edge = {}
Hsum = H.sum(0)
for ik in range(len(KS)):
    for ia, arm in enumerate(ARMS):
        y = Hsum[ik, ia, TI, iS, :].sum(0)[lo]; t = bc[lo]
        def nll(p):
            c, s, d = p
            if s <= 0.05 or c <= 0: return 1e30
            m = c * Phi((t + 5 + d) / s)
            return float((m - y * np.log(np.maximum(m, 1e-300))).sum())
        best = None
        for s0 in (0.30, 0.36, 0.42):
            for d0 in (-0.05, 0.0, 0.05, 0.10):
                r = minimize(nll, [y[-10:].mean(), s0, d0], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
                if best is None or r.fun < best.fun: best = r
        edge[(ik, ia)] = best.x
result["mirror_edge_fit"] = {f"k={KS[ik]}|{ARMS[ia]}": dict(s_ps=1000 * v[1], d_ps=1000 * v[2]) for (ik, ia), v in edge.items()}
print("mirror edge fits:", {k: (round(v['s_ps'], 1), round(v['d_ps'], 1)) for k, v in result["mirror_edge_fit"].items() if k.startswith("k=1.0")}, flush=True)


def dR_mirror(Hs):
    res_dR = np.zeros(Hs.shape[:-5] + (len(KS), len(ARMS), len(NAMES))); res_b = np.zeros_like(res_dR)
    for ik in range(len(KS)):
        for ia in range(len(ARMS)):
            _, s, d = edge[(ik, ia)]
            Aacc = Phi((5 + d - bc) / s) * Phi((bc + 5 + d) / s)
            hh = Hs[..., ik, ia, :, :, :]
            sig = hh[..., iS, :]
            c = sig[..., sbm].sum(-1) / Aacc[sbm].sum()
            h = hh[..., fw]
            sub = h[..., iS, :] - c[..., None] * Aacc[fw]
            ref = h[..., iG, :] + h[..., iPc, :]
            ts = tau(sub, CURVES["acc359"])
            res_dR[..., ik, ia, :] = 100 * (ts - tau(ref, CURVES["acc359"])) / TCFG
            res_b[..., ik, ia, :] = 100 * (ts - TCFG) / TCFG
    return res_dR, res_b


def bands(Hs):
    acc = Hs[..., iAP, :] + Hs[..., iPu, :]
    a = acc[..., TI, :].sum(-2); ref = a[..., sbm].sum(-1) / A_all[sbm].sum()
    return np.stack([a[..., (bc > l) & (bc < u)].sum(-1) / A_all[(bc > l) & (bc < u)].sum() / ref - 1
                     for l, u in ((1, 3), (3, 4.5), (4.5, 5))], -1)


S = H.shape[0]; flat = H.reshape(S, -1)
P0 = dR_const(H.sum(0)); M0 = dR_mirror(H.sum(0)); B0 = bands(H.sum(0))
rng = np.random.default_rng(20261019); PB, MB, BB = [], [], []
for ib in range(20):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(50)]).astype(float)
    Hb = (w @ flat).reshape((50,) + H.shape[1:])
    PB.append(dR_const(Hb)); m = dR_mirror(Hb); MB.append(np.stack(m, 1)); BB.append(bands(Hb))
    print("boot batch", ib, flush=True)
PB = np.concatenate(PB); MB = np.concatenate(MB); BB = np.concatenate(BB)
mean_t = lambda a: np.nanmean(a[..., TI], axis=-1)
A_out = {}
for ia, arm in enumerate(ARMS):
    y0 = mean_t(P0[:, ia]); yb = mean_t(PB[:, :, ia])                    # [k], [B, k]
    D21 = y0[2] - 2 * y0[1]; D21b = yb[:, 2] - 2 * yb[:, 1]
    D10 = y0[1] - 2 * y0[0]; D10b = yb[:, 1] - 2 * yb[:, 0]
    X = np.c_[mu, mu ** 2]
    c0 = np.linalg.lstsq(X, y0, rcond=None)[0]; cb = np.array([np.linalg.lstsq(X, v, rcond=None)[0] for v in yb])
    m0 = mean_t(M0[0][:, ia]); mb = mean_t(MB[:, 0, :, ia]); mbb0 = mean_t(M0[1][:, ia]); mbb = mean_t(MB[:, 1, :, ia])
    A_out[arm] = dict(dR=[[float(y0[i]), float(sd(yb[:, i]))] for i in range(3)],
                      D21=[float(D21), float(sd(D21b))], D10=[float(D10), float(sd(D10b))],
                      quad_a_per_0p01=[float(c0[0] / 100), float(sd(cb[:, 0]) / 100)],
                      quad_b_per_0p01sq=[float(c0[1] / 1e4), float(sd(cb[:, 1]) / 1e4)],
                      mirror_dR=[[float(m0[i]), float(sd(mb[:, i]))] for i in range(3)],
                      mirror_b=[[float(mbb0[i]), float(sd(mbb[:, i]))] for i in range(3)],
                      bands_k1=[[float(B0[1, ia, j] * 100), float(sd(BB[:, 1, ia, j]) * 100)] for j in range(3)],
                      bands_k2=[[float(B0[2, ia, j] * 100), float(sd(BB[:, 2, ia, j]) * 100)] for j in range(3)])
    e = A_out[arm]
    print(f"{arm:11s} dR " + " ".join(f"{v:+.2f}+-{s:.2f}" for v, s in e['dR']) +
          f" | D21 {e['D21'][0]:+.2f}+-{e['D21'][1]:.2f} D10 {e['D10'][0]:+.2f}+-{e['D10'][1]:.2f}"
          f" | mirror dR " + " ".join(f"{v:+.2f}+-{s:.2f}" for v, s in e['mirror_dR']) +
          f" | bands k1 " + " ".join(f"{v:+.1f}+-{s:.1f}" for v, s in e['bands_k1']), flush=True)
result["A_rate_and_B_mirror_and_D_bands"] = dict(ks=KS, mu=mu.tolist(), arms=A_out)
json.dump(result, open(OUT / "r18_const.json", "w"), indent=1)

# ---------------------------------------------------------------- C: published windows with the constant
PW = RES / "2026-10-09_pair_window/c0"
Zp = np.load(PW / "hists_pair_window.npz"); rows = json.loads((PW / "pair_window.json").read_text())["rows"]
tc = np.array([r["tau_cfg"] for r in rows]); bins = Zp["bins"]; bcp = 0.5 * (bins[:-1] + bins[1:])
assert np.allclose(bins, BINS)
C_out = {}
for arm in ("all4p2_t20", "adj2"):
    for rule in ("all", "closest"):
        HS_, HR_ = Zp[f"{arm}_{rule}_sig"][:, 1:, :], Zp[f"{arm}_{rule}_rnd"][:, 1:, :]   # tubes only [shard, 8, bin]
        def rd(hs, hr):
            o = {}
            o["b_delayed"] = 100 * (tau(hs[..., fw] - hr[..., fw], CURVES["sharp"]) - tc) / tc
            for nm, (l, u) in (("sb45", (-4.5, -2.0)), ("sb50", (-5.0, -2.0))):
                m = (bcp > l) & (bcp < u)
                if arm == "adj2":
                    c = hs[..., m].sum(-1) / A_all[m].sum()
                    o[f"b_const_{nm}"] = 100 * (tau(hs[..., fw] - c[..., None] * A_all[fw], CURVES["acc359"]) - tc) / tc
                else:
                    c = hs[..., m].mean(-1)
                    o[f"b_const_{nm}"] = 100 * (tau(hs[..., fw] - c[..., None], CURVES["sharp"]) - tc) / tc
                o[f"diff_{nm}"] = o[f"b_const_{nm}"] - o["b_delayed"]
            return o
        o0 = rd(HS_.sum(0), HR_.sum(0))
        r2 = np.random.default_rng(18); S2 = HS_.shape[0]; ob = {k: [] for k in o0}
        for _ in range(8):
            w = np.stack([np.bincount(r2.integers(0, S2, S2), minlength=S2) for _ in range(50)]).astype(float)
            o = rd(np.tensordot(w, HS_, axes=(1, 0)), np.tensordot(w, HR_, axes=(1, 0)))
            for k in o: ob[k].append(o[k])
        ob = {k: np.concatenate(v) for k, v in ob.items()}
        e = {}
        for k in o0:
            e[k] = dict(per_tube=o0[k].tolist(), mean=[float(o0[k].mean()), float(sd(ob[k].mean(-1)))],
                        mean_abs=[float(np.abs(o0[k]).mean()), float(sd(np.abs(ob[k]).mean(-1)))])
        C_out[f"{arm}|{rule}"] = e
        print(f"{arm:11s} {rule:8s} b_delayed mean {e['b_delayed']['mean'][0]:+.2f} |b| {e['b_delayed']['mean_abs'][0]:.2f} | "
              f"const45 mean {e['b_const_sb45']['mean'][0]:+.2f}+-{e['b_const_sb45']['mean'][1]:.2f} |b| {e['b_const_sb45']['mean_abs'][0]:.2f} | "
              f"diff45 {e['diff_sb45']['mean'][0]:+.2f}+-{e['diff_sb45']['mean'][1]:.2f} diff50 {e['diff_sb50']['mean'][0]:+.2f}+-{e['diff_sb50']['mean'][1]:.2f}", flush=True)
r = {}
for arm in ("all4p2_t20", "adj2"):
    for nm in ("b_delayed", "b_const_sb45"):
        a = np.array(C_out[f"{arm}|all"][nm]["per_tube"]); c = np.array(C_out[f"{arm}|closest"][nm]["per_tube"])
        r[f"{arm}|{nm}"] = float((a - c).mean())
C_out["selection_shift_mean_pp"] = r
print("selection shift (b_all - b_closest), mean over tubes:", {k: round(v, 2) for k, v in r.items()})
result["C_published_windows"] = C_out
json.dump(result, open(OUT / "r18_const.json", "w"), indent=1)
print("wrote", OUT / "r18_const.json")
