#!/usr/bin/env python3
"""Analysis of crt_sweep_c0_hists.npz: per timing resolution and rule, the per-tube signed bias of the windowed-moment
readout inverted with the sharp-edge moment (the production truncation correction) and with the a priori acceptance
sigma_e = hypot(353.2 ps, sigma_single / sqrt 2) (Sec. 2 of the paper: 359 ps at 90.9 ps). Shard bootstrap B = 1000
(seed 20261014), same resample for every resolution. Writes crt_sweep_c0_summary.json."""
import json, os
from pathlib import Path
import numpy as np
from scipy.special import erfc
D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-10-14_known_randoms"
Z = np.load(D / "crt_sweep_c0_hists.npz"); H = Z["hist"].astype(float)      # [shard, sigma, rule, voi+1, sig/rnd, bin]
SIG, CRT, RULES, NAMES, TCFG, BINS = Z["sigma_ps"], Z["crt_ps"], list(Z["rules"]), list(Z["voi_names"]), Z["tau_cfg"], Z["bins"]
TI = [i for i, n in enumerate(NAMES) if n.startswith("T")]
bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0)); tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def curve(se):
    A = np.ones_like(t) if se is None else Phi((5 - t) / se) * Phi((t + 5) / se)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


CUR = {"sharp": [curve(None)] * len(SIG), "acc": [curve(np.hypot(0.3532, s / 1000 / np.sqrt(2))) for s in SIG]}


def biases(Hs):
    h = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]
    m = (h * x).sum(-1) / h.sum(-1)                               # [sigma, rule, voi+1]
    out = {}
    for k, cs in CUR.items():
        tau = np.stack([np.interp(m[j], cs[j], tg, left=np.nan, right=np.nan) for j in range(len(SIG))])
        out[k] = 100 * (tau - TCFG) / TCFG
    return out


B0 = biases(H.sum(0))
rng = np.random.default_rng(20261014); S = H.shape[0]; flat = H.reshape(S, -1); BB = {k: [] for k in CUR}
for c in range(0, 1000, 100):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(100)]).astype(float)
    Hb = (w @ flat).reshape((100,) + H.shape[1:])
    for i in range(100):
        b = biases(Hb[i])
        for k in CUR: BB[k].append(b[k])
BB = {k: np.array(v) for k, v in BB.items()}
sd = lambda a: float(0.5 * (np.nanpercentile(a, 84) - np.nanpercentile(a, 16)))
rows = []
for j, (s, c) in enumerate(zip(SIG, CRT)):
    e = dict(sigma_single_ps=float(s), crt_ps=int(c), sigma_e_ps=float(1000 * np.hypot(0.3532, s / 1000 / np.sqrt(2))))
    for k in CUR:
        for ir, r in enumerate(RULES):
            v = B0[k][j, ir, TI]; vb = BB[k][:, j, ir][:, TI]
            e[f"{k}_{r}"] = dict(per_tube=v.round(3).tolist(), mean=float(v.mean()), mean_sd=sd(vb.mean(1)), mean_abs=float(np.abs(v).mean()),
                                 lo=float(v.min()), hi=float(v.max()))
            if j > 0:   # change against 303 ps, paired
                d0 = v - B0[k][0, ir, TI]; db = vb - BB[k][:, 0, ir][:, TI]
                e[f"{k}_{r}"]["change_vs_303_mean"] = float(d0.mean()); e[f"{k}_{r}"]["change_vs_303_sd"] = sd(db.mean(1))
        d = B0[k][j, 1, TI] - B0[k][j, 0, TI]; db = BB[k][:, j, 1][:, TI] - BB[k][:, j, 0][:, TI]
        e[f"{k}_rule_diff"] = dict(mean=float(d.mean()), mean_sd=sd(db.mean(1)), lo=float(d.min()), hi=float(d.max()))
    rows.append(e)
    print(f"CRT {c:3d} ps (sigma {s:5.1f}): closest sharp mean {e['sharp_closest']['mean']:+6.2f}+-{e['sharp_closest']['mean_sd']:.2f} "
          f"[{e['sharp_closest']['lo']:+.1f}, {e['sharp_closest']['hi']:+.1f}] change {e['sharp_closest'].get('change_vs_303_mean', 0):+.2f}+-{e['sharp_closest'].get('change_vs_303_sd', 0):.2f}; "
          f"all sharp {e['sharp_all']['mean']:+6.2f}; rule diff {e['sharp_rule_diff']['mean']:+.2f}+-{e['sharp_rule_diff']['mean_sd']:.2f}; "
          f"all acc mean|b| {e['acc_all']['mean_abs']:.2f} signed {e['acc_all']['mean']:+.2f}+-{e['acc_all']['mean_sd']:.2f}")
json.dump(dict(doc=__doc__, tubes=[NAMES[i] for i in TI], rows=rows), open(D / "crt_sweep_c0_summary.json", "w"), indent=1)
