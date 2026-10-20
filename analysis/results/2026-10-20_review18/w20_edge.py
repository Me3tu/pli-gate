#!/usr/bin/env python3
"""Review round 18 (third report on v2.17, 2026-10-20, item 3.2): true-triple acceptance edge of the 20 cm water
cylinder, by the method of Note S17 (2026-10-07_acceptance_factorization). No new simulation.

Data: the eight batches of 2026-10-19_review16/water20cm (water20_batch_200..207.h5, 128 shards), the clinical
builder of w20_extract.py (steinberger_extract_ap: 90.9 ps single-photon smear with the same seed and draw order,
adjacent pairs within 2 ns, phantom outline and LOR cut of the 20 cm cylinder, TOF annihilation point, per-event
path-corrected delay), all-prompts matching with the prompt search WIDENED to +-9 ns. Every candidate triple
records tau, s = t_p - tbar, whether its three singles share one eventID, and whether any of the three singles was
Compton-scattered in the phantom. The production selection is the subset |s| <= 5 (identical by construction to
the production all-prompts triples; checked by reproducing the w20_extract readouts). Randoms arm: search displaced
by +50 ns, built the same way.

Outputs: histograms H[shard, arm(0 signal, 1 delayed), class(0 true unscattered, 1 true scattered, 2 not true),
sel(0 |s|<=9, 1 |s|<=5), region(bg, VOI_5x15, VOI_10x20, VOI_15x25), bin] in w20_edge_hists.npz; and w20_edge.json with
  - edge of the true triples (all, unscattered, scattered), binomial ML fit of n(|s|<=5)/n(|s|<=9) per 50 ps bin of
    tau between 3 and 7 ns with Phi((5 + d - tau)/s), in VOI_15x25 and in all regions pooled
  - windowed mean of the selected true triples on (1, 5) ns against the mean predicted from all true candidates
    weighted by an acceptance model: a priori (350.1 ps, d = 0), randoms arm (w20_extract: 372.6 ps, +93.2 ps),
    and the fitted true-triple edge
  - lifetime of the production randoms-subtracted all-prompts spectrum of VOI_15x25 with the sharp, a priori,
    randoms-arm and fitted true-triple acceptances, against event identity (1.836 +- 0.008 ns) and the configured
    1.84 ns
Uncertainties: shard bootstrap B = 300 (seed 18), edge refits started at the full-sample estimate."""
import json, os, sys, time
from pathlib import Path
import numpy as np, h5py
from scipy.optimize import minimize
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
W20 = RES / "2026-10-19_review16/water20cm"
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa: E402

OUTLINE = {"type": "cylinder", "R_mm": 102.0, "halfz_mm": 160.0, "buffer_mm": 2.0, "lor_tol_mm": 10.0}
VOIS = [dict(name=n, shape="cylinder", cx_mm=0.0, cy_mm=0.0, cz_mm=0.0, r_mm=r, half_z_mm=hz)
        for n, r, hz in (("VOI_5x15", 25.0, 75.0), ("VOI_10x20", 50.0, 100.0), ("VOI_15x25", 75.0, 125.0))]
RNAMES = ["bg"] + [v["name"] for v in VOIS]
SIG_T_NS, W_WIDE = 0.0909, 9.0
EDGES = np.arange(-6.0, 9.0001, 0.05); NB = len(EDGES) - 1
bc = 0.5 * (EDGES[:-1] + EDGES[1:])
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))


def shard_hists(t_true, in5, inp, gx, gy, gz, eid, comp, rng):
    H = np.zeros((2, 3, 2, len(RNAMES), NB), np.int64)
    n = len(t_true)
    if n == 0: return H
    t = t_true + rng.normal(0.0, SIG_T_NS, size=n)          # identical draw to clinical_per_shard
    i5 = np.flatnonzero(in5.astype(bool)); ip = np.flatnonzero(inp.astype(bool))
    if len(i5) < 2 or len(ip) < 1: return H
    o5 = np.argsort(t[i5]); op = np.argsort(t[ip])            # unstable sort, as clinical_per_shard
    i5 = i5[o5]; ip = ip[op]; t5 = t[i5]; tp = t[ip]
    ok = np.diff(t5) <= X.W_PAIR_NS
    if not ok.any(): return H
    pa = i5[:-1][ok]; pb = i5[1:][ok]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    pann = X.tof_p_ann(p1, p2, t[pa], t[pb])
    geom = (X.lor_perp(pann, p1, p2) <= OUTLINE["lor_tol_mm"]) & X.in_phantom_outline(pann[:, 0], pann[:, 1], pann[:, 2], OUTLINE)
    pa, pb, pann, p1, p2 = pa[geom], pb[geom], pann[geom], p1[geom], p2[geom]
    if len(pa) == 0: return H
    tbar = 0.5 * (t[pa] + t[pb])
    for arm, shift in ((0, 0.0), (1, X.DELAY_NS)):
        sel, pick = X.all_prompt_indices(tbar + shift, tp, W_WIDE)
        if len(sel) == 0: continue
        q = ip[pick]
        pp = np.stack([gx[q], gy[q], gz[q]], -1)
        tau = X.per_event_tau(p1[sel], p2[sel], pp, t[pa[sel]], t[pb[sel]], t[q], pann[sel]) + shift
        s = t[q] - (tbar[sel] + shift)
        inside = np.abs(s) <= X.W_PROMPT_NS
        true = (eid[pa[sel]] == eid[pb[sel]]) & (eid[pa[sel]] == eid[q])
        scat = (comp[pa[sel]] > 0) | (comp[pb[sel]] > 0) | (comp[q] > 0)
        cls = np.where(true & ~scat, 0, np.where(true, 1, 2))
        b = np.searchsorted(EDGES, tau, side="right") - 1
        good = (b >= 0) & (b < NB)
        masks = [X.voi_bg_mask(pann[sel], VOIS)] + X.voi_masks(pann[sel], VOIS)
        for r, mr in enumerate(masks):
            for c in range(3):
                m = mr & good & (cls == c)
                H[arm, c, 0, r] += np.bincount(b[m], minlength=NB)
                H[arm, c, 1, r] += np.bincount(b[m & inside], minlength=NB)
    return H


def build():
    files = sorted((Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-10-19_review16/water20cm").glob("water20_batch_2*.h5"))
    if not files: raise SystemExit("water20_batch_2*.h5: raw data not included in the release, see analysis/README.md")
    rng = np.random.default_rng(int(90.9 * 1000) or 1)      # as w20_extract.run() for each rule
    HS = []; t0 = time.time()
    for fn in files:
        f = h5py.File(fn, "r"); g = f["singles"]
        for sh_id, s, e in X.find_shard_boundaries(g["shard_id"][:]):
            HS.append(shard_hists(g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e], g["gx"][s:e], g["gy"][s:e],
                                  g["gz"][s:e], g["eventID"][s:e], g["compton_phantom"][s:e], rng))
        f.close()
        print(f"{fn.name} done, {len(HS)} shards, {time.time() - t0:.0f} s", flush=True)
    HS = np.array(HS)
    np.savez_compressed(OUT / "w20_edge_hists.npz", H=HS, edges=EDGES, regions=RNAMES, files=[str(f) for f in files])
    return HS


def fit_edge(nsel, nall, x0=None):
    m = (bc > 3.0) & (bc < 7.0) & (nall > 0)
    k, n, t = nsel[m], nall[m], bc[m]
    def nll(p):
        s, d = p
        if s <= 0.05: return 1e30
        q = np.clip(Phi((5 + d - t) / s), 1e-12, 1 - 1e-12)
        return -float((k * np.log(q) + (n - k) * np.log(1 - q)).sum())
    starts = [x0] if x0 is not None else [[s0, d0] for s0 in (0.3, 0.36, 0.42) for d0 in (-0.05, 0.0, 0.05, 0.1)]
    best = min((minimize(nll, st, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6)) for st in starts), key=lambda r: r.fun)
    return best.x


fw = (bc > 1.0) & (bc < 5.0); x = bc[fw]
tg = np.linspace(0.3, 12.0, 6000); tt = np.linspace(1.0, 5.0, 8001)


from functools import lru_cache


@lru_cache(maxsize=64)
def curve(s=None, d=0.0):
    A = np.ones_like(tt) if s is None else Phi((5 + d - tt) / s) * Phi((tt + 5 - d) / s)
    E = np.exp(-tt[None, :] / tg[:, None]) * A[None, :]
    return (E * tt).sum(1) / E.sum(1)


def analyse(H, iv, models, x0s=None):
    """H summed over shards. iv: region index (or list). Returns edge fits, predicted-minus-selected means, readouts."""
    h = H[..., iv, :] if np.ndim(iv) == 0 else H[..., iv, :].sum(-2)
    out = {}
    for nm, cl in (("true", [0, 1]), ("true_unscattered", [0]), ("true_scattered", [1])):
        nall = h[0, cl, 0].sum(0); nsel = h[0, cl, 1].sum(0)
        s, d = fit_edge(nsel, nall, None if x0s is None else x0s[f"edge_{nm}"])
        out[f"edge_{nm}"] = (s, d)
        out[f"frac_{nm}_4p5_5p5"] = float(nall[(bc > 4.5) & (bc < 5.5)].sum())
    nall = h[0, [0, 1], 0].sum(0)[fw]; nsel = h[0, [0, 1], 1].sum(0)[fw]
    sel_mean = (nsel * x).sum() / nsel.sum()
    for k, (s, d) in models.items():
        A = Phi((5 + d - x) / s) * Phi((x + 5 - d) / s)
        out[f"pred_minus_sel_ps_{k}"] = 1000 * ((nall * A * x).sum() / (nall * A).sum() - sel_mean)
    s_t, d_t = out["edge_true"]
    out["pred_minus_sel_ps_fitted"] = 1000 * ((nall * Phi((5 + d_t - x) / s_t) * Phi((x + 5 - d_t) / s_t) * x).sum()
                                             / (nall * Phi((5 + d_t - x) / s_t) * Phi((x + 5 - d_t) / s_t)).sum() - sel_mean)
    out["frac_true_1to5_inside"] = float(nsel.sum() / nall.sum())
    # readouts of the production all-prompts randoms-subtracted spectrum (|s| <= 5, all classes)
    sig = h[0, :, 1].sum(0)[fw]; rnd = h[1, :, 1].sum(0)[fw]; sub = sig - rnd
    mean = (sub * x).sum() / sub.sum()
    cv = {"sharp": curve(), **{k: curve(s, d) for k, (s, d) in models.items()}, "fitted_true": curve(s_t, d_t)}
    for k, c in cv.items():
        out[f"tau_{k}"] = float(np.interp(mean, c, tg))
    tr = h[0, [0, 1], 1].sum(0)[fw]
    out["tau_true_selected_fitted"] = float(np.interp((tr * x).sum() / tr.sum(), cv["fitted_true"], tg))
    out["tau_true_selected_apriori"] = float(np.interp((tr * x).sum() / tr.sum(), cv["apriori"], tg))
    return out


def main():
    p = OUT / "w20_edge_hists.npz"
    HS = np.load(p)["H"] if p.exists() else build()
    print("histograms", HS.shape, flush=True)
    j = json.loads((W20 / "w20_extract.json").read_text())["all"]
    models = {"apriori": (0.35006, 0.0), "randoms_arm": (j["sigma_e_ps"] / 1000, j["edge_shift_ps"] / 1000)}
    regions = {"VOI_15x25": 3, "all_regions": [0, 1, 2, 3]}
    res = {"doc": __doc__, "models_ns": models, "event_identity_ns": [1.836, 0.008], "tau_cfg_ns": 1.84, "regions": {}}
    S = HS.shape[0]; rng = np.random.default_rng(18); B = 300
    for rn, iv in regions.items():
        o0 = analyse(HS.sum(0), iv, models)
        boots = []
        for _ in range(B):
            w = np.bincount(rng.integers(0, S, S), minlength=S).astype(float)
            boots.append(analyse(np.tensordot(w, HS, axes=(0, 0)), iv, models, x0s=o0))
        row = {}
        for k, v in o0.items():
            if isinstance(v, tuple):
                bs = np.array([b[k] for b in boots])
                row[k] = dict(s_ps=1000 * v[0], s_ps_sd=1000 * float(bs[:, 0].std(ddof=1)), d_ps=1000 * v[1], d_ps_sd=1000 * float(bs[:, 1].std(ddof=1)))
            else:
                bs = np.array([b[k] for b in boots])
                row[k] = [float(v), float(np.nanstd(bs, ddof=1))]
        res["regions"][rn] = row
        print(f"== {rn}")
        for k, v in row.items(): print(f"   {k}: {v}")
    # reproduction check of the production readouts (w20_extract, VOI_15x25, all prompts)
    res["check_w20_extract"] = {k: j["vois"]["VOI_15x25"][k] for k in ("sharp", "acc_measured", "acc359")}
    json.dump(res, open(OUT / "w20_edge.json", "w"), indent=1, default=float)
    print("wrote", OUT / "w20_edge.json")


if __name__ == "__main__":
    main()
