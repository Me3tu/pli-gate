#!/usr/bin/env python3
"""Review round 12 (2026-10-14): prompt-photon path correction, reconstructed annihilation point against true decay point.

Question (referee, v2.11): the delay subtracts the prompt flight time |p_prompt - p_ann|/c with p_ann the TOF point of
the annihilation pair, whereas the prompt gamma is emitted at the decay position.  How large is the error, and does it
change the recovered lifetimes, depending on configured lifetime or migration?

Variants, on the SAME accepted triples (selection, regions, and annihilation-time term unchanged):
  prod : tau = (t1 - |p1 - p_ann|/c) - (t_p - |p_p - p_ann|/c)            production (correct_travel_time / sim1_extract)
  src  : tau = (t1 - |p1 - p_ann|/c) - (t_p - |p_p - s_p|/c)              s_p = true emission point of the prompt photon
                                                                           (GATE source position of the prompt single)
  delta = tau_src - tau_prod = (|p_p - s_p| - |p_p - p_ann|)/c.

Part A, chain closure (12 tubes, event-identity triples, files of the archived closure, R_cut 75 mm, tau_min 1 ns,
         censored-exponential mean as analyze_per_tube_tau); per-tube paired bootstrap of the difference; closure trend
         refitted with the residual-scatter method of make_fig2_v3.py, but with the statistical sigmas of the
         analysis CSV: it gives -12.2+-7.8 ps/ns for the production delays, not the published -11.7+-6.0, so only the
         difference between the two variants (d_slope) is used. [Note added 2026-10-14 after the run; code unchanged.]
Part B, clinical chain on the corrected-yield zero-delay eight-tube stream C0 (production builder of sim1_extract,
         identical smear seeds as photon_conventions_c0.py, both matching rules, randoms arm +50 ns with the prompt
         single's own source position); sharp and a priori acceptance (359 ps) moment readouts; paired shard bootstrap.
         The prod histograms are checked bin for bin against photon_conventions_c0/photon_conventions.npz.
Writes prompt_path.json and prompt_path_c0.npz.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np, h5py
from scipy.special import erfc

REPO = Path(os.environ.get("PLI_REPO") or Path(__file__).resolve().parents[3])
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = RES / '2026-10-14_prompt_path'
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(RES / '2026-08-17_derivation_test'))
sys.path.insert(0, str(RES / '2026-04-29_2320_phase2b_validation_sims/step5_extraction'))
from reconstruction.correct_travel_time import correct as travel_time_correct  # noqa: E402
import analyze_per_tube_tau as AP  # noqa: E402
import sim1_extract as SX  # noqa: E402

C = 299.792458  # mm/ns
TAU_MIN = 1.0
DBINS = np.linspace(-0.5, 0.5, 1001)       # delta, ns (1 ps bins)
RBINS = np.linspace(0.0, 200.0, 401)       # |s - p_ann|, mm
ABINS = np.linspace(-150.0, 150.0, 601)    # along-LOR component, mm
PBINS = np.linspace(0.0, 30.0, 301)        # perpendicular component, mm


def geom_parts(s, p_ann, p1, p2):
    u = p2 - p1; u = u / np.maximum(np.linalg.norm(u, axis=-1, keepdims=True), 1e-9)
    v = s - p_ann
    along = (v * u).sum(-1)
    perp = np.linalg.norm(v - along[:, None] * u, axis=-1)
    return np.linalg.norm(v, axis=-1), along, perp


def stats_from_hist(h, edges):
    c = 0.5 * (edges[1:] + edges[:-1]); n = h.sum()
    if n == 0: return {}
    m = (h * c).sum() / n; sd = np.sqrt((h * (c - m) ** 2).sum() / n)
    cdf = np.cumsum(h) / n
    q = {f'p{int(p*100)}': float(np.interp(p, cdf, c)) for p in (0.01, 0.05, 0.5, 0.95, 0.99)}
    return dict(n=int(n), mean=float(m), sd=float(sd), **q)


# ============================================================ Part A: closure
def closure(sim, cfg, csv):
    base = Path(os.environ.get("PLI_RAW", "<missing>")) / '2026-08-13_overnight/O4_rerun' / sim
    det = base / f'triples_detector_{sim}_rcut75.h5'; tru = base / f'triples_truth_{sim}_rcut75.h5'
    tau, P = travel_time_correct(det)
    with h5py.File(det, 'r') as f:
        g = f['detector']
        pp = np.stack([g['p_x'][:], g['p_y'][:], g['p_z'][:]], -1)
        p1 = np.stack([g['a1_x'][:], g['a1_y'][:], g['a1_z'][:]], -1)
        p2 = np.stack([g['a2_x'][:], g['a2_y'][:], g['a2_z'][:]], -1)
    with h5py.File(tru, 'r') as f:
        t = f['truth']; s = np.stack([t['source_x'][:], t['source_y'][:], t['source_z'][:]], -1)
    delta = (np.linalg.norm(pp - s, axis=-1) - np.linalg.norm(pp - P, axis=-1)) / C
    tau_src = tau + delta
    r, along, perp = geom_parts(s, P, p1, p2)
    tubes = AP.TUBE_LAYOUTS[cfg]; asg = AP.assign_tube(s[:, 0], s[:, 1], tubes, AP.TUBE_R_MM_BY_CONFIG[cfg])
    ref = {}
    for line in Path(csv).read_text().splitlines()[1:]:
        x = line.split(','); ref[x[1]] = dict(tau_cfg=float(x[5]), tau_rec=float(x[6]), sigma=float(x[7]), n_used=int(x[9]))
    rng = np.random.default_rng(20261014)
    rows = []
    for i, (name, cx, cy) in enumerate(tubes):
        m = asg == i
        a, b = tau[m], tau_src[m]
        tp_, sp_, n_ = AP.mle_tau(a, TAU_MIN); ts_, ss_, ns_ = AP.mle_tau(b, TAU_MIN)
        n = int(m.sum()); B = 200; d = np.empty(B)
        for k in range(B):
            w = rng.poisson(1.0, n).astype(float)
            sa = (a > TAU_MIN) * w; sb = (b > TAU_MIN) * w
            d[k] = ((b * sb).sum() / sb.sum()) - ((a * sa).sum() / sa.sum())
        fw = m & (tau > TAU_MIN)
        dh = np.histogram(delta[m], DBINS)[0]; dfw = np.histogram(delta[fw], DBINS)[0]
        # independence of the error from the delay: mean delta in delay bins of the fit population
        dbin = {}
        for lo, hi in ((1, 2), (2, 3), (3, 5), (5, 10)):
            mm = m & (tau > lo) & (tau <= hi)
            dbin[f'{lo}-{hi}ns'] = dict(n=int(mm.sum()), mean_ps=float(1000 * delta[mm].mean()) if mm.any() else None,
                                        sem_ps=float(1000 * delta[mm].std() / np.sqrt(max(mm.sum(), 1))) if mm.any() else None)
        rows.append(dict(tube=f'T{i+1}', material=name, tau_cfg=ref[name]['tau_cfg'],
                         tau_prod=tp_, sigma_prod=sp_, n_prod=n_, tau_src=ts_, sigma_src=ss_, n_src=ns_,
                         reproduces_archived=bool(abs(tp_ - ref[name]['tau_rec']) < 1e-9 and n_ == ref[name]['n_used']),
                         d_ps=1000 * (ts_ - tp_), d_sd_ps=float(1000 * d.std()),
                         n_crossing_censor=int(((a > TAU_MIN) != (b > TAU_MIN)).sum()),
                         delta_all=stats_from_hist(dh, DBINS), delta_fitwindow=stats_from_hist(dfw, DBINS),
                         delta_by_delay=dbin,
                         disp_mm=stats_from_hist(np.histogram(r[fw], RBINS)[0], RBINS),
                         along_mm=stats_from_hist(np.histogram(along[fw], ABINS)[0], ABINS),
                         perp_mm=stats_from_hist(np.histogram(perp[fw], PBINS)[0], PBINS)))
        print(f"  {sim} {name:26s} prod {tp_:.5f} (archived {ref[name]['tau_rec']:.5f})  src {ts_:.5f}  "
              f"d {1000*(ts_-tp_):+.2f} +- {1000*d.std():.2f} ps", flush=True)
    return rows, dict(delta_all=np.histogram(delta, DBINS)[0], perp=np.histogram(perp[tau > TAU_MIN], PBINS)[0],
                      along=np.histogram(along[tau > TAU_MIN], ABINS)[0], disp=np.histogram(r[tau > TAU_MIN], RBINS)[0])


def trend(rows, key):
    cfg = np.array([r['tau_cfg'] for r in rows]); rec = np.array([r[key] for r in rows])
    sig = np.array([r['sigma_prod'] for r in rows]) * 1000
    rps = (rec - cfg) * 1000

    def chi2_dof(ss):
        w = 1 / (sig ** 2 + ss ** 2); mu = (w * rps).sum() / w.sum()
        return (w * (rps - mu) ** 2).sum() / (len(rps) - 1)
    lo, hi = 0.0, 60.0
    if chi2_dof(lo) <= 1: ss = 0.0
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if chi2_dof(mid) > 1: lo = mid
            else: hi = mid
        ss = 0.5 * (lo + hi)
    w = 1 / (sig ** 2 + ss ** 2); X = np.stack([np.ones_like(cfg), cfg], 1)
    cov = np.linalg.inv(X.T @ (w[:, None] * X)); beta = cov @ (X.T @ (w * rps))
    return dict(slope_ps_per_ns=float(beta[1]), slope_err=float(np.sqrt(cov[1, 1])), intercept_ps=float(beta[0]),
                sigma_sys_ps=float(ss), mean_abs_resid_pct_sc=float(np.mean([abs(r[key] / r['tau_cfg'] - 1) * 100 for r in rows[:8]])) if len(rows) == 12 else None)


# ============================================================ Part B: clinical C0
SRC = Path(os.environ.get("PLI_RAW", "<missing>")) / 'feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5'
CFG = json.loads((RES / '2026-08-17_derivation_test/sim1_voi_config.json').read_text())['sim1_tubes']
VOIS, OUTLINE = CFG['vois'], CFG['phantom_outline']
NB = len(SX.TAU_BINS) - 1
RULES = ['closest', 'all']


def one_shard(args):
    sh_id, a, b = args
    with h5py.File(SRC, 'r') as f:
        g = f['singles']
        t = g['time_ns'][a:b].astype(np.float64); in5 = g['in_511'][a:b].astype(bool); inp = g['in_prompt'][a:b].astype(bool)
        gx, gy, gz = g['gx'][a:b], g['gy'][a:b], g['gz'][a:b]
        sx, sy, sz = g['sx'][a:b], g['sy'][a:b], g['sz'][a:b]; eid = g['eventID'][a:b]
    H = np.zeros((2, len(RULES), len(VOIS) + 1, 2, NB), dtype=np.int32)
    DH = np.zeros((len(RULES), 2, len(DBINS) - 1), dtype=np.int64)     # [rule, true-triple(1)/all(0), bin], fit window, signal arm
    GH = np.zeros((3, len(RBINS) - 1), dtype=np.int64); GA = np.zeros(len(ABINS) - 1, dtype=np.int64); GP = np.zeros(len(PBINS) - 1, dtype=np.int64)
    for ir, rule in enumerate(RULES):
        rng = np.random.default_rng([90900, int(sh_id)])           # as photon_conventions_c0.py (per rule)
        ts = t + rng.normal(0.0, 0.0909, size=len(t))
        i5 = np.flatnonzero(in5); ip_ = np.flatnonzero(inp)
        if len(i5) < 2 or len(ip_) < 1: continue
        t5 = ts[i5]; tpr = ts[ip_]
        o5 = np.argsort(t5, kind='stable'); op = np.argsort(tpr, kind='stable')
        i5 = i5[o5]; t5 = t5[o5]; ip_ = ip_[op]; tpr = tpr[op]
        ok = np.diff(t5) <= SX.W_PAIR_NS
        if not ok.any(): continue
        pa = i5[:-1][ok]; pb = i5[1:][ok]
        P1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); P2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
        pann = SX.tof_p_ann(P1, P2, ts[pa], ts[pb])
        geom = (SX.lor_perp(pann, P1, P2) <= OUTLINE.get('lor_tol_mm', 10.0)) & SX.in_phantom_outline(pann[:, 0], pann[:, 1], pann[:, 2], OUTLINE)
        if not geom.any(): continue
        pa = pa[geom]; pb = pb[geom]; pann = pann[geom]
        tm = 0.5 * (ts[pa] + ts[pb])
        for arm, shift in ((0, 0.0), (1, SX.DELAY_NS)):
            if rule == 'all':
                sel, pick = SX.all_prompt_indices(tm + shift, tpr, SX.W_PROMPT_NS)
            else:
                has, pick = SX.closest_prompt_indices(tm + shift, tpr, SX.W_PROMPT_NS); sel = np.flatnonzero(has)
            if len(sel) == 0: continue
            a1 = pa[sel]; a2 = pb[sel]; pn = pann[sel]; q = ip_[pick]
            p1 = np.stack([gx[a1], gy[a1], gz[a1]], -1); p2 = np.stack([gx[a2], gy[a2], gz[a2]], -1)
            pq = np.stack([gx[q], gy[q], gz[q]], -1); sq = np.stack([sx[q], sy[q], sz[q]], -1)
            tau_p = SX.per_event_tau(p1, p2, pq, ts[a1], ts[a2], ts[q], pn) + shift
            delta = (np.linalg.norm(pq - sq, axis=-1) - np.linalg.norm(pq - pn, axis=-1)) / SX.C_MM_NS
            tau_s = (tau_p + delta).astype(np.float32)
            tau_p32 = tau_p.astype(np.float32)          # production histograms are filled from float32 (sim1_extract)
            masks = SX.voi_masks(pn, VOIS); bg = SX.voi_bg_mask(pn, VOIS)
            for iv, tt in enumerate((tau_p32, tau_s)):
                for i, mk in enumerate(masks):
                    if mk.any(): H[iv, ir, i + 1, arm] = np.histogram(tt[mk], bins=SX.TAU_BINS)[0]
                if bg.any(): H[iv, ir, 0, arm] = np.histogram(tt[bg], bins=SX.TAU_BINS)[0]
            if arm == 0:
                fw = (tau_p > 1.0) & (tau_p < 5.0)
                true = (eid[a1] == eid[a2]) & (eid[a1] == eid[q])
                DH[ir, 0] += np.histogram(delta[fw], DBINS)[0]; DH[ir, 1] += np.histogram(delta[fw & true], DBINS)[0]
                if rule == 'all':
                    r, al, pe = geom_parts(sq[fw & true], pn[fw & true], p1[fw & true], p2[fw & true])
                    GH[0] += np.histogram(r, RBINS)[0]; GA += np.histogram(al, ABINS)[0]; GP += np.histogram(pe, PBINS)[0]
    return sh_id, H, DH, GH, GA, GP


Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
TG = np.linspace(0.3, 12.0, 6000); TT = np.linspace(1.0, 5.0, 8001)


def moment_curve(sig_e):
    A = np.ones_like(TT) if sig_e is None else Phi((5 - TT) / sig_e) * Phi((TT + 5) / sig_e)
    E = np.exp(-TT[None, :] / TG[:, None]) * A[None, :]
    return (E * TT).sum(1) / E.sum(1)


CURVES = {'sharp': moment_curve(None), 'acc359': moment_curve(0.35897)}
BC = 0.5 * (SX.TAU_BINS[:-1] + SX.TAU_BINS[1:]); SEL = (BC > 1.0) & (BC < 5.0)


def readouts(Hs):  # Hs [..., voi+1, 2, bin]
    hs = Hs[..., 0, :][..., SEL] - Hs[..., 1, :][..., SEL]
    mean = (hs * BC[SEL]).sum(-1) / hs.sum(-1)
    return {k: np.interp(mean, m, TG, left=np.nan, right=np.nan) for k, m in CURVES.items()}


def clinical():
    with h5py.File(SRC, 'r') as f: sh = f['singles']['shard_id'][:]
    bounds = SX.find_shard_boundaries(sh); t0 = time.time(); res = {}
    with Pool(24) as pool:
        for i, out in enumerate(pool.imap_unordered(one_shard, bounds)):
            res[out[0]] = out[1:]
            if i % 16 == 0: print(f'  C0 shard {i+1}/{len(bounds)}  {time.time()-t0:.0f} s', flush=True)
    ids = sorted(res)
    H = np.stack([res[s][0] for s in ids]); DH = sum(res[s][1] for s in ids)
    GH = sum(res[s][2] for s in ids); GA = sum(res[s][3] for s in ids); GP = sum(res[s][4] for s in ids)
    np.savez_compressed(OUT / 'prompt_path_c0.npz', hist=H, delta_hist=DH, disp=GH, along=GA, perp=GP, shard_ids=np.array(ids),
                        bins=SX.TAU_BINS, dbins=DBINS, rbins=RBINS, abins=ABINS, pbins=PBINS)
    # exact reproduction of the production histograms
    Z = np.load(RES / '2026-09-30_feedback3/c0_sim1/photon_conventions_c0/photon_conventions.npz')
    assert list(Z['shard_ids']) == ids
    same = bool(np.array_equal(Z['hist'][:, 0], H[:, 0]))
    tcfg = np.array([v['tau_cfg_ns'] for v in VOIS]); names = [v['name'] for v in VOIS]
    T = readouts(H.sum(0))      # [var, rule, voi+1]
    rng = np.random.default_rng(14); S = H.shape[0]; B = 1000
    flat = H.reshape(S, -1).astype(float); boots = {k: [] for k in CURVES}
    for _ in range(0, B, 250):
        w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(250)]).astype(float)
        tb = readouts((w @ flat).reshape((250,) + H.shape[1:]))
        for k in CURVES: boots[k].append(tb[k])
    boots = {k: np.concatenate(v) for k, v in boots.items()}
    out = dict(prod_reproduces_production_histograms=same, tubes=[], summary={})
    for k in CURVES:
        for ir, rule in enumerate(RULES):
            bp = 100 * (T[k][0, ir, 1:] / tcfg - 1); bs = 100 * (T[k][1, ir, 1:] / tcfg - 1)
            dd = 1000 * (T[k][1, ir, 1:] - T[k][0, ir, 1:])
            dsd = 1000 * np.nanstd(boots[k][:, 1, ir, 1:] - boots[k][:, 0, ir, 1:], axis=0)
            mab = np.nanmean(np.abs(100 * (boots[k][:, 1, ir, 1:] / tcfg - 1)), 1) - np.nanmean(np.abs(100 * (boots[k][:, 0, ir, 1:] / tcfg - 1)), 1)
            out['summary'][f'{rule}_{k}'] = dict(mean_abs_bias_prod_pct=float(np.mean(np.abs(bp))), mean_abs_bias_src_pct=float(np.mean(np.abs(bs))),
                                                 d_mean_abs_bias_pp=float(np.mean(np.abs(bs)) - np.mean(np.abs(bp))), d_mean_abs_bias_sd=float(np.std(mab)),
                                                 mean_signed_prod_pct=float(bp.mean()), mean_signed_src_pct=float(bs.mean()),
                                                 per_tube_d_ps=[float(x) for x in dd], per_tube_d_sd_ps=[float(x) for x in dsd],
                                                 max_abs_d_over_sd=float(np.max(np.abs(dd) / dsd)))
        dB_p = 100 * (T[k][0, 1, 1:] - T[k][0, 0, 1:]) / tcfg; dB_s = 100 * (T[k][1, 1, 1:] - T[k][1, 0, 1:]) / tcfg
        ddb = np.nanmean(100 * ((boots[k][:, 1, 1, 1:] - boots[k][:, 1, 0, 1:]) - (boots[k][:, 0, 1, 1:] - boots[k][:, 0, 0, 1:])) / tcfg, 1)
        out['summary'][f'rule_difference_{k}'] = dict(mean_dB_prod_pp=float(dB_p.mean()), mean_dB_src_pp=float(dB_s.mean()),
                                                      change_pp=float(dB_s.mean() - dB_p.mean()), change_sd_pp=float(np.std(ddb)))
    for ir, rule in enumerate(RULES):
        out['summary'][f'delta_fitwindow_{rule}_all_ps'] = {k: (v * 1000 if k not in ('n',) else v) for k, v in stats_from_hist(DH[ir, 0], DBINS).items()}
        out['summary'][f'delta_fitwindow_{rule}_true_ps'] = {k: (v * 1000 if k not in ('n',) else v) for k, v in stats_from_hist(DH[ir, 1], DBINS).items()}
    out['summary']['true_triples_all_rule_disp_mm'] = stats_from_hist(GH[0], RBINS)
    out['summary']['true_triples_all_rule_along_mm'] = stats_from_hist(GA, ABINS)
    out['summary']['true_triples_all_rule_perp_mm'] = stats_from_hist(GP, PBINS)
    out['tubes'] = names; out['tau_cfg'] = [float(x) for x in tcfg]
    return out


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = dict(doc=__doc__)
    A1, h1 = closure('sim1_calibration', 'sim1', RES / '2026-08-13_overnight/O4_rerun/sim1_calibration/analysis/post_taumin1p0.csv')
    A2, h2 = closure('sim2_mercolli_te', 'sim2', RES / '2026-08-13_overnight/O4_rerun/sim2_mercolli_te/analysis/post_taumin1p0.csv')
    rows = A1 + A2
    res['closure'] = dict(rows=rows, trend_prod=trend(rows, 'tau_prod'), trend_src=trend(rows, 'tau_src'),
                          all_reproduce_archived=all(r['reproduces_archived'] for r in rows),
                          max_abs_d_ps=float(max(abs(r['d_ps']) for r in rows)),
                          max_abs_d_over_sd=float(max(abs(r['d_ps']) / r['d_sd_ps'] for r in rows)),
                          mean_d_ps=float(np.mean([r['d_ps'] for r in rows])),
                          d_slope=trend([dict(tau_cfg=r['tau_cfg'], tau_prod=r['tau_cfg'] + r['d_ps'] / 1000,
                                              sigma_prod=r['d_sd_ps'] / 1000) for r in rows], 'tau_prod'))
    sc = h1['delta_all'] + 0; res['closure']['delta_all_sc_ps'] = {k: (v * 1000 if k != 'n' else v) for k, v in stats_from_hist(sc, DBINS).items()}
    res['closure']['perp_fitwindow_sc_mm'] = stats_from_hist(h1['perp'], PBINS)
    res['closure']['along_fitwindow_sc_mm'] = stats_from_hist(h1['along'], ABINS)
    res['closure']['disp_fitwindow_sc_mm'] = stats_from_hist(h1['disp'], RBINS)
    res['closure']['delta_all_i_ps'] = {k: (v * 1000 if k != 'n' else v) for k, v in stats_from_hist(h2['delta_all'], DBINS).items()}
    print(json.dumps({k: v for k, v in res['closure'].items() if k != 'rows'}, indent=1), flush=True)
    res['clinical_c0'] = clinical()
    res['wall_s'] = time.time() - t0
    (OUT / 'prompt_path.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res['clinical_c0']['summary'], indent=1))


if __name__ == '__main__':
    main()
