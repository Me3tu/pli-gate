#!/usr/bin/env python3
"""Review round 9 (2026-10-09): annihilation-pair window and triple-candidate window of the published Quadra sorter.

Question (referee round 9): the round-8 builder test (2026-10-08_builder_sensitivity) varied adjacency at a fixed
2 ns pair window; the published Quadra sorter forms all triple combinations within a 20 ns coincidence window and
then requires 4.2 ns between the annihilation photons (Steinberger et al. 2024). Does the matching artifact, or the
all-prompts residual, change when all pairs within 4.2 ns are formed, and does a 20 ns triple-candidate window change
the triples inside the fit window?

Method: as builder_sensitivity.py (identical events, timing smear with the seed and draw order of sim1_extract.py,
geometric cuts, prompt windows, delayed window +50 ns, histograms, truncation inversion over [1, 5] ns). Arms:
  adj2        adjacent pairs, dt <= 2 ns, prompt search +-5 ns          (production, Algorithm S1)
  all2        all pairs,      dt <= 2 ns, prompt search +-5 ns          (round-8 test)
  all4p2      all pairs,      dt <= 4.2 ns, prompt search +-5 ns        (published annihilation-photon window)
  adj4p2      adjacent pairs, dt <= 4.2 ns, prompt search +-5 ns
  all4p2_t20  all pairs,      dt <= 4.2 ns, prompt search +-20 ns       (published 20 ns triple window; the estimator
                                                                         window [1, 5] ns is unchanged)
Both matching rules for every arm; per-shard histograms; paired shard bootstrap (B = 400) of every arm against adj2.
Writes pair_window.json and hists_pair_window.npz.
"""
import os, sys, json, time, argparse
from pathlib import Path
import numpy as np, h5py

D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / '2026-08-17_derivation_test'
sys.path.insert(0, str(D))
import sim1_extract as SX  # noqa: E402

TAU_BINS = SX.TAU_BINS
NB = len(TAU_BINS) - 1
ARMS = [('adj2', 'adjacent', 2.0, 5.0), ('all2', 'allpairs', 2.0, 5.0), ('all4p2', 'allpairs', 4.2, 5.0),
        ('adj4p2', 'adjacent', 4.2, 5.0), ('all4p2_t20', 'allpairs', 4.2, 20.0)]
RULES = ('closest', 'all')


def pair_lists(t_511, builder, w):
    if builder == 'adjacent':
        dt = np.diff(t_511); ok = dt <= w
        ia = np.flatnonzero(ok); return ia, ia + 1
    A, B = [], []
    k = 1
    while k < len(t_511):
        dt = t_511[k:] - t_511[:-k]; ok = dt <= w
        if not ok.any(): break
        ia = np.flatnonzero(ok); A.append(ia); B.append(ia + k); k += 1
    return (np.concatenate(A), np.concatenate(B)) if A else (np.zeros(0, int), np.zeros(0, int))


def shard_hists(t_true, in_511, in_prompt, gx, gy, gz, sigma_t_ns, rng, vois, outline):
    n = len(t_true)
    t_smear = t_true + rng.normal(0.0, sigma_t_ns, size=n)   # same draw as SX.clinical_per_shard
    in_511_b = in_511.astype(bool); in_prompt_b = in_prompt.astype(bool)
    idx_511 = np.flatnonzero(in_511_b); idx_prompt = np.flatnonzero(in_prompt_b)
    out = {}; counts = {}
    if len(idx_511) < 2 or len(idx_prompt) < 1: return out, counts
    t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
    o5 = np.argsort(t_511, kind='stable'); op = np.argsort(t_prompt, kind='stable')
    idx_511 = idx_511[o5]; t_511 = t_511[o5]; idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
    for arm, builder, wpair, wprompt in ARMS:
        ia, ib = pair_lists(t_511, builder, wpair)
        counts[arm] = int(len(ia))
        if len(ia) == 0: continue
        pa = idx_511[ia]; pb = idx_511[ib]
        p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
        p_ann = SX.tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
        perp = SX.lor_perp(p_ann, p1, p2)
        geom = (perp <= outline.get('lor_tol_mm', 10.0)) & SX.in_phantom_outline(p_ann[:, 0], p_ann[:, 1], p_ann[:, 2], outline)
        counts[arm + '_geom'] = int(geom.sum())
        pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
        if len(pa) == 0: continue
        t_pair = 0.5 * (t_smear[pa] + t_smear[pb])
        for rule in RULES:
            if rule == 'all':
                sel_s, pick_s = SX.all_prompt_indices(t_pair, t_prompt, wprompt)
                sel_r, pick_r = SX.all_prompt_indices(t_pair + SX.DELAY_NS, t_prompt, wprompt)
            else:
                has_s, pick_s = SX.closest_prompt_indices(t_pair, t_prompt, wprompt)
                has_r, pick_r = SX.closest_prompt_indices(t_pair + SX.DELAY_NS, t_prompt, wprompt)
                sel_s = np.flatnonzero(has_s); sel_r = np.flatnonzero(has_r)
            counts[f'{arm}_{rule}_sig_triples'] = counts.get(f'{arm}_{rule}_sig_triples', 0) + int(len(sel_s))
            H = []
            for sel, pick, shift in ((sel_s, pick_s, 0.0), (sel_r, pick_r, SX.DELAY_NS)):
                h = np.zeros((len(vois) + 1, NB))
                if len(sel):
                    a, b, pp_ = pa[sel], pb[sel], p_ann[sel]; ip = idx_prompt[pick]
                    tau = SX.per_event_tau(np.stack([gx[a], gy[a], gz[a]], -1), np.stack([gx[b], gy[b], gz[b]], -1),
                                           np.stack([gx[ip], gy[ip], gz[ip]], -1), t_smear[a], t_smear[b], t_smear[ip], pp_) + shift
                    masks = SX.voi_masks(pp_, vois)
                    for i, m in enumerate(masks):
                        if m.any(): h[i + 1] = np.histogram(tau[m], bins=TAU_BINS)[0]
                    bg = SX.voi_bg_mask(pp_, vois)
                    if bg.any(): h[0] = np.histogram(tau[bg], bins=TAU_BINS)[0]
                H.append(h)
            out[(arm, rule)] = (H[0], H[1])
    return out, counts


def readout(hs, hr):
    bc = 0.5 * (TAU_BINS[:-1] + TAU_BINS[1:]); sel = (bc > SX.TAU_MIN_NS) & (bc < SX.W_PROMPT_NS)
    s = hs[sel] - hr[sel]; n = s.sum()
    if n <= 0: return np.nan
    tau_obs = (bc[sel] * s).sum() / n - SX.TAU_MIN_NS
    return SX.truncation_correct(tau_obs, SX.W_PROMPT_NS - SX.TAU_MIN_NS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--singles', default=str(Path(os.environ.get("PLI_RAW", "<missing>")) / 'feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5'))
    ap.add_argument('--config', default=str(D / 'sim1_voi_config.json'))
    ap.add_argument('--phantom', default='sim1_tubes')
    ap.add_argument('--sigma-t-ps', type=float, default=90.9)
    ap.add_argument('--boot', type=int, default=400)
    ap.add_argument('--out', default=str(Path(__file__).parent))
    ap.add_argument('--max-shards', type=int, default=None)
    a = ap.parse_args()
    if not Path(a.singles).exists(): raise SystemExit(f"{a.singles}: raw data not included in the release, see analysis/README.md")
    Path(a.out).mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(a.config).read_text())[a.phantom]; vois = cfg['vois']; outline = cfg['phantom_outline']
    names = ['bg'] + [v['name'] for v in vois]; tcfg = [None] + [v.get('tau_cfg_ns') for v in vois]
    f = h5py.File(a.singles, 'r'); g = f['singles']
    bounds = SX.find_shard_boundaries(g['shard_id'][:])
    if a.max_shards: bounds = bounds[:a.max_shards]
    rng = np.random.default_rng(int(a.sigma_t_ps * 1000) or 1)
    keys = [(arm[0], r) for arm in ARMS for r in RULES]
    HS = {k: np.zeros((len(bounds), len(vois) + 1, NB)) for k in keys}
    HR = {k: np.zeros((len(bounds), len(vois) + 1, NB)) for k in keys}
    cnt = {}
    t0 = time.time()
    for j, (sh, s, e) in enumerate(bounds):
        out, c = shard_hists(g['time_ns'][s:e], g['in_511'][s:e], g['in_prompt'][s:e], g['gx'][s:e], g['gy'][s:e], g['gz'][s:e],
                             a.sigma_t_ps / 1000.0, rng, vois, outline)
        for k, (hs, hr) in out.items(): HS[k][j] = hs; HR[k][j] = hr
        for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
        if j % 16 == 0: print(f'shard {j}/{len(bounds)}  {time.time() - t0:.0f} s', flush=True)
    f.close()
    np.savez_compressed(Path(a.out) / 'hists_pair_window.npz', **{f'{k[0]}_{k[1]}_sig': HS[k] for k in keys},
                        **{f'{k[0]}_{k[1]}_rnd': HR[k] for k in keys}, bins=TAU_BINS)
    # fit-window identity of the 20 ns triple window (all4p2 vs all4p2_t20), summed over shards and regions
    bc = 0.5 * (TAU_BINS[:-1] + TAU_BINS[1:]); fw = (bc > SX.TAU_MIN_NS) & (bc < SX.W_PROMPT_NS)
    ident = {}
    for rule in RULES:
        for tag, H in (('sig', HS), ('rnd', HR)):
            x = H[('all4p2', rule)][:, :, fw].sum(); y = H[('all4p2_t20', rule)][:, :, fw].sum()
            ident[f'{rule}_{tag}_fitwindow_counts'] = [float(x), float(y), float(y - x)]

    def table(idx=None):
        res = {}
        for k in keys:
            hs = HS[k] if idx is None else HS[k][idx]; hr = HR[k] if idx is None else HR[k][idx]
            hs = hs.sum(0); hr = hr.sum(0)
            res[k] = np.array([readout(hs[v], hr[v]) for v in range(len(vois) + 1)])
        return res

    cen = table()
    brng = np.random.default_rng(20261009)
    boots = [table(brng.integers(0, len(bounds), len(bounds))) for _ in range(a.boot)]
    rows = []
    vidx = list(range(1, len(vois) + 1))
    for v in vidx:
        r = dict(voi=names[v], tau_cfg=tcfg[v])
        for arm in [x[0] for x in ARMS]:
            for rule in RULES:
                r[f'{arm}_{rule}_bias_pct'] = float(100 * (cen[(arm, rule)][v] / tcfg[v] - 1))
            r[f'{arm}_dB_pp'] = r[f'{arm}_all_bias_pct'] - r[f'{arm}_closest_bias_pct']
            if arm != 'adj2':
                dd = np.array([100 * ((x[(arm, 'all')][v] - x[(arm, 'closest')][v]) - (x[('adj2', 'all')][v] - x[('adj2', 'closest')][v])) / tcfg[v] for x in boots])
                r[f'{arm}_ddB_pp'] = r[f'{arm}_dB_pp'] - r['adj2_dB_pp']; r[f'{arm}_ddB_sd'] = float(np.nanstd(dd))
                for rule in RULES:
                    db = np.array([100 * (x[(arm, rule)][v] - x[('adj2', rule)][v]) / tcfg[v] for x in boots])
                    r[f'{arm}_d_{rule}_pp'] = r[f'{arm}_{rule}_bias_pct'] - r[f'adj2_{rule}_bias_pct']; r[f'{arm}_d_{rule}_sd'] = float(np.nanstd(db))
        rows.append(r)
    summ = {}
    for arm in [x[0] for x in ARMS]:
        for rule in RULES:
            summ[f'{arm}_{rule}_mean_abs_bias_pct'] = float(np.mean([abs(r[f'{arm}_{rule}_bias_pct']) for r in rows]))
            summ[f'{arm}_{rule}_mean_signed_bias_pct'] = float(np.mean([r[f'{arm}_{rule}_bias_pct'] for r in rows]))
        summ[f'{arm}_mean_dB_pp'] = float(np.mean([r[f'{arm}_dB_pp'] for r in rows]))
        if arm != 'adj2':
            summ[f'{arm}_mean_ddB_pp'] = float(np.mean([r[f'{arm}_ddB_pp'] for r in rows]))
            bm = np.array([np.mean([100 * ((x[(arm, 'all')][v] - x[(arm, 'closest')][v]) - (x[('adj2', 'all')][v] - x[('adj2', 'closest')][v])) / tcfg[v] for v in vidx]) for x in boots])
            summ[f'{arm}_mean_ddB_sd'] = float(np.nanstd(bm))
            for rule in RULES:
                bma = np.array([np.mean([abs(100 * (x[(arm, rule)][v] / tcfg[v] - 1)) for v in vidx]) - np.mean([abs(100 * (x[('adj2', rule)][v] / tcfg[v] - 1)) for v in vidx]) for x in boots])
                summ[f'{arm}_{rule}_d_mean_abs_bias_pp'] = summ[f'{arm}_{rule}_mean_abs_bias_pct'] - summ[f'adj2_{rule}_mean_abs_bias_pct']
                summ[f'{arm}_{rule}_d_mean_abs_bias_sd'] = float(np.nanstd(bma))
    summ['counts'] = cnt
    summ['t20_fitwindow_identity'] = ident
    res = dict(doc=__doc__, singles=a.singles, sigma_t_ps=a.sigma_t_ps, n_shards=len(bounds), boot=a.boot,
               arms=ARMS, rows=rows, summary=summ, wall_s=time.time() - t0)
    Path(a.out, 'pair_window.json').write_text(json.dumps(res, indent=1))
    for r in rows:
        print(f"{r['voi']:28s} " + '  '.join(f"{arm} c/a {r[f'{arm}_closest_bias_pct']:+6.2f}/{r[f'{arm}_all_bias_pct']:+6.2f}" for arm in [x[0] for x in ARMS]))
    print(json.dumps({k: v for k, v in summ.items() if k != 'counts'}, indent=1)); print(json.dumps(cnt))


if __name__ == '__main__':
    main()
