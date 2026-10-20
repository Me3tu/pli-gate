#!/usr/bin/env python3
"""Review round 8 (2026-10-08): pair-builder sensitivity of the matching artifact.

Question (referee): the artifact was measured with the adjacent-pair builder of Algorithm S1, which never
pairs singles i and i+2 when i+1 lies between them; the published Quadra sorter forms all combinations within
its coincidence window. Does the magnitude of the closest-prompt artifact change with the builder?

Method: identical events, identical timing smear (same RNG seed and draw order as the production extractor
sim1_extract.py), identical geometry cuts, prompt windows, delayed window (+50 ns), histogramming and
truncation inversion. Only the pair list differs:
  adjacent  : pairs (i, i+1) of time-sorted 511-window singles with dt <= 2 ns   (production, Algorithm S1)
  allpairs  : every pair (i, j), i < j, with t_j - t_i <= 2 ns                   (all combinations)
Both matching rules (closest, all) are run for both builders. Per-shard histograms are kept so that a paired
shard bootstrap gives the uncertainty of builder differences on the same events.
Writes builder_sensitivity.json and hists_builder.npz.
"""
import os, sys, json, time, argparse
from pathlib import Path
import numpy as np, h5py

D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / '2026-08-17_derivation_test'
sys.path.insert(0, str(D))
import sim1_extract as SX  # noqa: E402

TAU_BINS = SX.TAU_BINS
NB = len(TAU_BINS) - 1


def pair_lists(t_511, builder):
    if builder == 'adjacent':
        dt = np.diff(t_511); ok = dt <= SX.W_PAIR_NS
        ia = np.flatnonzero(ok); return ia, ia + 1
    A, B = [], []
    k = 1
    while k < len(t_511):
        dt = t_511[k:] - t_511[:-k]; ok = dt <= SX.W_PAIR_NS
        if not ok.any(): break
        ia = np.flatnonzero(ok); A.append(ia); B.append(ia + k); k += 1
    return np.concatenate(A), np.concatenate(B)


def shard_hists(t_true, in_511, in_prompt, gx, gy, gz, sigma_t_ns, rng, vois, outline):
    """Returns dict[(builder, rule)] -> (hs, hr) arrays [n_voi+1, NB] (index 0 = bg), plus pair counts."""
    n = len(t_true)
    t_smear = t_true + rng.normal(0.0, sigma_t_ns, size=n)   # same draw as SX.clinical_per_shard
    in_511_b = in_511.astype(bool); in_prompt_b = in_prompt.astype(bool)
    idx_511 = np.flatnonzero(in_511_b); idx_prompt = np.flatnonzero(in_prompt_b)
    out = {}; counts = {}
    if len(idx_511) < 2 or len(idx_prompt) < 1: return out, counts
    t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
    o5 = np.argsort(t_511, kind='stable'); op = np.argsort(t_prompt, kind='stable')
    idx_511 = idx_511[o5]; t_511 = t_511[o5]; idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
    for builder in ('adjacent', 'allpairs'):
        ia, ib = pair_lists(t_511, builder)
        counts[builder] = int(len(ia))
        if len(ia) == 0: continue
        pa = idx_511[ia]; pb = idx_511[ib]
        p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
        p_ann = SX.tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
        perp = SX.lor_perp(p_ann, p1, p2)
        geom = (perp <= outline.get('lor_tol_mm', 10.0)) & SX.in_phantom_outline(p_ann[:, 0], p_ann[:, 1], p_ann[:, 2], outline)
        counts[builder + '_geom'] = int(geom.sum())
        pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
        if len(pa) == 0: continue
        t_pair = 0.5 * (t_smear[pa] + t_smear[pb])
        for rule in ('closest', 'all'):
            if rule == 'all':
                sel_s, pick_s = SX.all_prompt_indices(t_pair, t_prompt, SX.W_PROMPT_NS)
                sel_r, pick_r = SX.all_prompt_indices(t_pair + SX.DELAY_NS, t_prompt, SX.W_PROMPT_NS)
            else:
                has_s, pick_s = SX.closest_prompt_indices(t_pair, t_prompt, SX.W_PROMPT_NS)
                has_r, pick_r = SX.closest_prompt_indices(t_pair + SX.DELAY_NS, t_prompt, SX.W_PROMPT_NS)
                sel_s = np.flatnonzero(has_s); sel_r = np.flatnonzero(has_r)
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
            out[(builder, rule)] = (H[0], H[1])
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
    cfg = json.loads(Path(a.config).read_text())[a.phantom]; vois = cfg['vois']; outline = cfg['phantom_outline']
    names = ['bg'] + [v['name'] for v in vois]; tcfg = [None] + [v.get('tau_cfg_ns') for v in vois]
    f = h5py.File(a.singles, 'r'); g = f['singles']
    bounds = SX.find_shard_boundaries(g['shard_id'][:])
    if a.max_shards: bounds = bounds[:a.max_shards]
    rng = np.random.default_rng(int(a.sigma_t_ps * 1000) or 1)
    keys = [(b, r) for b in ('adjacent', 'allpairs') for r in ('closest', 'all')]
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
    np.savez_compressed(Path(a.out) / 'hists_builder.npz', **{f'{b}_{r}_sig': HS[(b, r)] for b, r in keys},
                        **{f'{b}_{r}_rnd': HR[(b, r)] for b, r in keys}, bins=TAU_BINS)

    def table(idx=None):
        res = {}
        for k in keys:
            hs = HS[k] if idx is None else HS[k][idx]; hr = HR[k] if idx is None else HR[k][idx]
            hs = hs.sum(0); hr = hr.sum(0)
            res[k] = np.array([readout(hs[v], hr[v]) for v in range(len(vois) + 1)])
        return res

    cen = table()
    brng = np.random.default_rng(20261008)
    boots = [table(brng.integers(0, len(bounds), len(bounds))) for _ in range(a.boot)]
    rows = []
    for v in range(1, len(vois) + 1):
        r = dict(voi=names[v], tau_cfg=tcfg[v])
        for b in ('adjacent', 'allpairs'):
            for rule in ('closest', 'all'):
                r[f'{b}_{rule}_ns'] = float(cen[(b, rule)][v])
                r[f'{b}_{rule}_bias_pct'] = float(100 * (cen[(b, rule)][v] / tcfg[v] - 1))
            # rule artifact (signed and folded) per builder
            bc_, ba_ = r[f'{b}_closest_bias_pct'], r[f'{b}_all_bias_pct']
            r[f'{b}_dB_pp'] = ba_ - bc_; r[f'{b}_dAbsB_pp'] = abs(bc_) - abs(ba_)
            bb = np.array([100 * (x[(b, 'all')][v] - x[(b, 'closest')][v]) / tcfg[v] for x in boots])
            r[f'{b}_dB_sd'] = float(np.nanstd(bb))
        dd = np.array([100 * ((x[('allpairs', 'all')][v] - x[('allpairs', 'closest')][v]) - (x[('adjacent', 'all')][v] - x[('adjacent', 'closest')][v])) / tcfg[v] for x in boots])
        r['ddB_allpairs_minus_adjacent_pp'] = r['allpairs_dB_pp'] - r['adjacent_dB_pp']; r['ddB_sd'] = float(np.nanstd(dd))
        for rule in ('closest', 'all'):
            db = np.array([100 * (x[('allpairs', rule)][v] - x[('adjacent', rule)][v]) / tcfg[v] for x in boots])
            r[f'd_{rule}_allpairs_minus_adjacent_pp'] = r[f'allpairs_{rule}_bias_pct'] - r[f'adjacent_{rule}_bias_pct']
            r[f'd_{rule}_sd'] = float(np.nanstd(db))
        rows.append(r)
    summ = {}
    for b in ('adjacent', 'allpairs'):
        for rule in ('closest', 'all'):
            summ[f'{b}_{rule}_mean_abs_bias_pct'] = float(np.mean([abs(r[f'{b}_{rule}_bias_pct']) for r in rows]))
        summ[f'{b}_mean_dB_pp'] = float(np.mean([r[f'{b}_dB_pp'] for r in rows]))
    summ['mean_ddB_pp'] = float(np.mean([r['ddB_allpairs_minus_adjacent_pp'] for r in rows]))
    bm = np.array([np.mean([100 * ((x[('allpairs', 'all')][v] - x[('allpairs', 'closest')][v]) - (x[('adjacent', 'all')][v] - x[('adjacent', 'closest')][v])) / tcfg[v] for v in range(1, len(vois) + 1)]) for x in boots])
    summ['mean_ddB_sd'] = float(np.nanstd(bm))
    summ['pair_counts'] = cnt
    res = dict(doc=__doc__, singles=a.singles, sigma_t_ps=a.sigma_t_ps, n_shards=len(bounds), boot=a.boot,
               rows=rows, summary=summ, wall_s=time.time() - t0)
    Path(a.out, 'builder_sensitivity.json').write_text(json.dumps(res, indent=1))
    for r in rows:
        print(f"{r['voi']:28s} adj c/a {r['adjacent_closest_bias_pct']:+6.2f} {r['adjacent_all_bias_pct']:+6.2f}  "
              f"allp c/a {r['allpairs_closest_bias_pct']:+6.2f} {r['allpairs_all_bias_pct']:+6.2f}  "
              f"ddB {r['ddB_allpairs_minus_adjacent_pp']:+5.2f}+-{r['ddB_sd']:.2f}")
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
