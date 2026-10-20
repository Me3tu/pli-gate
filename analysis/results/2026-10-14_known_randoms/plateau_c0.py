#!/usr/bin/env python3
"""Review round 13: pedestal plateau numbers of the main text, Table S16 and Note S21, recomputed by a script (the
values in README section 3 were computed ad hoc with the f08_pedestal.py definition). For the corrected-yield arm C0
(2026-09-30_feedback3/c0_sim1/hists_c0_0ps.npz) and the discovery stream (2026-08-22_G3_consequences/hists_sim1.npz):
d = (all_sig - all_rnd) - (closest_sig - closest_rnd) per 50 ps bin; plateau = mean of d over 1.25-4.5 ns (bin
centers), per tube and averaged over the eight tubes; share = sum of d / sum of closest_rnd over the same bins, per
tube; slope = least-squares slope of the tube-averaged d against delay on 1.25-4.5 ns, shard bootstrap (1000
resamples, seed 20261014, half 16-84 width); first = mean of d in the bins with centers in 1.0-1.1 ns. Writes
plateau_c0.json."""
import json, os
import numpy as np
RES = os.environ.get('PLI_RES') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = {'c0': RES + '/2026-09-30_feedback3/c0_sim1/hists_c0_0ps.npz', 'discovery': RES + '/2026-08-22_G3_consequences/hists_sim1.npz'}
out = {}
for tag, path in SRC.items():
    Z = np.load(path, allow_pickle=True)
    names = [str(s) for s in Z['voi_names']]; idx = [i for i, n in enumerate(names) if n.startswith('T')]
    b = Z['tau_bins']; c = 0.5 * (b[:-1] + b[1:]); m = (c > 1.25) & (c < 4.5); mf = (c > 1.0) & (c < 1.1)
    CS, CR, AS, AR = (Z[k][:, idx, :].astype(float) for k in ('closest_sig', 'closest_rnd', 'all_sig', 'all_rnd'))

    def stats(sh):
        d = (AS[sh].sum(0) - AR[sh].sum(0)) - (CS[sh].sum(0) - CR[sh].sum(0))       # [tube, bin]
        dm = d.mean(0)
        return d, dm, np.polyfit(c[m], dm[m], 1)[0]
    nsh = CS.shape[0]; allsh = np.arange(nsh)
    d, dm, slope = stats(allsh)
    per = d[:, m].mean(1); share = 100 * d[:, m].sum(1) / CR[:, :, m].sum(0).sum(1)
    rng = np.random.default_rng(20261014)
    bs = np.array([stats(rng.integers(0, nsh, nsh))[2] for _ in range(1000)])
    lo, hi = np.percentile(bs, [16, 84])
    out[tag] = dict(plateau_mean=float(dm[m].mean()), plateau_tube_min=float(per.min()), plateau_tube_max=float(per.max()),
                    share_pct_min=float(share.min()), share_pct_max=float(share.max()), slope_per_ns=float(slope),
                    slope_sd=float(0.5 * (hi - lo)), first_0p1ns=float(dm[mf].mean()), tubes=[names[i] for i in idx],
                    plateau_per_tube=per.tolist(), share_pct_per_tube=share.tolist())
    print(tag, {k: (round(v, 3) if isinstance(v, float) else '') for k, v in out[tag].items() if not isinstance(v, list)})
json.dump(out, open('plateau_c0.json', 'w'), indent=1)
