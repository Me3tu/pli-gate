#!/usr/bin/env python3
"""Review round 13, verifier follow-up: within the fit window (1-5 ns), the mean delay of the all-prompts count excess
exc = rnd - AP - PA_unc (known_randoms_hists.npz, eight tubes pooled) against the mean delay of the genuine triples G
and against the window midpoint (3.0 ns, the mean of a flat pedestal). Poisson variances per bin (rnd + AP + PA_unc),
first-order propagation for the mean. Also the chi2/dof of the flat + G-shaped decomposition of excess_shape.py on the
fit window, as a check of that decomposition there. Writes excess_mean.json."""
import json, numpy as np
Z = np.load('known_randoms_hists.npz'); H = Z['hist'].astype(float); KS = list(Z['ks']); CL = list(Z['classes']); B = Z['bins']
c = 0.5 * (B[1:] + B[:-1]); names = list(Z['voi_names']); TI = [i for i, n in enumerate(names) if n.startswith('T')]
P = H.sum(0)[:, :, TI].sum(2)                              # [k, rule, cls, bin]
ir = list(Z['rules']).index('all'); g = CL.index
X = json.load(open('excess_shape.json'))
MU = json.load(open('known_randoms_summary.json'))['mu']
m = (c > 1.0) & (c < 5.0); mw = (c > -1.0) & (c < 8.0)
out = []
for ik, k in enumerate(KS):
    exc = P[ik, ir, g('rnd')] - P[ik, ir, g('AP')] - P[ik, ir, g('PA_unc')]
    var = P[ik, ir, g('rnd')] + P[ik, ir, g('AP')] + P[ik, ir, g('PA_unc')]
    G = P[ik, ir, g('G')]
    e, v, x = exc[m], var[m], c[m]
    S = e.sum(); me = (x * e).sum() / S; sme = np.sqrt((((x - me) ** 2) * v).sum()) / abs(S)
    mg = (x * G[m]).sum() / G[m].sum()
    # flat + G-shaped decomposition of excess_shape.py, evaluated on the fit window
    r = X[ik]; Gw = G[mw]; pred = r['flat_per_bin'] + r['G_shaped_counts'] * G / Gw.sum()
    chi2_fit = float((((exc[m] - pred[m]) ** 2) / np.maximum(var[m], 1)).sum()); dof_fit = int(m.sum())
    chi2_w = float((((exc[mw] - pred[mw]) ** 2) / np.maximum(var[mw], 1)).sum()); dof_w = int(mw.sum()) - 2
    row = dict(k=float(k), mu=MU[ik], exc_fit=float(S), exc_fit_sd=float(np.sqrt(v.sum())), exc_pct_G=float(100 * S / G[m].sum()),
               mean_exc_ns=float(me), mean_exc_sd=float(sme), mean_G_ns=float(mg), z_vs_G=float((me - mg) / sme),
               z_vs_flat=float((me - 3.0) / sme), model_fit_counts=float(pred[m].sum()),
               chi2dof_fitwindow=chi2_fit / dof_fit, chi2dof_wide=chi2_w / dof_w)
    out.append(row)
    print(f"mu {row['mu']:.4f}: exc {S:+7.0f}+-{row['exc_fit_sd']:.0f} ({row['exc_pct_G']:+.2f}% of G) mean {me:.2f}+-{sme:.2f} ns "
          f"(G {mg:.2f}; z_G {row['z_vs_G']:+.1f}, z_flat {row['z_vs_flat']:+.1f}) | model {row['model_fit_counts']:.0f} "
          f"| chi2/dof fit {row['chi2dof_fitwindow']:.2f} wide {row['chi2dof_wide']:.2f}")
# pooled over rungs (rungs share events, descriptive only)
w = np.array([1 / r['mean_exc_sd'] ** 2 for r in out]); d = np.array([r['mean_exc_ns'] - r['mean_G_ns'] for r in out])
pool = dict(diff_vs_G_ns=float((w * d).sum() / w.sum()), diff_sd_ns=float(1 / np.sqrt(w.sum())),
            note='inverse-variance mean over rungs; rungs share events, so the uncertainty is optimistic')
print('pooled mean(exc) - mean(G):', round(pool['diff_vs_G_ns'], 3), '+-', round(pool['diff_sd_ns'], 3))
json.dump(dict(rows=out, pooled=pool), open('excess_mean.json', 'w'), indent=1)
