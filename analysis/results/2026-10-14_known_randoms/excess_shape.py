#!/usr/bin/env python3
"""Shape of the all-prompts count excess of the randoms arm over the uncorrelated accidental content of the signal
window, exc = rnd - AP - PA_unc (known_randoms_hists.npz), pooled over the eight tubes, per rung: least-squares
decomposition on -1..8 ns into a flat term and a term shaped like the genuine-triple spectrum G (scaled per rung).
Diagnostic for the observation that the count excess grows with mu while the lifetime error dR does not."""
import json, numpy as np
Z = np.load('known_randoms_hists.npz'); H = Z['hist'].astype(float); KS = list(Z['ks']); CL = list(Z['classes']); B = Z['bins']
c = 0.5 * (B[1:] + B[:-1]); names = list(Z['voi_names']); TI = [i for i, n in enumerate(names) if n.startswith('T')]
P = H.sum(0)[:, :, TI].sum(2)                              # [k, rule, cls, bin]
ir = list(Z['rules']).index('all'); g = CL.index
out = []
for ik, k in enumerate(KS):
    exc = P[ik, ir, g('rnd')] - P[ik, ir, g('AP')] - P[ik, ir, g('PA_unc')]
    G = P[ik, ir, g('G')]; var = P[ik, ir, g('rnd')] + P[ik, ir, g('AP')] + P[ik, ir, g('PA_unc')]
    row = dict(k=k)
    for lab, (lo, hi) in {'fit': (1.0, 5.0), 'wide': (-1.0, 8.0), 'early': (-1.0, 1.0), 'late': (5.5, 8.0)}.items():
        m = (c > lo) & (c < hi)
        row[f'exc_{lab}'] = float(exc[m].sum()); row[f'exc_{lab}_sd'] = float(np.sqrt(var[m].sum()))
        row[f'G_{lab}'] = float(G[m].sum()); row[f'rnd_{lab}'] = float(P[ik, ir, g('rnd')][m].sum())
    m = (c > -1) & (c < 8); w = 1 / np.maximum(var[m], 1)
    X = np.c_[np.ones(m.sum()), G[m] / G[m].sum()]
    A = X.T @ (X * w[:, None]); bvec = X.T @ (w * exc[m]); coef = np.linalg.solve(A, bvec); cov = np.linalg.inv(A)
    row.update(flat_per_bin=float(coef[0]), flat_sd=float(np.sqrt(cov[0, 0])), G_shaped_counts=float(coef[1]), G_shaped_sd=float(np.sqrt(cov[1, 1])))
    out.append(row)
    print(f"k {k:.2f}: exc fit {row['exc_fit']:+8.0f}+-{row['exc_fit_sd']:.0f} ({100*row['exc_fit']/row['G_fit']:+.2f}% of G, {100*row['exc_fit']/row['rnd_fit']:+.2f}% of rnd) | "
          f"early(-1,1) {row['exc_early']:+8.0f}+-{row['exc_early_sd']:.0f} ({100*row['exc_early']/row['G_early']:+.3f}% of G) | late(5.5,8) {row['exc_late']:+6.0f}+-{row['exc_late_sd']:.0f} | "
          f"flat {coef[0]:+.2f}+-{np.sqrt(cov[0,0]):.2f}/bin, G-shaped {coef[1]:+.0f}+-{np.sqrt(cov[1,1]):.0f}")
json.dump(out, open('excess_shape.json', 'w'), indent=1)
