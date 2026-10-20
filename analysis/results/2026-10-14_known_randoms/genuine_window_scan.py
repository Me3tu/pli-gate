#!/usr/bin/env python3
"""Review round 14 (referee 3: origin of the -0.5 to -0.6% genuine-triple residual under all-prompts matching).
Hypothesis tested: short-delay events cross the fixed 1 ns lower edge of the fit window and lower the windowed mean.
(Correction after the v2.14 verification: the first version of this docstring attributed them to the fast components
smeared with sigma_e = 0.36 ns. sigma_e is the width of the prompt-window acceptance edge in the uncorrected pair time, not
the delay resolution. The zero-delay peak of the path-corrected delay has sigma = 0.17 ns (peak_width.py). Second
correction: the edge_excess below fits one exponential to the pooled tubes, which underpredicts the early bin of a sum of
exponentials. genuine_edge_v2.py redoes it per tube and attributes about 181 counts to para-positronium (the earlier
'< 0.01%' came from a convolution bug in peak_width.py v1).)
Test: the windowed-moment lifetime of the genuine triples G (event-identity label, known_randoms_hists.npz, rung
k = 1.00, i.e. the production rate) with the lower edge of the fit window moved from 1.0 to 2.0 ns (upper edge 5 ns
fixed, inverted with the matching truncated moment, sharp and with the a priori acceptance at the upper edge).
If the residual is leakage across the lower edge it vanishes as the edge moves out; if it is set at the upper edge or
by the o-Ps tail itself it does not. Shard bootstrap (B = 1000, seed 20261014), eight-tube mean.
Usage: genuine_window_scan.py [hists.npz] [tag]. Writes genuine_window_scan[_tag].json."""
import json, sys
from pathlib import Path
import numpy as np
from scipy.special import erfc

D = Path(__file__).resolve().parent
src = sys.argv[1] if len(sys.argv) > 1 else 'known_randoms_hists.npz'
tag = sys.argv[2] if len(sys.argv) > 2 else ''
Z = np.load(D / src)
H = Z['hist'].astype(float); KS = list(Z['ks']); RULES = list(Z['rules']); CLS = list(Z['classes'])
NAMES = list(Z['voi_names']); TCFG = Z['tau_cfg']; BINS = Z['bins']
TI = [i for i, n in enumerate(NAMES) if n.startswith('T')]
bc = 0.5 * (BINS[:-1] + BINS[1:]); ik = KS.index(1.0); iG = CLS.index('G')
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
tg = np.linspace(0.3, 12.0, 6000)
LOWS = [1.0, 1.25, 1.5, 1.75, 2.0]


def curve(lo, sig_e):
    t = np.linspace(lo, 5.0, 8001)
    A = np.ones_like(t) if sig_e is None else Phi((5 - t) / sig_e) * Phi((t + 5) / sig_e)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


CUR = {(lo, lab): curve(lo, s) for lo in LOWS for lab, s in (('sharp', None), ('acc359', 0.35897))}


def bias(G):                                         # G[..., voi+1, bin] -> {(lo, lab): [..., voi+1]} percent
    out = {}
    for lo in LOWS:
        m = (bc > lo) & (bc < 5.0); x = bc[m]; g = G[..., m]
        mean = (g * x).sum(-1) / g.sum(-1)
        for lab in ('sharp', 'acc359'):
            tau = np.interp(mean, CUR[(lo, lab)], tg)
            out[(lo, lab)] = 100 * (tau - TCFG) / TCFG
    return out


res = {}
rng = np.random.default_rng(20261014); S = H.shape[0]; B = 1000
for ir, rule in enumerate(RULES):
    G = H[:, ik, ir, :, iG, :]                      # [shard, voi+1, bin]
    b0 = bias(G.sum(0))
    bb = {k: [] for k in b0}
    for _ in range(B):
        w = np.bincount(rng.integers(0, S, S), minlength=S).astype(float)
        bx = bias(np.tensordot(w, G, axes=(0, 0)))
        for k in b0: bb[k].append(np.nanmean(bx[k][TI]))
    sd = lambda a: float(0.5 * (np.percentile(a, 84) - np.percentile(a, 16)))
    res[rule] = {f'{lo:.2f}_{lab}': dict(eight_mean=float(np.nanmean(b0[(lo, lab)][TI])), sd=sd(bb[(lo, lab)]),
                                          per_tube=[float(v) for v in b0[(lo, lab)][TI]],
                                          change_vs_1p00=float(np.nanmean(b0[(lo, lab)][TI]) - np.nanmean(b0[(1.0, lab)][TI])),
                                          change_sd=sd(np.array(bb[(lo, lab)]) - np.array(bb[(1.0, lab)])))
                 for (lo, lab) in b0}
    for lab in ('sharp', 'acc359'):
        print(rule, lab, '  '.join(f"lo {lo:.2f}: {res[rule][f'{lo:.2f}_{lab}']['eight_mean']:+.2f}+-{res[rule][f'{lo:.2f}_{lab}']['sd']:.2f} (chg {res[rule][f'{lo:.2f}_{lab}']['change_vs_1p00']:+.2f}+-{res[rule][f'{lo:.2f}_{lab}']['change_sd']:.2f})"
                                      for lo in LOWS))
# counts of G per tube below and above the edge, as the scale of the leak: fraction of G in 1.0-1.25 ns relative to an
# exponential extrapolation of 1.5-3.0 ns (eight tubes pooled, all-prompts)
GS = H[:, ik, RULES.index('all'), :, iG, :][:, TI].sum(1)          # [shard, bin]
m1 = (bc > 1.5) & (bc < 3.0); m0 = (bc > 1.0) & (bc < 1.25)
def edge(G):
    p = np.polyfit(bc[m1], np.log(G[m1]), 1); pred = np.exp(np.polyval(p, bc[m0])).sum(); return G[m0].sum(), pred, p[0]
obs, pred, sl = edge(GS.sum(0))
rb = np.random.default_rng(20261015); eb = []
for _ in range(1000):
    w = np.bincount(rb.integers(0, S, S), minlength=S).astype(float); o, q, _ = edge(w @ GS); eb.append(o / q - 1)
res['edge_excess'] = dict(observed_1p0_1p25=float(obs), extrapolated=float(pred), excess_frac=float(obs / pred - 1),
                          excess_sd=float(0.5 * (np.percentile(eb, 84) - np.percentile(eb, 16))), fit_slope_per_ns=float(sl))
print('G counts 1.00-1.25 ns: observed', int(obs), 'extrapolated from 1.5-3.0 ns', int(pred), f"excess {100 * (obs / pred - 1):+.1f}+-{100 * res['edge_excess']['excess_sd']:.1f}%")
res['doc'] = __doc__; res['source'] = src
json.dump(res, open(D / f'genuine_window_scan{("_" + tag) if tag else ""}.json', 'w'), indent=1)
