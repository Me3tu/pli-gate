#!/usr/bin/env python3
"""R2c analysis (REGISTERED_R2C.md): pairing check (P1), matched-yield arm difference against the registered
first-order values (P2), paired against independent uncertainty (P3). Readout as R2b: tail mean above 1 ns minus
1 ns on the true delays of the PLI_R2 records. Uncertainties: bootstrap over blocks of 1000 consecutive records
(records are independent annihilations), B = 2000, the same block draw for both arms (paired) or separate draws
(independent). Writes r2c_results.json."""
import gzip, json
from collections import Counter
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
TAU_MIN = 1.0
PRED = {'A': -82.06, 'B': -57.18}
TAU = {'A': {'adipose_moskal2023': 2.72, 'brain_tumor_moskal2024': 1.77},
       'B': {'adipose_avachat2024': 2.54, 'water_shibuya2020_n2': 1.9215}}
FALLBACK_TAU = 2.0


def load(arm):
    d, s, a, c = [], [], [], []
    with gzip.open(HERE / 'r2c_runs' / arm / 'pli_r2.txt.gz', 'rt') as fh:
        for ln in fh:
            p = ln.rstrip('\n').split(':')
            if len(p) < 5 or p[0] != 'PLI_R2':
                continue
            d.append(float(p[1])); s.append(p[2]); a.append(p[3]); c.append(p[4])
    return np.array(d), np.array(s, dtype=object), np.array(a, dtype=object), np.array(c, dtype=object)


def blocks(d, L=1000):
    t = d > TAU_MIN
    nb = len(d) // L
    n = np.add.reduceat(t[:nb * L].astype(float), np.arange(0, nb * L, L))
    sm = np.add.reduceat(np.where(t, d - TAU_MIN, 0.0)[:nb * L], np.arange(0, nb * L, L))
    return n, sm


out = {}
rng = np.random.default_rng(20261008)
for pair in ('A', 'B'):
    T = load(f'{pair}_transport'); S = load(f'{pair}_source')
    n = min(len(T[0]), len(S[0]))
    p1 = dict(n_transport=len(T[0]), n_source=len(S[0]),
              same_sequence=bool(len(T[0]) == len(S[0]) and all((T[k] == S[k]).all() for k in (1, 2, 3))))
    o = (T[3][:n] == 'oPs2g') & (S[0][:n] > 0)
    exp = np.array([TAU[pair].get(a, FALLBACK_TAU) / TAU[pair].get(s, FALLBACK_TAU) for a, s in zip(T[2][:n][o], S[1][:n][o])])
    rel = np.abs(T[0][:n][o] / S[0][:n][o] / exp - 1)
    p1.update(n_ops=int(o.sum()), max_rel_dev_ratio=float(rel.max()), n_rel_dev_gt_1e4=int((rel > 1e-4).sum()))
    tauT = T[0][T[0] > TAU_MIN].mean() - TAU_MIN; tauS = S[0][S[0] > TAU_MIN].mean() - TAU_MIN
    diff = 1e3 * (tauT - tauS)
    nT, sT = blocks(T[0][:n]); nS, sS = blocks(S[0][:n])
    nb = len(nT); B = 2000
    idx = rng.integers(0, nb, size=(B, nb))
    bT = sT[idx].sum(1) / nT[idx].sum(1); bS = sS[idx].sum(1) / nS[idx].sum(1)
    paired = 1e3 * (bT - bS)
    idx2 = rng.integers(0, nb, size=(B, nb))
    bS2 = sS[idx2].sum(1) / nS[idx2].sum(1)
    indep = 1e3 * (bT - bS2)
    se_p, se_i = float(paired.std(ddof=1)), float(indep.std(ddof=1))
    se_T, se_S = float(1e3 * bT.std(ddof=1)), float(1e3 * bS2.std(ddof=1))
    # first-order relation on this run's own transition matrix (equal yields: the yield cancels)
    core = Counter(S[1].tolist()).most_common(1)[0][0]
    m = S[1] == core
    tm = Counter(T[2][m].tolist()); tot = sum(tm.values())
    tau_m = {k: TAU[pair].get(k, FALLBACK_TAU) for k in tm}
    tail = {k: tm[k] / tot * np.exp(-TAU_MIN / tau_m[k]) for k in tm}
    fo = 1e3 * sum(tail[k] * (tau_m[k] - tau_m[core]) for k in tm) / sum(tail.values())
    out[pair] = dict(P1=p1, tau_transport_ns=float(tauT), tau_source_ns=float(tauS), diff_ps=float(diff),
                     se_paired_ps=se_p, se_independent_ps=se_i, se_transport_ps=se_T, se_source_ps=se_S,
                     first_order_registered_ps=PRED[pair], first_order_this_run_ps=float(fo),
                     migration_fraction=float((T[2][m] != core).mean()),
                     transition=dict((k, tm[k] / tot) for k in tm),
                     P2_met=bool(abs(diff - PRED[pair]) < 3.0), P3_met=bool(se_i / se_p >= 2.0),
                     n_tail_transport=int((T[0] > TAU_MIN).sum()), n_tail_source=int((S[0] > TAU_MIN).sum()))
    r = out[pair]
    print(f"pair {pair}: P1 {p1}")
    print(f"  tau T {tauT:.5f}  S {tauS:.5f}  diff {diff:+.2f} ps  paired se {se_p:.3f}  independent se {se_i:.3f} "
          f"(T {se_T:.3f}, S {se_S:.3f})  registered {PRED[pair]:+.2f}  first-order this run {fo:+.2f}  "
          f"migration {r['migration_fraction']:.4f}  P2 {r['P2_met']}  P3 {r['P3_met']}")
json.dump(out, open(HERE / 'r2c_results.json', 'w'), indent=1)
print('wrote r2c_results.json')
