#!/usr/bin/env python3
"""Round 6 (2026-10-20): first-arrival model against the measured closest-prompt bias, with the eight tubes of the
corrected-yield arm (C0) in place of the earlier-yield stream, so that the fourteen resolved regions are those of the
pedestal regression of the main text. Same model, estimator and rate definition as test_derivation.py of
2026-08-17_derivation_test:
    predicted relative bias = -r tau / (1 + r tau),  r = prompt-window singles per ns of the same acquisition,
    measured = 100 (tau_hat_trunc - tau) / tau, closest-prompt extraction, sigma_pct = 100 sigma_trunc / tau.
The disk, water-cylinder and NEMA points are copied from derivation_test.json unchanged.
Writes refutation_c0.json (schema of derivation_test.json)."""
import csv, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent
old = json.load(open(RES / '2026-08-17_derivation_test/derivation_test.json'))['points']
mu = json.load(open(RES / '2026-09-30_feedback3/c0_sim1/c0_compare.json'))['summary']['mu']['c0']
r = mu['n_prompt'] / mu['span_ns']
pts = []
for row in csv.DictReader(open(RES / '2026-09-30_feedback3/c0_sim1/out_closest/per_voi.csv')):
    if row['voi'] == 'bg' or not row['tau_cfg_ns']:
        continue
    tau = float(row['tau_cfg_ns']); th = float(row['clinical_tau_trunc_ns']); sg = float(row['clinical_sigma_trunc_ns'])
    pred = -100.0 * r * tau / (1.0 + r * tau); meas = 100.0 * (th - tau) / tau
    pts.append(dict(phantom='eight-tube C0', roi=row['voi'], tau=tau, tau_hat=th, sigma=sg, measured=meas,
                    n=float(row['clinical_n_signal']), r_per_ns=r, r_tau=r * tau, predicted=pred,
                    ratio=meas / pred, sigma_pct=100.0 * sg / tau))
pts += [d for d in old if d['phantom'] != 'Sim-1']
res = [d for d in pts if not d['roi'].startswith('sphere')]
rat = np.array([d['ratio'] for d in res])
out = dict(points=pts, n=len(pts), n_resolved=len(res), r_per_ns_c0=r,
           ratio_min_resolved=float(rat.min()), ratio_max_resolved=float(rat.max()),
           ratio_median_resolved=float(np.median(rat)),
           ratio_c0_tubes=[float(d['ratio']) for d in pts if d['phantom'] == 'eight-tube C0'])
json.dump(out, open(HERE / 'refutation_c0.json', 'w'), indent=1)
for d in pts:
    print(f"{d['phantom']:>18} {d['roi']:>28} tau {d['tau']:.4f} pred {d['predicted']:7.3f} meas {d['measured']:7.2f} "
          f"+- {d['sigma_pct']:5.2f} ratio {d['ratio']:6.1f}")
print('resolved', len(res), 'ratio', round(rat.min(), 2), 'to', round(rat.max(), 2), 'median', round(float(np.median(rat)), 2))
old_res = [d for d in old if not d['roi'].startswith('sphere')]
print('earlier-yield set: ratio', round(min(d['ratio'] for d in old_res), 2), 'to', round(max(d['ratio'] for d in old_res), 2))
