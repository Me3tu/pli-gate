#!/usr/bin/env python3
"""Review round 14 (referee 4): uncertainty of the slope of the signed pedestal regression with the C0 tubes
(signed_eq3_c0.json rows_c0 + the twelve non-tube regions of 2026-10-08_signed_eq3/signed_eq3.json, exactly the twenty
regions of Fig. 9c and Supplementary Fig. S1). Weighted least squares through the origin, meas = slope * pred:
 (a) free slope with the paired-bootstrap sigma of each region, standard error scaled by sqrt(chi2/dof),
 (b) slope with the constant scatter term s added in quadrature (s from signed_eq3_c0.json), its standard error,
 (c) the same for the fourteen resolved regions and for the discovery-stream tubes (published regression).
Regions within one phantom share events, so these standard errors treat the regions as independent. A jackknife over
the four phantoms is given as a check of that assumption. Writes slope_ci.json."""
import json, os
from pathlib import Path
import numpy as np
R = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
C0 = json.load(open(R / '2026-10-14_known_randoms/signed_eq3_c0.json'))
PUB = json.load(open(R / '2026-10-08_signed_eq3/signed_eq3.json'))
other = [r for r in PUB['rows'] if r['phantom'] != 'Sim-1']
tubes_c0 = [dict(phantom='Sim-1', roi=r['roi'], pred=r['pred'], meas=r['meas_signed'], sig=r['sig_signed']) for r in C0['rows_c0']]
tubes_disc = [dict(phantom='Sim-1', roi=r['roi'], pred=r['pred'], meas=r['meas_signed'], sig=r['sig_signed'])
              for r in PUB['rows'] if r['phantom'] == 'Sim-1']
oth = [dict(phantom=r['phantom'], roi=r['roi'], pred=r['pred'], meas=r['meas_signed'], sig=r['sig_signed']) for r in other]


def fit(rows, s=0.0):
    x = np.array([r['pred'] for r in rows]); y = np.array([r['meas'] for r in rows]); v = np.array([r['sig'] for r in rows]) ** 2 + s ** 2
    w = 1 / v; b = (w * x * y).sum() / (w * x * x).sum(); se = 1 / np.sqrt((w * x * x).sum())
    chi2 = (w * (y - b * x) ** 2).sum(); dof = len(x) - 1
    return b, se, chi2 / dof


out = {}
for lab, rows, s in (('c0_all20', tubes_c0 + oth, C0['signed']['all20']['scatter_pp']),
                     ('discovery_all20', tubes_disc + oth, PUB['signed']['all20']['scatter_pp'])):
    b, se, cd = fit(rows); bs, ses, cds = fit(rows, s)
    ph = sorted(set(r['phantom'] for r in rows)); jk = np.array([fit([r for r in rows if r['phantom'] != p], s)[0] for p in ph])
    sej = float(np.sqrt((len(ph) - 1) / len(ph) * ((jk - jk.mean()) ** 2).sum()))
    out[lab] = dict(n=len(rows), free_slope=b, free_se_scaled=se * np.sqrt(cd), chi2_dof_free=cd, scatter_pp=s,
                    slope_with_scatter=bs, se_with_scatter=ses, chi2_dof_with_scatter=cds,
                    phantom_jackknife_se_with_scatter=sej, jackknife_slopes=dict(zip(ph, jk.tolist())))
    print(lab, f"free {b:.3f} +- {se*np.sqrt(cd):.3f} (chi2/dof {cd:.2f}); with scatter {s:.3f}: {bs:.3f} +- {ses:.3f} "
               f"(chi2/dof {cds:.2f}); phantom jackknife se {sej:.3f} {np.round(jk, 3)}")
json.dump(out, open(R / '2026-10-14_known_randoms/slope_ci.json', 'w'), indent=1)
