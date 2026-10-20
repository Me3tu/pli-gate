"""Supplementary Fig. S4 (stem figS3_eq2_refutation): (a) Sim-1 tubes at common
accidental rate, measured |bias| vs configured tau, with the model's monotonic
shape r tau/(1+r tau) scaled to the mean measured bias (dashed). The 1.5 ps pair
ringed. (b) all regions against the model with the identity and x10, x100 guides.
Data: 2026-08-17_derivation_test/derivation_test.json."""
import json
import numpy as np
from common import *  # noqa

pts = json.load(open(RES / '2026-08-17_derivation_test/derivation_test.json'))['points']
sim1 = [d for d in pts if d['phantom'] == 'Sim-1']
PH = [('Sim-1', 'o', THIS, 'eight-tube phantom (earlier-yield stream)'), ('Steinberger disks', 's', THIS, 'reference disks'),
      ('Steinberger water', '^', THIS, 'water cylinder'), ('NEMA IQ', 'D', CTRL, 'NEMA background')]
# Review round 16: panel (a), the shape test at n = 8, had no power and is dropped. The model is rejected by the
# scale mismatch drawn here. The eight-tube points are the earlier-yield stream (Supplementary Note S8).
fig = mm_fig(110, 92)
axes = [ax_mm(fig, 16, 30, 88, 56)]
r = None
ax = axes[0]
lim = np.array([0.05, 60])
for f, ls in ((1, '--'), (10, ':'), (100, ':')):
    ax.plot(lim, f * lim, color=REF if f == 1 else CTRL, lw=0.9, ls=ls, zorder=1)
    if f > 1: ax.text(52.0 / f, 56, f'$\\times${f}', fontsize=8.0, color=CTRL, ha='left', va='top')
for ph, mk, col, lab in PH:
    sel = [d for d in pts if d['phantom'] == ph and not d['roi'].startswith('sphere')]
    if not sel: continue
    ax.errorbar([abs(d['predicted']) for d in sel], [abs(d['measured']) for d in sel], yerr=[d['sigma_pct'] for d in sel],
                fmt=mk, ms=4.2, color=col, ecolor=col, elinewidth=0.8, capsize=1.8, lw=0, zorder=3, label=lab)
sph = [d for d in pts if d['phantom'] == 'NEMA IQ' and d['roi'].startswith('sphere')]
ax.errorbar([abs(d['predicted']) for d in sph], [abs(d['measured']) for d in sph], yerr=[d['sigma_pct'] for d in sph],
            fmt='D', ms=4.2, color=CTRL, mfc='white', mec=CTRL, ecolor=CTRL, elinewidth=0.8, capsize=1.8, lw=0, zorder=3,
            label='NEMA spheres (unresolved)')
ax.plot([], [], color=REF, ls='--', lw=0.9, label='identity')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.08, 60); ax.set_ylim(0.3, 60)
ax.set_xlabel(r'predicted $|$bias$|$ $= r\tau/(1+r\tau)$ (%)'); ax.set_ylabel(r'measured $|$bias$|$ (%)')
hh, ll = ax.get_legend_handles_labels()
ORD = ['identity', 'eight-tube phantom (earlier-yield stream)', 'reference disks', 'water cylinder', 'NEMA background', 'NEMA spheres (unresolved)']
hh = [hh[ll.index(o)] for o in ORD]; ll = ORD
fig.legend(hh, ll, frameon=False, loc='lower right', bbox_to_anchor=(108 / 110, 0.6 / 92), ncol=2, fontsize=8.0,
           handlelength=1.4, columnspacing=1.0, borderaxespad=0.0)
save(fig, 'figS3_eq2_refutation')
