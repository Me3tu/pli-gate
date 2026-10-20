"""Fig. 7 (stem fig5_sensitivity), feedback round 3.
(a) expected signal: configured lifetime shift relative to 2 mmHg, in ps, with the full-range
    span and the ratio of the single-insert standard deviation to it.
(b) measurement: recovered minus configured lifetime per insert with 95% percentile intervals
    of an event-level bootstrap (tail events resampled within each insert, B = 20000), zero line,
    and the mean paired difference with its event-level bootstrap 95% interval (band).
Review round 8: the phantom was rerun with database v8 (plasma-based coefficient 2.5704e-5 ns^-1 mmHg^-1)
and the production binary (corrected I3 convention); per-insert values from the v8 triples
(2026-10-08_oxygen_v8/oxygen_boot_v8.json, 'reextraction'). Panel (a): solid, as simulated (plasma
coefficient); dashed, the 37 C water solubility (1.40 umol/L/mmHg) as sensitivity value, radius anchored at
the 15 mmHg reference as in the database."""
import json
import numpy as np
from common import *  # noqa

O = json.load(open(RES / '2026-10-08_oxygen_v8/oxygen_boot_v8.json'))
X = O['reextraction']; P = X['per_insert']
po2 = np.array([r['pO2'] for r in P]); cfg = np.array([r['tau_cfg'] for r in P])
d = np.array([r['diff_ps'] for r in P]); lo = np.array([r['boot95_ps'][0] for r in P]); hi = np.array([r['boot95_ps'][1] for r in P])
KAPPA = 0.0204e-3 * 1.26  # ns^-1 mmHg^-1, Eq. (2), database v8 (plasma solubility)
span = 1000 * (cfg[0] - cfg[-1]); ratio = X['mean_sigma_ps'] / span
print(f'span {span:.2f} ps, mean sigma {X["mean_sigma_ps"]:.1f} ps, ratio {ratio:.1f}; mean diff {X["mean_diff_ps"]:.1f} '
      f'95% {X["mean_diff_boot95_ps"]}; half-widths {X["halfwidth_range_ps"]}')

fig = mm_fig(150, 64)
ax = ax_mm(fig, 16, 13, 52, 44)
lam0 = 1 / cfg[0] - KAPPA * po2[0]
xx = np.linspace(0, 42, 200)
ax.plot(xx, 1000 * (1 / (lam0 + KAPPA * xx) - cfg[0]), color=REF, lw=1.2)
ax.plot(po2, 1000 * (cfg - cfg[0]), 'o', ms=4.6, color=REF, zorder=4)
# Sensitivity: the 37 C solubility of water (1.40 umol/L/mmHg); the radius is anchored at the 15 mmHg reference,
# so tau(15) is unchanged.
KAPPA_W = 0.0204e-3 * 1.406  # Battino et al at 37 C (review round 16: 1.406, giving 2.87e-5)
tau15 = 1 / (lam0 + KAPPA * 15.0)
lam_w = 1 / tau15 - KAPPA_W * 15.0
yw = 1000 * (1 / (lam_w + KAPPA_W * xx) - 1 / (lam_w + KAPPA_W * 2.0))
span_w = 1000 * (1 / (lam_w + KAPPA_W * 2.0) - 1 / (lam_w + KAPPA_W * 40.0))
print(f'water coefficient {KAPPA_W:.4e}: span {span_w:.2f} ps')
ax.plot(xx, yw, color=REF, lw=1.0, ls='--', zorder=3)
ax.axhline(0, color='#aaaaaa', lw=0.6)
ax.set_xlim(0, 42); ax.set_ylim(-5.0, 1.0)  # round 10: room below the dashed line for the key
ax.set_xlabel(r'oxygen tension $p\mathrm{O}_2$ (mmHg)')
ax.set_ylabel(r'configured shift $\Delta\tau_{\rm cfg}$ (ps)')
sd_lo, sd_hi = X['mean_bootsd_ps'], X['mean_sigma_ps']  # bootstrap s.d. and tau/sqrt(N), as in the text
sdtxt = (f'{sd_lo:.0f} ps ({sd_lo / span:.0f}$\\times$)' if (round(sd_lo) == round(sd_hi) and round(sd_lo / span) == round(sd_hi / span))
         else f'{sd_lo:.0f}\u2013{sd_hi:.0f} ps ({sd_lo / span:.0f}\u2013{sd_hi / span:.0f}$\\times$)')
key(ax, 0.04, 0.03, f'full range: {span:.2f} ps as simulated,\n{span_w:.1f} ps with water solubility (dashed)\nsingle-insert s.d. {sdtxt}',
    ha='left', va='bottom', fontsize=8.0)
st.panel_letter(ax, 'a')

ax = ax_mm(fig, 93, 13, 55, 44)
m, (mlo, mhi) = X['mean_diff_ps'], X['mean_diff_boot95_ps']
ax.axhline(0, color=REF, lw=0.9, zorder=2)
ax.axhspan(mlo, mhi, color=THIS, alpha=0.10, lw=0, zorder=1)
ax.axhline(m, color=THIS, lw=0.9, ls='--', zorder=2)
ax.errorbar(po2, d, yerr=[d - lo, hi - d], fmt='o', ms=4.8, color=THIS, elinewidth=1.0, capsize=2, lw=0, zorder=4)
ax.set_xlim(0, 42); ax.set_ylim(min(-420, 40 * np.floor(lo.min() / 40) - 40), max(340, 40 * np.ceil(hi.max() / 40) + 260))
ax.set_xlabel(r'oxygen tension $p\mathrm{O}_2$ (mmHg)')
ax.set_ylabel(r'$\hat\tau-\tau_{\rm cfg}$ (ps)')
key(ax, 0.03, 0.97, f'mean paired difference {sg(m, "+.0f")} ps\n95% interval [{sg(mlo, ".0f")}, {sg(mhi, "+.0f")}] ps',
    ha='left', va='top', fontsize=8.0, color=THIS)
st.panel_letter(ax, 'b')
save(fig, 'fig5_sensitivity')
