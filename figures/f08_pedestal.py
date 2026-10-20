"""Fig. 8 (stem fig9_pedestal_derivation), three panels, and Supplementary
Fig. S1 (stem figS2_pedestal_residuals, the former panel d).
(a) larger: hepatic tube spectra under both rules, shaded excess, fit window;
(b) difference between the rules: thin grey tubes, bold mean, a-priori
    acceptance roll-off; (c) measured vs predicted change in |bias| with the
    zero-parameter consistency line (conditional on the measured pedestal) and the scatter band.
Review round 6 (2026-10-08): panel (c) and Fig. S1 plot the SIGNED rule difference
b_all - b_closest, the quantity Eq. (4) predicts, with its paired-bootstrap uncertainty and the
scatter term of the signed regression (2026-10-08_signed_eq3/signed_eq3.json); the folded
|b_closest| - |b_all| differs from it only where a bias is positive (hepatic tube, two spheres).
Data: 2026-08-22_G3_consequences/hists_sim1.npz; signed_eq3.json."""
# Review round 9 (2026-10-09): the pedestal relation is Eq. (4) of the main text since the branching equation became Eq. (3).
import json
import numpy as np
from scipy.stats import norm
from common import *  # noqa

# Review round 13 (2026-10-14): eight-tube data from the corrected-yield arm C0 (hists_c0_0ps.npz, same dump script
# as hists_sim1.npz); panel (c) and Fig. S2 use the twenty-region signed regression with the C0 tubes in place of the
# discovery tubes (2026-10-14_known_randoms/signed_eq3_c0.json); the other twelve regions are unchanged.
Z = np.load(F3 / 'c0_sim1/hists_c0_0ps.npz', allow_pickle=False)  # v1.1.0: was allow_pickle=True; every array is numeric or <U27
names = [str(s) for s in Z['voi_names']]
tubes = [n for n in names if n.startswith('T')]
b = Z['tau_bins']; c = 0.5 * (b[:-1] + b[1:])
SG0 = json.load(open(RES / '2026-10-08_signed_eq3/signed_eq3.json'))
SG = json.load(open(RES / '2026-10-14_known_randoms/signed_eq3_c0.json'))
P = [dict(roi=r['roi'], phantom='Sim-1', pred=r['pred'], meas=r['meas_signed'], sigma_paired=r['sig_signed']) for r in SG['rows_c0']] + \
    [dict(roi=r['roi'], phantom=r['phantom'], pred=r['pred'], meas=r['meas_signed'], sigma_paired=r['sig_signed'])
     for r in SG0['rows'] if r['phantom'] != 'Sim-1']
SS = SG['signed']['all20']['scatter_pp']
print('signed regression: scatter', round(SS, 3), 'pp; slope', round(SG['signed']['all20']['slope'], 3),
      'chi2/region', round(SG['signed']['all20']['chi2_per_region'], 2))
PH_ORDER = ['Sim-1', 'Steinberger disks', 'Steinberger water', 'NEMA IQ']
PH_LABEL = {'Sim-1': 'eight tubes', 'Steinberger disks': 'reference disks',
            'Steinberger water': 'water cylinder, 20 mm', 'NEMA IQ': 'NEMA IQ'}
MRK = {'Sim-1': 'o', 'Steinberger disks': 's', 'Steinberger water': '^', 'NEMA IQ': 'D'}
COLP = {'Sim-1': THIS, 'Steinberger disks': THIS, 'Steinberger water': THIS, 'NEMA IQ': CTRL}


def arms(ki, shards=None):
    sl = slice(None) if shards is None else shards
    cs = Z['closest_sig'][sl, ki, :].sum(0); cr = Z['closest_rnd'][sl, ki, :].sum(0)
    as_ = Z['all_sig'][sl, ki, :].sum(0); ar = Z['all_rnd'][sl, ki, :].sum(0)
    return as_ - ar, cs - cr


idx = [names.index(t) for t in tubes]
nsh = Z['closest_sig'].shape[0]


def pooled(shards=None):
    d = np.zeros(len(c))
    for ki in idx:
        sa, sc = arms(ki, shards); d += sa - sc
    return d


d0 = pooled()
rng = np.random.default_rng(20260926)
boots = np.array([pooled(rng.integers(0, nsh, nsh)) for _ in range(300)])
lo16, hi84 = np.percentile(boots, [16, 84], axis=0); err = 0.5 * (hi84 - lo16)

fig = mm_fig(150, 126)
ax = ax_mm(fig, 18, 80, 130, 39)
ki = names.index('T4_hepatic_avachat2024'); sa, sc = arms(ki)
m = (c > -0.6) & (c < 5.6)
ax.axvspan(1.0, 5.0, color=st.BAND, alpha=0.05, lw=0, zorder=0)
ax.fill_between(c[m], np.clip(sc[m], 1, None), np.clip(sa[m], 1, None), where=(c[m] > 1.0) & (c[m] < 5.0),
                step='mid', facecolor='none', edgecolor=THIS, hatch='//////', lw=0, alpha=0.55, zorder=2)
ax.step(c[m], np.clip(sa[m], 1, None), where='mid', color=THIS, lw=1.2, zorder=3, label='all-prompts matching')
ax.step(c[m], np.clip(sc[m], 1, None), where='mid', color=THIS, lw=1.0, ls=(0, (3, 1.5)), zorder=3, label='closest-prompt matching')
ax.fill_between([], [], [], facecolor='none', edgecolor=THIS, hatch='//////', lw=0, alpha=0.55, label='excess subtracted under closest-prompt matching')
ax.set_yscale('log'); ax.set_xlim(-0.6, 5.6); ax.set_ylim(20, 8e4)
ax.set_xlabel('delay (ns)'); ax.set_ylabel('randoms-subtracted counts per 50 ps')
ax.text(1.05, 6.0e4, 'fit window', fontsize=8.0, ha='left', va='center', color='#666666')
ax.legend(frameon=False, fontsize=8.0, loc='upper right', bbox_to_anchor=(1.0, 0.98))
st.panel_letter(ax, 'a')
# Review round 17 (second report, Sec. 8): linear inset of the fit window, where the log scale hides the excess.
axi = ax_mm(fig, 51, 100.0, 33, 14.5)
mi = (c > 1.0) & (c < 5.0)
axi.fill_between(c[mi], sc[mi], sa[mi], step='mid', facecolor='none', edgecolor=THIS, hatch='//////', lw=0, alpha=0.55)
axi.step(c[mi], sa[mi], where='mid', color=THIS, lw=1.0)
axi.step(c[mi], sc[mi], where='mid', color=THIS, lw=0.9, ls=(0, (3, 1.5)))
top = 1.08 * sa[mi].max()
axi.set_xlim(1.0, 5.0); axi.set_ylim(0, top)
axi.set_xticks([1, 2, 3, 4, 5]); axi.tick_params(labelsize=7.0, length=2, pad=1.5)
axi.text(0.97, 0.93, 'linear scale', transform=axi.transAxes, ha='right', va='top', fontsize=7.0, color='#444444')
axi.set_yticks([0, 300, 600])
print('inset top', round(top, 0), 'counts at 1 ns', round(float(sa[(c > 1.0) & (c < 1.05)][0]), 0))

ax = ax_mm(fig, 23, 12, 50, 49)
ax.axvspan(1.0, 5.0, color=st.BAND, alpha=0.05, lw=0, zorder=0)
for t, ki in zip(tubes, idx):
    sa, sc = arms(ki); dd = sa - sc
    k = 10; n = (len(c) // k) * k
    cr_ = c[:n].reshape(-1, k).mean(1); dr = dd[:n].reshape(-1, k).mean(1)
    mm = (cr_ > 0.75) & (cr_ < 5.25)
    ax.plot(cr_[mm], dr[mm], '-', color=CTRL, lw=0.6, alpha=0.7, zorder=2)
mm = (c > 0.5) & (c < 5.5)
ax.errorbar(c[mm], d0[mm] / 8, yerr=err[mm] / 8, fmt='o', ms=2.4, color=THIS, ecolor=THIS, elinewidth=0.6, capsize=0,
            lw=0, zorder=4, label='mean of eight tubes, 50 ps bins')
ax.plot([], [], '-', color=CTRL, lw=0.6, label='individual tubes, 500 ps bins')
mp = (c > 1.25) & (c < 4.5); lvl = d0[mp].mean() / 8; SIG_E = 0.359
tt = np.linspace(0.5, 5.5, 400)
ax.plot(tt, lvl * norm.cdf((5.0 - tt) / SIG_E), color=REF, lw=1.1, ls='--', zorder=5, label='plateau $\\times$ window acceptance')
ax.axhline(0, color='#999999', lw=0.6, zorder=1)
ax.set_xlim(0.5, 5.5)
ax.set_xlabel('delay (ns)'); ax.set_ylabel('difference between rules\n(counts per 50 ps, per tube)')
ax.legend(frameon=False, fontsize=8.0, loc='upper left', borderaxespad=0.1)
ax.set_ylim(-20, 115)
st.panel_letter(ax, 'b')
print('plateau per tube', round(lvl, 1))

ax = ax_mm(fig, 97, 12, 51, 49)
RESOLVED = [p for p in P if not p['roi'].startswith('sphere_')]
lim = (0, 27); xx = np.array(lim)
ax.fill_between(xx, xx - SS, xx + SS, facecolor=st.BAND, alpha=0.10, lw=0, zorder=1)
ax.plot(lim, lim, color=REF, lw=1.0, ls='--', zorder=2)
ax.text(12.0, 15.4, r'Eq. (5) at measured $\phi$', rotation=45, rotation_mode='anchor', fontsize=8.0, ha='center', va='center', color=KEYC)
for ph in PH_ORDER:
    sel = [p for p in RESOLVED if p['phantom'] == ph]
    if not sel: continue
    ax.errorbar([p['pred'] for p in sel], [p['meas'] for p in sel], yerr=[p['sigma_paired'] for p in sel],
                fmt=MRK[ph], ms=4.2, color=COLP[ph], ecolor=COLP[ph], elinewidth=0.8, capsize=1.6, lw=0, zorder=4,
                label=PH_LABEL[ph])
ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect('equal', adjustable='box')
ax.set_xticks([0, 5, 10, 15, 20, 25]); ax.set_yticks([0, 5, 10, 15, 20, 25])
ax.set_xlabel(r'$\Delta b$ from Eq. (5) (pp)'); ax.set_ylabel(r'measured $\Delta b$ (pp)')
ax.legend(fontsize=8.0, loc='upper left', frameon=False, handletextpad=0.3, borderaxespad=0.1)
key(ax, 0.97, 0.03, '14 resolved of\n20 regions shown,\nregions within a\nphantom share events', ha='right', va='bottom', fontsize=8.0, color='#555555')  # round 10: clear of the low points
st.panel_letter(ax, 'c')
save(fig, 'fig9_pedestal_derivation')

# ---- Supplementary Fig. S2: normalized residuals, all twenty regions ----------
fig, ax = plt.subplots(figsize=(110 * MM, 112 * MM))
order = []
for ph in PH_ORDER: order += [p for p in P if p['phantom'] == ph]
yy = np.arange(len(order))[::-1]
ax.axvspan(-2, 2, color=st.BAND, alpha=0.05, lw=0, zorder=1)
ax.axvline(0, color=REF, lw=1.0, zorder=2)
for y, p in zip(yy, order):
    r = p['meas'] - p['pred']; s = np.hypot(p['sigma_paired'], SS)
    ax.errorbar([r / s], [y], xerr=[[1.0]], fmt=MRK[p['phantom']], ms=4.0, color=COLP[p['phantom']],
                ecolor=COLP[p['phantom']], elinewidth=0.8, capsize=1.6, zorder=4)
ax.set_yticks(yy); ax.set_yticklabels([NAME[p['roi']] for p in order], fontsize=8.0)
ax.tick_params(axis='y', length=0); ax.set_ylim(-0.8, len(order) - 0.2)
ax.set_xlim(-3, 3); ax.set_xticks([-2, -1, 0, 1, 2])
ax.set_xlabel(r'(measured $\Delta b$ $-$ Eq. (5)) / $\sigma$')
ax.spines['left'].set_visible(False)
n = 0
for ph in PH_ORDER[:-1]:
    n += sum(1 for p in order if p['phantom'] == ph)
    ax.axhline(len(order) - n - 0.5, color='#bbbbbb', lw=0.5, zorder=1)
fig.tight_layout()
save(fig, 'figS2_pedestal_residuals')
