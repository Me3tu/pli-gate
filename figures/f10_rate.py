"""Fig. 10 (stem fig10_rate_law): (a) rate replays of two event populations, twelve rungs each:
the discovery stream (filled, published power-law line, joint-bootstrap exponent band pivoted at the
rung centroid, per-rung shard-bootstrap errors) and, added in review round 5, the corrected-yield
zero-delay stream C0 (open, dashed OLS line, joint shard-bootstrap errors), same six tubes, same
method (2026-10-08_c0_ladder). (b) per-region rule difference per 0.01 in mu for every resolved region of the separately simulated
configurations (audit round 7b, v2.6 and earlier showed configuration means over the discovery-replay law).
Review round 6 (2026-10-08): the shaded bands are pointwise 95% intervals of the fitted power law from
the joint shard bootstrap with prefactor and exponent resampled together (2026-10-08_signed_eq3/
ladder_band.json, same resamples as the exponent intervals), one per stream, replacing the exponent
interval pivoted at the centroid. The discovery band belongs to the per-shard re-extraction, whose fit
(416 mu^0.997) lies within 1% of the drawn registered line (417 mu^0.999) over the plotted range.
Review round 13 (2026-10-14): the corrected-yield replay is drawn as the primary population (filled, solid) and the
discovery replay as the replication (open, dashed). Panels (c) and (d) add the known-content test of the randoms
correction on the corrected-yield replay (2026-10-14_known_randoms/known_randoms_summary.json, a priori acceptance
readout): (c) randoms-correction error tau(sig - rnd) - tau(G + PA_corr), (d) readout of the genuine triples alone."""
import json
import numpy as np
from common import *  # noqa
from matplotlib.lines import Line2D

d = json.load(open(MF / 'fig10_rate_law.json'))
L = d['ladder']; rows = L['rows']
# Code release v1.1.0: 2026-09-30_referee2/ladder_joint/ladder_boot.json was loaded here into (lo, hi) and never used,
# load removed.
mu = np.array([r['mu'] for r in rows]); did = np.array([r['did'] for r in rows])
C0, N0 = np.exp(L['ln_intercept']), L['n_exponent']
FH = 124
fig = mm_fig(150, FH)
ax = ax_mm(fig, 14, 70, 60, 47)
xx = np.logspace(np.log10(0.0035), np.log10(0.11), 60)
LB = json.load(open(RES / '2026-10-08_signed_eq3/ladder_band.json'))
gx = np.array(LB['mu_grid'])
# Review round 16: corrected-yield stream only (the discovery stream is reported in Supplementary Note S15).
for nm, al in (('corrected', 0.22),):
    bd = LB['streams'][nm]['folded']['band']
    assert np.allclose(LB['streams'][nm]['folded']['ci95'], LB['streams'][nm]['folded']['reference_ci95'])
    ax.fill_between(gx, bd['lo'], bd['hi'], color=st.BAND, alpha=al, lw=0, zorder=1)
rs = LB['streams']['discovery']['folded']
dev = np.abs(rs['prefactor'] * xx ** rs['exponent'] / (C0 * xx ** N0) - 1).max()
print(f'discovery re-extraction fit {rs["prefactor"]:.1f} mu^{rs["exponent"]:.4f}; max deviation from drawn line {100*dev:.2f}%')
# Review round 8: line and points both from the per-shard re-extraction (the analysis estimate). The original
# extraction fit (417 mu^0.999) lies within 1% of this line over the plotted range (printed above).
# Review round 5: the same replay on the corrected-yield stream C0 (six common tubes)
CL = RES / '2026-10-08_c0_ladder/ladder_boot_c0.json'
Lc = json.load(open(CL)); Sc = Lc['sets']['six_common']
muc = np.array(Lc['mu']); didc = np.array(Sc['mean_did']); sdc = np.array(Sc['mean_did_sd'])
ax.plot(xx, Sc['prefactor_ols'] * xx ** Sc['exponent_ols'], color=REF, lw=1.1, zorder=2)
ax.errorbar(muc, didc, yerr=sdc, fmt='o', color=THIS, ms=5, capsize=0, lw=0, elinewidth=0.9, zorder=4)
lc, hc = Sc['boot_ci95']
print(f"C0 replay: exponent {Sc['exponent_ols']:.3f} [{lc:.3f}, {hc:.3f}], GLS {Sc['gls_exponent']:.3f}+/-{Sc['gls_se']:.3f}, prefactor {Sc['prefactor_ols']:.0f}")
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.0035, 0.11); ax.set_ylim(1.2, 50)
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
ax.yaxis.set_major_locator(FixedLocator([2, 5, 10, 20, 40])); ax.yaxis.set_major_formatter(FixedFormatter(['2', '5', '10', '20', '40']))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_xlabel(r'mean accidental prompts per $\pm$5 ns window, $\mu$', fontsize=8.5)
ax.set_ylabel(r'selection shift $\Delta|b|$ (pp)')
key(ax, 0.04, 0.96, 'six tubes, zero delay,\ntwelve rate replays', ha='left', va='top', fontsize=8.6)
key(ax, 0.97, 0.04, f'exponent {Sc["exponent_ols"]:.3f}\n95% interval [{lc:.3f}, {hc:.3f}]', ha='right', va='bottom', fontsize=8.0)
st.panel_letter(ax, 'a')

ax = ax_mm(fig, 107, 70, 41, 47)
# Audit round 7b (2026-10-08): panel (b) shows the artifact per unit mu region by region, because a
# configuration-level ratio to the replay law depends on how regions are averaged (inverse-variance against
# plain means) and on the tube set of the law. Signed rule difference per region (signed_eq3.json for disks,
# water, NEMA, and c0_compare.json for the eight tubes at corrected yield and 388 ps), divided by the
# configuration's own mu, in percentage points per 0.01 accidental prompts per window. Resolved regions only.
SE = json.load(open(RES / '2026-10-08_signed_eq3/signed_eq3.json'))['rows']
CC = json.load(open(F3 / 'c0_sim1/c0_compare.json'))
FR = json.load(open(MF / 'fig11_realistic_forest.json'))
MU = {c['phantom'] + '|' + c['arm']: c['mu'] for c in d['configs']}
mu_c0, mu_disc, mu_c2 = (CC['summary']['mu'][k]['mu'] for k in ('c0', 'june0', 'c2'))
groups = []   # (label, list of (value, err), style)
groups.append(('eight tubes', [(r['c0']['did'] / mu_c0 / 100, r['c0']['sd_did'] / mu_c0 / 100) for r in CC['rows']], 'c0'))
groups.append(('disks', [(r['meas_signed'] / MU['Steinberger disks|0 ps'] / 100, r['sig_signed'] / MU['Steinberger disks|0 ps'] / 100)
                         for r in SE if r['phantom'] == 'Steinberger disks'], 'own'))
groups.append(('water cylinder, 20 mm', [(r['meas_signed'] / MU['Steinberger water|0 ps'] / 100, r['sig_signed'] / MU['Steinberger water|0 ps'] / 100)
                              for r in SE if r['phantom'] == 'Steinberger water'], 'own'))
groups.append(('NEMA background', [(r['meas_signed'] / MU['NEMA IQ|0 ps'] / 100, r['sig_signed'] / MU['NEMA IQ|0 ps'] / 100)
                           for r in SE if r['phantom'] == 'NEMA IQ' and r['roi'] == 'background_central'], 'own'))
groups.append(('eight tubes, 388 ps', [(r['c2']['did'] / mu_c2 / 100, r['c2']['sd_did'] / mu_c2 / 100) for r in CC['rows']], 'r388'))
nb = [r for r in FR['rows'] if r['phantom'] == 'NEMA IQ' and r['region'].startswith('background')][0]
groups.append(('NEMA bg, 388 ps', [(nb['did388'] / FR['mu_matched']['NEMA IQ 388 ps'] / 100, nb['sig388'] / FR['mu_matched']['NEMA IQ 388 ps'] / 100)], 'r388'))
law = Sc['prefactor_ols'] * 0.0292 ** Sc['exponent_ols'] / 2.92
# Referee read v2.7: horizontal dot plot with horizontal category labels (rotated labels were unreadable).
ax.axvline(law, color=REF, lw=0.9, ls='--', zorder=1)
n0 = 4
ypos = [len(groups) - 1 - i + (0.6 if i < n0 else 0.0) for i in range(len(groups))]
for i, (lab_, vals, kind) in enumerate(groups):
    v = np.array([x for x, _ in vals]); e = np.array([y for _, y in vals])
    ys = ypos[i] + (np.linspace(0.22, -0.22, len(v)) if len(v) > 1 else np.zeros(1))
    mk = '^' if kind == 'r388' else ('o' if kind == 'disc' else 's')
    kw = dict(color=THIS, mfc='white', mec=THIS) if kind == 'disc' else dict(color=THIS, mfc=THIS, mec=THIS)
    ax.errorbar(v, ys, xerr=e, fmt=mk, ms=3.8, elinewidth=0.8, capsize=0, zorder=4, **kw)
    print(lab_, kind, 'per 0.01 mu:', np.round(v, 2), '+/-', np.round(e, 2))
ax.set_yticks(ypos); ax.set_yticklabels([g[0] for g in groups], fontsize=8.0)
ax.tick_params(axis='y', length=0)
ax.set_ylim(-0.6, ypos[0] + 0.9); ax.set_xlim(0, 9.6)
ydiv = (ypos[n0 - 1] + ypos[n0]) / 2
ax.axhline(ydiv, color='#cccccc', lw=0.6, zorder=0)
ax.text(9.4, ypos[0] + 0.55, 'direct at 0 ps', ha='right', va='center', fontsize=8.0, color='#444444')
ax.text(9.4, ydiv - 0.32, 'direct at 388 ps', ha='right', va='center', fontsize=8.0, color='#444444')
ax.set_xlabel(r'$\Delta b$ per 0.01 in $\mu$ (pp)', fontsize=8.5)
key(ax, 4.5, 0.45, 'dashed: replay law,\nsix tubes', transform=ax.transData, ha='left', va='center', fontsize=8.0, color='#555555')
ax.set_xticks([0, 2, 4, 6, 8])
st.panel_letter(ax, 'b')
print('law per 0.01 at mu=0.0292:', round(law, 2))

# ---- (c), (d): known-content test of the randoms correction (review round 13) ----------------------------
KR = json.load(open(RES / '2026-10-14_known_randoms/known_randoms_summary.json'))
mk = np.array(KR['mu'])
def series(rule, key_):
    v = np.array([e[f'{key_}_eight'][0] for e in KR['per_rung'][rule]]); e_ = np.array([e[f'{key_}_eight'][1] for e in KR['per_rung'][rule]])
    return v, e_
for j, (key_, lab, ylim, panel) in enumerate([('dR_acc359', 'randoms-correction error (pp)', (-38, 5), 'c'),
                                             ('bG_acc359', 'genuine-triple residual (%)', (-6.0, 1.5), 'd')]):
    ax = ax_mm(fig, 14 + j * 74, 14, 60, 38)
    ax.axhline(0, color=REF, lw=0.8, ls='--', zorder=1)
    vc, ec = series('closest', key_); va, ea = series('all', key_)
    if panel == 'd':   # residuals against configured lifetimes carry the run-to-run factor 1.30 (v2.19)
        ec, ea = 1.30 * ec, 1.30 * ea
    ax.errorbar(mk, vc, yerr=ec, fmt='o', color=THIS, mfc='white', mec=THIS, ms=4.2, capsize=0, lw=0.8, elinewidth=0.8, zorder=4, label='closest-prompt')
    ax.errorbar(mk, va, yerr=ea, fmt='o', color=THIS, ms=4.2, capsize=0, lw=0.8, elinewidth=0.8, zorder=5, label='all-prompts')
    ax.set_xlim(0, 0.093); ax.set_ylim(*ylim)
    ax.set_xlabel(r'mean accidental prompts per $\pm$5 ns window, $\mu$', fontsize=8.5); ax.set_ylabel(lab)
    ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08])
    if panel == 'c':
        sl = KR['slopes']['closest']['dR_acc359']; sa = KR['slopes']['all']['dR_acc359']
        key(ax, 0.40, 0.84, f'slope per 0.01 in $\\mu$:\nclosest {sg(sl["slope_per_0p01"], "+.2f")}$\\pm${sl["slope_sd"]:.2f}\nall-prompts {sg(sa["slope_per_0p01"], "+.2f")}$\\pm${sa["slope_sd"]:.2f}',
            ha='left', va='top', fontsize=8.0)
        ax.legend(frameon=False, loc='lower left', fontsize=8.0, handletextpad=0.3, borderaxespad=0.2)
        print('dR slopes', sl, sa)
    else:
        print('bG all-prompts range', va.min(), va.max(), 'closest', vc.min(), vc.max())
    st.panel_letter(ax, panel)
save(fig, 'fig10_rate_law')
