"""Fig. 9 (stem fig6_clinical_mechanism_v2), feedback round 3.
(a) signed bias of the closest-prompt production readout against the simulated CRT: the eight
    tubes (small dots), their mean joined across the CRT configurations (thick line) and the
    tube range (light band), zero line.
(b) per tube at 303 ps on identical events, the two corrections as sequential shifts:
    closest-prompt (open) -> all-prompts (filled) -> all-prompts with the measured window
    acceptance (teal diamond), arrows along each row, mean |bias| of each arm in the panel.
Data: 2026-10-14_known_randoms/crt_sweep_c0_summary.json (C0, production extractor, review round 13) and
fig6data/crt_sweep_v3_all.csv (frozen stream, extractor of that stage, grey dashed mean), 2026-09-30_feedback3/c0_sim1/acceptance_tm_c0.json."""
import csv, json
import numpy as np
from common import *  # noqa
from matplotlib.lines import Line2D

rows = list(csv.DictReader(open(FS / 'fig6data/crt_sweep_v3_all.csv')))   # discovery stream (extractor of that stage)
# Review round 13 (2026-10-14): panel (a) is the corrected-yield arm C0 with the production extractor
# (2026-10-14_known_randoms/crt_sweep_c0_summary.json), including noiseless timing. The discovery-stream tube
# mean of the prespecified sweep is drawn as a grey dashed line for comparison.
CS = json.load(open(RES / '2026-10-14_known_randoms/crt_sweep_c0_summary.json'))
CRTS = [0, 50, 100, 228, 303]
TUBES = sorted(set(r['roi'] for r in rows))
# Review round 16: discovery-stream line removed from (a), categorical CRT axis labeled as such, per-tube error
# bars in (b) (bootstrap sd x 1.30, the run-to-run factor of the Methods, v2.19). Feedback round 3: panel (b) uses the corrected-yield zero-delay arm C0 (identical to the 388 ps arm
# except for the direct delay). The frozen stream (superseded I3 convention) gave 20.5 -> 6.0 -> 2.2 %.
A = json.load(open(F3 / 'c0_sim1/acceptance_tm_c0.json'))
A_frozen = json.load(open(R2 / 'acceptance_tm.json'))
print('frozen stream: mean |bias|', round(A_frozen['closest_sharp_mean_abs'], 2), round(A_frozen['all_sharp_mean_abs'], 2),
      round(A_frozen['all_acc_measured_mean_abs'], 2))

FH = 74.0
fig = mm_fig(150, FH)
ax = ax_mm(fig, 17, 22, 46, 44)
ax.axhline(0, color=REF, lw=0.8, ls='--')
byc = {r['crt_ps']: r for r in CS['rows']}
V = np.array([byc[c]['sharp_closest']['per_tube'] for c in CRTS])
VD = {c: np.mean([float(next(r for r in rows if r['roi'] == t and int(r['crt_ps']) == c)['bias_trunc_pct']) for t in TUBES])
      for c in (50, 100, 228, 303)}
xs = np.arange(len(CRTS))
ax.fill_between(xs, V.min(1), V.max(1), color=THIS, alpha=0.12, lw=0, zorder=1)
for i in range(len(CRTS)):
    ax.plot(i + np.linspace(-0.10, 0.10, V.shape[1]), V[i], 'o', ms=2.8, color=THIS, mec=THIS, alpha=0.55, zorder=3)
ax.plot(xs, V.mean(1), '-', color=THIS, lw=2.4, zorder=4, solid_capstyle='round')
print('discovery-stream tube means', {c: round(v, 2) for c, v in VD.items()})
for cr, v in zip(CRTS, V): print(cr, 'mean', round(v.mean(), 2), 'range', round(v.min(), 1), round(v.max(), 1))
ax.set_xticks(xs); ax.set_xticklabels(['none' if c == 0 else str(c) for c in CRTS])
ax.set_xlim(-0.35, len(CRTS) - 0.65); ax.set_ylim(-40, 4)
ax.set_xlabel('simulated CRT (ps), categories')
ax.set_ylabel(r'residual $(\hat\tau-\tau_{\rm cfg})/\tau_{\rm cfg}$ (%)')
key(ax, 0.5, 0.885, 'closest-prompt matching', ha='center', va='top', fontsize=8.0)
mlo, mhi = V.mean(1).min(), V.mean(1).max()
key(ax, 0.5, 0.785, f'tube mean {sg(mlo)}% to {sg(mhi)}%', ha='center', va='top', fontsize=8.0, color=THIS)
print('tube-mean range across CRT', round(mlo, 2), round(mhi, 2))
h = [Line2D([], [], marker='o', ms=2.8, color=THIS, alpha=0.55, ls='none'), Line2D([], [], color=THIS, lw=2.4),
     plt.Rectangle((0, 0), 1, 1, color=THIS, alpha=0.12, lw=0)]
ax.legend(h, ['tubes', 'mean', 'range'], frameon=False, loc='lower center', ncol=3, handlelength=1.4,
          handletextpad=0.4, columnspacing=1.0, borderaxespad=0.1, fontsize=8.0)
st.panel_letter(ax, 'a')

ax = ax_mm(fig, 101, 22, 47, 44)
cl = {r['tube']: r for r in A['closest']}; al = {r['tube']: r for r in A['all']}
order = sorted(cl, key=lambda t: cl[t]['sharp'][0])
yy = np.arange(len(order)).astype(float)
ax.axvline(0, color=REF, lw=0.8, zorder=1)
AP = dict(arrowstyle='-|>', color='#9a9a9a', lw=0.9, mutation_scale=7, shrinkA=4.5, shrinkB=4.5)
for y, t in zip(yy, order):
    xc, xa, xm = cl[t]['sharp'][0], al[t]['sharp'][0], al[t]['acc359'][0]   # a priori acceptance (Note S6)
    ax.annotate('', xy=(xa, y), xytext=(xc, y), arrowprops=AP, zorder=2)
    if abs(xm - xa) > 0.8:
        ax.annotate('', xy=(xm, y), xytext=(xa, y), arrowprops=dict(AP, ls=':'), zorder=2)
    ec, ea, em = (1.30 * cl[t]['sharp'][1], 1.30 * al[t]['sharp'][1], 1.30 * al[t]['acc359'][1])
    ax.errorbar([xc], [y], xerr=[ec], fmt='o', ms=5, color=THIS, mfc='white', mec=THIS, elinewidth=0.8, capsize=1.5, zorder=4)
    ax.errorbar([xa], [y], xerr=[ea], fmt='o', ms=5, color=THIS, elinewidth=0.8, capsize=1.5, zorder=4)
    ax.errorbar([xm], [y], xerr=[em], fmt='D', ms=4.6, color=ACC, elinewidth=0.8, capsize=1.5, zorder=5)
ax.set_yticks(yy); ax.set_yticklabels([NAME[t] for t in order], fontsize=8.0)
ax.set_ylim(-0.7, len(order) - 0.3); ax.set_xlim(-32, 6)
ax.set_xlabel(r'residual at 303 ps (%)')
ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
ytop = len(order) + 0.25
mc, ma, mm_ = A['closest_sharp_mean_abs'], A['all_sharp_mean_abs'], A['all_acc359_mean_abs']
key(ax, -0.16, 1.0, f'mean |residual|:  {mc:.1f}%  $\\rightarrow$  {ma:.1f}%  $\\rightarrow$  {mm_:.1f}%', ha='left', va='bottom', fontsize=8.0)
print('mean |bias| closest -> all -> acc:', round(mc, 2), round(ma, 2), round(mm_, 2))
h = [Line2D([], [], marker='o', ms=5, color=THIS, mfc='white', mec=THIS, ls='none'),
     Line2D([], [], marker='o', ms=5, color=THIS, ls='none'), Line2D([], [], marker='D', ms=4.6, color=ACC, ls='none')]
fig.legend(h, ['closest-prompt', 'all-prompts', 'all-prompts, a priori acceptance'], frameon=False,
           loc='lower right', bbox_to_anchor=(148 / 150, 0.8 / FH), ncol=3, handletextpad=0.2, columnspacing=0.9, borderaxespad=0.0)
st.panel_letter(ax, 'b', dx_pt=-92)
save(fig, 'fig6_clinical_mechanism_v2')
