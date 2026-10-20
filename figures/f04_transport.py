"""Fig. 4 (stem fig_r2_transport_novelty): (a) radial cross-section, core and
shell, migration across the boundary; (b) Delta tau-hat (vertex minus source arm)
against the lifetime step: null control at step 0, the two material pairs, and
the first-order estimator-specific consistency check (same-run shares and yields)."""
import json
import numpy as np
from common import *  # noqa
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

R = json.load(open(FS / '2026-08-13_referee_response/r2_results.json'))
PAIRS = [('A', 2.72, 1.77), ('B', 2.54, 1.9215)]
# Audit round 7 (2026-10-08): first-order prediction recomputed from the as-run accounting instead of the
# frozen PRED of findings_R2 (2.6368, 2.4824), which weighted the shell by survival alone and omitted the
# annihilations outside core and shell. Shift = sum_m w_m (tau_m - tau_core), w_m = tail share (delay > 1 ns)
# of material m among the transport-arm annihilations of core-origin positrons, from the measured transition
# matrix and the (P_oPs, tau) with which the runs were made. As-run yields: P_oPs = I3 (o-Ps-intensity
# convention) for the cores and the pair-A shell; the pair-B water shell ran before the 2026-08-14 I3-convention
# correction at 0.75*0.27; structural materials (LSO, PTFE) and world air carried the fallback (0.75*0.35, 2.0 ns).
# These yields reproduce the observed transport/source o-Ps ratios of all three pairs within 0.05% (printed).
TAU_MIN = 1.0
FALLBACK = (0.75 * 0.35, 2.0)
ASRUN = {'A': {'adipose_moskal2023': (0.291, 2.72), 'brain_tumor_moskal2024': (0.30, 1.77)},
         'B': {'adipose_avachat2024': (0.186, 2.54), 'water_shibuya2020_n2': (0.75 * 0.27, 1.9215)},
         'null': {'adipose_moskal2023': (0.291, 2.0), 'brain_tumor_moskal2024': (0.30, 2.0)}}


def first_order(t, params):
    core = t['core_material']; tm = t['transition_matrix'][core]
    comp = {m: params.get(m, FALLBACK) for m in tm}
    tail = {m: tm[m] * comp[m][0] * np.exp(-TAU_MIN / comp[m][1]) for m in tm}
    T = sum(tail.values())
    shift = sum(tail[m] * (comp[m][1] - comp[core][1]) for m in tm) / T
    ops_ratio = sum(tm[m] * comp[m][0] for m in tm) / comp[core][0]
    return shift, ops_ratio, {m: tail[m] / T for m in tm}


def ops_obs(t):
    c = t['channel_counts']; return (c['oPs2g'] + c['oPs3g']) / t['n_annihilations']


for tag, key_t, key_s in (('null', 'r2a_transport', 'r2a_source'), ('A', 'r2b_A_transport', 'r2b_A_source'),
                          ('B', 'r2b_B_transport', 'r2b_B_source')):
    sh, opr, w = first_order(R[key_t], ASRUN[tag])
    print(tag, 'o-Ps ratio predicted', round(opr, 5), 'observed', round(ops_obs(R[key_t]) / ops_obs(R[key_s]), 5),
          'w', {k: round(v, 4) for k, v in w.items() if v > 1e-4})

fig = mm_fig(150, 66)

# (a) radial cross-section ------------------------------------------------
ax = ax_mm(fig, 8, 14, 52, 45)
ax.add_patch(Rectangle((0, 0), 10, 1, facecolor=THIS, alpha=0.35, lw=0))
ax.add_patch(Rectangle((10, 0), 50, 1, facecolor=CTRL, alpha=0.18, lw=0))
ax.axvline(10, color=REF, lw=1.0)
ax.text(5, 0.84, 'core', ha='center', va='center', fontsize=8.0)
ax.text(5, 0.70, r'$\tau_{\rm core}$', ha='center', va='center', fontsize=8.0)
ax.text(35, 0.84, 'shell', ha='center', va='center', fontsize=8.0)
ax.text(35, 0.70, r'$\tau_{\rm shell}$', ha='center', va='center', fontsize=8.0)
ax.text(5, 0.08, 'source\nhere', ha='center', va='bottom', fontsize=8.0, color='#333333')
rng = np.random.default_rng(4)
for yy in (0.28, 0.40, 0.52):
    x0 = rng.uniform(5.5, 8.5); L = rng.uniform(5, 9)
    ax.annotate('', xy=(x0 + L, yy), xytext=(x0, yy), arrowprops=dict(arrowstyle='-|>', color=REF, lw=0.7, mutation_scale=7))
ax.text(21, 0.40, '10.5% of positrons\nannihilate outside\nthe core (9.7%\nin the shell)', ha='left', va='center', fontsize=8.0)
ax.set_xlim(0, 60); ax.set_ylim(0, 1)
ax.set_xticks([0, 10, 20, 40, 60]); ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.set_xlabel('radius (mm)')
st.panel_letter(ax, 'a', dx_pt=-8)

# (b) Delta tau vs lifetime step -----------------------------------------
ax = ax_mm(fig, 84, 14, 63, 45)
ax.axhline(0, color='#aaaaaa', lw=0.7, zorder=1)
ct, cs = R['r2a_transport'], R['r2a_source']
d0 = 1000 * (ct['tau_tail_ns'] - cs['tau_tail_ns']); e0 = 1000 * np.hypot(ct['tau_err_ns'], cs['tau_err_ns'])
ax.errorbar([0.0], [d0], yerr=[e0], fmt='o', ms=5, color=CTRL, elinewidth=1.0, capsize=2, zorder=5)
for tag, tc, ts in PAIRS:
    t, s = R[f'r2b_{tag}_transport'], R[f'r2b_{tag}_source']
    step = tc - ts
    d = 1000 * (t['tau_tail_ns'] - s['tau_tail_ns']); e = 1000 * np.hypot(t['tau_err_ns'], s['tau_err_ns'])
    p = 1000 * first_order(t, ASRUN[tag])[0]
    ax.errorbar([step], [d], yerr=[e], fmt='o', ms=5, color=THIS, elinewidth=1.0, capsize=2, zorder=5)
    ax.plot([step + 0.06], [p], 'D', ms=5, mfc='white', mec=REF, mew=1.0, zorder=6)
    ax.text(step - 0.04, d, f'core {tc:.2f},\nshell {ts:.2f} ns', fontsize=8.0, ha='right', va='center', color='#333333')
    print(tag, 'step', round(step, 3), 'delta ps', round(d, 1), '+/-', round(e, 1), 'pred ps', round(p, 1))
h = [Line2D([], [], marker='o', ms=5, color=THIS, ls='none'), Line2D([], [], marker='D', ms=5, mfc='white', mec=REF, ls='none'),
     Line2D([], [], marker='o', ms=5, color=CTRL, ls='none')]
labs = ['measured, material pair', 'first-order consistency check\nfrom the same-run shares and yields', 'null control (equal lifetimes)']
# Review round 8: registered matched-yield arms (equal o-Ps yields everywhere, three-photon channel off), paired event
# for event; 2026-10-08_r2c_matched_yield/r2c_results.json
# Code release v1.1.0: this input is required (it was optional, and without it the figure silently lost the
# matched-yield points and changed its y-range).
R2C = RES / '2026-10-08_r2c_matched_yield/r2c_results.json'
C = json.load(open(R2C))
for tag, tc, ts in PAIRS:
    r = C[tag]
    ax.errorbar([tc - ts + 0.12], [r['diff_ps']], yerr=[r['se_paired_ps']], fmt='s', ms=4.6, mfc='white', mec=THIS, color=THIS,
                mew=1.0, elinewidth=1.0, capsize=2, zorder=6)
    ax.plot([tc - ts + 0.15, tc - ts + 0.19], [r['first_order_registered_ps']] * 2, '-', lw=1.3, color=REF, zorder=7)
    print(tag, 'matched-yield paired', round(r['diff_ps'], 2), '+/-', round(r['se_paired_ps'], 2))
h.insert(1, Line2D([], [], marker='s', ms=4.6, mfc='white', mec=THIS, ls='none'))
labs.insert(1, 'equal o-Ps yields, paired arms\n(dash: registered first order)')
ax.legend(h, labs, frameon=False, fontsize=8.0, loc='upper right', borderaxespad=0.1)
ax.set_xlim(-0.1, 1.2); ax.set_ylim(-110, 112); ax.set_yticks([-100, -80, -60, -40, -20, 0, 20])
ax.set_xlabel(r'lifetime step $\tau_{\rm core}-\tau_{\rm shell}$ (ns)')
ax.set_ylabel(r'$\hat\tau_{\rm stop}-\hat\tau_{\rm source}$ (ps)')  # review round 13: stopping-point terminology
st.panel_letter(ax, 'b')
save(fig, 'fig_r2_transport_novelty')
