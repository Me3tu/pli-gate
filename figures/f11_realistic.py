"""Fig. 11 (stem fig11_realistic_forest), feedback round 3 redesign: the two findings of
Sec. 3.9 are separated.
(a) early-time model: per resolved region, the rule artifact Delta|b| = |b_closest| - |b_all|
    with direct annihilation at 0 ps (square) and at 388 ps (triangle), joined. Truncated-mean
    readout, rate-matched arms, paired-bootstrap errors.
(b) the same regions as a retention ratio Delta|b|(388 ps) / Delta|b|(0 ps) on its own axis,
    with the registered band 0.56-0.74 (the unit mismatch of drawing a ratio band on the pp axis
    is gone). The six unresolved NEMA spheres are listed in Supplementary Table S14.
(c) estimator dependence at 388 ps: rule shift 100 (tau_all - tau_closest)/tau_cfg under a
    three-component estimator against the truncated mean on identical events. Count-level joint
    likelihood (filled, primary) and the prespecified pseudo-likelihood with clipping (open),
    NEMA background (grey diamond, pseudo-likelihood only, no paired bootstrap). 1:1 and 2:1.
Review round 6 (2026-10-08): panel (c) names the two estimator families on its axes and marks the
pooled fixed-effect values on each axis (truncated mean on x, joint and prespecified on y), so the
doubling reads directly. Joint and prespecified points of a tube are drawn 0.2 pp apart.
Data: fig11_realistic_forest.json, fig11_component_shift.json (frozen), comp_joint_ratio.json."""
import json
import numpy as np
from common import *  # noqa
from matplotlib.lines import Line2D

D = json.load(open(MF / 'fig11_realistic_forest.json'))
K = json.load(open(FS / 'fig11_component_shift.json'))
J = json.load(open(F3 / 'estimator/comp_joint_ratio.json'))
PUB = {r['region']: r for r in K['rows']}
JNT = {r['region']: r for r in J['rows']}
R = [dict(r) for r in D['rows'] if r['resolved']]
BAND = D['band']
# Feedback round 3: the frozen zero-delay eight-tube stream predates the I3-convention correction
# (o-Ps yield x3/4 in every tube), so the eight-tube zero-delay arm is taken from C0, the corrected
# rerun identical to C2 except for the direct delay (2026-09-30_feedback3/c0_sim1/c0_compare.json).
# Code release v1.1.0: this input is required (it was optional, with a fallback to the superseded frozen stream).
C0F = F3 / 'c0_sim1/c0_compare.json'
C0 = json.load(open(C0F))
C0R = {r['region']: r for r in C0['rows']}
for r in R:
    if r['phantom'] == 'Sim-1':
        c = C0R[r['region']]
        assert abs(c['c2']['did'] - r['did388']) < 1e-6, r['region']   # same 388 ps arm
        r['did0'], r['sig0'] = c['c0']['did'], c['c0']['sd_did']
        r['retention'], r['retention_err'] = c['retention_vs_c0'], c['retention_vs_c0_se']
        r['ret_june'], r['ret_june_err'] = c['retention_vs_june0'], c['retention_vs_june0_se']
print('eight-tube zero-delay arm from C0; pooled retention', C0['summary']['retention_ivw_c0'])

# pooled numbers quoted in the text (asserted in the PDF): all three estimators refitted on the
# same 400 paired shard resamples (comp_joint_ratio.json). The frozen analysis gave 9.31 / 1.95 / 4.77.
TM_FE, PUB_FE, JNT_FE = J['hetero_tm']['fe'], J['hetero_pub']['fe'], J['hetero_joint']['fe']
RAT_PUB, RAT_JNT = J['ratio_pub_over_tm']['point'], J['ratio_joint_over_tm']['point']
for r in R:
    if r['phantom'] == 'Sim-1':
        assert abs(JNT[r['region']]['point']['pub'] - PUB[r['region']]['shift_pp']) < 1e-6, r['region']  # same point fits
        assert abs(JNT[r['region']]['point']['tm'] - r['did388']) < 1e-6, r['region']  # both TM biases negative
print(f'TM {TM_FE:.2f}  PUB {PUB_FE:.2f} ({RAT_PUB:.2f}x)  JOINT {JNT_FE:.2f} ({RAT_JNT:.2f}x, CI {J["ratio_joint_over_tm"]["ci95"]})')
print('retention IVW', D['retention_ivw'], D['retention_ivw_err'])

FH = 85
# Review round 16: panel (a) of v2.15 (rule artifact per region) duplicated Table S14 and is dropped. Retention is
# shown against the corrected-yield zero-delay arm only (the discovery stream is in Note S15). The rule-shift panel
# shows the joint likelihood only (the clipped estimator is in Note S11), without fixed-effect pooled markers.
fig = mm_fig(150, FH)
yy = np.arange(len(R))[::-1].astype(float)
yy[-1] -= 0.5                                    # NEMA background set apart from the eight tubes
axr = ax_mm(fig, 36, 14, 38, 60)
axr.axvspan(BAND[0], BAND[1], color=st.BAND, alpha=0.10, lw=0, zorder=1)
for y, r in zip(yy, R):
    axr.errorbar([r['retention']], [y], xerr=[r['retention_err']], fmt='o', ms=4.4, color=THIS, elinewidth=0.8, capsize=1.6, zorder=4)
axr.set_yticks(yy); axr.set_yticklabels([NAME[r['region']] for r in R], fontsize=8.0)
axr.set_ylim(yy[-1] - 0.7, yy[0] + 0.7)
axr.set_xlim(0.15, 0.85); axr.set_xticks([0.2, 0.4, 0.6, 0.8])
axr.tick_params(axis='y', length=0); axr.spines['left'].set_visible(False)
axr.set_xlabel(r'retention $\Delta|b|_{388}/\Delta|b|_{0}$')
key(axr, 0.65, 1.005, 'registered band', ha='center', va='bottom', transform=axr.get_xaxis_transform(), fontsize=7.6, color='#555555')
axr.axhline(yy[-2] - 0.75, color='#dddddd', lw=0.6, zorder=0)
st.panel_letter(axr, 'a', dx_pt=-92)

ax = ax_mm(fig, 100, 14, 46, 60)
xl, yl = (0, 8), (0, 22)
ax.plot([0, 8], [0, 8], color=REF, lw=0.9, ls='--', zorder=1)
ax.plot([0, 8], [0, 16], color=CTRL, lw=0.8, ls=':', zorder=1)
key(ax, 7.9, 7.4, '1:1', transform=ax.transData, ha='right', va='top', fontsize=8.0, color='#555555')
key(ax, 7.9, 16.4, '2:1', transform=ax.transData, ha='right', va='bottom', fontsize=8.0, color='#555555')
for r in R:
    if r['phantom'] != 'Sim-1':
        continue
    j = JNT[r['region']]; x = j['point']['tm']
    ax.errorbar([x], [j['point']['joint']], xerr=[j['sig_shift']['tm']], yerr=[j['sig_shift']['joint']], fmt='o', ms=4.8, color=THIS,
                elinewidth=0.8, capsize=1.6, zorder=4)
ax.set_xlim(*xl); ax.set_ylim(*yl)
ax.set_xlabel('selection shift, truncated mean (pp)'); ax.set_ylabel('selection shift, joint likelihood (pp)')
ax.set_yticks([0, 5, 10, 15, 20])
key(ax, 0.03, 0.97, 'eight tubes, 388 ps', ha='left', va='top', fontsize=8.0)
st.panel_letter(ax, 'b')
save(fig, 'fig11_realistic_forest')
