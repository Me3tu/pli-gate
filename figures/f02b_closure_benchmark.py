"""Fig. 2 (stem fig2_closure_benchmark), review round 16: Secs. 3.1 and 3.2 of v2.15 merged into one figure.
(a) signed closure residual per calibration tube (panel (a) of the former fig2_recipe_validation), error bars
    multiplied by 1.30, the measured run-to-run factor of the three-seed study (Methods). The exploratory trend panel is dropped.
(b) the 124I benchmark per sample (panel (b) of the former fig12_replication): published lifetime with its published
    uncertainty, this work under all-prompts and closest-prompt matching. The statistical part of each error bar of
    this work is multiplied by 1.30 and combined with half the endpoint spread, and z is recomputed with it.
Data: fig2_recipe_validation.json, 2026-10-08_replication_330/repl_sig0.3306.json."""
import json
import math
import numpy as np
from common import *  # noqa
from matplotlib.lines import Line2D

INFL = 1.30   # measured run-to-run factor (Methods), v2.19 (was 1.8)
d = json.load(open(MF / 'fig2_recipe_validation.json'))
T = d['tubes']
NM = {'water_shibuya2020_n2': 'water (N$_2$-saturated)', 'adipose_avachat2024': 'adipose',
      'adipose_moskal2023': 'mediastinal adipose', 'hepatic_avachat2024': 'hepatic',
      'muscle_avachat2024': 'muscle', 'brain_healthy_moskal2024': 'brain (healthy)',
      'brain_tumor_moskal2024': 'brain (tumor)', 'myxoma_moskal2023': 'cardiac myxoma',
      'XAD4_dry': 'XAD-4 (dry)', 'XAD4_wet': 'XAD-4 (wet)', 'XAD4_gelatin': 'XAD-4 (gelatin)',
      'XAD4_water': 'water (demineralized)'}
MK = {'Sc-44': 'o', 'I-124': 's'}
R = json.load(open(RES / '2026-10-08_replication_330/repl_sig0.3306.json'))
TUBES = ['T1_dry', 'T2_wet', 'T3_gelatin', 'T4_water']
NICE2 = {'T1_dry': 'dry\nXAD-4', 'T2_wet': 'wet\nXAD-4', 'T3_gelatin': 'XAD-4 in\ngelatin', 'T4_water': 'demineralized\nwater'}
get = lambda t, r: next(x for x in R['forest'] if x['tube'] == t and x['rule'] == r)

fig = mm_fig(150, 80)
ax = ax_mm(fig, 34, 14, 38, 54)
ax.axvline(0, color=REF, lw=0.8, zorder=2)
rows, y = [], 0
for iso in ('Sc-44', 'I-124'):
    for t in [t for t in T if t['isotope'] == iso]:
        rows.append((y, t)); y -= 1
    y -= 1.1
for yy, t in rows:
    ax.errorbar(t['residual_pct'], yy, xerr=INFL * 100 * t['sigma_ns'] / t['tau_cfg'], fmt=MK[t['isotope']], ms=4.6,
                color=THIS, ecolor=THIS, capsize=1.8, elinewidth=0.9, lw=0, zorder=4)
ax.set_yticks([r[0] for r in rows]); ax.set_yticklabels([NM[r[1]['material']] for r in rows], fontsize=8.0)
ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
ax.set_xlim(-4.6, 3.0); ax.set_xticks([-4, -2, 0, 2])
ax.set_xlabel(r'$(\hat\tau-\tau_{\rm cfg})/\tau_{\rm cfg}$ (%)')
ax.set_ylim(rows[-1][0] - 0.7, rows[0][0] + 1.1)
for i0, lab in ((0, r'$^{44}$Sc tubes'), (8, r'$^{124}$I tubes')):
    ax.text(-4.6, rows[i0][0] + 0.75, lab, ha='left', va='center', fontsize=8.0, color='#444444', fontweight='bold')
st.panel_letter(ax, 'a', dx_pt=-92)
mae = {iso: np.mean([abs(t['residual_pct']) for t in T if t['isotope'] == iso]) for iso in ('Sc-44', 'I-124')}
print('closure mean |residual|', {k: round(v, 3) for k, v in mae.items()})

ax = ax_mm(fig, 90, 14, 28, 54)
yy = np.arange(len(TUBES))[::-1]
zs = {}
for y0, t in zip(yy, TUBES):
    a = get(t, 'all'); c = get(t, 'closest')
    for r in (a, c):
        r['tot_i'] = math.sqrt((INFL * r['stat']) ** 2 + r['syst_half'] ** 2)
        r['z_i'] = (r['ours'] - r['pub']) / math.sqrt(r['tot_i'] ** 2 + r['pub_err'] ** 2)
        zs[(t, r['rule'])] = round(r['z_i'], 2)
    ax.plot([a['pub'], a['ours']], [y0 + 0.13, y0 + 0.13], '-', color='#bbbbbb', lw=1.2, zorder=2)
    ax.plot([a['pub'], c['ours']], [y0 - 0.13, y0 - 0.13], '-', color='#bbbbbb', lw=1.2, zorder=2)
    ax.errorbar([a['pub']], [y0], xerr=[a['pub_err']], fmt='s', ms=5, color=REF, elinewidth=1.4, capsize=2.6, zorder=5)
    ax.errorbar([a['ours']], [y0 + 0.13], xerr=[a['tot_i']], fmt='o', ms=4.6, color=THIS, elinewidth=0.9, capsize=1.8, zorder=6)
    ax.errorbar([c['ours']], [y0 - 0.13], xerr=[c['tot_i']], fmt='o', ms=4.6, color=THIS, mfc='white', mec=THIS,
                elinewidth=0.9, capsize=1.8, zorder=6)
    tr = ax.get_yaxis_transform()
    for xcol, r in ((1.36, a), (1.76, c)):
        pct = 100 * (r['ours'] - r['pub']) / r['pub']
        ax.text(xcol, y0 + 0.07, sg(pct) + '%', transform=tr, fontsize=8.0, ha='right', va='bottom', color=THIS, clip_on=False)
        ax.text(xcol, y0 + 0.02, '$z$ ' + sg(r['z_i']), transform=tr, fontsize=8.0, ha='right', va='top', color='#777777', clip_on=False)
print('benchmark z (inflated)', zs)
tr = ax.get_yaxis_transform()
ax.text(1.36, len(TUBES) - 0.42, 'all', transform=tr, fontsize=8.0, ha='right', va='bottom', color='#444444', clip_on=False)
ax.text(1.76, len(TUBES) - 0.42, 'closest', transform=tr, fontsize=8.0, ha='right', va='bottom', color='#444444', clip_on=False)
h = [Line2D([], [], marker='s', ms=5, color=REF, ls='none'),
     Line2D([], [], marker='o', ms=4.6, color=THIS, ls='none'),
     Line2D([], [], marker='o', ms=4.6, color=THIS, mfc='white', mec=THIS, ls='none')]
fig.legend(h, ['published', 'all-prompts', 'closest-prompt'], frameon=False,
           loc='upper left', bbox_to_anchor=(80 / 150, 79.5 / 80), ncol=3, borderaxespad=0.0, handletextpad=0.3,
           columnspacing=0.8, fontsize=8.0, handlelength=1.4)
ax.set_yticks(yy); ax.set_yticklabels([NICE2[t] for t in TUBES], fontsize=8.0, linespacing=1.0)
ax.set_ylim(-0.5, len(TUBES) - 0.25)
ax.set_xlim(1.68, 2.62); ax.set_xticks([1.8, 2.2, 2.6])
ax.set_xlabel('o-Ps lifetime (ns)')
ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
st.panel_letter(ax, 'b', dx_pt=-64)
save(fig, 'fig2_closure_benchmark')
json.dump(dict(closure_mae=mae, benchmark_z_inflated={f'{k[0]}_{k[1]}': v for k, v in zs.items()}),
          open(OUT / 'fig2_closure_benchmark.json', 'w'), indent=1)
