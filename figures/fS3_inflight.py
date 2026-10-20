"""Supplementary Fig. S8 (stem figS4_inflight), feedback round 3: effect of suppressing
in-flight annihilation on the headline observables.
(a) physics level, 1e6 44Sc positrons at a water/adipose interface: annihilation-delay spectrum
    of the shipped list and of the hybrid list (stock in-flight annihilation kept, custom
    processes at rest). Lower strip: hybrid/shipped per 0.5 ns with the expectation 1 - f_inflight.
(b) photon energies of the in-flight annihilations (hybrid list), 44Sc and 124I, per photon, with
    the tight annihilation window 481-541 keV and the prompt window 1100-1200 keV.
(c) full GATE chain, eight-tube phantom at 388 ps, hybrid minus shipped per tube, in pp of the
    configured lifetime: truth-grouped lifetime, all-prompts and closest-prompt truncated-mean
    lifetimes, and the rule artifact Delta|b|. Independent-arm errors, pooled values below.
(d) relative change of the count observables against the ranges registered before the run
    (REGISTERED_INFLIGHT.md): 511-window singles, prompt-window singles, prompt competition mu,
    genuine triples (mean over tubes, range as bar), and the physics-level delayed yield.
Data: 2026-09-30_feedback3/inflight/{physics_level.json, physics_level_delay_hists.npz,
physics_level_energy_hists.npz, compare_inflight.json}."""
import json
import os
import numpy as np
from common import *  # noqa
from matplotlib.lines import Line2D

I = F3 / 'inflight'
P = json.load(open(I / 'physics_level.json'))
C = json.load(open(os.environ.get('PLI_INFLIGHT_JSON', str(I / 'compare_inflight.json'))))
DH = np.load(I / 'physics_level_delay_hists.npz'); EH = np.load(I / 'physics_level_energy_hists.npz')
N0 = 1_000_000
ph, pp = P['sc44_hybrid'], P['sc44_pli']; cmp_ = P['sc44_compare']
f_if = ph['inflight_frac']

FH = 132.0
fig = mm_fig(150, FH)
# ---------------- (a) delay spectrum ----------------
e = DH['edges']; c = 0.5 * (e[1:] + e[:-1]); w = np.diff(e)
sel = (e[:-1] >= 0.05) & (e[1:] <= 8.0)
ax = ax_mm(fig, 16, 94, 54, 31)
ax.step(c[sel], DH['sc44_pli'][sel] / N0, where='mid', color=REF, lw=1.0, label='shipped list')
ax.step(c[sel], DH['sc44_hybrid'][sel] / N0, where='mid', color=THIS, lw=1.0, ls='--', label='in-flight kept')
ax.set_yscale('log'); ax.set_xlim(0, 8); ax.set_ylim(1e-4, 0.3)
ax.set_ylabel('positrons per 50 ps'); ax.tick_params(axis='x', labelbottom=False)
ax.legend(frameon=False, loc='upper right', handlelength=1.8, borderaxespad=0.1)
st.panel_letter(ax, 'a')
ax2 = ax_mm(fig, 16, 77, 54, 15, sharex=ax)
eb = np.arange(0.5, 8.001, 0.5)
hp = np.array([DH['sc44_pli'][(c >= lo) & (c < hi)].sum() for lo, hi in zip(eb[:-1], eb[1:])])
hh = np.array([DH['sc44_hybrid'][(c >= lo) & (c < hi)].sum() for lo, hi in zip(eb[:-1], eb[1:])])
r = hh / hp; rs = r * np.sqrt(1 / hh + 1 / hp); cb_ = 0.5 * (eb[1:] + eb[:-1])
ax2.axhline(1.0, color='#bbbbbb', lw=0.7)
ax2.axhline(1 - f_if, color=THIS, lw=0.9, ls=':')
ax2.errorbar(cb_, r, yerr=rs, fmt='o', ms=2.6, color=THIS, elinewidth=0.7, capsize=0, lw=0)
ax2.set_ylim(0.90, 1.05); ax2.set_yticks([0.9, 0.95, 1.0])
ax2.set_xlabel('annihilation delay (ns)'); ax2.set_ylabel('ratio', fontsize=8.0)
key(ax2, 0.98, 0.06, f'1 $-$ in-flight fraction = {1 - f_if:.3f}', ha='right', va='bottom', fontsize=8.0, color=THIS)
print(f"44Sc in-flight {100*f_if:.2f}% KE {ph['inflight_meanKE_keV']:.0f} keV; tail ratio {cmp_['tail_yield_ratio_hybrid_over_pli']:.4f}"
      f"+-{cmp_['tail_yield_ratio_se']:.4f}; tail mean diff {cmp_['tail_mean_diff_ps']:+.1f}+-{cmp_['tail_mean_diff_se_ps']:.1f} ps")

# ---------------- (b) in-flight photon energies ----------------
ax = ax_mm(fig, 92, 77, 56, 48)
ee = EH['edges']; ec = 0.5 * (ee[1:] + ee[:-1])
for iso, col, ls, lab in (('sc44', THIS, '-', r'$^{44}$Sc'), ('i124', CTRL, '--', r'$^{124}$I')):
    h = EH[f'{iso}_inflight_photons'].astype(float)
    ax.step(ec, h / h.sum() / 0.010, where='mid', color=col, lw=1.0, ls=ls, label=lab)
ax.axvspan(481, 541, color=REF, alpha=0.10, lw=0); ax.axvspan(1100, 1200, color=REF, alpha=0.10, lw=0)
ax.text(511, 1.02, '481–541', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8.0, color='#555555')
ax.text(1150, 1.02, '1100–1200', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8.0, color='#555555')
ax.set_yscale('log'); ax.set_xlim(0, 2000); ax.set_ylim(1e-3, 5); ax.set_xticks([0, 500, 1000, 1500])  # referee read v2.7: 2000 label clipped
ax.set_xlabel('photon energy (keV)'); ax.set_ylabel('in-flight photons per MeV\n(normalized)')
ax.legend(frameon=False, loc='upper right', handlelength=1.8, borderaxespad=0.1)
key(ax, 0.03, 0.04, f"both photons in 481–541 keV:\n{100*ph['both_in_w481_541_inflight']:.2f}% of in-flight pairs", ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.0))
st.panel_letter(ax, 'b', dx_pt=-30)

# ---------------- (c) GATE per tube ----------------
T = C['tubes']; PO = C['pooled']
ax = ax_mm(fig, 36, 22, 52, 43)
yy = np.arange(len(T))[::-1].astype(float)
Q = [('oracle_pct', 'oracle_pct_se', 'o', REF, REF, 'event-identity lifetime'),
     ('tm_all', 'tm_all_se', 'o', THIS, THIS, 'all-prompts lifetime'),
     ('tm_closest', 'tm_closest_se', 'o', THIS, 'white', 'closest-prompt lifetime'),
     ('did_tm', 'did_tm_se', 'D', CTRL, CTRL, r'selection shift $\Delta|b|$')]
off = np.array([0.27, 0.09, -0.09, -0.27])
ax.axvline(0, color='#bbbbbb', lw=0.8, zorder=1)
for (k, ks, m, ec_, fc, lab), dy in zip(Q, off):
    v = np.array([t['d'][k] for t in T]); s = np.array([t['d'][ks] for t in T])
    ax.errorbar(v, yy + dy, xerr=s, fmt=m, ms=3.6, color=ec_, mfc=fc, mec=ec_, elinewidth=0.7, capsize=0, lw=0, zorder=3, label=lab)
    pm = PO[k]
    ax.errorbar([pm['mean']], [-1.3 + dy], xerr=[pm['se']], fmt=m, ms=4.2, color=ec_, mfc=fc, mec=ec_, elinewidth=1.0, capsize=1.6, lw=0, zorder=4)
ax.axhline(-0.65, color='#dddddd', lw=0.6)
ax.set_yticks(list(yy) + [-1.3]); ax.set_yticklabels([NAME[t['region']] for t in T] + ['pooled'], fontsize=8.0)
ax.set_ylim(-1.8, len(T) - 0.45)
lim = max(3.0, float(np.ceil(np.nanmax([abs(t['d'][k]) + t['d'][ks] for t in T for (k, ks, *_ ) in Q]))))
ax.set_xlim(-lim, lim)
ax.set_xlabel(r'in-flight kept $-$ shipped (pp of $\tau_{\rm cfg}$)')
ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, frameon=False, loc='lower left', bbox_to_anchor=(2 / 150, 0.3 / FH), ncol=2, fontsize=8.0,
           handletextpad=0.2, columnspacing=0.8, borderaxespad=0.0)
st.panel_letter(ax, 'c', dx_pt=-92)
for k, *_ in Q:
    print(f"POOLED {k:11s} {PO[k]['mean']:+.3f} +- {PO[k]['se']:.3f}  z {PO[k]['z']:+.2f}  chi2 {PO[k]['chi2']:.1f}/{PO[k]['df']}")

# ---------------- (d) count observables against registered ranges ----------------
ax = ax_mm(fig, 118, 22, 30, 43)
SR = C['singles']['rel']; ts = np.array([t['d']['true_sig_rel'] for t in T])
rows = [('511-window singles', 100 * SR['n_511'], None, (-3.0, -1.0)),
        ('prompt-window singles', 100 * SR['n_prompt'], None, (-0.5, 0.5)),
        (r'accidental-prompt rate $\mu$', 100 * SR['mu'], None, (-0.5, 0.5)),
        ('signal triples (event identity)', 100 * ts.mean(), (100 * ts.min(), 100 * ts.max()), (-3.5, -1.5)),
        ('delayed yield (physics)', 100 * (cmp_['tail_yield_ratio_hybrid_over_pli'] - 1), 100 * cmp_['tail_yield_ratio_se'], None)]
yd = np.arange(len(rows))[::-1].astype(float)
ax.axvline(0, color='#bbbbbb', lw=0.8, zorder=1)
for y, (lab, v, err, band) in zip(yd, rows):
    if band is not None:
        ax.fill_betweenx([y - 0.32, y + 0.32], band[0], band[1], color=st.BAND, alpha=0.14, lw=0, zorder=1)
    if isinstance(err, tuple):
        ax.plot(err, [y, y], color=THIS, lw=1.0, zorder=2)
        ax.plot([v], [y], 'o', ms=4.2, color=THIS, zorder=3)
    elif err is not None:
        ax.errorbar([v], [y], xerr=[err], fmt='o', ms=4.2, color=THIS, elinewidth=0.8, capsize=1.4, lw=0, zorder=3)
        ax.plot([-100 * f_if] * 2, [y - 0.3, y + 0.3], color=REF, lw=1.4, zorder=4)
    else:
        ax.plot([v], [y], 'o', ms=4.2, color=THIS, zorder=3)
    print(f"{lab:26s} {v:+.3f}%  band {band}")
ax.set_yticks(yd); ax.set_yticklabels([r[0] for r in rows], fontsize=8.0)
ax.set_xlim(-4.5, 1.5); ax.set_xticks([-4, -2, 0])
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xlabel('relative change (%)')
ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
fig.legend([plt.Rectangle((0, 0), 1, 1, color=st.BAND, alpha=0.14, lw=0), Line2D([0], [0], marker='|', markersize=10, markeredgewidth=1.6, markeredgecolor=REF, color=REF, linestyle='None')],
           ['registered range', r'expected: $-$in-flight fraction'], frameon=False, loc='lower right', bbox_to_anchor=(148 / 150, 0.3 / FH),
           ncol=2, fontsize=8.0, handletextpad=0.3, columnspacing=0.8, borderaxespad=0.0)
st.panel_letter(ax, 'd', dx_pt=-100)
print('scores', C.get('scores'))
save(fig, 'figS4_inflight')
