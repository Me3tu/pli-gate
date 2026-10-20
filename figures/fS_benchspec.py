"""Supplementary figure (stem figS_benchspec), review round 16: the benchmark spectrum moved from main Fig. 3(a).
Dry XAD-4, all-prompts, randoms-subtracted spectrum over the fit range with the three-component model as in the fit
objective (total thin, so that the data points stay visible), o-Ps and summed short components. Lower panel:
Poisson-normalized residual for every bin of the fit range. Data: 2026-10-08_replication_330/repl_sig0.3306.json."""
import json
import numpy as np
from common import *  # noqa

R = json.load(open(RES / '2026-10-08_replication_330/repl_sig0.3306.json'))
assert abs(R['sigma_e_ns'] - 0.3306) < 1e-9
S = R['spectrum']
fig = mm_fig(110, 84)
bc = np.array(S['bin_centres_ns']); y = np.array(S['data']); mt = np.array(S['model_total'])
mo = np.array(S['model_ops']); ms_ = np.array(S['model_short']); pull = np.array(S['pull']); var = np.array(S['var'])
ax = ax_mm(fig, 18, 36, 86, 44)
pos = y > 0
ax.plot(bc, mt, lw=0.7, color=REF, zorder=3, label='three-component fit')
ax.plot(bc, mo, lw=0.8, ls='--', color=THIS, zorder=2, label='o-Ps component')
ax.plot(bc, ms_, lw=0.8, ls=':', color=THIS, zorder=2, label='p-Ps + direct')
ax.errorbar(bc[pos], y[pos], yerr=np.sqrt(var[pos]), fmt='o', ms=2.4, color=CTRL, mfc='white', mec=CTRL, mew=0.7, ecolor=CTRL,
            elinewidth=0.5, capsize=0, lw=0, zorder=5, label='data (randoms subtracted)')
ax.set_yscale('log'); ax.set_xlim(-1.0, 4.5); ax.set_ylim(30, y.max() * 6)
ax.set_ylabel('counts per 50 ps'); ax.tick_params(axis='x', labelbottom=False)
ax.legend(frameon=False, loc='upper right', fontsize=7.8, handlelength=1.8, borderaxespad=0.3)
key(ax, 0.03, 0.97, 'dry XAD-4,\nall-prompts', ha='left', va='top', fontsize=8.0)
ax2 = ax_mm(fig, 18, 14, 86, 18, sharex=ax)
ax2.axhspan(-2, 2, color=st.BAND, alpha=0.06, lw=0, zorder=0)
ax2.axhline(0, color=REF, lw=0.7)
ax2.plot(bc, pull, 'o', ms=2.0, color=CTRL, mfc=CTRL, mec=CTRL, zorder=3)
ax2.set_ylim(-4, 4); ax2.set_yticks([-3, 0, 3])
ax2.set_xlabel('delay (ns)'); ax2.set_ylabel('pull', fontsize=8.0)
print(f'chi2 {S["chi2"]:.1f}/{S["ndf"]}, max |pull| {S["max_abs_pull"]:.2f} at {S["t_max_abs_pull"]} ns')
save(fig, 'figS_benchspec')
