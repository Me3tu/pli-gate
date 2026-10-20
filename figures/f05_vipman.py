"""Fig. 5 (stem fig3_vipman_axial), post-correction O4 triples (the binary of
every other result), no truth in any selection step.
(a) MLEM activity (4 mm, the two 4 mm slices inside the 8 mm slice averaged).
(b) per-voxel tau-hat at 8 mm (Table 4 estimator: mean(tau-1 | tau>1)), scale
    1.7-2.6 ns, voxels with N<30 masked.
(c) per-voxel counting-statistics scale tau-hat/sqrt(N) (ps).
(d) per tissue, paired: the full chain at 8 mm (filled) against the same estimator and
    aggregation applied to truth-tagged annihilations without any detector or reconstruction
    step (open). The vertical distance is the contribution of detection and reconstruction.
    The 4 and 6 mm values are in Supplementary Table S2 (feedback round 3).
The 8 mm tau-hat map is read from a cache written from the same triples
(2026-09-30_feedback3/figdata/vipman_tau8_cache.npz) when present."""
import json
import numpy as np
try:
    import h5py
except ImportError:
    h5py = None
from common import *  # noqa
from matplotlib.lines import Line2D

O4 = RES / '2026-08-13_overnight/O4_rerun/vipman'
RECON = RES / '2026-05-05_2122_phase3q_production/recon_phase3q'
CACHE = F3 / 'figdata/vipman_tau8_cache.npz'
n8 = (70, 35, 10)
if CACHE.exists():
    _c = np.load(CACHE); s, N = _c['s'], _c['N']
else:
    with h5py.File(O4 / 'triples_detector_phase3q.h5', 'r') as f:
        P = f['derived']['P_annih_tof_mm'][:].astype(np.float32)
        tau = f['derived']['tau_corrected_ns'][:].astype(np.float32)
    EXT0 = np.array([-280.0, -140.0, -40.0])
    i8 = np.floor((P - EXT0) / 8.0).astype(np.int64)
    ok = np.all((i8 >= 0) & (i8 < np.array(n8)), 1) & (tau > 1.0)
    fl = np.ravel_multi_index(i8[ok].T, n8)
    s = np.bincount(fl, weights=(tau[ok] - 1.0).astype(np.float64), minlength=np.prod(n8)).reshape(n8)
    N = np.bincount(fl, minlength=np.prod(n8)).reshape(n8)
valid = N >= 30
tmap = np.where(valid, s / np.maximum(N, 1), np.nan)
smap = np.where(valid, tmap / np.sqrt(np.maximum(N, 1)), np.nan)
sup4 = np.load(RECON / 'support_mask.npy')
Z8, Z4 = 5, [10, 11]
EXT = (-280, 280, -140, 140)
# Review round 17 (second report, Sec. 8): panel (a) shows the MLEM reconstruction of the corrected-yield acquisition
# (2026-10-19_review16/vipman_recon_O4, simple_recon single, 4 mm, 10 iterations, the settings of recon_phase3q),
# central two 4 mm slices. The support mask is identical.
# Code release v1.1.0: the slice is read in float32 from 2026-10-20_review19/vipman_slice/vipman_O4_slice.npz
# (written by make_vipman_slice.py in that folder from vipman_recon_O4/activity_x.npy). The figure in the submitted
# paper was drawn from an 8-bit quantized transfer of the same slice (2026-10-19_review16/vipman_O4_slice.npz).
# The recon_phase3q activity, previously loaded here and then overwritten by this slice, is no longer read.
_V = np.load(RES / '2026-10-20_review19/vipman_slice/vipman_O4_slice.npz')
assert (_V['sup2'] == sup4[:, :, Z4]).all()
act_z = np.where(_V['sup2'].any(axis=2), _V['act_z'], np.nan).T
sup_z = sup4[:, :, Z4].any(axis=2).astype(float).T
# body support at 8 mm (any 4 mm support voxel in the 2x2x2 block). Lifetime and
# uncertainty are shown inside the reconstructed body support only (TOF-blurred
# placement puts a thin, low-count rim of events outside it). Voxels inside the
# body without an estimate (N<30) would be drawn grey, outside the body white.
sup8 = sup4.reshape(70, 2, 35, 2, 10, 2).any(axis=(1, 3, 5))
out_rim = valid[:, :, Z8] & ~sup8[:, :, Z8]
print('valid voxels outside body support (not shown)', int(out_rim.sum()), 'median N there', float(np.median(N[:, :, Z8][out_rim])),
      'median N inside', float(np.median(N[:, :, Z8][valid[:, :, Z8] & sup8[:, :, Z8]])))
tau_z = np.where(sup8[:, :, Z8], tmap[:, :, Z8], np.nan).T; sig_z = 1000 * np.where(sup8[:, :, Z8], smap[:, :, Z8], np.nan).T
nodata_z = np.where(sup8[:, :, Z8] & ~valid[:, :, Z8], 1.0, np.nan).T
from matplotlib.colors import ListedColormap
GREYMASK = ListedColormap(['#c8c8c8'])
print('in-body voxels without estimate in slice', int(np.nansum(nodata_z)))
rows_, cols_ = np.where(sup4[:, :, Z4].any(axis=2))
x0, x1 = (rows_.min() - 2) * 4 - 280, (rows_.max() + 3) * 4 - 280
y0, y1 = (cols_.min() - 2) * 4 - 140, (cols_.max() + 3) * 4 - 140
print('valid voxels in slice', int(valid[:, :, Z8].sum()), 'tau range', np.nanmin(tau_z), np.nanmax(tau_z),
      'sigma ps median', np.nanmedian(sig_z))

fig = mm_fig(150, 96)
ar = (x1 - x0) / (y1 - y0); wm = 53.0; hm = wm / ar; ytop = 87.0
ax_a = ax_mm(fig, 14, ytop - hm, wm, hm); cax_a = ax_mm(fig, 14 + wm + 1.5, ytop - hm, 2.0, hm)
ax_b = ax_mm(fig, 84, ytop - hm, wm, hm); cax_b = ax_mm(fig, 84 + wm + 1.5, ytop - hm, 2.0, hm)
ax_c = ax_mm(fig, 14, 12, wm, hm); cax_c = ax_mm(fig, 14 + wm + 1.5, 12, 2.0, hm)
ax_d = ax_mm(fig, 95, 14, 52, 31)
print('map height mm', round(hm, 1))
cm_a = plt.get_cmap('magma').copy(); cm_a.set_bad('#1a1a1a')
im = ax_a.imshow(act_z, cmap=cm_a, vmin=0, vmax=float(np.nanpercentile(act_z, 99)), extent=EXT, origin='lower',
                 aspect='equal', interpolation='nearest')
ax_a.contour(sup_z, levels=[0.5], colors='white', linewidths=0.6, extent=EXT, origin='lower')
cb = fig.colorbar(im, cax=cax_a); cb.ax.set_title('activity\n(a.u.)', fontsize=8.0, pad=2); cb.ax.tick_params(labelsize=8.0)
cm_b = plt.get_cmap('cividis').copy(); cm_b.set_bad('#ffffff')
im = ax_b.imshow(tau_z, cmap=cm_b, vmin=1.7, vmax=2.6, extent=EXT, origin='lower', aspect='equal', interpolation='nearest')
ax_b.imshow(nodata_z, cmap=GREYMASK, vmin=0, vmax=1, extent=EXT, origin='lower', aspect='equal', interpolation='nearest')
ax_b.contour(sup_z, levels=[0.5], colors='black', linewidths=0.6, extent=EXT, origin='lower')
cb = fig.colorbar(im, cax=cax_b, ticks=[1.8, 2.0, 2.2, 2.4, 2.6]); cb.ax.set_yticklabels(['1.8', '2.0', '2.2', '2.4', '2.6'])
cb.ax.set_title(r'$\hat\tau$ (ns)', fontsize=8.0, pad=3); cb.ax.tick_params(labelsize=8.0)
cm_c = plt.get_cmap('Greys').copy(); cm_c.set_bad('#ffffff')
im = ax_c.imshow(sig_z, cmap=cm_c, vmin=0, vmax=float(np.nanpercentile(sig_z, 98)), extent=EXT, origin='lower',
                 aspect='equal', interpolation='nearest')
ax_c.imshow(nodata_z, cmap=GREYMASK, vmin=0, vmax=1, extent=EXT, origin='lower', aspect='equal', interpolation='nearest')
ax_c.contour(sup_z, levels=[0.5], colors=THIS, linewidths=0.6, extent=EXT, origin='lower')
cb = fig.colorbar(im, cax=cax_c); cb.ax.set_title(r'$\hat\tau/\sqrt{N}$ (ps)', fontsize=8.0, pad=3); cb.ax.tick_params(labelsize=8.0)
for ax, L in ((ax_a, 'a'), (ax_b, 'b'), (ax_c, 'c')):
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xlabel('x (mm)'); ax.set_yticks([-100, 0, 100])
    if ax is ax_b:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('y (mm)')
    st.panel_letter(ax, L)
ax_b.set_ylabel(''); ax_b.set_yticklabels([])

Dc = json.load(open(R2 / 'vipman_decomp.json')); Or = json.load(open(R2 / 'vipman_oracle_mix.json'))
CFG = {"Lung": 1.85, "Adipose": 2.55, "Heart": 2.00, "Muscle": 2.03, "Body": 2.00, "Blood": 1.80, "Cartilage": 2.10, "SpineBone": 1.60}
tis = list(CFG)
xs = np.arange(len(tis))
chain = np.array([Dc['B_truthfree_8mm'][t]['resid'] for t in tis])
truth = np.array([100 * (Or['E2'][t] - CFG[t]) / CFG[t] for t in tis])
ax_d.axhline(0, color=REF, lw=0.8, zorder=1)
for x, c_, t_ in zip(xs, chain, truth):
    ax_d.plot([x, x], [t_, c_], '-', color='#9a9a9a', lw=1.3, zorder=2, solid_capstyle='butt')
ax_d.plot(xs, truth, 'o', ms=5.4, mfc='white', mec=THIS, mew=1.1, zorder=4)
ax_d.plot(xs, chain, 'o', ms=5.0, color=THIS, zorder=5)
print('chain - truth (pp):', dict(zip(tis, np.round(chain - truth, 2))), 'max', round(float(np.max(np.abs(chain - truth))), 2))
LOW = {'SpineBone': 'spine bone'}
ax_d.set_xticks(xs); ax_d.set_xticklabels([LOW.get(t, t.lower()) for t in tis], rotation=30, ha='right', fontsize=8.0)
ax_d.set_xlim(-0.6, len(tis) - 0.4)
ax_d.set_ylabel(r'$(\hat\tau-\tau_{\rm cfg})/\tau_{\rm cfg}$ (%)'); ax_d.set_ylim(-14, 40); ax_d.set_yticks([-10, 0, 10, 20, 30])
h = [Line2D([], [], marker='o', ms=5.0, color=THIS, ls='none'), Line2D([], [], marker='o', ms=5.4, mfc='white', mec=THIS, ls='none')]
dmax = float(np.max(np.abs(chain - truth)))
key(ax_d, 0.985, 0.04, f'|chain $-$ records| $\\leq$ {dmax:.1f} pp', ha='right', va='bottom', fontsize=8.0)
key(ax_d, xs[-1] + 0.25, truth[-1] - 5.0, 'already at\ntruth level', transform=ax_d.transData, ha='right', va='top', fontsize=8.0)
ax_d.legend(h, ['full chain, 8 mm voxels', 'annihilation records,\nsame aggregation'], frameon=False, loc='upper left',
            fontsize=8.0, handletextpad=0.3)
st.panel_letter(ax_d, 'd')
save(fig, 'fig3_vipman_axial')
