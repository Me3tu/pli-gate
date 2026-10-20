"""Fig. 6 (stem fig4_brain_tumor), post-correction.
(a) MLEM activity, central 4 mm slice, activity-derived ROIs solid, configured
    insert positions faint dashed. (b) per-voxel Delta tau = tau-hat - <tau_bg>
    on the central 8 mm slice, diverging scale centred at 0. (c) contrast
    tau_ROI - tau_bg vs diameter, configured contrast dashed, 8 mm marked as not
    localized. Data: 2026-08-16_fig4_rebuild, W3 activity ROIs, cluster geometry."""
import json
import numpy as np
import scipy.ndimage as ndi
import matplotlib.patches as mpatches
from common import *  # noqa
from matplotlib.lines import Line2D

SESS = RES / '2026-08-16_fig4_rebuild'
# Code release v1.1.0: 2026-08-13_overnight/W3_detection_sigma/w3_activity_rois.json was loaded here but never used. Load removed.
meta = json.load(open(SESS / 'tau_voxel_map_8mm.json')); TAU_BG = meta['tau_bg_global']
TRUTH = [dict(r=4, cx=-66), dict(r=8, cx=-49), dict(r=12, cx=-24), dict(r=16, cx=9), dict(r=20, cx=50)]
# Review round 17: regions redefined on the corrected-yield reconstruction (2026-10-19_review16/brain_corrected_rois.py,
# maxima of 2026-10-08_lesion_peaks), lifetimes from the same corrected-yield triples.
BR = json.load(open(RES / '2026-10-19_review16/brain_corrected_rois.json'))['corrected_regions']
ROIC = {int(r): (v['centre'][0], v['centre'][1], v['R_ROI']) for r, v in BR['rois'].items()}
act = np.load(SESS / 'recon_4mm/activity_x.npy'); sup = np.load(SESS / 'recon_4mm/support_mask.npy')
z4 = act.shape[2] // 2
act_disp = np.where(sup[:, :, z4], ndi.gaussian_filter(act, 1.0)[:, :, z4], np.nan).T
NX4, NY4 = act.shape[:2]; EXT4 = (-NX4 * 2.0, NX4 * 2.0, -NY4 * 2.0, NY4 * 2.0)
tau_s = np.load(SESS / 'tau_map_smoothed.npy'); mv = np.load(SESS / 'mask_valid.npy')
z8 = tau_s.shape[2] // 2; NX8, NY8 = tau_s.shape[:2]; EXT8 = (-NX8 * 4.0, NX8 * 4.0, -NY8 * 4.0, NY8 * 4.0)
dtau = np.where(mv[:, :, z8], tau_s[:, :, z8] - TAU_BG, np.nan).T

FW, FH = 150.0, 64.0
fig = plt.figure(figsize=(FW * MM, FH * MM))
def mm_axes(x0, y0, w, h):
    return fig.add_axes([x0 / FW, y0 / FH, w / FW, h / FH])
ax_a = mm_axes(13.0, 17.0, 33.0, 33.0); cax_a = mm_axes(47.5, 17.0, 2.0, 33.0)
ax_b = mm_axes(61.0, 17.0, 33.0, 33.0); cax_b = mm_axes(95.5, 17.0, 2.0, 33.0)
ax_c = mm_axes(123.5, 17.0, 25.5, 33.0)
cm_a = plt.get_cmap('magma').copy(); cm_a.set_bad('#1a1a1a')
im = ax_a.imshow(act_disp, cmap=cm_a, vmin=0, vmax=float(np.nanpercentile(act_disp, 99)), extent=EXT4, origin='lower',
                 aspect='equal', interpolation='nearest')
cb = fig.colorbar(im, cax=cax_a); cb.ax.set_title('activity\n(a.u.)', fontsize=8.0, pad=3); cb.ax.tick_params(labelsize=8.0)
for t in TRUTH:
    ax_a.add_patch(mpatches.Circle((t['cx'], 0), t['r'], fill=False, ec='#b0b0b0', lw=0.45, ls=(0, (1.5, 2.0)), alpha=0.8))
for r, (cx, cy, R) in ROIC.items():
    ax_a.add_patch(mpatches.Circle((cx, cy), R, fill=False, ec='#6FD3FF', lw=1.1))
cm_b = plt.get_cmap('RdBu_r').copy(); cm_b.set_bad('#d9d9d9')
im = ax_b.imshow(dtau, cmap=cm_b, vmin=-0.8, vmax=0.8, extent=EXT8, origin='lower', aspect='equal', interpolation='nearest')
cb = fig.colorbar(im, cax=cax_b, ticks=[-0.5, 0.0, 0.5]); cb.ax.set_title(r'$\Delta\hat\tau$ (ns)', fontsize=8.0, pad=3); cb.ax.tick_params(labelsize=8.0)
for r, (cx, cy, R) in ROIC.items():
    ax_b.add_patch(mpatches.Circle((cx, cy), R, fill=False, ec=REF, lw=0.8))
for ax in (ax_a, ax_b):
    ax.set_xlim(-110, 110); ax.set_ylim(-110, 110); ax.set_xticks([-100, 0, 100]); ax.set_yticks([-100, 0, 100])
    ax.set_xlabel('x (mm)')
ax_a.set_ylabel('y (mm)'); ax_b.set_yticklabels([])
h = [Line2D([], [], color='#6FD3FF', lw=1.4), Line2D([], [], color='#8a8a8a', lw=0.6, ls=(0, (1.5, 2.0)))]
fig.legend(h, ['analysis region (from activity)', 'configured position (reference only)'], frameon=False, fontsize=8.0,
           loc='lower left', bbox_to_anchor=(8.0 / FW, 0.0), ncol=2, handlelength=1.6, columnspacing=1.4)
st.panel_letter(ax_a, 'a'); st.panel_letter(ax_b, 'b', dx_pt=-12)

rois = [dict(r_inferred=int(k), **v) for k, v in sorted(BR['rois'].items(), key=lambda kv: int(kv[0]))]; bg = BR['background']
d = np.array([2 * r['r_inferred'] for r in rois])
dt = np.array([r['tau'] - bg['tau'] for r in rois])
e = np.array([np.hypot(r['tau'] / np.sqrt(r['N']), bg['tau'] / np.sqrt(bg['N'])) for r in rois])
ax_c.axhline(0, color='#aaaaaa', lw=0.7, zorder=1)
ax_c.axhline(1.77 - 2.72, color=REF, lw=1.0, ls='--')
ax_c.errorbar(d, dt, yerr=e, fmt='o', ms=5, color=THIS, elinewidth=1.0, capsize=2, zorder=4)
ax_c.axvspan(2, 12.5, color=CTRL, alpha=0.12, lw=0, zorder=0)
ax_c.plot([8], [0], 'x', ms=7, mew=1.6, color=CTRL, mec=CTRL, zorder=5)
ax_c.text(7.25, -0.30, 'not\nlocal-\nized', ha='center', va='center', fontsize=8.0, color='#555555')
ax_c.annotate('', xy=(43, 0.2), xytext=(14, 0.2), arrowprops=dict(arrowstyle='-', color='#777777', lw=0.8))
ax_c.text(28.5, 0.23, 'localized here:\n$\\geq$16 mm', ha='center', va='bottom', fontsize=8.0, color=KEYC)  # review round 8: this configuration only
ax_c.text(43, 1.77 - 2.72 - 0.03, 'configured contrast', ha='right', va='top', fontsize=8.0, color='#444444')
ax_c.set_xlim(2, 44); ax_c.set_ylim(-1.1, 0.35); ax_c.set_xticks([8, 16, 24, 32, 40]); ax_c.tick_params(axis='x', labelsize=8.0); ax_c.set_yticks([-1.0, -0.5, 0.0])
ax_c.set_xlabel('insert diameter (mm)'); ax_c.set_ylabel(r'$\hat\tau_{\rm region}-\hat\tau_{\rm bg}$ (ns)')
st.panel_letter(ax_c, 'c')
for r, x, y, ee in zip(rois, d, dt, e): print('diam', x, 'contrast ns', round(y, 3), '+/-', round(ee, 3), 'nsigma', r['nsigma_analytic'])
save(fig, 'fig4_brain_tumor')
