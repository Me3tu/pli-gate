#!/usr/bin/env python3
"""Review round 5, item 12: the activity peak search of
2026-05-07_2140_phase3r3_batch1_clinical/run_roi_analysis_and_render.py (Part A), copied
unchanged, applied to (i) the corrected-yield activity reconstruction shown in Fig. 6a
(2026-08-16_fig4_rebuild/recon_4mm) and reported next to (ii) the May 10^9-primary maxima
that define the published regions (per_roi_tau_FINAL.csv). Adds the coverage of each configured
insert sphere by its region sphere (radius sqrt(r^2 + 7^2) mm, centre in the z = 0 plane)."""
import json, csv, os
from pathlib import Path
import numpy as np
import scipy.ndimage as ndi
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = RES / "2026-10-08_lesion_peaks"
VOX = 4.0; EXT_X = (-280, 280); EXT_Y = (-140, 140)
TRUTH = [(-66, 0, 4), (-49, 0, 8), (-24, 0, 12), (9, 0, 16), (50, 0, 20)]
RECON_FWHM = 7.0

def lens(r1, r2, d):
    if d >= r1 + r2: return 0.0
    if d <= abs(r2 - r1): return 4 / 3 * np.pi * min(r1, r2) ** 3
    return np.pi * (r1 + r2 - d) ** 2 * (d * d + 2 * d * (r1 + r2) - 3 * (r1 - r2) ** 2) / (12 * d)

def search(activity, support):
    nx, ny, nz = activity.shape; z0 = nz // 2
    act_smooth = ndi.gaussian_filter(activity, sigma=1.0)
    sl = act_smooth[:, :, z0]; sup = support[:, :, z0]
    ig, jg = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    xg = (ig + 0.5) * VOX + EXT_X[0]; yg = (jg + 0.5) * VOX + EXT_Y[0]
    rxy = np.hypot(xg, yg)
    bg = sup & (rxy > 60) & (rxy < 90)
    mean, std = float(sl[bg].mean()), float(sl[bg].std())
    thr = mean + 3 * std
    pk = (ndi.maximum_filter(sl, size=5) == sl) & (sl > thr) & sup
    peaks = sorted(((float(sl[i, j]), i, j) for i, j in np.argwhere(pk)), reverse=True)
    keep = []
    for v, i, j in peaks:
        if any(np.hypot((i - i2) * VOX, (j - j2) * VOX) < 12.0 for _, i2, j2 in keep): continue
        keep.append((v, i, j))
    keep = keep[:5]
    rec = [dict(x_mm=(i + 0.5) * VOX + EXT_X[0], y_mm=(j + 0.5) * VOX + EXT_Y[0], val=v) for v, i, j in keep]
    # value and nearest local maximum around the r = 4 mm insert
    i4 = int(round((-66 - EXT_X[0]) / VOX - 0.5)); j4 = int(round((0 - EXT_Y[0]) / VOX - 0.5))
    near = (np.hypot(xg + 66, yg) <= 8.0) & sup
    return dict(bg_mean=mean, bg_std=std, threshold=thr, n_bg=int(bg.sum()), maxima=rec,
                r4_insert_max_value_within_8mm=float(sl[near].max()))

def score(maxima):
    mx = sorted(maxima, key=lambda p: p["x_mm"])
    truth = TRUTH if len(mx) == 5 else TRUTH[1:]   # assignment by x-order as in the published code
    out = []
    for p, (tx, ty, tr) in zip(mx, truth):
        R = float(np.hypot(tr, RECON_FWHM)); d = float(np.hypot(p["x_mm"] - tx, p["y_mm"] - ty))
        out.append(dict(r_mm=tr, centre=[p["x_mm"], p["y_mm"]], offset_mm=d, R_roi_mm=R,
                        truth_covered=lens(tr, R, d) / (4 / 3 * np.pi * tr ** 3)))
    d4 = min(float(np.hypot(p["x_mm"] + 66, p["y_mm"])) for p in maxima)
    return out, d4

res = {}
act = np.load(RES / "2026-08-16_fig4_rebuild/recon_4mm/activity_x.npy")
sup = np.load(RES / "2026-08-16_fig4_rebuild/recon_4mm/support_mask.npy")
A = search(act, sup); s, d4 = score(A["maxima"])
res["corrected_yield_recon"] = dict(A, regions=s, nearest_maximum_to_r4_mm=d4, shape=list(act.shape))
rows = list(csv.DictReader(open(RES / "2026-05-08_0010_phase3s3_brain_tumor_1B_KEPT_SMALLFILES/brain_tumor_recon/roi_analysis/per_roi_tau_FINAL.csv")))
may = [dict(x_mm=float(r["x_mm"]), y_mm=float(r["y_mm"]), val=float("nan")) for r in rows if r["rank"] != "BG"]
s2, d42 = score(may)
res["may_1e9_published_regions"] = dict(maxima=may, regions=s2, nearest_maximum_to_r4_mm=d42)
json.dump(res, open(OUT / "peaks_result.json", "w"), indent=1)
for k, v in res.items():
    print("==", k, {kk: v[kk] for kk in ("threshold", "bg_mean", "bg_std", "r4_insert_max_value_within_8mm") if kk in v})
    for r in v["regions"]:
        print(f"   r={r['r_mm']:2d}  centre=({r['centre'][0]:+.0f},{r['centre'][1]:+.0f})  offset={r['offset_mm']:.1f} mm  R={r['R_roi_mm']:.2f}  covered={100*r['truth_covered']:.1f}%")
    print("   nearest maximum to r=4 insert:", round(v["nearest_maximum_to_r4_mm"], 1), "mm")
