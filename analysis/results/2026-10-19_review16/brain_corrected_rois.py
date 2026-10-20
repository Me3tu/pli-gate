#!/usr/bin/env python3
"""Review round 17 (second report, Sec. 8 and first-round point 4.8): the activity-guided lesion regions of the
brain-tumor phantom redefined on the corrected-yield activity reconstruction (2026-08-16_fig4_rebuild/recon_4mm,
maxima of 2026-10-08_lesion_peaks/peaks_result.json, the published peak search copied unchanged), with the
lifetimes read from the same corrected-yield triples (O4 rerun) and the definition of
2026-08-13_overnight/W3_detection_sigma/june_activity_rois.py imported unchanged:
  sphere R = sqrt(r^2 + 7^2) mm at the activity maximum (z = 0), r the configured radius assigned by order along
  the row; background annulus 80 < r_xy < 95 mm, |z| < 20 mm; tau_hat = mean(tau | tau > 1 ns) - 1 ns;
  sigma = tau/sqrt(N) (analytic) and s/sqrt(N) (empirical); n_sigma = |dtau| / hypot(sigma_ROI, sigma_bg).
Also reports the published-region values on the same triples for comparison, the offset of each maximum from the
configured insert centre, and the partial-volume shift tau_ROI - tau_insert (configured).
Writes brain_corrected_rois.json."""
import json, os, sys
from pathlib import Path
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-13_overnight/W3_detection_sigma"))
import june_activity_rois as J  # noqa: the frozen definition
PK = json.load(open(RES / "2026-10-08_lesion_peaks/peaks_result.json"))["corrected_yield_recon"]["regions"]
TRUTH = {4: (-66, 0), 8: (-49, 0), 12: (-24, 0), 16: (9, 0), 20: (50, 0)}
new = [(i + 1, float(r["centre"][0]), float(r["centre"][1]), int(r["r_mm"])) for i, r in
       enumerate(sorted(PK, key=lambda r: -r["r_mm"]))]
if not Path(J.O4_H5).exists(): raise SystemExit(f"{J.O4_H5}: raw data not included in the release, see analysis/README.md")
out = {}
for label, rois in (("corrected_regions", new), ("published_regions", list(J.ROIS))):
    J.ROIS = rois
    res = J.run(J.O4_H5, label)
    for rank, cx, cy, ri in rois:
        tx, ty = TRUTH[ri]
        res["rois"][ri]["centre"] = [cx, cy]
        res["rois"][ri]["offset_mm"] = float(np.hypot(cx - tx, cy - ty))
    out[label] = res
json.dump(out, open(Path(__file__).resolve().parent / "brain_corrected_rois.json", "w"), indent=1, default=float)
print("wrote brain_corrected_rois.json")
