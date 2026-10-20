#!/usr/bin/env python3
"""Activity-ROI detection statistic on the JUNE (pre-correction) brain-tumor
triples, then the June-vs-O4 comparison on one axis.

DEFINITION, transcribed from the published pipeline
(2026-05-07_2140_phase3r3_batch1_clinical/run_roi_analysis_and_render.py):

    ROI       sphere, R_ROI = sqrt(r_inferred^2 + RECON_FWHM^2), RECON_FWHM = 7.0 mm,
              centred on the activity-map peak at z = 0
    background cylindrical annulus 80 < r_xy < 95 mm, |z| < 20 mm
    tau_hat   mean(tau | tau > 1 ns) - 1 ns
    sigma     tau_hat / sqrt(N)                       (analytic, manuscript form)
    n_sigma   |tau_ROI - tau_bg| / hypot(sigma_ROI, sigma_bg)

The ROI CENTRES are held fixed at the published activity peaks. Re-running the
MLEM activity reconstruction on each dataset would move the centres and mix a
localisation change into a lifetime comparison. Holding them fixed isolates the
lifetime term. Activity is I3-independent, so the centres are not expected to
move -- but that is a separate measurement, not an input to this one.

sigma_empirical = s/sqrt(N) is reported alongside sigma_analytic, as at W3.
"""
from __future__ import annotations
import os
import sys
import json
from pathlib import Path

import numpy as np

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from reconstruction.correct_travel_time import correct        # noqa: E402

W3 = RES / "2026-08-13_overnight/W3_detection_sigma"
JUNE_H5 = Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/W3_detection_sigma/june/triples_detector_brain_tumor_june_rcut75.h5"
O4_H5 = (Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/O4_rerun/"
                "brain_tumor/triples_detector_brain_tumor_rcut75.h5")

TAU_MIN = 1.0
RECON_FWHM = 7.0
# rank -> (peak_x, peak_y, r_inferred). Published activity peaks
ROIS = [(1, 54.0, -2.0, 20), (2, 10.0, -2.0, 16), (3, -26.0, -2.0, 12), (4, -50.0, -2.0, 8)]


def est(t):
    s = t[(t > TAU_MIN) & np.isfinite(t)]
    n = int(s.size)
    if n < 20:
        return None
    d = s - TAU_MIN
    th = float(d.mean())
    sd = float(d.std(ddof=1))
    return dict(N=n, tau=th, s=sd, sig_a=th / np.sqrt(n), sig_e=sd / np.sqrt(n))


def run(h5: Path, label: str):
    print(f"[{label}] correct_travel_time on {h5.name}", flush=True)
    tau, P = correct(h5)
    r_xy = np.hypot(P[:, 0], P[:, 1])
    bg = est(tau[(r_xy > 80) & (r_xy < 95) & (np.abs(P[:, 2]) < 20.0)])
    out = dict(n_triples=int(len(tau)), background=bg, rois={})
    print(f"[{label}] N_triples={len(tau):,}  bg tau={bg['tau']:.4f} N={bg['N']:,} "
          f"s/tau={bg['s']/bg['tau']:.3f}", flush=True)
    for rank, cx, cy, ri in ROIS:
        R = float(np.sqrt(ri ** 2 + RECON_FWHM ** 2))
        d2 = (P[:, 0] - cx) ** 2 + (P[:, 1] - cy) ** 2 + P[:, 2] ** 2
        a = est(tau[d2 <= R * R])
        dt = abs(a["tau"] - bg["tau"])
        sa = float(np.hypot(a["sig_a"], bg["sig_a"]))
        se = float(np.hypot(a["sig_e"], bg["sig_e"]))
        out["rois"][ri] = dict(R_ROI=round(R, 3), **{k: (round(v, 6) if isinstance(v, float) else v)
                                                     for k, v in a.items()},
                               nsigma_analytic=round(dt / sa, 4),
                               nsigma_empirical=round(dt / se, 4),
                               ratio_emp_over_ana=round(se / sa, 4))
        print(f"[{label}]  r={ri:2d} N={a['N']:>7,} tau={a['tau']:.4f} "
              f"nsig_a={dt/sa:7.3f} nsig_e={dt/se:7.3f}", flush=True)
    return out


def main():
    if not Path(JUNE_H5).exists(): raise SystemExit(f"{JUNE_H5}: raw data not included in the release, see analysis/README.md")
    res = {"june_precorrection": run(JUNE_H5, "june"),
           "o4_postcorrection": run(O4_H5, "o4")}

    print("\n=== June (pre) vs O4 (post), one axis, same lineage, same code path ===")
    print(f"{'r':>3} {'N_june':>8} {'N_O4':>8} {'ratio':>6} {'tau_june':>8} {'tau_O4':>7} "
          f"{'d tau %':>8} {'nsig_june':>9} {'nsig_O4':>8}")
    j, o = res["june_precorrection"], res["o4_postcorrection"]
    delta = {}
    for ri in (20, 16, 12, 8):
        a, b = j["rois"][ri], o["rois"][ri]
        dp = 100 * (b["tau"] - a["tau"]) / a["tau"]
        delta[ri] = dict(N_ratio=round(b["N"] / a["N"], 4), dtau_pct=round(dp, 3),
                         nsigma_june=a["nsigma_analytic"], nsigma_o4=b["nsigma_analytic"])
        print(f"{ri:>3} {a['N']:>8,} {b['N']:>8,} {b['N']/a['N']:6.3f} {a['tau']:8.4f} "
              f"{b['tau']:7.4f} {dp:+8.2f} {a['nsigma_analytic']:9.2f} {b['nsigma_analytic']:8.2f}")
    bj, bo = j["background"], o["background"]
    delta["background"] = dict(N_ratio=round(bo["N"] / bj["N"], 4),
                               dtau_pct=round(100 * (bo["tau"] - bj["tau"]) / bj["tau"], 3))
    print(f"{'bg':>3} {bj['N']:>8,} {bo['N']:>8,} {bo['N']/bj['N']:6.3f} {bj['tau']:8.4f} "
          f"{bo['tau']:7.4f} {100*(bo['tau']-bj['tau'])/bj['tau']:+8.2f}")
    print(f"\ntriples ratio O4/June: {o['n_triples']/j['n_triples']:.4f} "
          f"({j['n_triples']:,} -> {o['n_triples']:,})")

    res["comparison"] = delta
    res["definition"] = (
        "sphere R=sqrt(r_inf^2+7.0^2) mm at the published activity peak (z=0); "
        "bg annulus 80<r_xy<95 mm, |z|<20 mm; tau_hat=mean(tau|tau>1)-1; "
        "sigma=tau/sqrt(N); nsigma=|dtau|/hypot(sig_ROI,sig_bg)")
    res["axis"] = ("June truetime (pre-correction) vs O4 (post-correction). SAME macro "
                   "lineage and SAME code path. Carries the O1 3-gamma compilation fix "
                   "and the I3 convention correction TOGETHER; not split here.")
    res["r_min_mm"] = 8
    res["r_min_basis"] = ("smallest insert the activity localiser finds; r=4 mm has no "
                          "activity peak, so it has no ROI under this definition and the "
                          "floor is set by localisation, not by lifetime contrast")
    json.dump(res, open(W3 / "w3_activity_rois_june_vs_o4.json", "w"), indent=1)
    print(f"\nwrote {W3/'w3_activity_rois_june_vs_o4.json'}")


if __name__ == "__main__":
    main()
