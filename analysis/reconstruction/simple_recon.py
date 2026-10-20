"""
SIMPLE reconstruction (Huang et al. Commun Phys 2025) on parallelproj.

Reconstructs mean positron lifetime m_j = w_j / x_j via two MLEM passes.
Supports parameter sweeps over voxel size and iteration count.

Uses ONLY detector fields. No truth leakage.
"""
from __future__ import annotations
import argparse
import json
import os, h5py, numpy as np
from datetime import datetime
from pathlib import Path
from scipy.ndimage import gaussian_filter
import time

os.environ.setdefault("PARALLELPROJ_C_LIB", "/usr/local/lib/libparallelproj_c.so")
os.environ.setdefault("PARALLELPROJ_CUDA_LIB", "/usr/local/lib/libparallelproj_cuda.so")

import parallelproj

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- constants ----------
C_MM_NS = 299.792458
GATE_TIME_SCALE_NS = 1e9
EPS = 1e-9
W_FLOOR = 1e-6  # positivity floor for w after each iteration

# ---------- paths ----------
_SCRIPT_DIR = Path(__file__).resolve().parent
DET = _SCRIPT_DIR / "data" / "triples_detector_vipman_500M.h5"
OUT_BASE = _SCRIPT_DIR / "results" / "simple"

# ---------- TOF parameters (shared across all runs) ----------
# Siemens Biograph Vision Quadra coincidence timing resolution (CRT)
# FWHM = 247 ps. parallelproj.TOFParameters.sigma_tof expects σ (not FWHM):
# per parallelproj 1.4.0, the library's internal Gaussian uses σ = FWHM /
# (2·sqrt(2·ln2)) = FWHM / 2.355. Convert time→LOR-distance via c·T/2
# (the /2 accounts for TOF localization being half the path-length
# uncertainty from two-photon arrival times).
#   sigma_tof_mm = c · (FWHM_CRT / 2.355) / 2
CRT_FWHM_NS = 0.247
SIGMA_TOF_MM = CRT_FWHM_NS * C_MM_NS / 2 / 2.355   # ~15.72 mm
TOFBIN_WIDTH_MM = 13.0
N_TOF_BINS = 29
MAX_TOF_BIN = N_TOF_BINS // 2        # 14

# ---------- Gaussian post-smoothing ----------
# Fixed physical FWHM = 8 mm, sigma in mm, converted to voxels per run
SMOOTH_FWHM_MM = 8.0
SMOOTH_SIGMA_MM = SMOOTH_FWHM_MM / 2.3548  # ~3.397 mm

# ---------- physical extent (fixed across voxel sizes) ----------
PHYS_EXTENT_MM = (560.0, 280.0, 80.0)  # x, y, z


# ================================================================
# Data loading
# ================================================================

def load_data(det_path: Path, tau_lo: float, tau_hi: float) -> dict:
    """Load triples, apply τ window, compute TOF bins. Returns dict of arrays."""
    print(f"Loading triples from {det_path}...")
    t0 = time.time()
    with h5py.File(det_path, "r") as f:
        d = f["detector"]
        a1_x = d["a1_x"][:].astype(np.float32)
        a1_y = d["a1_y"][:].astype(np.float32)
        a1_z = d["a1_z"][:].astype(np.float32)
        a2_x = d["a2_x"][:].astype(np.float32)
        a2_y = d["a2_y"][:].astype(np.float32)
        a2_z = d["a2_z"][:].astype(np.float32)
        a1_time = d["a1_time"][:] * GATE_TIME_SCALE_NS
        a2_time = d["a2_time"][:] * GATE_TIME_SCALE_NS
        tau_all = f["derived"]["tau_corrected_ns"][:]

    N_total = len(tau_all)
    print(f"  {N_total:,} events loaded in {time.time()-t0:.1f}s")

    mask = (tau_all > tau_lo) & (tau_all < tau_hi)
    N_filt = mask.sum()
    print(f"  τ window [{tau_lo}, {tau_hi}] ns: {N_filt:,} / {N_total:,} "
          f"({100*mask.mean():.1f}%)")

    tau = tau_all[mask].astype(np.float32)
    p1 = np.stack([a1_x[mask], a1_y[mask], a1_z[mask]], axis=-1)
    p2 = np.stack([a2_x[mask], a2_y[mask], a2_z[mask]], axis=-1)

    tof_offset_mm = 0.5 * C_MM_NS * (a1_time[mask] - a2_time[mask])
    tof_bins = np.clip(
        np.round(tof_offset_mm / TOFBIN_WIDTH_MM).astype(np.int16),
        -MAX_TOF_BIN, MAX_TOF_BIN,
    )
    print(f"  TOF bins: min={tof_bins.min()}, max={tof_bins.max()}, "
          f"std={tof_bins.std():.2f}")
    print(f"  τ stats: mean={tau.mean():.3f}, median={np.median(tau):.3f}, "
          f"τ>0: {(tau>0).sum():,}/{len(tau):,} ({100*(tau>0).mean():.1f}%)")

    return dict(p1=p1, p2=p2, tof_bins=tof_bins, tau=tau)


# ================================================================
# Image geometry
# ================================================================

def compute_img_geometry(voxel_mm: float) -> tuple:
    """Return (img_shape, voxel_size) keeping physical extent fixed."""
    img_shape = tuple(int(round(ext / voxel_mm)) for ext in PHYS_EXTENT_MM)
    voxel_size = (voxel_mm, voxel_mm, voxel_mm)
    return img_shape, voxel_size


def load_mumap_to_recon_grid(mumap_mhd: Path,
                              img_shape: tuple,
                              voxel_mm: float,
                              phys_extent_mm: tuple) -> np.ndarray:
    """Phase 3d2: load a MetaImage µmap and resample onto the recon grid.

    Trilinear interpolation at recon voxel centres. Voxels outside the
    µmap bounding box → 0 (air/void: the recon FOV is larger than the
    cropped VIPman atlas, so the outside-atlas region is physical
    scanner-bore air. Documented in Methods.).

    Returns: float32 (nx, ny, nz) µmap in cm⁻¹, parallelproj-ordered,
    voxel centres aligned with the recon grid (centred at world origin).
    """
    from scipy.ndimage import map_coordinates
    with open(mumap_mhd) as f:
        h = {k.strip(): v.strip()
             for line in f if "=" in line
             for k, _, v in [line.partition("=")]}
    nx_mu, ny_mu, nz_mu = map(int, h["DimSize"].split())
    voxel_mu = float(h["ElementSize"].split()[0])
    offset_mu = np.array([float(x) for x in h["Offset"].split()],
                         dtype=np.float64)
    raw_path = mumap_mhd.parent / h["ElementDataFile"]
    arr = np.fromfile(raw_path, dtype=np.float32).reshape(nz_mu, ny_mu, nx_mu)

    half = np.array(phys_extent_mm, dtype=np.float64) / 2.0
    rxs = (np.arange(img_shape[0]) + 0.5) * voxel_mm - half[0]
    rys = (np.arange(img_shape[1]) + 0.5) * voxel_mm - half[1]
    rzs = (np.arange(img_shape[2]) + 0.5) * voxel_mm - half[2]

    # Recon centres → µmap voxel-index coords (centre of voxel-i is i+0.5
    # in mm, reverse: idx = (mm - origin)/voxel - 0.5 to land on voxel
    # CENTRE in numpy index convention).
    mux = (rxs - offset_mu[0]) / voxel_mu - 0.5
    muy = (rys - offset_mu[1]) / voxel_mu - 0.5
    muz = (rzs - offset_mu[2]) / voxel_mu - 0.5

    zz, yy, xx = np.meshgrid(muz, muy, mux, indexing="ij")
    sampled = map_coordinates(arr,
                               np.array([zz.ravel(), yy.ravel(), xx.ravel()]),
                               order=1, mode="constant", cval=0.0)
    sampled = sampled.reshape(img_shape[2], img_shape[1], img_shape[0])
    return sampled.transpose(2, 1, 0).astype(np.float32)


# ================================================================
# Single SIMPLE run
# ================================================================

def run_simple(
    data: dict,
    voxel_mm: float,
    n_iter: int,
    out_dir: Path,
    run_label: str,
    preflight: bool = False,
    sensitivity_path: Path | None = None,
    sigma_tof_mm_override: float | None = None,
    mumap_path: Path | None = None,
) -> dict:
    """Run full SIMPLE pipeline for one configuration.

    Returns dict with summary statistics.
    """
    sigma_tof_mm = sigma_tof_mm_override if sigma_tof_mm_override is not None \
        else SIGMA_TOF_MM
    out_dir.mkdir(parents=True, exist_ok=True)
    t_run_start = time.time()

    p1, p2 = data["p1"], data["p2"]
    tof_bins, tau = data["tof_bins"], data["tau"]
    img_shape, voxel_size = compute_img_geometry(voxel_mm)

    print(f"\n{'='*60}")
    print(f"RUN {run_label}: {voxel_mm}mm, {n_iter} iter, "
          f"img={img_shape}, N={len(tau):,}")
    print(f"{'='*60}")

    # ---- FOV coverage check ----
    img_half = np.array([s * voxel_mm / 2 for s in img_shape])
    with h5py.File(DET, "r") as f:
        # Reuse the τ mask implicitly: P_annih aligns with detector order
        # For coverage check, sample first 500k events (fast enough)
        P_sample = f["derived"]["P_annih_tof_mm"][:500_000]
    print(f"  FOV half-extents: x={img_half[0]:.0f}, y={img_half[1]:.0f}, "
          f"z={img_half[2]:.0f} mm")
    for i, ax in enumerate("xyz"):
        p1_pct = np.percentile(P_sample[:, i], 1)
        p99_pct = np.percentile(P_sample[:, i], 99)
        span = p99_pct - p1_pct
        cov = span / (2 * img_half[i]) * 100
        tag = " ***" if cov > 90 else ""
        print(f"  {ax}: 1-99%=[{p1_pct:.0f},{p99_pct:.0f}] span={span:.0f}mm "
              f"cov={cov:.0f}%{tag}")
    del P_sample

    # ---- build projector ----
    proj = parallelproj.ListmodePETProjector(
        event_start_coordinates=p1,
        event_end_coordinates=p2,
        img_shape=img_shape,
        voxel_size=voxel_size,
        img_origin=None,
    )
    tof_params = parallelproj.TOFParameters(
        num_tofbins=N_TOF_BINS,
        tofbin_width=TOFBIN_WIDTH_MM,
        sigma_tof=sigma_tof_mm,
        num_sigmas=3.0,
        tofcenter_offset=0.0,
    )
    if sigma_tof_mm_override is not None:
        print(f"  TOF σ override: {sigma_tof_mm:.4f} mm "
              f"(default {SIGMA_TOF_MM:.4f} mm)")
    proj.tof_parameters = tof_params
    proj.event_tofbins = tof_bins
    proj.tof = True

    # ---- Phase 3d2: optional attenuation correction (--mumap) ---------------
    # Compute per-event line integral of µ along each LOR using a 2nd
    # ListmodePETProjector instance with TOF disabled, then take exp(-∫ µ ds)
    # as the attenuation factor. Conversion: parallelproj returns line
    # integral in voxel_size units (mm here). µmap is cm⁻¹, so divide by 10.
    if mumap_path is not None:
        print(f"\n  loading µmap from {mumap_path}")
        mu_recon = load_mumap_to_recon_grid(mumap_path, img_shape,
                                              voxel_mm, PHYS_EXTENT_MM)
        print(f"  µmap recon-grid: shape={mu_recon.shape}  "
              f"min={mu_recon.min():.4f}  max={mu_recon.max():.4f}  "
              f"mean={mu_recon.mean():.4f} cm⁻¹")
        proj_atten = parallelproj.ListmodePETProjector(
            event_start_coordinates=p1,
            event_end_coordinates=p2,
            img_shape=img_shape,
            voxel_size=voxel_size,
            img_origin=None,
        )
        mu_integral_dimensionless = proj_atten(mu_recon) / 10.0
        atten_LOR = np.exp(-mu_integral_dimensionless).astype(np.float32)
        # Sanity report (Phase 3d2 check before validation).
        a_p5, a_p50, a_p95 = np.percentile(atten_LOR, [5, 50, 95])
        print(f"  attenuation per event: min={atten_LOR.min():.4f}  "
              f"p5={a_p5:.4f}  median={a_p50:.4f}  p95={a_p95:.4f}  "
              f"max={atten_LOR.max():.4f}  mean={atten_LOR.mean():.4f}")
        # 10-bin histogram (atten ∈ [0, 1])
        h_counts, h_edges = np.histogram(atten_LOR, bins=10, range=(0.0, 1.0))
        print(f"  atten_LOR 10-bin histogram (bin edges [0..1]):")
        for i in range(10):
            print(f"    [{h_edges[i]:.2f}, {h_edges[i+1]:.2f}]:  "
                  f"{int(h_counts[i]):>10,} ({h_counts[i]/len(atten_LOR)*100:5.2f} %)")
    else:
        atten_LOR = None
    # -------------------------------------------------------------------------

    if sensitivity_path is not None:
        sensitivity = np.load(sensitivity_path).astype(np.float32)
        if sensitivity.shape != tuple(img_shape):
            raise ValueError(
                f"sensitivity shape {sensitivity.shape} != img_shape "
                f"{tuple(img_shape)}"
            )
        print(f"  sensitivity loaded: {sensitivity_path}  "
              f"min={sensitivity.min():.3f} max={sensitivity.max():.3f} "
              f"mean={sensitivity.mean():.3f}")
    else:
        sensitivity = np.ones(img_shape, dtype=np.float32)

    # ---- optional pre-flight (first run only) ----
    if preflight:
        print("\n  Pre-flight: 1-iter activity check...")
        x_pre = np.ones(img_shape, dtype=np.float32)
        y_bar = proj(x_pre)
        x_pre = x_pre * proj.adjoint(1.0 / np.maximum(y_bar, EPS))
        roi_mid = tuple(slice(s//4, 3*s//4) for s in img_shape)
        roi = x_pre[roi_mid]
        contrast = roi.max() / max(roi.mean(), EPS)
        print(f"  Pre-flight ROI contrast: {contrast:.1f} "
              f"({'OK' if contrast > 1.5 else 'LOW — check TOF sign'})")

    # ---- MLEM: activity x ----
    print(f"\n  MLEM start: {datetime.now().isoformat()}")
    print(f"\n  Activity MLEM ({n_iter} iter):")
    x = np.ones(img_shape, dtype=np.float32)
    x_iter_snapshots = []
    for it in range(n_iter):
        t0 = time.time()
        y_bar = proj(x)
        if atten_LOR is None:
            x = x * proj.adjoint(1.0 / np.maximum(y_bar, EPS)) / np.maximum(sensitivity, EPS)
        else:
            x = x * proj.adjoint(1.0 / atten_LOR / np.maximum(y_bar, EPS)) / np.maximum(sensitivity, EPS)
        x_iter_snapshots.append(x.copy())
        print(f"    iter {it+1:2d}/{n_iter}: {time.time()-t0:.1f}s  "
              f"sum={x.sum():.3e}  max={x.max():.3e}")

    # ---- MLEM: lifetime-weighted w (raw signed τ, floor-clamped) ----
    print(f"\n  Lifetime-weighted MLEM ({n_iter} iter, raw signed τ):")
    w = np.ones(img_shape, dtype=np.float32) * max(float(tau.mean()), 0.1)

    # Phase 3n: µmap-based support is REQUIRED. Replaces the previous
    # contrast-dependent ``x > 0.1 * x.max()`` mask which collapsed onto
    # bright voxels in high-contrast phantoms and pulled in noise-floor
    # in low-contrast patient recons. The µmap threshold 0.005 cm⁻¹
    # separates air (~0.0001) from any tissue. The sensitivity-based
    # fallback explored in CHECKPOINT 1B was rejected as degenerate
    # for production sensitivity files (uniform across FOV).
    if mumap_path is None:
        raise ValueError(
            "--mumap is required for the support mask. "
            "Use reconstruction/build_phantom_mumap.py to synthesise one "
            "for phantom geometries (Sim 1, Sim 2). "
            "For VIPman recons, use the harmonized atlas µmap."
        )
    support = mu_recon > 0.005
    iter_log_w = []

    for it in range(n_iter):
        t0 = time.time()
        y_bar = proj(w)
        ratio = tau / np.maximum(y_bar, EPS)
        if atten_LOR is not None:
            ratio = ratio / atten_LOR
        w = w * proj.adjoint(ratio) / np.maximum(sensitivity, EPS)
        w = np.maximum(w, W_FLOOR)
        elapsed = time.time() - t0

        m_check = w / np.maximum(x, 1e-3)
        m_vals = m_check[support]
        p10, p50, p90 = np.percentile(m_vals, [10, 50, 90])
        iter_log_w.append((float(p10), float(p50), float(p90)))
        print(f"    iter {it+1:2d}/{n_iter}: {elapsed:.1f}s  "
              f"sum_w={w.sum():.3e}  "
              f"m[10/50/90%]=[{p10:.3f}, {p50:.3f}, {p90:.3f}] ns")
    print(f"\n  MLEM end: {datetime.now().isoformat()}")

    # ---- per-iter activity percentiles within final support (retrospective) ----
    iter_log_x = []
    for x_k in x_iter_snapshots:
        vals = x_k[support]
        p10, p50, p90 = np.percentile(vals, [10, 50, 90])
        iter_log_x.append((float(p10), float(p50), float(p90)))
    iter_log = iter_log_w  # preserve legacy key for sweep main()

    # ---- write iter_log.json sidecar ----
    iter_log_payload = {
        "x_pass": [
            {"iter": i + 1, "p10": p10, "p50": p50, "p90": p90}
            for i, (p10, p50, p90) in enumerate(iter_log_x)
        ],
        "w_pass": [
            {"iter": i + 1, "p10": p10, "p50": p50, "p90": p90}
            for i, (p10, p50, p90) in enumerate(iter_log_w)
        ],
        "support_mask_threshold": "mu_recon > 0.005 cm-1 (Phase 3n; --mumap required)",
    }
    with open(out_dir / "iter_log.json", "w") as f:
        json.dump(iter_log_payload, f, indent=2)

    # ---- compute raw m ----
    m_raw = w / np.maximum(x, 1e-3)
    m_raw_masked = np.where(support, m_raw, 0.0)

    # ---- Gaussian post-smoothing (fixed physical FWHM) ----
    sigma_vox = SMOOTH_SIGMA_MM / voxel_mm
    print(f"\n  Gaussian smoothing: FWHM={SMOOTH_FWHM_MM}mm, "
          f"σ={SMOOTH_SIGMA_MM:.2f}mm = {sigma_vox:.3f} voxels")
    x_smooth = gaussian_filter(x, sigma=sigma_vox)
    w_smooth = gaussian_filter(w, sigma=sigma_vox)
    m_smooth = w_smooth / np.maximum(x_smooth, 1e-3)
    m_smooth_masked = np.where(support, m_smooth, 0.0)

    # ---- statistics ----
    wall_time = time.time() - t_run_start

    def stats_block(label, arr):
        vals = arr[support]
        print(f"\n  {label} (within support, {support.sum():,} voxels):")
        print(f"    min={vals.min():.3f}  10th={np.percentile(vals,10):.3f}  "
              f"median={np.median(vals):.3f}  "
              f"90th={np.percentile(vals,90):.3f}  max={vals.max():.3f} ns")
        print(f"    mean={vals.mean():.3f}  std={vals.std():.3f} ns")
        return dict(
            min=float(vals.min()), p10=float(np.percentile(vals, 10)),
            median=float(np.median(vals)), p90=float(np.percentile(vals, 90)),
            max=float(vals.max()), mean=float(vals.mean()), std=float(vals.std()),
        )

    raw_stats = stats_block("m (raw)", m_raw_masked)
    smooth_stats = stats_block("m (smoothed)", m_smooth_masked)

    # ---- secondary "unmasked" percentiles (Phase 3n: now consolidated
    # with the primary support mask. Both bugged the same way under
    # the old contrast-dependent thresholds. Under Option B they share
    # the same anatomically-meaningful mask) ----
    m_all = m_raw[support]
    ms_all = m_smooth[support]
    print(f"\n  m within support ({len(m_all):,} voxels):")
    print(f"    raw:      median={np.median(m_all):.3f}  "
          f"[10,90]=[{np.percentile(m_all,10):.3f}, {np.percentile(m_all,90):.3f}]")
    print(f"    smoothed: median={np.median(ms_all):.3f}  "
          f"[10,90]=[{np.percentile(ms_all,10):.3f}, {np.percentile(ms_all,90):.3f}]")

    print(f"\n  Wall time: {wall_time:.1f}s")

    # ---- save arrays ----
    np.save(out_dir / "activity_x.npy", x)
    np.save(out_dir / "lifetime_weighted_w.npy", w)
    np.save(out_dir / "mean_lifetime_m.npy", m_raw_masked)
    np.save(out_dir / "m_smoothed.npy", m_smooth_masked)
    np.save(out_dir / "support_mask.npy", support)

    # ---- per-run orthoslice figure ----
    mid = [s // 2 for s in img_shape]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    vmax_act = np.percentile(x, 99.5)
    for c, (slc, lbl) in enumerate([
        (x[:, :, mid[2]].T, "axial"), (x[:, mid[1], :].T, "coronal"),
        (x[mid[0], :, :].T, "sagittal"),
    ]):
        ax = axes[0, c]
        im = ax.imshow(slc, origin="lower", cmap="inferno", vmin=0, vmax=vmax_act)
        ax.set_title(f"activity — {lbl}", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.8)

    for c, (slc, lbl) in enumerate([
        (m_smooth_masked[:, :, mid[2]].T, "axial"),
        (m_smooth_masked[:, mid[1], :].T, "coronal"),
        (m_smooth_masked[mid[0], :, :].T, "sagittal"),
    ]):
        ax = axes[1, c]
        im = ax.imshow(slc, origin="lower", cmap="viridis", vmin=0.0, vmax=2.0)
        ax.set_title(f"m smoothed — {lbl}", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.8, label="ns")

    fig.suptitle(f"Run {run_label}: {voxel_mm}mm, {n_iter} iter "
                 f"(smoothed FWHM={SMOOTH_FWHM_MM}mm)",
                 fontsize=12, weight="bold")
    plt.tight_layout()
    plt.savefig(out_dir / "fig_orthoslices.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ---- per-run m histogram ----
    fig, ax = plt.subplots(figsize=(7, 4))
    m_support = m_smooth_masked[support]
    ax.hist(m_support, bins=150, range=(0, 3), color="#3c6e9d",
            edgecolor="black", linewidth=0.3, alpha=0.8)
    ax.set_yscale("log")
    for ref, lbl, col in [(0.8, "water", "#2196F3"), (1.0, "muscle", "#FF9800"),
                           (1.7, "adipose", "#F44336")]:
        ax.axvline(ref, color=col, linestyle="--", alpha=0.7, linewidth=1.5, label=lbl)
    ax.legend(fontsize=9)
    ax.set_xlabel("mean lifetime m (ns)")
    ax.set_ylabel("voxel count")
    ax.set_title(f"Run {run_label}: m histogram (smoothed, within support)")
    ax.set_xlim(0, 3)
    plt.tight_layout()
    plt.savefig(out_dir / "fig_m_hist.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"  Figures saved to {out_dir}/")

    return dict(
        run_label=run_label, voxel_mm=voxel_mm, n_iter=n_iter,
        img_shape=img_shape, wall_time=wall_time,
        raw_stats=raw_stats, smooth_stats=smooth_stats,
        iter_log=iter_log, iter_log_x=iter_log_x, iter_log_w=iter_log_w,
        out_dir=str(out_dir),
    )


# ================================================================
# Main: load once, run sweep
# ================================================================

def run_single(args: argparse.Namespace) -> None:
    """Single-run mode (used by run_C_geosens.py wrapper)."""
    tau_lo, tau_hi = args.tau_range
    det_path = Path(args.triples_override) if args.triples_override else DET
    data = load_data(det_path, tau_lo, tau_hi)
    out_dir = Path(args.outdir)
    summary = run_simple(
        data=data,
        voxel_mm=args.voxel,
        n_iter=args.iters,
        out_dir=out_dir,
        run_label=args.label,
        preflight=args.preflight,
        sensitivity_path=args.sensitivity,
        sigma_tof_mm_override=args.tof_sigma_mm,
        mumap_path=args.mumap,
    )
    summary["tau_window"] = (tau_lo, tau_hi)
    print(f"\n  Summary written; activity+m saved to {out_dir}")
    return summary


def main():
    print("=" * 60)
    print("SIMPLE reconstruction — Phase 1.5 parameter sweep")
    print("=" * 60)
    print(f"Gaussian smoothing: FWHM={SMOOTH_FWHM_MM}mm (physical), "
          f"σ={SMOOTH_SIGMA_MM:.2f}mm")

    # Define sweep configurations
    runs = [
        ("A", 4.0, 3,  -1.0, 15.0),
        ("B", 4.0, 10, -1.0, 15.0),
        ("C", 8.0, 3,  -1.0, 15.0),
        ("Cp", 8.0, 3,  0.0, 10.0),   # C': tighter τ window
        ("D", 8.0, 10, -1.0, 15.0),
    ]

    # Load data for each unique τ window
    data_cache = {}
    for label, voxel_mm, n_iter, tau_lo, tau_hi in runs:
        key = (tau_lo, tau_hi)
        if key not in data_cache:
            data_cache[key] = load_data(DET, tau_lo, tau_hi)

    # Run sweep
    results = []
    for i, (label, voxel_mm, n_iter, tau_lo, tau_hi) in enumerate(runs):
        out_dir = OUT_BASE / f"run_{label}_{int(voxel_mm)}mm_{n_iter}iter"
        data = data_cache[(tau_lo, tau_hi)]
        summary = run_simple(
            data, voxel_mm, n_iter, out_dir,
            run_label=label,
            preflight=(i == 0),
        )
        summary["tau_window"] = (tau_lo, tau_hi)
        results.append(summary)

    # ---- summary table ----
    print(f"\n{'='*60}")
    print("SWEEP SUMMARY")
    print(f"{'='*60}")
    print(f"{'Run':>4s}  {'vox':>4s}  {'iter':>4s}  {'τ win':>10s}  "
          f"{'m_sm med':>8s}  {'m_sm [10,90]':>18s}  {'time':>6s}")
    print("-" * 70)
    for r in results:
        s = r["smooth_stats"]
        tw = r["tau_window"]
        print(f"  {r['run_label']:>3s}  {r['voxel_mm']:4.0f}  {r['n_iter']:4d}  "
              f"[{tw[0]:4.1f},{tw[1]:4.1f}]  "
              f"{s['median']:8.3f}  "
              f"[{s['p10']:.3f}, {s['p90']:.3f}]  "
              f"{r['wall_time']:5.1f}s")

    # ---- C vs C' comparison ----
    c_res = next(r for r in results if r["run_label"] == "C")
    cp_res = next(r for r in results if r["run_label"] == "Cp")
    c_med = c_res["smooth_stats"]["median"]
    cp_med = cp_res["smooth_stats"]["median"]
    diff_pct = abs(c_med - cp_med) / c_med * 100
    print(f"\n  C vs C' (τ window effect):")
    print(f"    C  [-1,15]: m_smooth median = {c_med:.3f} ns")
    print(f"    C' [ 0,10]: m_smooth median = {cp_med:.3f} ns")
    print(f"    Difference: {diff_pct:.1f}%")
    if diff_pct < 5:
        print(f"    → Negative-τ events are noise-dominated; "
              f"window can be tightened.")
    else:
        print(f"    → Negative-τ events carry real information; "
              f"keep wide window.")

    print(f"\nAll done. Outputs in {OUT_BASE}/")


def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("single", help="Run a single SIMPLE configuration.")
    s.add_argument("--label", default="C")
    s.add_argument("--voxel", type=float, default=8.0)
    s.add_argument("--iters", type=int, default=3)
    s.add_argument("--tau-range", type=float, nargs=2, default=[-1.0, 15.0],
                   metavar=("LO", "HI"))
    s.add_argument("--outdir", required=True)
    s.add_argument("--sensitivity", type=Path, default=None,
                   help="Path to precomputed sensitivity .npy. "
                        "If omitted, uniform sensitivity is used.")
    s.add_argument("--mumap", type=Path, default=None,
                   help="Path to harmonized µmap MHD for attenuation "
                        "correction. If omitted, uniform attenuation = 1 "
                        "(bit-identical to current behaviour).")
    s.add_argument("--preflight", action="store_true")
    s.add_argument("--triples-override", type=str, default=None,
                   help="Override the triples HDF5 path (same schema as DET). "
                        "Used for synthetic forward-projection tests.")
    s.add_argument("--tof-sigma-mm", type=float, default=None,
                   help="Override SIGMA_TOF_MM passed to parallelproj "
                        "TOFParameters. Used for TOF σ ablation sweeps.")
    sub.add_parser("sweep", help="Run the full 5-run sweep (default).")
    args = ap.parse_args()
    if args.cmd == "single":
        run_single(args)
    else:
        main()


if __name__ == "__main__":
    _cli()
