"""Phase 3r-2 PART B — publication-quality matplotlib styling baseline.

Targets PMB / EJNMMI Physics / IEEE-TMI house-style envelope:
- Vector PDF + 300-DPI PNG
- Serif typography via STIX (built into matplotlib; no system-font dep)
- Color-blind-safe Wong palette (Wong 2011, Nature Methods)
- Tight axis layout, minimal grid, clean ticks
- Single-column 88 mm / double-column 180 mm targets

Helper functions:
    setup_publication_rcparams()    — apply global rcParams
    save_figure(fig, path, ext)     — save vector + raster with metadata
    tight_layout_with_margins(fig)  — multi-panel-safe tight layout
    figure_dims(width, height)      — convert mm → inches for figsize
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

# ── Color palettes (color-blind-safe) ─────────────────────────────────────────

# Wong 2011 (Nature Methods 8, 441) — 7-color color-blind-safe palette.
WONG = {
    "black":         "#000000",
    "orange":        "#E69F00",
    "sky_blue":      "#56B4E9",
    "bluish_green":  "#009E73",
    "yellow":        "#F0E442",
    "blue":          "#0072B2",
    "vermillion":    "#D55E00",
    "reddish_purple":"#CC79A7",
}

# Paul Tol "muted" palette — alternative; cleaner for many-category plots
TOL_MUTED = {
    "indigo":   "#332288",
    "cyan":     "#88CCEE",
    "teal":     "#44AA99",
    "green":    "#117733",
    "olive":    "#999933",
    "sand":     "#DDCC77",
    "rose":     "#CC6677",
    "wine":     "#882255",
    "purple":   "#AA4499",
}

# Default palette choice for prototypes: Wong (better contrast on bar charts).
PALETTE = WONG


# ── Figure dimensions ─────────────────────────────────────────────────────────

# PMB / EJNMMI Phys / IEEE-TMI envelope (mm widths)
COLUMN_WIDTH_MM = {
    "single": 88.0,    # PMB single-col
    "1.5":    120.0,   # 1.5-col where allowed
    "double": 180.0,   # PMB double-col / IEEE 2-col
}


def figure_dims(width_mm: float, height_mm: float) -> tuple[float, float]:
    """Convert mm to inches for plt.figure(figsize=...)."""
    MM_PER_INCH = 25.4
    return (width_mm / MM_PER_INCH, height_mm / MM_PER_INCH)


# ── rcParams ──────────────────────────────────────────────────────────────────

def setup_publication_rcparams() -> None:
    """Apply publication-quality rcParams globally.

    Targets:
    - Serif font (STIX, ships with matplotlib; no system-font dep)
    - Math font matched to text (mathtext.fontset='stix')
    - 9 pt body / 10 pt title (within 8-11 pt journal envelope)
    - Vector-friendly defaults (PDF metadata, embedded fonts)
    """
    matplotlib.rcParams.update({
        # Typography
        "font.family":        "serif",
        "font.serif":         ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset":   "stix",
        "font.size":          9.0,        # body
        "axes.titlesize":     10.0,       # subplot title
        "axes.labelsize":     9.0,
        "xtick.labelsize":    8.0,
        "ytick.labelsize":    8.0,
        "legend.fontsize":    8.0,
        "figure.titlesize":   11.0,       # suptitle

        # Axes
        "axes.linewidth":     0.8,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.edgecolor":     "#222222",
        "axes.labelcolor":    "#222222",

        # Ticks
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        "xtick.major.size":   3.0,
        "ytick.major.size":   3.0,
        "xtick.major.width":  0.7,
        "ytick.major.width":  0.7,
        "xtick.color":        "#222222",
        "ytick.color":        "#222222",

        # Grid (off by default; opt-in per axis)
        "axes.grid":          False,
        "grid.linewidth":     0.4,
        "grid.color":         "#cccccc",
        "grid.alpha":         0.6,

        # Lines + scatter
        "lines.linewidth":    1.2,
        "lines.markersize":   4.0,
        "patch.linewidth":    0.5,
        "patch.edgecolor":    "#222222",

        # Legend
        "legend.frameon":     False,
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.borderpad":   0.3,
        "legend.columnspacing": 1.0,

        # Save defaults (override in save_figure)
        "savefig.dpi":        300,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype":       42,         # embed TrueType (Acrobat-friendly)
        "ps.fonttype":        42,
    })


def save_figure(fig, out_stem: Path, formats=("pdf", "png"),
                metadata: dict | None = None) -> list[Path]:
    """Save figure to multiple formats. Returns list of written paths."""
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    md = {"Creator": "matplotlib", "Producer": "phase3r2_baseline"}
    if metadata:
        md.update(metadata)
    for ext in formats:
        p = out_stem.with_suffix(f".{ext}")
        if ext == "pdf":
            fig.savefig(p, format="pdf", metadata={k: str(v) for k, v in md.items()})
        elif ext == "png":
            fig.savefig(p, format="png", dpi=300)
        else:
            fig.savefig(p, format=ext)
        written.append(p)
    return written


def tight_layout_with_margins(fig, *, top=0.96, bottom=0.10,
                                left=0.10, right=0.98) -> None:
    """tight_layout that handles edge cases for multi-panel figures
    by clamping margins to safe envelope first.
    """
    fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right,
                          hspace=0.30, wspace=0.30)


# ── Demo / self-test ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    setup_publication_rcparams()
    print("Wong palette:", list(PALETTE.keys()))
    print("Single-col 88 mm × 70 mm = inches:", figure_dims(88, 70))
    print("Double-col 180 mm × 110 mm = inches:", figure_dims(180, 110))
    print("\nrcParams applied; helper functions ready.")
