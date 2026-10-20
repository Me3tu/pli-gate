"""Locked EJMP figure style.

v2, 2026-08-28 — the polish pass. Typography, geometry, palette and mark
conventions are set ONCE here so twelve figures change together and no
generator needs a per-figure hack.

NOTHING in this module changes a number. The palette keys keep their names and
their SEMANTIC assignments; only the hex values are tuned for print and screen:

    green  -> teal  #0F7173   ours / signal / 388 ps
    verm   -> coral #C8502D   0 ps / theirs-contrast
    blue   -> slate #3D5A80   component-fit / canonical
    orange -> gold  #E09F3E   auxiliary
    purple -> plum  #7B4B94   auxiliary
    grey           #8A8988    unresolved / de-emphasised

Import and call apply(); save with save_fig(fig, name).
Panel letters: panel_letter(ax, "a") puts a bold 11 pt letter OUTSIDE the axes,
top-left, so a panel's own descriptive title keeps its wording underneath it.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Semantic palette. Keys are historical (Wong); values are the v2 tuning.
WONG = {
    "black":  "#333333",
    "orange": "#E09F3E",   # gold, auxiliary
    "sky":    "#6FA8C7",   # auxiliary light
    "green":  "#0F7173",   # teal   — ours / signal / 388 ps
    "yellow": "#E6C34A",
    "blue":   "#3D5A80",   # slate  — component fit / canonical
    "verm":   "#C8502D",   # coral  — 0 ps / theirs-contrast
    "purple": "#7B4B94",   # plum, auxiliary
    "grey":   "#8A8988",   # unresolved / de-emphasised
}
GREY = WONG["grey"]
RULE = "#333333"           # fit / identity lines
REF = "#999999"            # reference / registered lines
BAND = "#000000"           # shaded bands, used at low alpha
BAND_ALPHA = 0.06

CYCLE = [WONG[k] for k in ("blue", "verm", "green", "orange",
                           "sky", "purple", "yellow", "black")]

MM = 1 / 25.4
FULL_W = 180 * MM   # double-column


def apply():
    plt.rcParams.update({
        # typography
        "font.family": "serif",
        "font.serif": ["STIXGeneral"],
        "mathtext.fontset": "stix",
        "font.size": 8.8,
        "axes.labelsize": 8.8,
        "axes.titlesize": 8.8,
        "legend.fontsize": 7.6,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        # geometry: left and bottom spines only, ticks outward
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        # grid, off by default; generators switch it on where it aids reading
        "grid.linewidth": 0.4,
        "grid.alpha": 0.22,
        "grid.color": "#000000",
        "axes.grid": False,
        "axes.axisbelow": True,
        # marks
        "axes.prop_cycle": plt.cycler(color=CYCLE),
        "lines.linewidth": 1.2,
        "lines.markersize": 5.0,
        "lines.markeredgewidth": 0.8,
        "lines.markeredgecolor": "white",
        "errorbar.capsize": 2.0,
        # legends
        "legend.frameon": False,
        "legend.handletextpad": 0.55,
        "legend.borderaxespad": 0.35,
        # output
        "figure.dpi": 110,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })


def panel_letter(ax, letter, dx_pt=-22.0, dy_pt=9.0):
    """Bold panel letter OUTSIDE the axes, top-left.

    The offset is in POINTS, not axes fractions, so the gap between the letter
    and a panel's descriptive title is identical whether the panel is wide or
    narrow. The title keeps its exact wording and sits to the right of the
    letter via ax.set_title(..., loc="left").
    """
    from matplotlib.transforms import ScaledTranslation
    t = ax.transAxes + ScaledTranslation(dx_pt / 72.0, dy_pt / 72.0,
                                         ax.figure.dpi_scale_trans)
    return ax.text(0.0, 1.0, f"({letter})", transform=t, fontsize=11,
                   fontweight="bold", va="bottom", ha="left", color=WONG["black"])


def band(ax, orientation="h", lo=None, hi=None, **kw):
    """A shaded band at the house alpha with a hairline edge."""
    f = ax.axhspan if orientation == "h" else ax.axvspan
    kw.setdefault("color", BAND)
    kw.setdefault("alpha", BAND_ALPHA)
    kw.setdefault("lw", 0.0)
    p = f(lo, hi, zorder=kw.pop("zorder", 1), **kw)
    return p


def save_fig(fig, stem):
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=300)
    print(f"saved {stem}.pdf/.png")
