"""Shared setup for the figure generators (referee round 2, 2026-09-30, and feedback round 3,
2026-09-30 evening).

Unified visual language:
  THIS  filled blue      this work / all-prompts / primary analysis
  OPEN  open blue        closest-prompt matching / alternative or prespecified analysis
  REF   black            reference or configured value
  CTRL  grey             control / unresolved / secondary
  ACC   teal diamond     acceptance-aware readout (third arm where needed)
Timing models are separated by marker shape (0 ps squares, 388 ps triangles),
never by color alone, so every figure reads in grayscale.

Round 3: figures are drawn at the IOP double-column width (W = 150 mm), and every text
element is 8 pt or larger at that size (IOP: 8 to 12 pt at final size, 15 cm wide). A panel may
carry ONE key number where the conclusion rests on it. Every such number is read
from the same JSON as the text that quotes it and is asserted in the built PDF.
No interpretive sentences are drawn into a figure.
"""
import os
from pathlib import Path
# Code release v1.1.0: inputs default to the shipped source_data/ next to this file (dated result-folder
# names kept), outputs to out/. Both can be redirected with PLI_RES and PLI_FIGOUT.
HERE = Path(__file__).resolve().parent
RES = Path(os.environ.get('PLI_RES', str(HERE / 'source_data')))
FS = RES / '2026-09-24_jp_a3_figures/figsrc'
MF = FS / '2026-05-11_paper_PMB/manuscript_ejmp/figures'
R2 = RES / '2026-09-30_referee2'
F3 = RES / '2026-09-30_feedback3'
import ejmp_style as st  # noqa: E402  (shipped next to this file, was results/2026-08-11_ejmp_figure_regen/ejmp_style.py)
st.apply()
import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

THIS = st.WONG['blue']
REF = '#1a1a1a'
CTRL = st.WONG['grey']
ACC = st.WONG['green']
MM = st.MM
W = 150 * MM          # IOP double-column width
SMALL = 8.0           # smallest text in any figure, at printed size
KEYC = '#333333'      # color of in-panel key numbers and short labels
plt.rcParams.update({'font.size': 9.0, 'axes.labelsize': 9.0, 'axes.titlesize': 9.0, 'legend.fontsize': 8.0,
                     'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'savefig.pad_inches': 0.0,
                     'savefig.bbox': None})
OUT = Path(os.environ.get('PLI_FIGOUT', str(HERE / 'out')))
OUT.mkdir(parents=True, exist_ok=True)

NAME = {'T1_water_shibuya2020_n2': 'water (N$_2$)', 'T2_adipose_avachat2024': 'adipose',
        'T3_adipose_moskal2023': 'mediastinal adipose', 'T4_hepatic_avachat2024': 'hepatic',
        'T5_muscle_avachat2024': 'muscle', 'T6_brain_healthy_moskal2024': 'brain (healthy)',
        'T7_brain_tumor_moskal2024': 'brain (tumor)', 'T8_myxoma_moskal2023': 'cardiac myxoma',
        'quartz': 'quartz disk', 'polycarbonate': 'polycarbonate disk',
        'VOI_5x15': 'water VOI 5 mm$\\times$15 cm', 'VOI_10x20': 'water VOI 10 mm$\\times$20 cm',
        'VOI_15x25': 'water VOI 15 mm$\\times$25 cm',
        'sphere_37mm': 'sphere 37 mm', 'sphere_28mm': 'sphere 28 mm', 'sphere_22mm': 'sphere 22 mm',
        'sphere_17mm': 'sphere 17 mm', 'sphere_13mm': 'sphere 13 mm', 'sphere_10mm': 'sphere 10 mm',
        'background_central': 'NEMA background'}
NAME_SHORT = dict(NAME, T3_adipose_moskal2023='adipose (mediastinal)', T1_water_shibuya2020_n2='water (N$_2$)')


def sg(v, f='+.1f'):
    """Signed number with a typographic minus."""
    return format(v, f).replace('-', '\u2212')


def key(ax, x, y, s, transform=None, **kw):
    """A short key number or label inside a panel (axes coordinates by default)."""
    kw.setdefault('fontsize', 8.0); kw.setdefault('color', KEYC)
    return ax.text(x, y, s, transform=ax.transAxes if transform is None else transform, **kw)


def mm_fig(w_mm, h_mm):
    """A figure of exact size. Axes are placed in millimetres with ax_mm."""
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM))
    fig._mm = (w_mm, h_mm)
    return fig


def ax_mm(fig, x0, y0, w, h, **kw):
    fw, fh = fig._mm
    return fig.add_axes([x0 / fw, y0 / fh, w / fw, h / fh], **kw)


def save(fig, stem):
    fig.savefig(OUT / f'{stem}.pdf')
    fig.savefig(OUT / f'{stem}.png', dpi=200)
    print('saved', stem)
