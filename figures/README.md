# Figures

Scripts and source data that regenerate all 14 figures (6 main, 8 supplementary) of

> D. Dhar and J. Peter, "Material-dependent positronium lifetimes in GATE and prompt-selection bias" (submitted to Physics in Medicine and Biology).

Part of pli-gate v1.1.0. The figures are drawn from processed result files (JSON, NPY, NPZ, CSV) shipped in `source_data/`. The simulations and reconstructions that produced these files are not rerun here.

## Quick start

```bash
cd figures
python3 -m pip install -r requirements.txt     # pip wheels, see "Environment" below
./make_figures.sh                              # writes out/
./make_figures.sh /some/other/outdir           # or any output folder
PYTHON=/path/to/python3 ./make_figures.sh      # choose the interpreter
```

`make_figures.sh` runs the 12 Python generators, then compiles `fig1_methods.tex` with pdflatex in a temporary folder. It takes about one minute (most of it the 300 shard-bootstrap resamples of `f08_pedestal.py`). It writes `<stem>.pdf` and `<stem>.png` (200 dpi) for every matplotlib figure, `fig2_closure_benchmark.json` (the numbers quoted with Fig. 2), `fig1_methods.pdf`, and one log per generator in `<outdir>/logs/` with the key numbers each script prints. It always reads `source_data/` next to the script: it sets `PLI_RES` and `PLI_FIGOUT` and unsets `PLI_INFLIGHT_JSON`.

To run a single generator: `cd figures && python3 f07_oxygen.py`. Without environment variables the scripts read `source_data/` and write `out/` (both next to `common.py`). `PLI_RES` and `PLI_FIGOUT` redirect input and output.

Check the files with `sha256sum -c SHA256SUMS` (from inside `figures/`).

## Contents

| file | purpose |
|---|---|
| `common.py` | shared paths, style, palette and helpers, imported by every generator |
| `ejmp_style.py` | house style (fonts, sizes, colour palette), unchanged from `reconstruction/results/2026-08-11_ejmp_figure_regen/` |
| `f02b_closure_benchmark.py` ... `fS_benchspec.py` | the 12 generators (table below) |
| `fig1_methods.tex` | Fig. 1, standalone TikZ |
| `source_data/` | the 36 data files the generators read, plus `make_vipman_slice.py` (the generator of the Fig. S6(a) input) |
| `make_figures.sh` | runs everything |
| `requirements.txt` | Python package versions |
| `EDITS.md` | every change to the scripts relative to the versions that produced the paper's figures |
| `SHA256SUMS` | sha256 of every other file in this folder |

## Figure, script and inputs

Paths of inputs are relative to `source_data/`. Roles: plotted (reaches the figure), print only (feeds a printed comparison only), assert only (feeds a consistency assertion only). The docstrings of the scripts use the figure numbers of earlier manuscript versions. The output stem identifies the figure.

| paper | output stem | generator | inputs | role |
|---|---|---|---|---|
| Fig. 1 | `fig1_methods` | `fig1_methods.tex` (pdflatex, TikZ) | none | |
| Fig. 2 | `fig2_closure_benchmark` (+ `.json`) | `f02b_closure_benchmark.py` | `2026-09-24_jp_a3_figures/figsrc/2026-05-11_paper_PMB/manuscript_ejmp/figures/fig2_recipe_validation.json` | plotted |
| | | | `2026-10-08_replication_330/repl_sig0.3306.json` | plotted |
| Fig. 3 | `fig6_clinical_mechanism_v2` | `f09_mechanism.py` | `2026-10-14_known_randoms/crt_sweep_c0_summary.json` | plotted |
| | | | `2026-09-30_feedback3/c0_sim1/acceptance_tm_c0.json` | plotted |
| | | | `2026-09-24_jp_a3_figures/figsrc/fig6data/crt_sweep_v3_all.csv` | print only |
| | | | `2026-09-30_referee2/acceptance_tm.json` | print only |
| Fig. 4 | `fig9_pedestal_derivation` | `f08_pedestal.py` | `2026-09-30_feedback3/c0_sim1/hists_c0_0ps.npz` | plotted |
| | | | `2026-10-08_signed_eq3/signed_eq3.json` | plotted |
| | | | `2026-10-14_known_randoms/signed_eq3_c0.json` | plotted |
| Fig. 5 | `fig11_realistic_forest` | `f11_realistic.py` | `2026-09-24_jp_a3_figures/figsrc/2026-05-11_paper_PMB/manuscript_ejmp/figures/fig11_realistic_forest.json` | plotted |
| | | | `2026-09-30_feedback3/estimator/comp_joint_ratio.json` | plotted |
| | | | `2026-09-30_feedback3/c0_sim1/c0_compare.json` | plotted |
| | | | `2026-09-24_jp_a3_figures/figsrc/fig11_component_shift.json` | assert only |
| Fig. 6 | `fig10_rate_law` | `f10_rate.py` | `2026-09-24_jp_a3_figures/figsrc/2026-05-11_paper_PMB/manuscript_ejmp/figures/fig10_rate_law.json` | plotted |
| | | | `2026-10-08_signed_eq3/ladder_band.json` | plotted |
| | | | `2026-10-08_c0_ladder/ladder_boot_c0.json` | plotted |
| | | | `2026-10-08_signed_eq3/signed_eq3.json` | plotted |
| | | | `2026-09-30_feedback3/c0_sim1/c0_compare.json` | plotted |
| | | | `2026-09-24_jp_a3_figures/figsrc/2026-05-11_paper_PMB/manuscript_ejmp/figures/fig11_realistic_forest.json` | plotted |
| | | | `2026-10-14_known_randoms/known_randoms_summary.json` | plotted |
| Fig. S1 | `figS_benchspec` | `fS_benchspec.py` | `2026-10-08_replication_330/repl_sig0.3306.json` | plotted |
| Fig. S2 | `fig_r2_transport_novelty` | `f04_transport.py` | `2026-09-24_jp_a3_figures/figsrc/2026-08-13_referee_response/r2_results.json` | plotted |
| | | | `2026-10-08_r2c_matched_yield/r2c_results.json` | plotted |
| Fig. S3 | `figS2_pedestal_residuals` | `f08_pedestal.py` (second output of the Fig. 4 run) | as Fig. 4 | plotted |
| Fig. S4 | `fig5_sensitivity` | `f07_oxygen.py` | `2026-10-08_oxygen_v8/oxygen_boot_v8.json` | plotted |
| Fig. S5 | `fig4_brain_tumor` | `f06_brain.py` | `2026-08-16_fig4_rebuild/tau_voxel_map_8mm.json` | plotted |
| | | | `2026-10-19_review16/brain_corrected_rois.json` | plotted (block `corrected_regions` only) |
| | | | `2026-08-16_fig4_rebuild/recon_4mm/activity_x.npy` | plotted |
| | | | `2026-08-16_fig4_rebuild/recon_4mm/support_mask.npy` | plotted |
| | | | `2026-08-16_fig4_rebuild/tau_map_smoothed.npy` | plotted |
| | | | `2026-08-16_fig4_rebuild/mask_valid.npy` | plotted |
| Fig. S6 | `fig3_vipman_axial` | `f05_vipman.py` | `2026-09-30_feedback3/figdata/vipman_tau8_cache.npz` | plotted |
| | | | `2026-05-05_2122_phase3q_production/recon_phase3q/support_mask.npy` | plotted |
| | | | `2026-10-20_review19/vipman_slice/vipman_O4_slice.npz` | plotted (panel a) |
| | | | `2026-09-30_referee2/vipman_decomp.json` | plotted |
| | | | `2026-09-30_referee2/vipman_oracle_mix.json` | plotted |
| Fig. S7 | `figS3_eq2_refutation` | `fS1_refutation.py` | `2026-08-17_derivation_test/derivation_test.json` | plotted |
| Fig. S8 | `figS4_inflight` | `fS3_inflight.py` | `2026-09-30_feedback3/inflight/physics_level.json` | plotted |
| | | | `2026-09-30_feedback3/inflight/compare_inflight.json` | plotted |
| | | | `2026-09-30_feedback3/inflight/physics_level_delay_hists.npz` | plotted |
| | | | `2026-09-30_feedback3/inflight/physics_level_energy_hists.npz` | plotted |

## Provenance of `source_data/`

The folder keeps the dated result-folder names of the project results tree (`reconstruction/results/` of the simulation project), so every file has the same relative path as its original. 34 of the 36 data files are byte copies of those originals (sha256 checked against the inventory of the paper build). The other two:

- `2026-10-19_review16/brain_corrected_rois.json` is the results-tree version (2,971 bytes). The copy used for the paper build held only the block `corrected_regions` without the fields `s` and `ratio_emp_over_ana`. The values `f06_brain.py` reads are identical in both (equal canonical-JSON hash of `corrected_regions` with those two fields removed).
- `2026-10-20_review19/vipman_slice/vipman_O4_slice.npz` (36,201 bytes, written with `numpy.savez_compressed`) was made for this release by `make_vipman_slice.py` in the same folder from `2026-10-19_review16/vipman_recon_O4/activity_x.npy` (sha256 `b6c9de5b...`, not shipped, 784 kB) and `support_mask.npy` of the same reconstruction. It holds `act_z` (float32, 140 x 70, the mean of z slices 10 and 11, unmasked) and `sup2` (bool, 140 x 70 x 2, the support mask in those slices, 14,094 voxels). An earlier draft of the paper drew Fig. S6(a) from an 8-bit quantized transfer of the same slice (`2026-10-19_review16/vipman_O4_slice.npz`, not shipped): equal support mask, and inside the support every value equals the float32 value rounded to a step of 7.24 a.u. and capped at 1846.2 a.u. (6 voxels capped, exact maximum 1937.8 a.u.). Outside the support the 8-bit file held zeros, the float32 file the unmasked MLEM values. `f05_vipman.py` draws only the support.

`2026-08-16_fig4_rebuild/tau_voxel_map_8mm.json` contains a provenance string with an absolute path on the project workstation (`/home/dhritman/...`), the location of the triples file the map was computed from. That file is not shipped.

## Environment

The figure PDFs in the submitted paper were produced with Python 3.11.15, numpy 2.4.4, scipy 1.17.1 and matplotlib 3.10.9 from pip wheels (STIXGeneral fonts bundled with matplotlib, TrueType embedding). Fig. 1 was compiled with pdfTeX 1.40.25 (TeX Live 2023, pgf 3.1.10) on a system without `stix.sty`, so `fig1_methods.tex` fell back to `mathptmx` (Nimbus Roman).

Text placement in matplotlib depends on the FreeType library matplotlib is built against. The pip wheels of matplotlib bundle FreeType 2.6.1, as used for the paper. Linux distribution packages link the system FreeType instead (for example 2.14.1 on Arch Linux), which places text labels up to a few tenths of a point differently and, through `tight_layout`, can shift axes slightly. Data, marks and every printed number are unaffected. For renders identical to the paper, install matplotlib from pip as in `requirements.txt`.

`fig1_methods.tex` uses `stix.sty` when it is installed and `mathptmx` otherwise. With `stix.sty` installed, Fig. 1 is set in STIX instead of Nimbus Roman (same content). `make_figures.sh` prints a note when this happens.

## Reproduction test (2026-10-20)

`make_figures.sh` was run on a second machine (Arch Linux workstation) from this folder alone, with no `PLI_*` variable set, and every output was compared with the figure PDFs included in the built paper. Both sides were rendered at 150 dpi with pdftoppm (poppler 26.01 and 24.02, which gave identical renders for identical PDFs). Two Python environments were used on the same machine:

- run B: Python 3.14.2, numpy 2.3.5, scipy 1.17.0, matplotlib 3.10.8 pip wheel (FreeType 2.6.1)
- run A: Python 3.14.2, numpy 2.3.5, scipy 1.17.0, matplotlib 3.10.8 Arch package (FreeType 2.14.1)

TeX: pdfTeX 1.40.29 (TeX Live 2026) with `stix.sty` installed in the user tree.

| paper | stem | run B (pip wheel) | run A (distribution package) |
|---|---|---|---|
| Fig. 1 | `fig1_methods` | differs: STIX fonts (stix.sty installed). Compiled with stix.sty hidden: same fonts, same text, all 181 word boxes within 0.005 pt of the paper, render not pixel-identical (pdfTeX 1.40.29 writes glyph positions as TJ kerns, 1.40.25 as Td moves) | as run B |
| Fig. 2 | `fig2_closure_benchmark` | pixel-identical | text placement, 2,461 of 419,078 px |
| Fig. 3 | `fig6_clinical_mechanism_v2` | pixel-identical | text placement, 4,497 of 388,068 px |
| Fig. 4 | `fig9_pedestal_derivation` | pixel-identical | text placement, 6,501 of 660,070 px |
| Fig. 5 | `fig11_realistic_forest` | pixel-identical | text placement, 4,596 of 444,772 px |
| Fig. 6 | `fig10_rate_law` | pixel-identical | text placement, 11,481 of 649,438 px |
| Fig. S1 | `figS_benchspec` | pixel-identical | text placement and legend, 5,963 of 323,050 px |
| Fig. S2 | `fig_r2_transport_novelty` | pixel-identical | text placement, 3,472 of 345,540 px |
| Fig. S3 | `figS2_pedestal_residuals` | pixel-identical | text placement and axes position (tight_layout), 23,165 of 430,300 px |
| Fig. S4 | `fig5_sensitivity` | pixel-identical | text placement, 4,373 of 334,908 px |
| Fig. S5 | `fig4_brain_tumor` | pixel-identical (see note) | text placement, 1,256 of 334,908 px |
| Fig. S6 | `fig3_vipman_axial` | panel (a) differs from the earlier 8-bit drawing by design, panels (b) to (d) pixel-identical (see note). The submitted figure is drawn from the shipped float32 slice | panel (a) as in run B, plus text placement, 2,127 of 502,362 px against run B |
| Fig. S7 | `figS3_eq2_refutation` | pixel-identical | text placement, 2,117 of 353,600 px |
| Fig. S8 | `figS4_inflight` | pixel-identical | text placement, 7,354 of 691,080 px |

Pixel-identical means an identical 150-dpi render. In run B the PDFs of Figs. 2 to 6, S1 to S4, S7 and S8 are also byte-identical to the paper's apart from the metadata fields CreationDate, Creator and Producer. In both runs `fig2_closure_benchmark.json` is byte-identical to the paper build, and the printed output of every generator is identical between runs A and B and the paper environment (also with the unedited scripts).

Notes:

- Fig. S5, run B: the PDF differs from the paper's only in the coordinates of the colour-bar ticks of panel (a), in the fifth decimal of a point (the colour-scale top is the 99th percentile of the Gaussian-smoothed activity, which differs in its last digits between scipy 1.17.0 and 1.17.1). The render is identical.
- Fig. S6: panel (a) now uses the float32 slice described above instead of the 8-bit transfer. Against the paper, 14,169 of the 46,810 pixels of the panel (a) image change, each by at most 3 of 255 in any colour channel, and 6 pixels change at one tick mark of the panel (a) colour bar, because the colour-scale top (99th percentile inside the support) is 1711.44 a.u. instead of 1708.63 a.u. (+0.16%). No other pixel of the figure changes. The run B PDF is byte-identical (apart from metadata) to the one made with the same input in the paper environment.
- Run A against run B: every difference in the PDF drawing streams is in text-positioning operators, except in Fig. S3 (axes geometry follows the text extents through `tight_layout`) and Fig. S1 (legend handle positions follow the label extents).

## Further notes

- `fS3_inflight.py` reads `compare_inflight.json` from the path in `PLI_INFLIGHT_JSON` if that variable is set (a development override). `make_figures.sh` unsets it.
- `f05_vipman.py` reads the 8 mm lifetime map from `2026-09-30_feedback3/figdata/vipman_tau8_cache.npz`. Only if that file is missing does it try to recompute the map from `2026-08-13_overnight/O4_rerun/vipman/triples_detector_phase3q.h5` (3.2 GB, not shipped, needs h5py).
- `f09_mechanism.py` prints comparisons with the discovery stream from two inputs that do not reach the figure (`crt_sweep_v3_all.csv`, `2026-09-30_referee2/acceptance_tm.json`). `f11_realistic.py` uses `fig11_component_shift.json` only in an assertion that guards consistency with the published component fits.
- Figure scripts of earlier manuscript versions and figures not in the paper are not included.
