# Changes to the figure scripts for code release v1.1.0

Notation: `<PROJECT>` is the project tree on the analysis workstation, `<DATA>` its data disk, `<HOME>` the home folder. The shipped files use the environment variables instead.

The scripts in this folder were taken from the figure working copy that produced the figure PDFs of the submitted paper (container folder `figdev16/`, current versions without `.v2xx` suffix, and `figdev16/tikz/fig1_methods.tex`). This file lists every change made for the release. No change alters a plotted value, with the one deliberate exception of the Fig. S6(a) input (item 3), which replaces an 8-bit quantized slice by the float32 slice it was quantized from.

## Files

| file | sha256 of the source version | sha256 in this release | changed |
|---|---|---|---|
| `common.py` | `1e465d0684aa475c2a289a8e79c5d541f358fc9a20c524ed2b1d4909f18c7608` | `6e62d0195610ef25b4ad147de82ed8770deb1c4c177a02f46c2f6340f7b3281a` | yes (1, 6) |
| `ejmp_style.py` | `29f832d52b3fdd7d786ca0d095b0654f885f3918cacc90de4538412889e703c4` | `73ffc6589b12764fa6cbe4f1c3bf44cc3cea4e0bbb21d740c70319490506f728` | moved, yes (1, 6) |
| `f02b_closure_benchmark.py` | `0ae67a7fd6aa24f79e580c30f6205ef3895aeca1e5218eaf6c64380926568526` | `b714d774f434d96cdbfb945b436ef6505fd682e3ffcfc760cf7c4390519eeb0a` | yes (6) |
| `f04_transport.py` | `873d018c285949ce1ae5380ae5d511205c8e07caf0cbad83ca9bb6a7b1a5c955` | `80e3d4f5e1097454bc98668a07185732d02552cb60a87bd3e909ecf332efd70d` | yes (2, 6) |
| `f05_vipman.py` | `5db0dd5c639ce706d0aece0ed9b07ab1274670a634eb96d2c449785a6c4cb1fa` | `a59446f418839d65df210a060dc40b81f6158a503d60bfed0dc2333e0627feb9` | yes (3, 4, 6) |
| `f06_brain.py` | `48f4414aef45b2a2968c15d29047a57785a0ffe793b586faf306a81ea07fb52f` | `90a7f7c624df04a38743d66842c2e1d549502e3f53ba471c580c3ad79a24f9bf` | yes (4, 6) |
| `f07_oxygen.py` | `7560b36858c60ef5494967f317b6aef718e4dcf32eaf5f3c00ff6bb300cfa3ee` | `fac7c7530a1e52728773333f8b2bf9b013dd1a59589d2a61e52b82d89660f2e0` | yes (6) |
| `f08_pedestal.py` | `5f04948eb9286e4541ad9b05acdd973a65c98edc407c8843b48a98ea6153762f` | `2c5754e71e66ffbde260c99230ffd6f6d5ab91391e8821e036751a0092dbb377` | yes (5, 6, 7) |
| `f09_mechanism.py` | `dbdeebfb0e51e01fcf72be129315c741eb70a4cb341bcd39f51e18241222dffb` | `368e7e83001c2cc3a7667fb0f4c97d14184f79b120a82dc285d9ee6cc99c058d` | yes (6, 7) |
| `f10_rate.py` | `b415031686d37e20a3e7f6513f43cc21f8fb36924d5d9ad14abdd971ef3770bf` | `4f77934b18046257d2e9738f76d1f5761ef6f0d6228941b48861f9887a47b3df` | yes (4, 6, 7) |
| `f11_realistic.py` | `06f6831e3da2480295582d78b3d03bf8e069428621a459d1c80e1fa6f86d69ec` | `f15a7ec458c66cf6ba208d4eeac76dbd12f0eea516148289eec540a4a04c2824` | yes (2, 6, 7) |
| `fS1_refutation.py` | `dd0ef0ab536bf1c6d8c54e03b01303cf6510c6944820ef8f9dae4710204f731d` | `02b5e487a025dac1028538df04e1f2fa3200eee8f198f5ea429a1d6b869367cc` | yes (6, 7) |
| `fS3_inflight.py` | `52a5e9bb15fa43d04826600802d31bd48afc8f26ccd4014f580dd50e9b1cf429` | `a11c4263dfc061355a7a4f460546fd574664f68062021b1435482d18370d9fbc` | yes (6, 7) |
| `fS_benchspec.py` | `83d8861547d17723e9f04228f4a2a003384484ceee178f701e40c05a476f7780` | `eb73b38a62eff332ae0f14eb16f38ec47b523bb143ab88d11ea71efce35e5b8e` | yes (6) |
| `fig1_methods.tex` | `b1d2a77cc62e531b28621079cf89f658396d6d5985b87ff79e221fbbf2767615` | same | no |

`ejmp_style.py` is copied, with the punctuation of its comments changed (6), from `reconstruction/results/2026-08-11_ejmp_figure_regen/ejmp_style.py`. `make_figures.sh`, `requirements.txt`, `README.md`, `SHA256SUMS` and this file are new.

## Changes

1. `common.py`, environment
   - `PLI_RES` defaults to `source_data/` next to `common.py` (was the absolute workstation path `<PROJECT>/reconstruction/results`).
   - `PLI_FIGOUT` defaults to `out/` next to `common.py` (was `<PLI_RES>/2026-09-30_referee2/figs/out`, which wrote into the input tree).
   - `ejmp_style.py` is imported from the folder of `common.py`. The `sys.path.insert(0, RES / '2026-08-11_ejmp_figure_regen')` line and the then unused `import sys` are removed.
   - New module-level name `HERE` (the folder of `common.py`). No script defines or uses a name `HERE`.
2. Optional inputs made required
   - `f04_transport.py`: `2026-10-08_r2c_matched_yield/r2c_results.json` was read only `if R2C.exists()`. Without it the figure silently lost the matched-yield points and its y-range changed (upper limit 58 instead of 112). The block now runs unconditionally (de-indented, otherwise identical) and the y-limit is the constant 112.
   - `f11_realistic.py`: `2026-09-30_feedback3/c0_sim1/c0_compare.json` was read only if present, with a fallback (and a printed warning) to the superseded frozen zero-delay stream. The block now runs unconditionally (de-indented, otherwise identical) and the fallback branch is removed.
3. `f05_vipman.py`, Fig. S7 panel (a) input (Fig. S6 in v1.1.0)
   - The slice is read from `2026-10-20_review19/vipman_slice/vipman_O4_slice.npz` (float32, written by `make_vipman_slice.py` in the same folder from `2026-10-19_review16/vipman_recon_O4/activity_x.npy`, z slices 10 and 11 averaged) instead of `2026-10-19_review16/vipman_O4_slice.npz`, which had been decoded from an 8-bit quantized transfer of the same slice (step 7.24 a.u., values above 1846.2 a.u. clipped). The support mask `sup2` is identical in both files. The comment describing the input is updated. Effect on the figure: see README.md, section "Reproduction test".
4. Dead inputs removed (each was loaded and never reached the figure)
   - `f05_vipman.py`: `2026-05-05_2122_phase3q_production/recon_phase3q/activity_x.npy` was loaded into `act` and used only for an `act_z` that was overwritten three lines later. Both statements are removed. `support_mask.npy` of the same folder is still read (assert and outlines).
   - `f06_brain.py`: `2026-08-13_overnight/W3_detection_sigma/w3_activity_rois.json` was loaded into `W3`, never used. The load is replaced by a comment.
   - `f10_rate.py`: `2026-09-30_referee2/ladder_joint/ladder_boot.json` was loaded into `B` and unpacked into `lo, hi`, never used. Both statements are replaced by a comment.
5. `f08_pedestal.py`: `np.load(..., allow_pickle=True)` is now `allow_pickle=False`. Every array in `hists_c0_0ps.npz` is numeric or a fixed-width unicode string (`<U27`), so nothing in the file needs unpickling.

6. Punctuation in comments and docstrings (`common.py`, `ejmp_style.py`, `f02b_closure_benchmark.py`, `f04_transport.py`, `f05_vipman.py`, `f06_brain.py`, `f07_oxygen.py`, `f08_pedestal.py`, `f09_mechanism.py`, `f10_rate.py`, `f11_realistic.py`, `fS1_refutation.py`, `fS3_inflight.py`, `fS_benchspec.py`, `make_figures.sh`, `source_data/2026-10-20_review19/vipman_slice/make_vipman_slice.py`): em dashes and semicolons were replaced by colons, commas, parentheses or periods. No code line changed.

7. Supplement numbering (v1.1.1, 2026-10-20)
   - Comments and docstrings that name supplementary notes, tables and figures follow the numbering of the final supplement (Notes S1 to S8, Tables S1 to S9, Figures S1 to S8): `f08_pedestal.py`, `f09_mechanism.py`, `f10_rate.py`, `f11_realistic.py`, `fS1_refutation.py`, `fS3_inflight.py`, the comments of `make_figures.sh` and the docstring of `source_data/2026-10-20_review19/vipman_slice/make_vipman_slice.py`. No code line changed in these files except the one below.
   - `fS1_refutation.py`: the legend label of the eight-tube points reads "earlier-yield stream" instead of "discovery stream" (two string literals). This is the only change to a figure.

Not changed, documented in README.md instead: the optional override `PLI_INFLIGHT_JSON` in `fS3_inflight.py` (unset by `make_figures.sh`), the `h5py` fallback in `f05_vipman.py` that is used only if the shipped 8 mm cache is absent, the print-only inputs of `f09_mechanism.py` and the assert-only input of `f11_realistic.py`, and the docstrings (some still name the former inputs).

## Unified diff (source version against release version)

The diff shows the release text of v1.1.0 (items 1 to 5). Item 7 (v1.1.1) later changed some of the comment lines shown, as listed under item 7.

```diff
--- container figdev16/common.py
+++ release figures/common.py
@@ -17,15 +17,16 @@
 No interpretive sentences are drawn into a figure.
 """
 import os
-import sys
 from pathlib import Path
-RES = Path(os.environ.get('PLI_RES', '<PROJECT>/reconstruction/results'))
+# Code release v1.1.0: inputs default to the shipped source_data/ next to this file (dated result-folder
+# names kept), outputs to out/. Both can be redirected with PLI_RES and PLI_FIGOUT.
+HERE = Path(__file__).resolve().parent
+RES = Path(os.environ.get('PLI_RES', str(HERE / 'source_data')))
 FS = RES / '2026-09-24_jp_a3_figures/figsrc'
 MF = FS / '2026-05-11_paper_PMB/manuscript_ejmp/figures'
 R2 = RES / '2026-09-30_referee2'
 F3 = RES / '2026-09-30_feedback3'
-sys.path.insert(0, str(RES / '2026-08-11_ejmp_figure_regen'))
-import ejmp_style as st  # noqa: E402
+import ejmp_style as st  # noqa: E402  (shipped next to this file, was results/2026-08-11_ejmp_figure_regen/ejmp_style.py)
 st.apply()
 import matplotlib  # noqa: E402
 import matplotlib.pyplot as plt  # noqa: E402
@@ -41,7 +42,7 @@
 plt.rcParams.update({'font.size': 9.0, 'axes.labelsize': 9.0, 'axes.titlesize': 9.0, 'legend.fontsize': 8.0,
                      'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'savefig.pad_inches': 0.0,
                      'savefig.bbox': None})
-OUT = Path(os.environ.get('PLI_FIGOUT', str(R2 / 'figs' / 'out')))
+OUT = Path(os.environ.get('PLI_FIGOUT', str(HERE / 'out')))
 OUT.mkdir(parents=True, exist_ok=True)
 
 NAME = {'T1_water_shibuya2020_n2': 'water (N$_2$)', 'T2_adipose_avachat2024': 'adipose',
--- container figdev16/f04_transport.py
+++ release figures/f04_transport.py
@@ -88,19 +88,20 @@
 labs = ['measured, material pair', 'first-order consistency check\nfrom the same-run shares and yields', 'null control (equal lifetimes)']
 # Review round 8: registered matched-yield arms (equal o-Ps yields everywhere, three-photon channel off), paired event
 # for event, 2026-10-08_r2c_matched_yield/r2c_results.json
+# Code release v1.1.0: this input is required (it was optional, and without it the figure silently lost the
+# matched-yield points and changed its y-range).
 R2C = RES / '2026-10-08_r2c_matched_yield/r2c_results.json'
-if R2C.exists():
-    C = json.load(open(R2C))
-    for tag, tc, ts in PAIRS:
-        r = C[tag]
-        ax.errorbar([tc - ts + 0.12], [r['diff_ps']], yerr=[r['se_paired_ps']], fmt='s', ms=4.6, mfc='white', mec=THIS, color=THIS,
-                    mew=1.0, elinewidth=1.0, capsize=2, zorder=6)
-        ax.plot([tc - ts + 0.15, tc - ts + 0.19], [r['first_order_registered_ps']] * 2, '-', lw=1.3, color=REF, zorder=7)
-        print(tag, 'matched-yield paired', round(r['diff_ps'], 2), '+/-', round(r['se_paired_ps'], 2))
-    h.insert(1, Line2D([], [], marker='s', ms=4.6, mfc='white', mec=THIS, ls='none'))
-    labs.insert(1, 'equal o-Ps yields, paired arms\n(dash: registered first order)')
+C = json.load(open(R2C))
+for tag, tc, ts in PAIRS:
+    r = C[tag]
+    ax.errorbar([tc - ts + 0.12], [r['diff_ps']], yerr=[r['se_paired_ps']], fmt='s', ms=4.6, mfc='white', mec=THIS, color=THIS,
+                mew=1.0, elinewidth=1.0, capsize=2, zorder=6)
+    ax.plot([tc - ts + 0.15, tc - ts + 0.19], [r['first_order_registered_ps']] * 2, '-', lw=1.3, color=REF, zorder=7)
+    print(tag, 'matched-yield paired', round(r['diff_ps'], 2), '+/-', round(r['se_paired_ps'], 2))
+h.insert(1, Line2D([], [], marker='s', ms=4.6, mfc='white', mec=THIS, ls='none'))
+labs.insert(1, 'equal o-Ps yields, paired arms\n(dash: registered first order)')
 ax.legend(h, labs, frameon=False, fontsize=8.0, loc='upper right', borderaxespad=0.1)
-ax.set_xlim(-0.1, 1.2); ax.set_ylim(-110, 112 if R2C.exists() else 58); ax.set_yticks([-100, -80, -60, -40, -20, 0, 20])
+ax.set_xlim(-0.1, 1.2); ax.set_ylim(-110, 112); ax.set_yticks([-100, -80, -60, -40, -20, 0, 20])
 ax.set_xlabel(r'lifetime step $\tau_{\rm core}-\tau_{\rm shell}$ (ns)')
 ax.set_ylabel(r'$\hat\tau_{\rm stop}-\hat\tau_{\rm source}$ (ps)')  # review round 13: stopping-point terminology
 st.panel_letter(ax, 'b')
--- container figdev16/f05_vipman.py
+++ release figures/f05_vipman.py
@@ -38,14 +38,17 @@
 valid = N >= 30
 tmap = np.where(valid, s / np.maximum(N, 1), np.nan)
 smap = np.where(valid, tmap / np.sqrt(np.maximum(N, 1)), np.nan)
-act = np.load(RECON / 'activity_x.npy'); sup4 = np.load(RECON / 'support_mask.npy')
+sup4 = np.load(RECON / 'support_mask.npy')
 Z8, Z4 = 5, [10, 11]
 EXT = (-280, 280, -140, 140)
-act_z = np.where(sup4[:, :, Z4].any(axis=2), act[:, :, Z4].mean(axis=2), np.nan).T
 # Review round 17 (second report, Sec. 8): panel (a) shows the MLEM reconstruction of the corrected-yield acquisition
 # (2026-10-19_review16/vipman_recon_O4, simple_recon single, 4 mm, 10 iterations, the settings of recon_phase3q),
-# central two 4 mm slices, transferred as an 8-bit slice (vipman_O4_slice.npz); the support mask is identical.
-_V = np.load(RES / '2026-10-19_review16/vipman_O4_slice.npz')
+# central two 4 mm slices. The support mask is identical.
+# Code release v1.1.0: the slice is read in float32 from 2026-10-20_review19/vipman_slice/vipman_O4_slice.npz
+# (written by make_vipman_slice.py in that folder from vipman_recon_O4/activity_x.npy). The figure in the submitted
+# paper was drawn from an 8-bit quantized transfer of the same slice (2026-10-19_review16/vipman_O4_slice.npz).
+# The recon_phase3q activity, previously loaded here and then overwritten by this slice, is no longer read.
+_V = np.load(RES / '2026-10-20_review19/vipman_slice/vipman_O4_slice.npz')
 assert (_V['sup2'] == sup4[:, :, Z4]).all()
 act_z = np.where(_V['sup2'].any(axis=2), _V['act_z'], np.nan).T
 sup_z = sup4[:, :, Z4].any(axis=2).astype(float).T
--- container figdev16/f06_brain.py
+++ release figures/f06_brain.py
@@ -12,7 +12,7 @@
 from matplotlib.lines import Line2D
 
 SESS = RES / '2026-08-16_fig4_rebuild'
-W3 = json.load(open(RES / '2026-08-13_overnight/W3_detection_sigma/w3_activity_rois.json'))
+# Code release v1.1.0: 2026-08-13_overnight/W3_detection_sigma/w3_activity_rois.json was loaded here but never used. Load removed.
 meta = json.load(open(SESS / 'tau_voxel_map_8mm.json')); TAU_BG = meta['tau_bg_global']
 TRUTH = [dict(r=4, cx=-66), dict(r=8, cx=-49), dict(r=12, cx=-24), dict(r=16, cx=9), dict(r=20, cx=50)]
 # Review round 17: regions redefined on the corrected-yield reconstruction (2026-10-19_review16/brain_corrected_rois.py,
--- container figdev16/f08_pedestal.py
+++ release figures/f08_pedestal.py
@@ -18,7 +18,7 @@
 # Review round 13 (2026-10-14): eight-tube data from the corrected-yield arm C0 (hists_c0_0ps.npz, same dump script
 # as hists_sim1.npz). Panel (c) and Fig. S2 use the twenty-region signed regression with the C0 tubes in place of the
 # discovery tubes (2026-10-14_known_randoms/signed_eq3_c0.json). The other twelve regions are unchanged.
-Z = np.load(F3 / 'c0_sim1/hists_c0_0ps.npz', allow_pickle=True)
+Z = np.load(F3 / 'c0_sim1/hists_c0_0ps.npz', allow_pickle=False)  # v1.1.0: was allow_pickle=True. Every array is numeric or <U27
 names = [str(s) for s in Z['voi_names']]
 tubes = [n for n in names if n.startswith('T')]
 b = Z['tau_bins']; c = 0.5 * (b[:-1] + b[1:])
--- container figdev16/f10_rate.py
+++ release figures/f10_rate.py
@@ -20,10 +20,10 @@
 
 d = json.load(open(MF / 'fig10_rate_law.json'))
 L = d['ladder']; rows = L['rows']
-B = json.load(open(R2 / 'ladder_joint/ladder_boot.json'))
+# Code release v1.1.0: 2026-09-30_referee2/ladder_joint/ladder_boot.json was loaded here into (lo, hi) and never used,
+# load removed.
 mu = np.array([r['mu'] for r in rows]); did = np.array([r['did'] for r in rows])
 C0, N0 = np.exp(L['ln_intercept']), L['n_exponent']
-lo, hi = B['boot_ci95']
 FH = 124
 fig = mm_fig(150, FH)
 ax = ax_mm(fig, 14, 70, 60, 47)
--- container figdev16/f11_realistic.py
+++ release figures/f11_realistic.py
@@ -29,20 +29,18 @@
 # Feedback round 3: the frozen zero-delay eight-tube stream predates the I3-convention correction
 # (o-Ps yield x3/4 in every tube), so the eight-tube zero-delay arm is taken from C0, the corrected
 # rerun identical to C2 except for the direct delay (2026-09-30_feedback3/c0_sim1/c0_compare.json).
+# Code release v1.1.0: this input is required (it was optional, with a fallback to the superseded frozen stream).
 C0F = F3 / 'c0_sim1/c0_compare.json'
-C0 = json.load(open(C0F)) if C0F.exists() else None
-if C0 is not None:
-    C0R = {r['region']: r for r in C0['rows']}
-    for r in R:
-        if r['phantom'] == 'Sim-1':
-            c = C0R[r['region']]
-            assert abs(c['c2']['did'] - r['did388']) < 1e-6, r['region']   # same 388 ps arm
-            r['did0'], r['sig0'] = c['c0']['did'], c['c0']['sd_did']
-            r['retention'], r['retention_err'] = c['retention_vs_c0'], c['retention_vs_c0_se']
-            r['ret_june'], r['ret_june_err'] = c['retention_vs_june0'], c['retention_vs_june0_se']
-    print('eight-tube zero-delay arm from C0; pooled retention', C0['summary']['retention_ivw_c0'])
-else:
-    print('WARNING: c0_compare.json absent; eight-tube zero-delay arm is the frozen stream')
+C0 = json.load(open(C0F))
+C0R = {r['region']: r for r in C0['rows']}
+for r in R:
+    if r['phantom'] == 'Sim-1':
+        c = C0R[r['region']]
+        assert abs(c['c2']['did'] - r['did388']) < 1e-6, r['region']   # same 388 ps arm
+        r['did0'], r['sig0'] = c['c0']['did'], c['c0']['sd_did']
+        r['retention'], r['retention_err'] = c['retention_vs_c0'], c['retention_vs_c0_se']
+        r['ret_june'], r['ret_june_err'] = c['retention_vs_june0'], c['retention_vs_june0_se']
+print('eight-tube zero-delay arm from C0; pooled retention', C0['summary']['retention_ivw_c0'])
 
 # pooled numbers quoted in the text (asserted in the PDF): all three estimators refitted on the
 # same 400 paired shard resamples (comp_joint_ratio.json). The frozen analysis gave 9.31 / 1.95 / 4.77.
```
