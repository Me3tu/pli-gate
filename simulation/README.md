# Simulation: macros, runners and run scripts

This folder holds the GATE macros, the parallel runner and the run scripts
behind every simulation in the paper. Each file is a copy of the production
file with workstation paths replaced by variables. `EDITS.md` lists every
change, file by file, and `RUNS.json` is the run table below in
machine-readable form.

The scanner geometry and two further inputs are not part of the release (see
"Inputs you must supply"). Nothing here runs without them, apart from the
annihilation-vertex unit test.

## Layout

| folder | contents |
|---|---|
| `runner/` | `run_gate_parallel.sh` (splits a macro into K seeded shards, runs them under GNU parallel, merges with `hadd`), `run_gate_parallel_bin.sh` (the same runner with a binary override, used by the in-flight arm), `run_bounded.sh` (per-shard memory cap through a systemd user scope), `env_example.sh` (puts Geant4, ROOT and Gate on the path), `paths.sh` (sets and checks the variables below) |
| `geometry/` | project material files for GATE (`materials_pli.db`, `materials_pli_phase2b.db`) and the phantom geometries that the run macros include |
| `sources/` | source sub-macros (124I tubes, uniform 44Sc, brain-tumour inserts) |
| `macros/` | one run macro per phantom, plus `make_core_shell.py`, which writes the two core-shell macros |
| `runs/` | the run script of each production run |
| `lu176/` | analytic 176Lu background generator and the singles injector used by the VIP-Man run |

The project material databases for the positronium physics are in
`../materials/` (see `../materials/PROVENANCE.md`).

## Inputs you must supply

These inputs are derived from J. Peter's Musiré simulation package or from
the VIP-Man phantom, and their licence for redistribution is not cleared.
They are available from the authors on request (dhritman@dhar.dev). Every
macro refers to them only through the aliases below, so a run that lacks one
stops at once with `<NAME> is not defined as a shell variable`.

| alias / variable | what it must point to | needed by |
|---|---|---|
| `PLI_SCANNER_MACRO` | the Siemens Biograph Vision Quadra geometry macro (`quadra_detector.mac`): a `cylindricalPET` system of 38 rsectors x 4 blocks x 20 x 80 LSO crystals of 3.2 x 3.2 x 20 mm, whose crystal-level sensitive volume is named `LSO` (the digitizer commands in the run macros address `/gate/digitizerMgr/LSO/...`). It uses the materials `Air`, `PTFE` and `LSO` of `geometry/materials_pli.db` | every macro in `macros/` except `unit_vertex_*` |
| `PLI_GATE_MATERIALS_DB` | the GATE material database `gate-materials.db` | `nema_iq`, `vipman_slab`, `stock_vs_custom_*`, `unit_vertex_stock` |
| `PLI_VIPMAN_DIR` | a directory holding `vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd` with its `.raw`, `vipman-atlas-gate-materials.dat` and `vipman-atlas-gate-activities-f18-fdg.dat` (and `vipman-atlas-lungs-sim-mumap-harmonized.mhd` for the reconstruction step of `runs/vipman_batches.py`) | `vipman_slab`, `stock_vs_custom_*` |

`geometry/materials_pli.db` and `geometry/materials_pli_phase2b.db` are the
project's own files and are shipped. Their element table and the `Vacuum`,
`Water` and `LSO` entries have the same values as GATE's own
`GateMaterials.db` (LGPL-3.0). Nothing in them is specific to Musiré.

## Variables

| variable | meaning | default |
|---|---|---|
| `PLI_SIM` | this folder | set by `runner/paths.sh` |
| `PLI_PKG` | the release root | set by `runner/paths.sh` |
| `PLI_OUT` | output root. Data go to `$PLI_OUT/<data folder>` and logs to `$PLI_OUT/results/<dated folder>`, the folder names of the production runs | required |
| `PLI_ANALYSIS` | the analysis part of the release, where the run scripts find the singles converter, the concatenator and the triple extractors (`results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py`, `results/2026-08-16_T6_nema/concat_singles_u16.py`, `reconstruction/*.py`) | `$PLI_PKG/analysis` |
| `PLI_ENV_SCRIPT` | script that puts Geant4, ROOT and Gate on the path | `runner/env_example.sh` |
| `PLI_STACK`, `GEANT4_INSTALL`, `ROOT_INSTALL`, `GATE_INSTALL` | install locations read by `env_example.sh` | `$HOME/pli-stack/install/{geant4,root,gate}` (the layout of `docs/BOOTSTRAP.md`) |
| `PLI_VENV` | Python virtual environment to activate before the conversion steps | none |
| `PLI_SCANNER_MACRO`, `PLI_GATE_MATERIALS_DB`, `PLI_VIPMAN_DIR` | the inputs above | required where listed |
| `PLI_MATERIALS_JSON` | positronium database, a runtime switch of the physics | the run scripts use `$PLI_PKG/materials/pli_materials.json` unless it is set |
| `GATE_INFLIGHT_BIN_DIR` | `bin/` of the in-flight variant build (`../integration/PATCH_NOTES.md`) | required by `runs/inflight_arm.sh` |
| `PLI_BOUNDED_LOGDIR` | scratch folder of `run_bounded.sh` | `${TMPDIR:-/tmp}/pli-bounded` |
| `PLI_PYTHON`, `PLI_EXTRA_PYTHONPATH`, `PLI_SENSITIVITY_NPY`, `PARALLELPROJ_C_LIB`, `PARALLELPROJ_CUDA_LIB` | Python interpreter, extra module path, sensitivity image and parallelproj libraries of the VIP-Man reconstruction step | the running interpreter, none, `$PLI_ANALYSIS/results/sensitivity_4mm_fullfidelity_norm.npy` (not shipped: the scanner sensitivity image of the VIP-Man reconstruction, computed from the scanner geometry, available from the authors), `/usr/local/lib/libparallelproj_{c,cuda}.so` |

Inside the macros the paths are GATE aliases (`{PLI_SIM}`,
`{PLI_SCANNER_MACRO}`, `{PLI_GATE_MATERIALS_DB}`, `{PLI_VIPMAN_DIR}`,
`{PLI_OUT}`), read from the environment by `/control/getEnv` lines at the top
of each macro. To run a macro by hand, export the variables it names and call
`Gate macros/<name>.mac`.

## Running

```bash
export PLI_OUT=/scratch/pli-out
export PLI_SCANNER_MACRO=/path/to/quadra_detector.mac     # from the authors
export PLI_ENV_SCRIPT=/path/to/your/env.sh                 # or set PLI_STACK
bash simulation/runs/c2_eight_tube_388ps.sh
```

Every run script sources `runner/paths.sh`, changes to this folder, sets the
runtime switches of its run and calls the runner once per batch. The runner
gives shard t of batch b the seed 1,000,000 + 1000 b + t (t = 1..K), patches
the output names and the per-shard time or primary count, and merges the
shards with `hadd`. It needs GNU `parallel`, ROOT's `hadd`, Python with
`uproot`, and a systemd user session for the memory caps of `run_bounded.sh`.
Thread counts (12, 16, 30) and caps (4 to 5 GB per shard) are those of the
16-core, 78 GB workstation the runs were made on. The Python drivers
(`closure_o4_driver.py`, `vipman_batches.py`, `seed_repeatability.py`) run
the same runner with the batch layout of the table.

Two waits that only coordinated runs on one machine are disabled
(`c0_eight_tube.sh` waited for the in-flight arm, `core_shell_matched_yield.sh`
for the oxygen run). `closure_o4_driver.py` is the version of 2026-08-16. Its
product registry refuses the `sim1_calibration` recipe because that recipe
keeps no singles, while the Sim-1 closure itself ran on 2026-08-14 with an
earlier version of the same driver. Run that phantom with
`singles=True` added to its recipe, or call the runner with the layout of row
4a.

## Run table

Binaries: e37be8b4 is the production binary of `../ENVIRONMENT_MANIFEST.json`.
3522d567 is the archived build of 2026-08-14 14:54, used from then until the
production build of 2026-08-19 20:19, and kept without a source snapshot.
cc552a68 is the in-flight variant (`../integration/variants/`). Runs made
before 2026-08-14 14:54 used builds that were not retained and cannot be
rebuilt bit for bit. "Logged" means the run log records the binary hash, "by
date" that it follows from the run date and the install time.

Materials databases: v8 and v7 are `../materials/pli_materials.json` and
`../materials/previous/pli_materials.json`. They differ at run time only in
the five `TumorTissue_pO2_*` entries of the oxygen phantom (at most 0.11 ps),
so the v7 runs below reproduce with v8. The database hashes in the table below
are those of the files the runs used (`../materials/PROVENANCE.md` gives them
beside the hashes of the shipped copies). The database 2a26b68e of the O4 reruns
and the earlier ones are not retained. GATE material file: runs before
2026-08-22 loaded `materials_pli.db` with the old XAD-4 densities (1.02 and
1.08 g/cm3), which matters only for the Mercolli tubes. The shipped file has
the corrected 0.67 g/cm3.

| # | phantom | paper role | binary | macro | run script | batches x shards | seeds | switches | materials database |
|---|---|---|---|---|---|---|---|---|---|
| 1 | eight tubes, 44Sc, 8 x 1.77 MBq | C0: primary eight-tube stream, corrected yield, direct annihilation instantaneous | e37be8b4 (logged) | `macros/eight_tube_sc44.mac` | runs/c0_eight_tube.sh | time-driven, batches 0-9 x 12 shards + batch 10 x 8 = 128 shards, 0.137933 s per shard (about 2.5e8 decays) | 1,000,001-1,009,012 and 1,010,001-1,010,008 | PLI_DIRECT_TAU_PS unset | ../materials/previous/pli_materials.json (version 7, f6d5a53c) (run 2026-09-30) |
| 2 | eight tubes, 44Sc | C2: direct annihilation at 388 ps | e37be8b4 | `macros/eight_tube_sc44.mac` | runs/c2_eight_tube_388ps.sh | as C0 | as C0 (paired) | PLI_DIRECT_TAU_PS=388 | ../materials/previous/pli_materials.json (version 7, f6d5a53c) |
| 3 | eight tubes, 44Sc | in-flight annihilation sensitivity arm | cc552a68 | `macros/eight_tube_sc44.mac` | runs/inflight_arm.sh (runner/run_gate_parallel_bin.sh) | as C2 | as C2 | PLI_DIRECT_TAU_PS=388, PLI_INFLIGHT_STOCK=1 | ../materials/previous/pli_materials.json (version 7, f6d5a53c) |
| 4a | 12-tube calibration closure, Sim-1 (8 x 44Sc tubes) | calibration closure | 3522d567 | `macros/closure_sim1_sc44.mac` | runs/closure_o4_driver.py (phantom sim1_calibration) | primaries mode, 10 batches x 16 shards x 1,562,500 = 2.5e8 primaries. The macro keeps setTimeStop 1.0 s, so each shard spreads its primaries evenly over 1 s | 1,000,000 + 1000 x batch + shard (shard 1..K), batches 0-9 | none | database 2a26b68e of the O4 reruns (not retained) |
| 4b | 12-tube calibration closure, Sim-2 (4 x 124I, Mercolli XAD-4 tubes) | calibration closure | 3522d567 | `macros/mercolli_xad4_i124.mac + geometry/phantom_mercolli_xad4.mac + sources/iodine124_tubes_poisson.mac` | runs/closure_o4_driver.py (phantom sim2_mercolli_te) | time mode, 4 batches x 16 shards, 2.520161 s per batch | 1,000,000 + 1000 x batch + shard (shard 1..K), batches 0-3 | none | database 2a26b68e of the O4 reruns (not retained), with the GATE material file from before the 2026-08-22 XAD-4 density correction (not retained) |
| 5 | Mercolli XAD-4 tubes, 124I | replication of the Mercolli benchmark at 388 ps | e37be8b4 (by date) | as 4b | runs/mercolli_replication.sh | 30 batches x 30 shards, 30 s per batch (900 s) | 1,000,000 + 1000 x batch + shard (shard 1..K), batches 0-29 | PLI_DIRECT_TAU_PS=388 | ../materials/previous/pli_materials.json (version 7, f6d5a53c), densities corrected. geometry/phantom_mercolli_xad4.mac was edited on 2026-08-30, after this run, and its header says only a density comment changed |
| 6 | Steinberger disks (quartz glass, polycarbonate) | reference-material comparison (Table 5) | June build (not retained) and 3522d567 | `macros/steinberger_disks.mac` | June launch script not shipped. O4 rerun: runs/closure_o4_driver.py (steinberger_disks) | O4: 12 batches x 16 shards x 5 s | 1,000,000 + 1000 x batch + shard (shard 1..K) | none | June run: database not retained. O4 rerun: database 2a26b68e of the O4 reruns (not retained). Which acquisition feeds Table 5 is set in the analysis part |
| 7 | water cylinder, 20 mm | rate, retention and pedestal checks | build of 2026-08-12 (not retained) and 3522d567 | `macros/water_cylinder_20mm.mac` | runs/closure_o4_driver.py (steinberger_water). The 2026-08-12 driver is not shipped | 8 batches x 16 shards, 2.93 s per batch (23.44 s), 7 MBq | 1,000,000 + 1000 x batch + shard (shard 1..K) | none | database 2a26b68e of the O4 reruns (not retained) |
| 8 | water cylinder, 20 cm x 30 cm | Table 5 | e37be8b4 (by date) | `macros/water_cylinder_20cm.mac` | runs/water_20cm.sh 200 201 202 203 204 205 206 207 | batches 200-207 x 16 shards, 2.93 s per batch | 1,200,001-1,207,016 | none | ../materials/pli_materials.json (version 8, 23f607ea) (water entries are the same in v7 and v8) |
| 9a | NEMA IEC image-quality phantom | NEMA IQ, direct annihilation instantaneous | build of 2026-08-11 (not retained) | macros/nema_iq.mac (needs PLI_GATE_MATERIALS_DB) | runs/nema_iq_0ps.sh (batches 11-15 were rerun with a script that is not shipped) | 16 batches x 16 shards, 2.2 s per batch | 1,000,000 + 1000 x batch + shard (shard 1..K), batches 0-15 | none | database of 2026-08-11 (not retained) |
| 9b | NEMA IEC image-quality phantom | NEMA IQ, direct annihilation at 388 ps | e37be8b4 (by date) | `macros/nema_iq.mac` | runs/nema_iq_388ps.sh | as 9a | as 9a | PLI_DIRECT_TAU_PS=388 | ../materials/previous/pli_materials.json (version 7, f6d5a53c) |
| 10 | oxygen (hypoxia) gradient, 5 inserts | oxygen sensitivity | e37be8b4 (logged) | `macros/oxygen_gradient.mac + geometry/phantom_hypoxia_gradient.mac + sources/scandium44_uniform_poisson.mac` | runs/oxygen_v8.sh | 4 batches x 16 shards, 62.5 s per batch (2.5e8) | 1,000,000 + 1000 x batch + shard (shard 1..K), batches 0-3 | none | ../materials/pli_materials.json (version 8, 23f607ea) (hash logged by the run) |
| 11 | brain tumour, 5 inserts | brain-tumour detectability | 3522d567 | `macros/brain_tumor.mac + geometry/phantom_brain_tumor_visibility.mac + sources/brain_tumor_sources_poisson.mac` | runs/closure_o4_driver.py (brain_tumor) | 16 batches x 16 shards, 61.2 s per batch (1e9) | 1,000,000 + 1000 x batch + shard (shard 1..K) | none | database 2a26b68e of the O4 reruns (not retained), and the May 1e9 acquisition that defines the regions, and its 176Lu-injected variant, ran on a May build that was not retained |
| 12 | VIP-Man thoracic slab, 1.4e9, 176Lu background injected | anthropomorphic example | 3522d567 | macros/vipman_slab.mac (needs PLI_GATE_MATERIALS_DB and PLI_VIPMAN_DIR) | runs/vipman_batches.py, with PLI_MATERIALS_JSON=$PLI_PKG/materials/pli_materials_vipman.json | 32 batches x 16 shards, 0.836562 s per batch | 1,000,000 + 1000 x batch + shard (shard 1..K), 176Lu generator 12,345 + batch, injection 42 + batch | PLI_MATERIALS_JSON (atlas database) | atlas database before v7 (not retained) |
| 13 | 176Lu singles injection (no GATE) | LSO background for #12 | n/a | lu176/lu176_bg_generator.py, lu176/lso_background_inject.py | called by runs/vipman_batches.py |  | see #12 |  |  |
| 14 | core-shell transport pairs (null control, pair A, pair B) | lifetime transport test | build of 2026-08-13 (not retained) | macros/core_shell_pairA.mac, macros/core_shell_pairB.mac (from macros/make_core_shell.py) | runs/core_shell.sh | 8 arms x 16 shards, 1e7 primaries per arm, batch 0 | 1,000,001-1,000,016 in every arm | PLI_LIFETIME_MODE=transport/source/fixed, PLI_FIXED_TAU_NS=2.0, PLI_LOG_XMAT=1 | database of 2026-08-13 (not retained). Null arm: ../materials/variants/r2a_equal_tau_materials.json, which this release cannot load |
| 15 | core-shell pairs A and B, matched yield | R2c matched-yield arms | e37be8b4 (logged) | as #14 | runs/core_shell_matched_yield.sh | 4 arms x 16 shards x 1e7 primaries | as #14 | as #14 plus PLI_DISABLE_3GAMMA=1 | ../materials/variants/r2c_equal_yield_materials.json (25e3441c, made by ../materials/variants/make_r2c_materials.py from v8) |
| 16 | Mercolli XAD-4 tubes, wet XAD-4 against a lookup table | lookup-arm comparison (May) | May build (not retained) | macros/mercolli_xad4_i124.mac with setTimeStop 0.787550 s | May orchestrator (not shipped) |  |  | PLI_MATERIALS_JSON=../materials/variants/lookup_baseline.json | ../materials/variants/lookup_baseline.json (80da91ba), which this release cannot load |
| 17 | discovery eight-tube stream (the earlier-yield stream of the paper) | replication only | May build (not retained) | not shipped (May archive macro sim1_te_poisson.mac) | not shipped |  |  |  | May database (not retained) |
| 18 | 12-tube closure, Sim-1 half | three-seed repeatability | build of 2026-08-12 (not retained) | `macros/closure_sim1_sc44.mac` | runs/seed_repeatability.py | 8 batches x 16 shards x 750,000 = 9.6e7 primaries per seed (the script docstring's 4 x 16 x 4.425 s = 250 M is wrong) | bases 2,000,000 / 3,000,000 / 4,000,000 + 1000 x batch + shard | none | database of 2026-08-12 (not retained) |
| 19 | VIP-Man slab, stock against custom annihilation, 1e6 primaries | Table S7 run time | build of 2026-08-11 (not retained) | macros/stock_vs_custom_pli.mac, macros/stock_vs_custom_stock.mac (need PLI_GATE_MATERIALS_DB and PLI_VIPMAN_DIR) | by hand: Gate <macro> | one process each | 424242 (in the macros) |  | not retained |
| 20 | annihilation-vertex unit test | Note S3 | e37be8b4 (by date) | macros/unit_vertex_custom.mac, macros/unit_vertex_stock.mac (the stock one needs PLI_GATE_MATERIALS_DB) | by hand: Gate <macro> (the original run_both.sh is not shipped) | 1e6 primaries each | in the macros | PLI_OUT for the phase-space output | ../materials/previous/pli_materials.json (version 7, f6d5a53c) |

Seeds: every production run uses the base 1,000,000, so seeds repeat across
studies (C0 batch 0 shard 1 and NEMA batch 0 shard 1 both use 1,000,001) and
are distinct within a study. The bases 2,000,000 to 4,000,000 belong only to
the seed study (#18).

The header comments of some macros are older than their last use. The
eight-tube macro still says "Primaries: 250M" and "1 s window", while the
C0, C2 and in-flight runs drove it in time mode at 0.137933 s per shard. It
differs from `closure_sim1_sc44.mac` only in the output name and in that time
setting. Such comments were left as they were, because only paths were edited.

## Registration records

The pre-run registrations of C0, C2, the in-flight arm, the density correction
and the matched-yield arms (`REGISTERED_C0.md`, `REGISTERED_C2.md`,
`REGISTERED_INFLIGHT.md`, `REGISTERED_density.md`, `REGISTERED_R2C.md`) are
kept with the other registration records in `../registration/`.
