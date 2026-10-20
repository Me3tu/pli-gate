# Edits made when the production files were copied into this release

Each file below is a copy of a production file from the project tree of the
workstation (`<PROJECT>`, with data on `<DATA>`). Only paths were changed,
apart from the exceptions listed next. The diff of every file follows, with
no context lines, so the unchanged rest of each file is the production text.

Applied to every file:

* Code files (`.sh`, `.py`, `.mac`, `.cc`, `CMakeLists.txt`) begin with two
  comment lines, `SPDX-License-Identifier: LGPL-3.0-only` and
  `Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)`, after the `#!` line when there is one.
  They are not repeated in the diffs. Data files (`.db`, `.json`) are unchanged.
* A macro that names a path through an alias begins, after those two lines,
  with a two-line comment and one `/control/getEnv <NAME>` line per alias it
  uses. These lines are part of the diffs.

Edits that are not path replacements:

* `runner/run_gate_parallel.sh`, `runner/run_gate_parallel_bin.sh`: `usage()` prints
  lines 4-34 of the file instead of 2-32, because of the two header lines.
* `runner/env_example.sh`: the two lines that put the standalone
  `libPLICustomPhysics` on `LD_LIBRARY_PATH`, and the message announcing it,
  are removed. Gate does not link that library (inventory: no
  `libPLICustomPhysics` in `readelf -d` of the production binary).
* `runs/c0_eight_tube.sh`: the wait for the in-flight arm and the call of its
  analysis are commented out (`#[v1.1.0 disabled]`). C0 does not depend on them.
* `runs/core_shell_matched_yield.sh`: the wait for the oxygen run is commented out.
* `runs/oxygen_v8.sh`: the `sed` that wrote the run macro from an archived macro
  is commented out, and `MAC` points at the shipped result of that `sed`,
  `macros/oxygen_gradient.mac`.
* Several run scripts create their log or output folder (`mkdir -p`), which
  existed on the workstation before the run. `closure_o4_driver.py` does the
  same for its output folder.
* `runs/closure_o4_driver.py`: the preflight check expands the GATE aliases from
  the environment and checks the paths named by `setMaterialDatabase`, `setImage`,
  `setRangeToMaterialFile`, `readTable` and `readFile`, instead of every
  `/home/...` string. `/control/execute` targets are followed as before.
* `macros/make_core_shell.py`: the macro template gains the header and
  `/control/getEnv` lines, and the macros are written as
  `core_shell_pair{A,B}.mac`. Its output equals the shipped macros byte for
  byte. Run from the current database it writes an `r2a_equal_tau_materials.json`
  that differs from the file in `../materials/variants/`, which was made
  from the database of 2026-08-13.
* `runs/c2_eight_tube_388ps.sh`, `runs/seed_repeatability.py`,
  `sources/iodine124_tubes_poisson.mac`: the wording of a comment (last section).
* `../integration/variants/pli_inflight_test/CMakeLists.txt`: the variant
  source is made at configure time by applying `../gate-inflight-stock.patch` to
  the shipped `PLIPhysicsList.cc`, and the bundled nlohmann/json is put on the
  include path (the production build used the system copy of the same version).

New files: `README.md`, `RUNS.json`, `EDITS.md`, `runner/paths.sh`.

Kept as they were: header comments that describe an earlier use of a file, and
comments that name production paths relative to the project.

## Files

| shipped file | production source | source sha256 | shipped sha256 |
|---|---|---|---|
| `simulation/runner/run_gate_parallel.sh` | `<PROJECT>/scripts/run_gate_parallel_musire.sh` | `b0b659fe137833bf` | `880d1b5f7af9c67c` |
| `simulation/runner/run_gate_parallel_bin.sh` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/run_gate_parallel_musire_bin.sh` | `13b0015c783bfa30` | `87dc6036a990ab06` |
| `simulation/runner/run_bounded.sh` | `<PROJECT>/scripts/run_bounded.sh` | `8fd9c31dd8de9ca0` | `ef0d57cd408e8615` |
| `simulation/runner/env_example.sh` | `<PROJECT>/scripts/setup/setup_environment_st.sh` | `dd3b1586bbbcd481` | `d74e884b9901a950` |
| `simulation/geometry/materials_pli.db` | `<PROJECT>/gate/geometry/materials_pli.db` | `5b6337c6b55ac908` | `5b6337c6b55ac908` |
| `simulation/geometry/materials_pli_phase2b.db` | `<PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db` | `3a7dbbe8c09d953f` | `3a7dbbe8c09d953f` |
| `simulation/geometry/phantom_mercolli_xad4.mac` | `<PROJECT>/gate/geometry/phantom_mercolli_xad4.mac` | `6a314f23ac07469e` | `38ee823a9b7f3c24` |
| `simulation/geometry/phantom_brain_tumor_visibility.mac` | `<PROJECT>/gate/geometry/phantom_brain_tumor_visibility.mac` | `b0aad9594e4bdd86` | `291633c45c7c9d09` |
| `simulation/geometry/phantom_hypoxia_gradient.mac` | `<PROJECT>/gate/geometry/phantom_hypoxia_gradient.mac` | `9f398596f29621ef` | `2e2f6461a7cf823c` |
| `simulation/sources/iodine124_tubes_poisson.mac` | `<PROJECT>/gate/macros/poisson_paired/iodine124_tubes_poisson.mac` | `9391ec7e8acb076a` | `6edcc998d5e13edc` |
| `simulation/sources/scandium44_uniform_poisson.mac` | `<PROJECT>/gate/macros/poisson_paired/scandium44_uniform_poisson.mac` | `ff6fe5fa0844e843` | `39d46d68da2dc206` |
| `simulation/sources/brain_tumor_sources_poisson.mac` | `<PROJECT>/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/brain_tumor_sources_poisson.mac` | `c29860a704bad160` | `c92cc2f64566c87a` |
| `simulation/macros/eight_tube_sc44.mac` | `<PROJECT>/reconstruction/results/2026-08-24_C2_sim1_388/sim1_sc44_388_timedriven.mac` | `95dc73617ca451e7` | `259dddefd1b4e2ae` |
| `simulation/macros/closure_sim1_sc44.mac` | `<PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/sim1_sc44_validation_phantom.mac` | `6cf0a901382f8132` | `d0e21c14be075324` |
| `simulation/macros/mercolli_xad4_i124.mac` | `<PROJECT>/reconstruction/results/2026-08-13_overnight/macros/mercolli_te_poisson.mac` | `4c6f82ca81f1f431` | `1b61749b479b37fb` |
| `simulation/macros/steinberger_disks.mac` | `<PROJECT>/reconstruction/results/2026-06-19_steinberger_validation/macros/steinberger_disk_phantom.mac` | `a59f39c1193b0d39` | `6c71c6bda80428e2` |
| `simulation/macros/water_cylinder_20mm.mac` | `<PROJECT>/reconstruction/results/2026-08-12_steinberger_allprompts/water_7MBq.mac` | `1a86c49750d3646b` | `dc00a7cb43127fdd` |
| `simulation/macros/water_cylinder_20cm.mac` | `<PROJECT>/reconstruction/results/2026-10-19_review16/water20cm/water20cm_7MBq.mac` | `f639f27780932a1a` | `d23c1afd4ea2358a` |
| `simulation/macros/nema_iq.mac` | `<PROJECT>/reconstruction/results/2026-04-24_1039_nema_iq/nema_iq_200M.mac` | `30815b06fd5df909` | `87f902a548cfcf23` |
| `simulation/macros/oxygen_gradient.mac` | `<PROJECT>/reconstruction/results/2026-10-08_oxygen_v8/hypoxia_te_v8_poisson.mac` | `ead4a55cafd4ba0c` | `5dc63eda7e1199dd` |
| `simulation/macros/brain_tumor.mac` | `<PROJECT>/reconstruction/results/2026-08-13_overnight/macros/brain_tumor_te_poisson_1B_truetime.mac` | `a05d70670e1e69bd` | `4ae9029d08fb2c93` |
| `simulation/macros/vipman_slab.mac` | `<PROJECT>/macros/vipman_sc44_1400M_poisson.mac` | `1fab5f8fb9427488` | `31ab5e5cd8d2e06c` |
| `simulation/macros/core_shell_pairA.mac` | `<PROJECT>/reconstruction/results/2026-08-13_referee_response/r2_pairA.mac` | `1597ba57479b14c7` | `c2cce7a58b09642f` |
| `simulation/macros/core_shell_pairB.mac` | `<PROJECT>/reconstruction/results/2026-08-13_referee_response/r2_pairB.mac` | `a66dacbc6afd3923` | `c253feaedbe50fa9` |
| `simulation/macros/stock_vs_custom_pli.mac` | `<PROJECT>/reconstruction/results/2026-08-11_ejmp_figure_regen/runtime_overhead/timing_pli_eplus.mac` | `0589589a0a263529` | `0d733a77536cbf70` |
| `simulation/macros/stock_vs_custom_stock.mac` | `<PROJECT>/reconstruction/results/2026-08-11_ejmp_figure_regen/runtime_overhead/timing_stock_eplus.mac` | `2fcc8023bcf97b90` | `8dac189973a80e15` |
| `simulation/macros/unit_vertex_custom.mac` | `<PROJECT>/reconstruction/results/2026-08-30_S5_annihilation_vertex/run_custom_1M.mac` | `232d499d8cc4ece5` | `b2ae3af1dcda0ee2` |
| `simulation/macros/unit_vertex_stock.mac` | `<PROJECT>/reconstruction/results/2026-08-30_S5_annihilation_vertex/run_stock_1M.mac` | `ebdbdb12e23e86d6` | `0b61b6a86289e81b` |
| `simulation/macros/make_core_shell.py` | `<PROJECT>/reconstruction/results/2026-08-13_referee_response/make_r2_inputs.py` | `74e759c63a886eec` | `18cff9a4af64319a` |
| `simulation/runs/c0_eight_tube.sh` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/c0_sim1/run_c0_chain.sh` | `13b4936f1f52b9ee` | `41f55d4d84087243` |
| `simulation/runs/c2_eight_tube_388ps.sh` | `<PROJECT>/reconstruction/results/2026-08-24_C2_sim1_388/run_c2.sh` | `5acc565addec17c1` | `7332f4c92830a5d1` |
| `simulation/runs/inflight_arm.sh` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/run_hybrid_c2.sh` | `2edfacbc58a9b87b` | `6a68184c71ad88ca` |
| `simulation/runs/mercolli_replication.sh` | `<PROJECT>/reconstruction/results/2026-08-22_G1_exposure/run_exposure.sh` | `84d9d7f41c420e71` | `8c53f2ca8e73146a` |
| `simulation/runs/water_20cm.sh` | `<PROJECT>/reconstruction/results/2026-10-19_review16/water20cm/driver.sh` | `3b517ef173d60622` | `67645bea52f8dc00` |
| `simulation/runs/nema_iq_0ps.sh` | `<PROJECT>/reconstruction/results/2026-08-11_nema_iec_250M/launch_nema_250M.sh` | `58235f5e6588d198` | `ca5ef4c851008cc0` |
| `simulation/runs/nema_iq_388ps.sh` | `<PROJECT>/reconstruction/results/2026-08-24_C1_nema388/run_c1.sh` | `08aa0967f22a96e0` | `dc691950d51305cb` |
| `simulation/runs/oxygen_v8.sh` | `<PROJECT>/reconstruction/results/2026-10-08_oxygen_v8/run_oxygen_v8.sh` | `348a2a17e9778bf6` | `545189834750f6d0` |
| `simulation/runs/core_shell.sh` | `<PROJECT>/reconstruction/results/2026-08-13_referee_response/run_r2.sh` | `98f4fa569048f3d6` | `074e779074f40351` |
| `simulation/runs/core_shell_matched_yield.sh` | `<PROJECT>/reconstruction/results/2026-10-08_r2c_matched_yield/run_r2c.sh` | `5f4d5e9422250e20` | `f51b281f971988df` |
| `simulation/runs/closure_o4_driver.py` | `<PROJECT>/reconstruction/results/2026-08-13_overnight/o4_batched_driver.py` | `817ae78e04ef76ec` | `31936e276a082e01` |
| `simulation/runs/vipman_batches.py` | `<PROJECT>/reconstruction/utils/run_phase3q_per_batch.py` | `f270dad8a77b8a64` | `c406812b74ee433f` |
| `simulation/runs/seed_repeatability.py` | `<PROJECT>/reconstruction/results/2026-08-12_seed_repeatability/run_seed_repeatability.py` | `0f073af591492faf` | `65277c1e4ede3e2b` |
| `simulation/lu176/lu176_bg_generator.py` | `<PROJECT>/reconstruction/utils/lu176_bg_generator.py` | `03f9953d89514546` | `d2e105d46d8471f3` |
| `simulation/lu176/lso_background_inject.py` | `<PROJECT>/reconstruction/utils/lso_background_inject.py` | `a1458848c092429c` | `b7fc5e0f868ddba5` |
| `materials/variants/make_r2c_materials.py` | `<PROJECT>/reconstruction/results/2026-10-08_r2c_matched_yield/make_r2c_materials.py` | `1b9eab6ea03ea340` | `ca168d31b898c335` |
| `integration/variants/build_gate_inflight.sh` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/build_gate_inflight.sh` | `2ffa3e7f34d6704b` | `3e19511e51b456e8` |
| `integration/variants/pli_inflight_test/pli_inflight_test.cc` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/pli_inflight_test.cc` | `67f3b81124c59ae4` | `c1f790e4a25b9324` |
| `integration/variants/pli_inflight_test/CMakeLists.txt` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/CMakeLists.txt` | `24b4507af0eaed0f` | `10bccf4cb6222ec7` |
| `integration/variants/pli_inflight_test/run_physics_level.sh` | `<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/run_physics_level.sh` | `7f03ffeca753bf58` | `1b612ba8566befba` |

## Diffs

### `simulation/runner/run_gate_parallel.sh`

```diff
--- original
+++ shipped
@@ -57 +57 @@
-    sed -n '2,32p' "${BASH_SOURCE[0]}"
+    sed -n '4,34p' "${BASH_SOURCE[0]}"
@@ -125 +125 @@
-if [[ -f "$PROJECT_ROOT/scripts/setup/setup_environment_st.sh" ]]; then
+if [[ -f "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}" ]]; then
@@ -127 +127 @@
-    source "$PROJECT_ROOT/scripts/setup/setup_environment_st.sh"
+    source "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}"
@@ -142 +142 @@
-        VIPMAN_JSON="$PROJECT_ROOT/materials/pli_materials_vipman.json"
+        VIPMAN_JSON="$PROJECT_ROOT/../materials/pli_materials_vipman.json"
@@ -284 +284 @@
-        "$PROJECT_ROOT/scripts/run_bounded.sh" \
+        "$PROJECT_ROOT/runner/run_bounded.sh" \
@@ -360 +360 @@
-if [[ -f "$PROJECT_ROOT/venv/bin/activate" ]]; then
+if [[ -f "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate" ]]; then
@@ -362 +362 @@
-    source "$PROJECT_ROOT/venv/bin/activate"
+    source "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate"
```

### `simulation/runner/run_gate_parallel_bin.sh`

```diff
--- original
+++ shipped
@@ -43 +43 @@
-PROJECT_ROOT="<PROJECT>"  # feedback-3 copy lives outside scripts/
+PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # v1.1.0: this copy lives in runner/
@@ -57 +57 @@
-    sed -n '2,32p' "${BASH_SOURCE[0]}"
+    sed -n '4,34p' "${BASH_SOURCE[0]}"
@@ -125 +125 @@
-if [[ -f "$PROJECT_ROOT/scripts/setup/setup_environment_st.sh" ]]; then
+if [[ -f "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}" ]]; then
@@ -127 +127 @@
-    source "$PROJECT_ROOT/scripts/setup/setup_environment_st.sh"
+    source "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}"
@@ -144 +144 @@
-        VIPMAN_JSON="$PROJECT_ROOT/materials/pli_materials_vipman.json"
+        VIPMAN_JSON="$PROJECT_ROOT/../materials/pli_materials_vipman.json"
@@ -286 +286 @@
-        "$PROJECT_ROOT/scripts/run_bounded.sh" \
+        "$PROJECT_ROOT/runner/run_bounded.sh" \
@@ -362 +362 @@
-if [[ -f "$PROJECT_ROOT/venv/bin/activate" ]]; then
+if [[ -f "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate" ]]; then
@@ -364 +364 @@
-    source "$PROJECT_ROOT/venv/bin/activate"
+    source "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate"
```

### `simulation/runner/run_bounded.sh`

```diff
--- original
+++ shipped
@@ -29 +29 @@
-LOGDIR="$HOME/research/pli-quadra-simulation/logs/bounded"
+LOGDIR="${PLI_BOUNDED_LOGDIR:-${TMPDIR:-/tmp}/pli-bounded}"
```

### `simulation/runner/env_example.sh`

```diff
--- original
+++ shipped
@@ -4 +4,3 @@
-export PLI_PROJECT_ROOT=${HOME}/research/pli-quadra-simulation
+# v1.1.0: install locations come from the environment; the defaults follow the
+# one-prefix layout of docs/BOOTSTRAP.md ($PLI_STACK/{geant4,root,gate}).
+export PLI_STACK=${PLI_STACK:-${HOME}/pli-stack/install}
@@ -7,3 +9,3 @@
-export GEANT4_INSTALL=${PLI_PROJECT_ROOT}/software/install/geant4-st
-export ROOT_INSTALL=${PLI_PROJECT_ROOT}/software/install/root
-export GATE_INSTALL=${PLI_PROJECT_ROOT}/software/install/gate-st
+export GEANT4_INSTALL=${GEANT4_INSTALL:-${PLI_STACK}/geant4}
+export ROOT_INSTALL=${ROOT_INSTALL:-${PLI_STACK}/root}
+export GATE_INSTALL=${GATE_INSTALL:-${PLI_STACK}/gate}
@@ -26,4 +27,0 @@
-# Add PLI Custom Physics library
-export PLI_PHYSICS_INSTALL=${PLI_PROJECT_ROOT}/software/install/pli_physics
-export LD_LIBRARY_PATH=${PLI_PHYSICS_INSTALL}/lib:${LD_LIBRARY_PATH}
-
@@ -31 +28,0 @@
-echo "✓ PLI Physics library added"
```

### `simulation/geometry/materials_pli.db`

Byte-identical copy.

### `simulation/geometry/materials_pli_phase2b.db`

Byte-identical copy.

### `simulation/geometry/phantom_mercolli_xad4.mac`

No edit apart from the header lines.

### `simulation/geometry/phantom_brain_tumor_visibility.mac`

No edit apart from the header lines.

### `simulation/geometry/phantom_hypoxia_gradient.mac`

No edit apart from the header lines.

### `simulation/sources/iodine124_tubes_poisson.mac`

No edit apart from the header lines.

### `simulation/sources/scandium44_uniform_poisson.mac`

No edit apart from the header lines.

### `simulation/sources/brain_tumor_sources_poisson.mac`

No edit apart from the header lines.

### `simulation/macros/eight_tube_sc44.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -36,2 +40,2 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
-/gate/geometry/setMaterialDatabase <PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli_phase2b.db
@@ -40 +44 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/closure_sim1_sc44.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -36,2 +40,2 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
-/gate/geometry/setMaterialDatabase <PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli_phase2b.db
@@ -40 +44 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/mercolli_xad4_i124.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -33 +37 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -36 +40 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -39 +43 @@
-/control/execute <PROJECT>/gate/geometry/phantom_mercolli_xad4.mac
+/control/execute {PLI_SIM}/geometry/phantom_mercolli_xad4.mac
@@ -102 +106 @@
-/control/execute <PROJECT>/gate/macros/poisson_paired/iodine124_tubes_poisson.mac
+/control/execute {PLI_SIM}/sources/iodine124_tubes_poisson.mac
```

### `simulation/macros/steinberger_disks.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -30 +34 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -35 +39 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/water_cylinder_20mm.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -33 +37 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -35 +39 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/water_cylinder_20cm.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -41 +45 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -43 +47 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/nema_iq.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_GATE_MATERIALS_DB
+/control/getEnv PLI_SCANNER_MACRO
@@ -39 +43 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/phantom/gate-materials.db
+/gate/geometry/setMaterialDatabase {PLI_GATE_MATERIALS_DB}
@@ -42 +46 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -250 +254 @@
-/gate/output/root/setFileName reconstruction/results/2026-04-24_1039_nema_iq/gate_output/sim_nema_iq_200M
+/gate/output/root/setFileName sim_nema_iq_200M
```

### `simulation/macros/oxygen_gradient.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -29 +33 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -32 +36 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -35 +39 @@
-/control/execute <PROJECT>/gate/geometry/phantom_hypoxia_gradient.mac
+/control/execute {PLI_SIM}/geometry/phantom_hypoxia_gradient.mac
@@ -71 +75 @@
-/control/execute <PROJECT>/gate/macros/poisson_paired/scandium44_uniform_poisson.mac
+/control/execute {PLI_SIM}/sources/scandium44_uniform_poisson.mac
```

### `simulation/macros/brain_tumor.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -28,2 +32,2 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
-/gate/geometry/setMaterialDatabase <PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli_phase2b.db
@@ -31,2 +35,2 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
-/control/execute <PROJECT>/gate/geometry/phantom_brain_tumor_visibility.mac
+/control/execute {PLI_SCANNER_MACRO}
+/control/execute {PLI_SIM}/geometry/phantom_brain_tumor_visibility.mac
@@ -60 +64 @@
-/control/execute <PROJECT>/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/brain_tumor_sources_poisson.mac
+/control/execute {PLI_SIM}/sources/brain_tumor_sources_poisson.mac
```

### `simulation/macros/vipman_slab.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,5 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_GATE_MATERIALS_DB
+/control/getEnv PLI_SCANNER_MACRO
+/control/getEnv PLI_VIPMAN_DIR
@@ -27 +32 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/phantom/gate-materials.db
+/gate/geometry/setMaterialDatabase {PLI_GATE_MATERIALS_DB}
@@ -30 +35 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -35,2 +40,2 @@
-/gate/vipman/geometry/setImage                      <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
-/gate/vipman/geometry/setRangeToMaterialFile        <PROJECT>/phantom/vipman-atlas-gate-materials.dat
+/gate/vipman/geometry/setImage                      {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/vipman/geometry/setRangeToMaterialFile        {PLI_VIPMAN_DIR}/vipman-atlas-gate-materials.dat
@@ -68,2 +73,2 @@
-/gate/source/vox/imageReader/rangeTranslator/readTable <PROJECT>/phantom/vipman-atlas-gate-activities-f18-fdg.dat
-/gate/source/vox/imageReader/readFile               <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/source/vox/imageReader/rangeTranslator/readTable {PLI_VIPMAN_DIR}/vipman-atlas-gate-activities-f18-fdg.dat
+/gate/source/vox/imageReader/readFile               {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
```

### `simulation/macros/core_shell_pairA.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -11,2 +15,2 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
-/gate/geometry/setMaterialDatabase <PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli_phase2b.db
@@ -15 +19 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/core_shell_pairB.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,4 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
@@ -11,2 +15,2 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
-/gate/geometry/setMaterialDatabase <PROJECT>/reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros/materials_pli_phase2b.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli_phase2b.db
@@ -15 +19 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
```

### `simulation/macros/stock_vs_custom_pli.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,5 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_GATE_MATERIALS_DB
+/control/getEnv PLI_SCANNER_MACRO
+/control/getEnv PLI_VIPMAN_DIR
@@ -29 +34 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/phantom/gate-materials.db
+/gate/geometry/setMaterialDatabase {PLI_GATE_MATERIALS_DB}
@@ -32 +37 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -37,2 +42,2 @@
-/gate/vipman/geometry/setImage                      <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
-/gate/vipman/geometry/setRangeToMaterialFile        <PROJECT>/phantom/vipman-atlas-gate-materials.dat
+/gate/vipman/geometry/setImage                      {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/vipman/geometry/setRangeToMaterialFile        {PLI_VIPMAN_DIR}/vipman-atlas-gate-materials.dat
@@ -70,2 +75,2 @@
-/gate/source/vox/imageReader/rangeTranslator/readTable <PROJECT>/phantom/vipman-atlas-gate-activities-f18-fdg.dat
-/gate/source/vox/imageReader/readFile               <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/source/vox/imageReader/rangeTranslator/readTable {PLI_VIPMAN_DIR}/vipman-atlas-gate-activities-f18-fdg.dat
+/gate/source/vox/imageReader/readFile               {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
```

### `simulation/macros/stock_vs_custom_stock.mac`

```diff
--- original
+++ shipped
@@ -0,0 +1,5 @@
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_GATE_MATERIALS_DB
+/control/getEnv PLI_SCANNER_MACRO
+/control/getEnv PLI_VIPMAN_DIR
@@ -29 +34 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/phantom/gate-materials.db
+/gate/geometry/setMaterialDatabase {PLI_GATE_MATERIALS_DB}
@@ -32 +37 @@
-/control/execute <PROJECT>/gate/geometry/quadra_detector.mac
+/control/execute {PLI_SCANNER_MACRO}
@@ -37,2 +42,2 @@
-/gate/vipman/geometry/setImage                      <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
-/gate/vipman/geometry/setRangeToMaterialFile        <PROJECT>/phantom/vipman-atlas-gate-materials.dat
+/gate/vipman/geometry/setImage                      {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/vipman/geometry/setRangeToMaterialFile        {PLI_VIPMAN_DIR}/vipman-atlas-gate-materials.dat
@@ -70,2 +75,2 @@
-/gate/source/vox/imageReader/rangeTranslator/readTable <PROJECT>/phantom/vipman-atlas-gate-activities-f18-fdg.dat
-/gate/source/vox/imageReader/readFile               <PROJECT>/phantom/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
+/gate/source/vox/imageReader/rangeTranslator/readTable {PLI_VIPMAN_DIR}/vipman-atlas-gate-activities-f18-fdg.dat
+/gate/source/vox/imageReader/readFile               {PLI_VIPMAN_DIR}/vipman-atlas-lungs-32mm-1mm-cropped-fixed.mhd
```

### `simulation/macros/unit_vertex_custom.mac`

```diff
--- original
+++ shipped
@@ -1 +1,5 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials_pli.db
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_OUT
+/gate/geometry/setMaterialDatabase {PLI_SIM}/geometry/materials_pli.db
@@ -22 +26 @@
-/gate/actor/PS/save <PROJECT>/reconstruction/results/2026-08-30_S5_annihilation_vertex/ps_custom_1M.root
+/gate/actor/PS/save {PLI_OUT}/ps_custom_1M.root
```

### `simulation/macros/unit_vertex_stock.mac`

```diff
--- original
+++ shipped
@@ -1 +1,5 @@
-/gate/geometry/setMaterialDatabase <PROJECT>/gate/geometry/materials.db
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_GATE_MATERIALS_DB
+/control/getEnv PLI_OUT
+/gate/geometry/setMaterialDatabase {PLI_GATE_MATERIALS_DB}
@@ -22 +26 @@
-/gate/actor/PS/save <PROJECT>/reconstruction/results/2026-08-30_S5_annihilation_vertex/ps_stock_1M.root
+/gate/actor/PS/save {PLI_OUT}/ps_stock_1M.root
```

### `simulation/macros/make_core_shell.py`

```diff
--- original
+++ shipped
@@ -34,4 +34,4 @@
-ROOT = Path("<PROJECT>")
-HERE = ROOT / "reconstruction/results/2026-08-13_referee_response"
-PH2B = ROOT / ("reconstruction/results/2026-04-29_2320_phase2b_validation_sims/"
-               "macros/materials_pli_phase2b.db")
+# v1.1.0: paths relative to the release; GATE aliases resolved at run time
+ROOT = Path(__file__).resolve().parents[2]          # release root
+HERE = Path(__file__).resolve().parent              # simulation/macros/
+PH2B = "{PLI_SIM}/geometry/materials_pli_phase2b.db"
@@ -47 +47,7 @@
-MACRO = """# ==========================================================================
+MACRO = """# SPDX-License-Identifier: LGPL-3.0-only
+# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
+# Paths (release v1.1.0): the aliases below are read from environment variables,
+# set by simulation/runner/paths.sh or by hand. See simulation/README.md.
+/control/getEnv PLI_SIM
+/control/getEnv PLI_SCANNER_MACRO
+# ==========================================================================
@@ -159,2 +165,2 @@
-    quadra = ROOT / "gate/geometry/quadra_detector.mac"
-    mat_main = ROOT / "gate/geometry/materials_pli.db"
+    quadra = "{PLI_SCANNER_MACRO}"
+    mat_main = "{PLI_SIM}/geometry/materials_pli.db"
@@ -168 +174 @@
-        out = HERE / f"r2_pair{pair}.mac"
+        out = HERE / f"core_shell_pair{pair}.mac"
```

### `simulation/runs/c0_eight_tube.sh`

```diff
--- original
+++ shipped
@@ -6,2 +6,2 @@
-R=<PROJECT>
-F3=$R/reconstruction/results/2026-09-30_feedback3
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+F3=$PLI_OUT/results/2026-09-30_feedback3
@@ -9,3 +9,3 @@
-OUT=<DATA>/feedback3_c0/c0_sim1_0ps
-MACRO=$R/reconstruction/results/2026-08-24_C2_sim1_388/sim1_sc44_388_timedriven.mac
-S=$R/reconstruction/results/2026-08-17_derivation_test
+OUT=$PLI_OUT/feedback3_c0/c0_sim1_0ps
+MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
+S=$PLI_ANALYSIS/results/2026-08-17_derivation_test
@@ -15,7 +15,8 @@
-say "chain started; waiting for HYB_DONE"
-while ! grep -q "HYB_DONE\|HYB_FAILED" "$I/run.log"; do sleep 60; done
-grep -q HYB_FAILED "$I/run.log" && { say "hybrid FAILED; C0 not started"; exit 1; }
-say "hybrid done; in-flight analysis"
-bash "$I/run_inflight_analysis.sh"; say "in-flight analysis rc=$?"
-cd "$R"; source scripts/setup/setup_environment_st.sh >/dev/null 2>&1; source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+# v1.1.0: the wait for the in-flight arm and its analysis are disabled; C0 does not depend on them.
+#[v1.1.0 disabled] say "chain started; waiting for HYB_DONE"
+#[v1.1.0 disabled] while ! grep -q "HYB_DONE\|HYB_FAILED" "$I/run.log"; do sleep 60; done
+#[v1.1.0 disabled] grep -q HYB_FAILED "$I/run.log" && { say "hybrid FAILED; C0 not started"; exit 1; }
+#[v1.1.0 disabled] say "hybrid done; in-flight analysis"
+#[v1.1.0 disabled] bash "$I/run_inflight_analysis.sh"; say "in-flight analysis rc=$?"
+cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -30 +31 @@
-  scripts/run_gate_parallel_musire.sh --macro "$MACRO" --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5 \
+  runner/run_gate_parallel.sh --macro "$MACRO" --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5 \
@@ -33,2 +34,2 @@
-  grep -q "gate-inflight" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo C0_FAILED >> "$L/run.log"; exit 1; }
-  python "$R/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
+  grep -q "${GATE_INFLIGHT_BIN_DIR:-gate-inflight}" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo C0_FAILED >> "$L/run.log"; exit 1; }
+  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
@@ -39 +40 @@
-  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
+  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
@@ -42 +43 @@
-python "$R/reconstruction/results/2026-08-16_T6_nema/concat_singles_u16.py" --in-dir "$OUT" --out "$OUT/singles_c0_sim1_0ps.h5" >> "$L/concat.log" 2>&1
+python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" --in-dir "$OUT" --out "$OUT/singles_c0_sim1_0ps.h5" >> "$L/concat.log" 2>&1
@@ -52 +53 @@
-python "$R/reconstruction/results/2026-08-19_V4_stats/dump_shard_hists.py" --singles "$M" --config "$S/sim1_voi_config.json" \
+python "$PLI_ANALYSIS/results/2026-08-19_V4_stats/dump_shard_hists.py" --singles "$M" --config "$S/sim1_voi_config.json" \
```

### `simulation/runs/c2_eight_tube_388ps.sh`

```diff
--- original
+++ shipped
@@ -15,4 +15,4 @@
-R=<PROJECT>
-L=$R/reconstruction/results/2026-08-24_C2_sim1_388
-OUT=<DATA>/c2_sim1_388
-MACRO=$R/reconstruction/results/2026-08-24_C2_sim1_388/sim1_sc44_388_timedriven.mac
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+L=$PLI_OUT/results/2026-08-24_C2_sim1_388
+OUT=$PLI_OUT/c2_sim1_388
+MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
@@ -29,2 +29,2 @@
-cd "$R"; source scripts/setup/setup_environment_st.sh >/dev/null 2>&1; source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -33 +33 @@
-mkdir -p "$OUT"
+mkdir -p "$OUT" "$L"
@@ -40 +40 @@
-  scripts/run_gate_parallel_musire.sh \
+  runner/run_gate_parallel.sh \
@@ -45 +45 @@
-  python "$R/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
+  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
@@ -67 +67 @@
-  say "  batch $bd/$((NBATCH-1)) done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
+  say "  batch $bd/$((NBATCH-1)) done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
@@ -70 +70 @@
-python "$R/reconstruction/results/2026-08-16_T6_nema/concat_singles_u16.py" \
+python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
```

### `simulation/runs/inflight_arm.sh`

```diff
--- original
+++ shipped
@@ -4,4 +4,4 @@
-R=<PROJECT>
-L=<PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight
-OUT=<DATA>/feedback3_inflight/c2_hybrid
-MACRO=$R/reconstruction/results/2026-08-24_C2_sim1_388/sim1_sc44_388_timedriven.mac
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+L=$PLI_OUT/results/2026-09-30_feedback3/inflight
+OUT=$PLI_OUT/feedback3_inflight/c2_hybrid
+MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
@@ -9,2 +9,2 @@
-cd "$R"; source scripts/setup/setup_environment_st.sh >/dev/null 2>&1; source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -13 +13 @@
-export GATE_BIN_DIR=$R/software/install/gate-inflight/bin
+export GATE_BIN_DIR=${GATE_INFLIGHT_BIN_DIR:?set GATE_INFLIGHT_BIN_DIR to the bin/ directory of the in-flight variant build (integration/PATCH_NOTES.md)}
@@ -15 +15 @@
-mkdir -p "$OUT"
+mkdir -p "$OUT" "$L"
@@ -23 +23 @@
-  $L/run_gate_parallel_musire_bin.sh       --macro "$MACRO"       --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5       --batch-id $b --label c2_hybrid --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
+  "$PLI_SIM/runner/run_gate_parallel_bin.sh"       --macro "$MACRO"       --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5       --batch-id $b --label c2_hybrid --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
@@ -25,2 +25,2 @@
-  grep -q "gate-inflight" "$L/gate.log" || { say "ABORT: runner did not use the hybrid binary"; echo HYB_FAILED >> "$L/run.log"; exit 1; }
-  python "$R/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py"       --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b       --keep-511 435 870 --keep-prompt 870 1300       --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
+  grep -q "$GATE_BIN_DIR" "$L/gate.log" || { say "ABORT: runner did not use the hybrid binary"; echo HYB_FAILED >> "$L/run.log"; exit 1; }
+  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py"       --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b       --keep-511 435 870 --keep-prompt 870 1300       --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
@@ -29 +29 @@
-  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
+  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
@@ -32 +32 @@
-python "$R/reconstruction/results/2026-08-16_T6_nema/concat_singles_u16.py"     --in-dir "$OUT" --out "$OUT/singles_c2_hybrid.h5" >> "$L/concat.log" 2>&1
+python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py"     --in-dir "$OUT" --out "$OUT/singles_c2_hybrid.h5" >> "$L/concat.log" 2>&1
```

### `simulation/runs/mercolli_replication.sh`

```diff
--- original
+++ shipped
@@ -5,5 +5,5 @@
-R=<PROJECT>
-L=$R/reconstruction/results/2026-08-22_G1_exposure
-OUT=<DATA>/g1_exposure   # 30 batches + concat need ~55 GB; / has only 47
-cd "$R"; source scripts/setup/setup_environment_st.sh >/dev/null 2>&1; source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+L=$PLI_OUT/results/2026-08-22_G1_exposure
+OUT=$PLI_OUT/g1_exposure   # 30 batches + concat need ~55 GB; / has only 47
+cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -15 +15 @@
-mkdir -p "$OUT"
+mkdir -p "$OUT" "$L"
@@ -20,2 +20,2 @@
-  scripts/run_gate_parallel_musire.sh \
-      --macro "$R/reconstruction/results/2026-08-13_overnight/macros/mercolli_te_poisson.mac" \
+  runner/run_gate_parallel.sh \
+      --macro "$PLI_SIM/macros/mercolli_xad4_i124.mac" \
@@ -34 +34 @@
-  python "$R/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
+  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
@@ -54 +54 @@
-  say "  batch $bd/29 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1)  free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
+  say "  batch $bd/29 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1)  free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
@@ -57 +57 @@
-python "$R/reconstruction/results/2026-08-16_T6_nema/concat_singles_u16.py" \
+python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
```

### `simulation/runs/water_20cm.sh`

```diff
--- original
+++ shipped
@@ -5,5 +5,7 @@
-cd <PROJECT>
-source scripts/setup/setup_environment_st.sh >/dev/null 2>&1
-export PLI_MATERIALS_JSON=$PWD/materials/pli_materials.json
-OUT=$PWD/reconstruction/results/2026-10-19_review16/water20cm
-CONV=$PWD/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+cd "$PLI_SIM"
+source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
+OUT=$PLI_OUT/results/2026-10-19_review16/water20cm
+CONV=$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py
+mkdir -p "$OUT"
@@ -12 +14 @@
-  scripts/run_gate_parallel_musire.sh --macro $OUT/water20cm_7MBq.mac \
+  runner/run_gate_parallel.sh --macro $PLI_SIM/macros/water_cylinder_20cm.mac \
@@ -18 +20 @@
-  ( source venv/bin/activate; python3 $CONV --root $OUT/sim/batch_$b/sim-merged.root \
+  ( [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"; python3 $CONV --root $OUT/sim/batch_$b/sim-merged.root \
```

### `simulation/runs/nema_iq_0ps.sh`

```diff
--- original
+++ shipped
@@ -7,4 +7,6 @@
-cd <PROJECT>
-source scripts/setup/setup_environment_st.sh >/dev/null 2>&1 || true
-export PLI_MATERIALS_JSON=$PWD/materials/pli_materials.json
-OUT=$PWD/reconstruction/results/2026-08-11_nema_iec_250M
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+cd "$PLI_SIM"
+source "$PLI_ENV_SCRIPT" >/dev/null 2>&1 || true
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
+OUT=$PLI_OUT/results/2026-08-11_nema_iec_250M
+mkdir -p "$OUT"
@@ -12,2 +14,2 @@
-  scripts/run_gate_parallel_musire.sh \
-    --macro reconstruction/results/2026-04-24_1039_nema_iq/nema_iq_200M.mac \
+  runner/run_gate_parallel.sh \
+    --macro macros/nema_iq.mac \
```

### `simulation/runs/nema_iq_388ps.sh`

```diff
--- original
+++ shipped
@@ -12,5 +12,5 @@
-R=<PROJECT>
-L=$R/reconstruction/results/2026-08-24_C1_nema388
-OUT=<DATA>/c1_nema388
-cd "$R"; source scripts/setup/setup_environment_st.sh >/dev/null 2>&1; source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+L=$PLI_OUT/results/2026-08-24_C1_nema388
+OUT=$PLI_OUT/c1_nema388
+cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -19 +19 @@
-mkdir -p "$OUT"
+mkdir -p "$OUT" "$L"
@@ -25,2 +25,2 @@
-  scripts/run_gate_parallel_musire.sh \
-      --macro "$R/reconstruction/results/2026-04-24_1039_nema_iq/nema_iq_200M.mac" \
+  runner/run_gate_parallel.sh \
+      --macro "$PLI_SIM/macros/nema_iq.mac" \
@@ -30 +30 @@
-  python "$R/reconstruction/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
+  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
@@ -46 +46 @@
-  say "  batch $bd/15 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
+  say "  batch $bd/15 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
@@ -49 +49 @@
-python "$R/reconstruction/results/2026-08-16_T6_nema/concat_singles_u16.py" \
+python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
```

### `simulation/runs/oxygen_v8.sh`

```diff
--- original
+++ shipped
@@ -11,4 +11,4 @@
-R=<PROJECT>
-L=$R/reconstruction/results/2026-10-08_oxygen_v8
-OUT=<DATA>/oxygen_v8
-A=<DATA>/pli-quadra-archive-2026-08-12/2026-05-06_2258_phase3s_2_paired_reruns
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+L=$PLI_OUT/results/2026-10-08_oxygen_v8
+OUT=$PLI_OUT/oxygen_v8
+A=$PLI_ANALYSIS/pli-quadra-archive-2026-08-12/2026-05-06_2258_phase3s_2_paired_reruns
@@ -17,4 +17,4 @@
-cd "$R"
-source scripts/setup/setup_environment_st.sh >/dev/null 2>&1
-source venv/bin/activate
-export PLI_MATERIALS_JSON=$R/materials/pli_materials.json
+cd "$PLI_SIM"
+source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
+[ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
@@ -22,7 +22,8 @@
-export PYTHONPATH=/usr/lib/python3.14/site-packages:${PYTHONPATH:-}
-export PARALLELPROJ_C_LIB=/usr/local/lib/libparallelproj_c.so PARALLELPROJ_CUDA_LIB=/usr/local/lib/libparallelproj_cuda.so
-say "start; materials $(sha256sum $PLI_MATERIALS_JSON | cut -c1-12); Gate $(command -v Gate) $(sha256sum $(command -v Gate) | cut -c1-12); free=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')GB"
-MAC=$L/hypoxia_te_v8_poisson.mac
-sed -e 's|^/gate/output/root/setFileName .*|/gate/output/root/setFileName hypoxia_te_v8|' \
-    -e 's|^# Phase 3s-2 — Hypoxia te arm Poisson master mac (auto-gen)|# Review round 8 oxygen rerun (database v8); copy of the Phase 3s-2 hypoxia te master mac|' \
-    "$A/hypoxia_te_v2/hypoxia_te_poisson.mac" > "$MAC"
+export PYTHONPATH=${PLI_EXTRA_PYTHONPATH:+$PLI_EXTRA_PYTHONPATH:}${PYTHONPATH:-}
+export PARALLELPROJ_C_LIB=${PARALLELPROJ_C_LIB:-/usr/local/lib/libparallelproj_c.so} PARALLELPROJ_CUDA_LIB=${PARALLELPROJ_CUDA_LIB:-/usr/local/lib/libparallelproj_cuda.so}
+say "start; materials $(sha256sum $PLI_MATERIALS_JSON | cut -c1-12); Gate $(command -v Gate) $(sha256sum $(command -v Gate) | cut -c1-12); free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
+# v1.1.0: the macro this sed generated from the archived Phase 3s-2 macro is shipped as macros/oxygen_gradient.mac
+MAC=$PLI_SIM/macros/oxygen_gradient.mac
+#[v1.1.0 disabled] sed -e 's|^/gate/output/root/setFileName .*|/gate/output/root/setFileName hypoxia_te_v8|' \
+#[v1.1.0 disabled]     -e 's|^# Phase 3s-2 — Hypoxia te arm Poisson master mac (auto-gen)|# Review round 8 oxygen rerun (database v8); copy of the Phase 3s-2 hypoxia te master mac|' \
+#[v1.1.0 disabled]     "$A/hypoxia_te_v2/hypoxia_te_poisson.mac" > "$MAC"
@@ -34 +35 @@
-  scripts/run_gate_parallel_musire.sh --macro "$MAC" --total-time-sec 62.5 --threads 16 --bounded-max-gb 4 \
+  runner/run_gate_parallel.sh --macro "$MAC" --total-time-sec 62.5 --threads 16 --bounded-max-gb 4 \
@@ -37 +38 @@
-  grep -q "gate-inflight" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
+  grep -q "${GATE_INFLIGHT_BIN_DIR:-gate-inflight}" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
@@ -40 +41 @@
-$R/software/install/root/bin/hadd -f "$OUT/sim/sim-merged.root" $OUT/sim/batch_0*/sim-merged.root >> "$L/hadd.log" 2>&1 \
+${ROOT_INSTALL:-$ROOTSYS}/bin/hadd -f "$OUT/sim/sim-merged.root" $OUT/sim/batch_0*/sim-merged.root >> "$L/hadd.log" 2>&1 \
@@ -45 +46 @@
-python reconstruction/extract_triples_chunked.py --input-format root --root-dir "$OUT/sim" --out-dir "$OUT/extracted" \
+python "$PLI_ANALYSIS/reconstruction/extract_triples_chunked.py" --input-format root --root-dir "$OUT/sim" --out-dir "$OUT/extracted" \
@@ -47 +48 @@
-python reconstruction/filter_triples_tof_pann.py --triples-detector "$OUT/extracted/triples_detector_hypoxia_te_v8.h5" \
+python "$PLI_ANALYSIS/reconstruction/filter_triples_tof_pann.py" --triples-detector "$OUT/extracted/triples_detector_hypoxia_te_v8.h5" \
```

### `simulation/runs/core_shell.sh`

```diff
--- original
+++ shipped
@@ -13,2 +13,3 @@
-ROOT=<PROJECT>
-HERE=$ROOT/reconstruction/results/2026-08-13_referee_response
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+HERE=$PLI_OUT/results/2026-08-13_referee_response
+mkdir -p "$HERE"
@@ -16,2 +17,2 @@
-cd $ROOT
-source scripts/setup/setup_environment_st.sh >/dev/null 2>&1
+cd "$PLI_SIM"
+source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
@@ -30 +31 @@
-    scripts/run_gate_parallel_musire.sh \
+    runner/run_gate_parallel.sh \
@@ -45,2 +46,2 @@
-PROD_JSON=$ROOT/materials/pli_materials.json
-EQ_JSON=$HERE/r2a_equal_tau_materials.json
+PROD_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
+EQ_JSON=$PLI_PKG/materials/variants/r2a_equal_tau_materials.json
@@ -49,2 +50,2 @@
-run_arm r2a_transport $HERE/r2_pairA.mac transport $EQ_JSON
-run_arm r2a_source    $HERE/r2_pairA.mac source    $EQ_JSON
+run_arm r2a_transport $PLI_SIM/macros/core_shell_pairA.mac transport $EQ_JSON
+run_arm r2a_source    $PLI_SIM/macros/core_shell_pairA.mac source    $EQ_JSON
@@ -53,3 +54,3 @@
-run_arm r2b_A_transport $HERE/r2_pairA.mac transport $PROD_JSON
-run_arm r2b_A_source    $HERE/r2_pairA.mac source    $PROD_JSON
-run_arm r2b_A_fixed     $HERE/r2_pairA.mac fixed     $PROD_JSON
+run_arm r2b_A_transport $PLI_SIM/macros/core_shell_pairA.mac transport $PROD_JSON
+run_arm r2b_A_source    $PLI_SIM/macros/core_shell_pairA.mac source    $PROD_JSON
+run_arm r2b_A_fixed     $PLI_SIM/macros/core_shell_pairA.mac fixed     $PROD_JSON
@@ -58,3 +59,3 @@
-run_arm r2b_B_transport $HERE/r2_pairB.mac transport $PROD_JSON
-run_arm r2b_B_source    $HERE/r2_pairB.mac source    $PROD_JSON
-run_arm r2b_B_fixed     $HERE/r2_pairB.mac fixed     $PROD_JSON
+run_arm r2b_B_transport $PLI_SIM/macros/core_shell_pairB.mac transport $PROD_JSON
+run_arm r2b_B_source    $PLI_SIM/macros/core_shell_pairB.mac source    $PROD_JSON
+run_arm r2b_B_fixed     $PLI_SIM/macros/core_shell_pairB.mac fixed     $PROD_JSON
```

### `simulation/runs/core_shell_matched_yield.sh`

```diff
--- original
+++ shipped
@@ -8,5 +8,5 @@
-ROOT=<PROJECT>
-HERE=$ROOT/reconstruction/results/2026-10-08_r2c_matched_yield
-R2=$ROOT/reconstruction/results/2026-08-13_referee_response
-OXY=$ROOT/reconstruction/results/2026-10-08_oxygen_v8/run.log
-OUT=<DATA>/r2c_runs
+source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
+HERE=$PLI_OUT/results/2026-10-08_r2c_matched_yield
+R2=$PLI_SIM/macros
+OXY=$PLI_OUT/results/2026-10-08_oxygen_v8/run.log
+OUT=$PLI_OUT/r2c_runs
@@ -15,7 +15,8 @@
-mkdir -p $OUT
-say "waiting for the oxygen GATE stage"
-while ! grep -q -E "GATE done|OXY_FAILED|OXY_DONE" $OXY; do sleep 60; done
-say "oxygen GATE stage finished; starting R2c"
-cd $ROOT
-source scripts/setup/setup_environment_st.sh >/dev/null 2>&1
-JSON=$HERE/r2c_equal_yield_materials.json
+mkdir -p $OUT "$HERE"
+# v1.1.0: the wait for the oxygen run (one GATE runner per machine) is disabled.
+#[v1.1.0 disabled] say "waiting for the oxygen GATE stage"
+#[v1.1.0 disabled] while ! grep -q -E "GATE done|OXY_FAILED|OXY_DONE" $OXY; do sleep 60; done
+#[v1.1.0 disabled] say "oxygen GATE stage finished; starting R2c"
+cd "$PLI_SIM"
+source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
+JSON=$PLI_PKG/materials/variants/r2c_equal_yield_materials.json
@@ -33 +34 @@
-    scripts/run_gate_parallel_musire.sh --macro $macro --total-primaries 10000000 --threads 16 --bounded-max-gb 5 \
+    runner/run_gate_parallel.sh --macro $macro --total-primaries 10000000 --threads 16 --bounded-max-gb 5 \
@@ -43,4 +44,4 @@
-run_arm A_transport $R2/r2_pairA.mac transport
-run_arm A_source    $R2/r2_pairA.mac source
-run_arm B_transport $R2/r2_pairB.mac transport
-run_arm B_source    $R2/r2_pairB.mac source
+run_arm A_transport $R2/core_shell_pairA.mac transport
+run_arm A_source    $R2/core_shell_pairA.mac source
+run_arm B_transport $R2/core_shell_pairB.mac transport
+run_arm B_source    $R2/core_shell_pairB.mac source
```

### `simulation/runs/closure_o4_driver.py`

```diff
--- original
+++ shipped
@@ -53,2 +53,5 @@
-REPO = Path("<PROJECT>")
-OUT = REPO / "reconstruction/results/2026-08-13_overnight"
+# v1.1.0: paths relative to the release and from the environment (simulation/README.md)
+REPO = Path(__file__).resolve().parents[1]          # simulation/, working dir of the runner
+ANALYSIS = Path(os.environ.get("PLI_ANALYSIS", REPO.parent / "analysis"))
+OUT = Path(os.environ.get("PLI_OUT") or sys.exit("set PLI_OUT (simulation/README.md)")) / "results/2026-08-13_overnight"
+OUT.mkdir(parents=True, exist_ok=True)
@@ -56 +59 @@
-RUNNER = REPO / "scripts/run_gate_parallel_musire.sh"
+RUNNER = REPO / "runner/run_gate_parallel.sh"
@@ -66,2 +69 @@
-P2B = REPO / "reconstruction/results/2026-04-29_2320_phase2b_validation_sims/macros"
-STB = REPO / "reconstruction/results/2026-06-19_steinberger_validation/macros"
+P2B = REPO / "macros"
@@ -76 +78 @@
-         macro=P2B / "sim2_i124_xad4_phantom.mac",
+         macro=P2B / "sim2_i124_xad4_phantom.mac",  # not shipped (superseded arm)
@@ -90 +92 @@
-         macro=OUT / "macros/mercolli_te_poisson.mac",
+         macro=REPO / "macros/mercolli_xad4_i124.mac",
@@ -95 +97 @@
-         macro=P2B / "sim1_sc44_validation_phantom.mac",
+         macro=REPO / "macros/closure_sim1_sc44.mac",
@@ -101 +103 @@
-         macro=STB / "steinberger_disk_phantom.mac",
+         macro=REPO / "macros/steinberger_disks.mac",
@@ -114 +116 @@
-         macro=REPO / "reconstruction/results/2026-08-12_steinberger_allprompts/water_7MBq.mac",
+         macro=REPO / "macros/water_cylinder_20mm.mac",
@@ -119 +121 @@
-         macro=OUT / "macros/brain_tumor_te_poisson_1B_truetime.mac",
+         macro=REPO / "macros/brain_tumor.mac",
@@ -146 +148 @@
-    st = os.statvfs("/home")
+    st = os.statvfs(OUT)
@@ -299 +301,10 @@
-RE_ABSPATH = __import__("re").compile(r"/home/[^\s\"']+")
+# v1.1.0: macros name their inputs through GATE aliases ({PLI_SIM}, ...) taken from the
+# environment, so expand those and check every database, image and table path.
+RE_ABSPATH = __import__("re").compile(
+    r"^/\S*(?:setMaterialDatabase|setImage|setRangeToMaterialFile|readTable|readFile)\s+(\S+)",
+    __import__("re").M)
+RE_ALIAS = __import__("re").compile(r"\{(\w+)\}")
+
+
+def _expand(txt: str) -> str:
+    return RE_ALIAS.sub(lambda m: os.environ.get(m.group(1), m.group(0)), txt)
@@ -323 +334 @@
-        txt = mac.read_text()
+        txt = _expand(mac.read_text())
@@ -401 +412 @@
-            subprocess.call([sys.executable, f"reconstruction/{ph['extractor']}",
+            subprocess.call([sys.executable, str(ANALYSIS / "reconstruction" / ph["extractor"]),
@@ -423,2 +434,2 @@
-                     "reconstruction/results/2026-06-10_brain_tumor_1B_truetime/"
-                     "scripts/gate_root_to_singles_h5.py",
+                     str(ANALYSIS / "results/2026-06-10_brain_tumor_1B_truetime/"
+                         "scripts/gate_root_to_singles_h5.py"),
@@ -517 +528 @@
-        rc = subprocess.call([sys.executable, "reconstruction/filter_triples_tof_pann.py",
+        rc = subprocess.call([sys.executable, str(ANALYSIS / "reconstruction/filter_triples_tof_pann.py"),
```

### `simulation/runs/vipman_batches.py`

```diff
--- original
+++ shipped
@@ -56 +56,3 @@
-REPO = Path("<PROJECT>").resolve()
+# v1.1.0: paths relative to the release and from the environment (simulation/README.md)
+REPO = Path(__file__).resolve().parents[1]          # simulation/
+ANALYSIS = Path(os.environ.get("PLI_ANALYSIS", REPO.parent / "analysis"))
@@ -66,4 +68,4 @@
-MACRO = REPO / "macros/vipman_sc44_1400M_poisson.mac"
-RUNNER = REPO / "scripts/run_gate_parallel_musire.sh"
-SENSITIVITY = REPO / "reconstruction/results/sensitivity_4mm_fullfidelity_norm.npy"
-MUMAP = REPO / "phantom/vipman-atlas-lungs-sim-mumap-harmonized.mhd"
+MACRO = REPO / "macros/vipman_slab.mac"
+RUNNER = REPO / "runner/run_gate_parallel.sh"
+SENSITIVITY = Path(os.environ.get("PLI_SENSITIVITY_NPY", ANALYSIS / "results/sensitivity_4mm_fullfidelity_norm.npy"))
+MUMAP = Path(os.environ.get("PLI_VIPMAN_DIR", "PLI_VIPMAN_DIR-not-set")) / "vipman-atlas-lungs-sim-mumap-harmonized.mhd"
@@ -75,2 +77,2 @@
-VENV_PY = REPO / "venv/bin/python"
-PY_PATH = "/usr/lib/python3.14/site-packages"
+VENV_PY = Path(os.environ.get("PLI_PYTHON", sys.executable))
+PY_PATH = os.environ.get("PLI_EXTRA_PYTHONPATH", "")
@@ -169 +171 @@
-            "reconstruction/utils/lu176_bg_generator.py",
+            str(REPO / "lu176/lu176_bg_generator.py"),
@@ -190 +192 @@
-    sys.path.insert(0, str(REPO / "reconstruction/utils"))
+    sys.path.insert(0, str(REPO / "lu176"))
@@ -207 +209 @@
-            "reconstruction/extract_triples_chunked.py",
+            str(ANALYSIS / "reconstruction/extract_triples_chunked.py"),
@@ -383 +385 @@
-    sys.path.insert(0, str(REPO / "reconstruction"))
+    sys.path.insert(0, str(ANALYSIS / "reconstruction"))
@@ -409 +411 @@
-            "reconstruction/simple_recon.py", "single",
+            str(ANALYSIS / "reconstruction/simple_recon.py"), "single",
@@ -419,2 +421,2 @@
-    env["PARALLELPROJ_C_LIB"] = "/usr/local/lib/libparallelproj_c.so"
-    env["PARALLELPROJ_CUDA_LIB"] = "/usr/local/lib/libparallelproj_cuda.so"
+    env["PARALLELPROJ_C_LIB"] = os.environ.get("PARALLELPROJ_C_LIB", "/usr/local/lib/libparallelproj_c.so")
+    env["PARALLELPROJ_CUDA_LIB"] = os.environ.get("PARALLELPROJ_CUDA_LIB", "/usr/local/lib/libparallelproj_cuda.so")
```

### `simulation/runs/seed_repeatability.py`

```diff
--- original
+++ shipped
@@ -27,5 +27,7 @@
-ROOT = Path("<PROJECT>")
-HERE = ROOT / "reconstruction/results/2026-08-12_seed_repeatability"
-MACRO = ROOT / ("reconstruction/results/2026-04-29_2320_phase2b_validation_sims/"
-                "macros/sim1_sc44_validation_phantom.mac")
-RUNNER = ROOT / "scripts/run_gate_parallel_musire.sh"
+# v1.1.0: paths relative to the release and from the environment (simulation/README.md)
+import os
+ROOT = Path(__file__).resolve().parents[1]          # simulation/
+ANALYSIS = Path(os.environ.get("PLI_ANALYSIS", ROOT.parent / "analysis"))
+HERE = Path(os.environ.get("PLI_OUT") or sys.exit("set PLI_OUT (simulation/README.md)")) / "results/2026-08-12_seed_repeatability"
+MACRO = ROOT / "macros/closure_sim1_sc44.mac"
+RUNNER = ROOT / "runner/run_gate_parallel.sh"
@@ -54 +56 @@
-    dst = ROOT / "scripts" / f"runner_seed{base_seed}.sh"
+    dst = ROOT / "runner" / f"runner_seed{base_seed}.sh"
@@ -136 +138 @@
-            sh([sys.executable, ROOT / "reconstruction/extract_triples_chunked.py",
+            sh([sys.executable, ANALYSIS / "reconstruction/extract_triples_chunked.py",
@@ -149 +151 @@
-        sh([sys.executable, ROOT / "reconstruction/filter_triples_tof_pann.py",
+        sh([sys.executable, ANALYSIS / "reconstruction/filter_triples_tof_pann.py",
```

### `simulation/lu176/lu176_bg_generator.py`

No edit apart from the header lines.

### `simulation/lu176/lso_background_inject.py`

No edit apart from the header lines.

### `materials/variants/make_r2c_materials.py`

```diff
--- original
+++ shipped
@@ -7 +7 @@
-REPO = Path('<PROJECT>')
+REPO = Path(__file__).resolve().parents[2]          # v1.1.0: release root
```

### `integration/variants/build_gate_inflight.sh`

```diff
--- original
+++ shipped
@@ -2,5 +2,6 @@
-R=<PROJECT>
-source $R/scripts/setup/setup_environment_st.sh >/dev/null 2>&1
-SRC=$R/software/source/gate-9.4.1-inflight
-BLD=$R/software/build/gate-inflight
-INS=$R/software/install/gate-inflight
+# v1.1.0: locations from the environment; preparing SRC is described in integration/PATCH_NOTES.md
+PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
+source "${PLI_ENV_SCRIPT:-$PKG/simulation/runner/env_example.sh}" >/dev/null 2>&1
+SRC=${GATE_INFLIGHT_SRC:?set GATE_INFLIGHT_SRC to a GATE 9.4.1 tree prepared with install_into_gate.sh and gate-inflight-stock.patch}
+BLD=${GATE_INFLIGHT_BUILD:-$PWD/build-gate-inflight}
+INS=${GATE_INFLIGHT_INSTALL:-$PWD/install-gate-inflight}
@@ -11 +12 @@
-  -DGeant4_DIR=$R/software/install/geant4-st/lib/cmake/Geant4 -DROOT_DIR=$R/software/install/root/cmake \
+  -DGeant4_DIR=$GEANT4_INSTALL/lib/cmake/Geant4 -DROOT_DIR=$ROOT_INSTALL/cmake \
```

### `integration/variants/pli_inflight_test/pli_inflight_test.cc`

No edit apart from the header lines.

### `integration/variants/pli_inflight_test/CMakeLists.txt`

```diff
--- original
+++ shipped
@@ -9 +9,10 @@
-set(PLI <PROJECT>/custom_physics)
+# v1.1.0: relative to this folder. The variant source is the shipped PLIPhysicsList.cc
+# with ../gate-inflight-stock.patch applied, made at configure time in the build tree.
+set(PLI ${CMAKE_CURRENT_SOURCE_DIR}/../../../custom_physics)
+set(PLI_VARIANT_SRC ${CMAKE_CURRENT_BINARY_DIR}/PLIPhysicsList_inflight.cc)
+configure_file(${PLI}/src/PLIPhysicsList.cc ${PLI_VARIANT_SRC} COPYONLY)
+execute_process(COMMAND patch --silent ${PLI_VARIANT_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/../gate-inflight-stock.patch
+                RESULT_VARIABLE PLI_PATCH_RC)
+if(NOT PLI_PATCH_RC EQUAL 0)
+  message(FATAL_ERROR "gate-inflight-stock.patch does not apply to ${PLI}/src/PLIPhysicsList.cc")
+endif()
@@ -13 +22 @@
-  <PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/PLIPhysicsList_inflight.cc
+  ${PLI_VARIANT_SRC}
@@ -16 +25 @@
-target_include_directories(pli_inflight_test PRIVATE ${PLI}/include)
+target_include_directories(pli_inflight_test PRIVATE ${PLI}/include ${PLI}/third_party)
```

### `integration/variants/pli_inflight_test/run_physics_level.sh`

```diff
--- original
+++ shipped
@@ -2,3 +2,7 @@
-source <PROJECT>/software/install/geant4-st/bin/geant4.sh
-export PLI_MATERIALS_JSON=<PROJECT>/materials/pli_materials.json
-cd <PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight
+# v1.1.0: expects the CMake build of this test in ./build (cmake -S . -B build);
+# output goes to $PLI_OUT/feedback3_inflight/physics.
+HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
+source "${GEANT4_INSTALL:?set GEANT4_INSTALL}/bin/geant4.sh"
+export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$HERE/../../../materials/pli_materials.json}
+PHYS=${PLI_OUT:?set PLI_OUT}/feedback3_inflight/physics; mkdir -p "$PHYS"
+cd "$HERE"
@@ -8 +12 @@
-    ( nice -n 10 ./build/pli_inflight_test $m $iso 1000000 1.0 $seed <DATA>/feedback3_inflight/physics/${iso}_${m}.csv > <DATA>/feedback3_inflight/physics/${iso}_${m}.log 2>&1; echo "$(date +%T) $iso $m rc=$?" >> <PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/physics_level.log ) &
+    ( nice -n 10 ./build/pli_inflight_test $m $iso 1000000 1.0 $seed $PHYS/${iso}_${m}.csv > $PHYS/${iso}_${m}.log 2>&1; echo "$(date +%T) $iso $m rc=$?" >> $PHYS/physics_level.log ) &
@@ -12 +16 @@
-echo "$(date +%T) ALL_DONE" >> <PROJECT>/reconstruction/results/2026-09-30_feedback3/inflight/physics_level.log
+echo "$(date +%T) ALL_DONE" >> $PHYS/physics_level.log
```

## Comment wording (release v1.1.0)

In the lines below, words of the project's internal working notes were replaced by neutral wording ("decision", "the authors'", "planned", "specified"). Only comments, docstrings and one free-text field changed. No code line, macro command or value changed.

| file | line in the release file | reads now |
|---|---|---|
| `runs/c2_eight_tube_388ps.sh` | 8 | # CONVERTER WINDOW (decision, 2026-08-24): keep the FULL digitizer band with no |
| `runs/seed_repeatability.py` | 7 | from the planned 4 because 4 x (32 min GATE + convert/extract) |
| `runs/seed_repeatability.py` | 8 | exceeds the ~2 h budget; stated here rather than silently |
| `sources/iodine124_tubes_poisson.mac` | 29 | # Activity per tube (from Mercolli 2025a Table): |
