# PLI-GATE — positronium lifetime physics for GATE

Release v1.1.0, the code release of the paper "Material-dependent positronium
lifetimes in GATE and prompt-selection bias" (D. Dhar and J. Peter, submitted
to Physics in Medicine and Biology).

## What this is

A physics extension for [GATE](https://opengatecollaboration.org/) 9.4.1 that
simulates positronium lifetime imaging, together with the simulation inputs,
analysis code, figure scripts and registration records of the paper.

A positron stopping in matter either annihilates directly or first forms
positronium. Ortho-positronium survives a few nanoseconds before being quenched
by the surrounding electrons, and how long it survives depends on the size of
the free volume it occupies and on the oxygen dissolved there. That delay is
measurable, which is what makes it an imaging signal. This package gives GATE
the physics to produce it: per-material positronium formation, Tao-Eldrup
lifetimes with oxygen quenching, and the Ore-Powell three-photon channel, all
driven by a materials database rather than by constants compiled into the code.

## What the release contains

| folder or file | contents | reproduces |
|---|---|---|
| `custom_physics/` | the positronium physics (five classes compiled into GATE, plus an optional standalone library and two C++ tests) | the physics of the production binary (`ENVIRONMENT_MANIFEST.json`) |
| `integration/` | `install_into_gate.sh`, the 11-file GATE patch, `PATCH_NOTES.md`, and `variants/` with the in-flight annihilation variant | the GATE build of every run, and the in-flight sensitivity arm |
| `materials/` | the positronium materials database: version 8 (current), version 7 in `previous/`, the run-specific variants in `variants/`, and `PROVENANCE.md` | the lifetime and yield of every material in every run |
| `simulation/` | GATE macros, the parallel runner and the retained run scripts of the simulations, with the run table (binary, macro, batches, seeds, switches, database) | the simulations of the paper, given the scanner geometry (below). `simulation/RUNS.json` marks the few runs whose macro or launch script was not retained |
| `analysis/` | analysis code and intermediate data | the numbers of the paper |
| `figures/` | figure scripts and their source data | the figures of the paper |
| `registration/` | pre-run registration records | the registered predictions and their dates |
| `examples/quickstart/` | a four-second run with no scanner | that the physics returns the lifetime it is given |
| `tests/` | the test suite | checks of the physics and of the shipped values |
| `docs/` | `PHYSICS.md`, `VALIDATION.md`, `BOOTSTRAP.md` | |
| `ENVIRONMENT_MANIFEST.json` | build, binaries, source and database hashes, runtime switches | |
| `CHANGELOG.md`, `CITATION.cff`, `.zenodo.json`, `LICENSE`, `COPYING` | | |

**Simulation (`simulation/`).** Each production macro and run script is
shipped with its workstation paths replaced by variables, and
`simulation/EDITS.md` lists every change. `simulation/README.md` maps each
paper result to its phantom, binary (the production binary e37be8b4, the
archived build 3522d567, the in-flight variant cc552a68, or an earlier build
that was not retained), macro, run script, batch layout, seeds, runtime
switches and materials database. The scanner geometry and two further inputs
are not included (see "Third-party inputs not included").

**Analysis (`analysis/`).** The scripts that turn the simulated singles into
the numbers of the paper, with the intermediate data they read. The run
scripts of `simulation/` call the singles converter and the triple extractors
of this part. `analysis/README.md` maps each result to its script and data.

**Figures (`figures/`).** One script per figure of the paper and of the
supplement, with the source data each one reads. `figures/README.md` lists
which script makes which figure and how to run them.

**Registration (`registration/`).** The records in which predictions were
registered before the runs that tested them, including those of the C0, C2,
in-flight and matched-yield runs of `simulation/`. `registration/README.md`
explains the records and how their dates are established.

## Requirements

| component | version | note |
|---|---|---|
| GATE | 9.4.1 | upstream commit `2d03a89bfec37c21ec69d196b4876d5e3bf5ed83` |
| Geant4 | 11.3.2 | see the note below |
| ROOT | 6.36.04 | needed by GATE, not by this physics |
| CMake | ≥ 3.16 | built with 4.2.1 |
| C++ | C++17 | built with GCC 15.2.1. The GATE patch includes the fixes GCC 15 needs |
| Python | ≥ 3.9 | standard library only, for the tests and the quickstart |

nlohmann/json is **bundled** (`custom_physics/third_party/`), so no system
package is needed.

The simulation runner also needs GNU `parallel`, ROOT's `hadd`, Python with
`uproot`, and a systemd user session (`simulation/README.md`).

If you do not already have Geant4 and ROOT, `docs/BOOTSTRAP.md` gives the
complete from-source path with the exact options used to validate this
package: about 21 minutes and 5.4 GB on a 32-core machine, nothing installed
system-wide.

**On the Geant4 version.** GATE 9.4.1 is validated against Geant4 11.2.1 and
prints a warning for any other version. That is a warning, not a requirement:
`CMakeLists.txt:64` calls `FIND_PACKAGE(Geant4 REQUIRED)` with no version pin,
and the version check at lines 101-102 only prints a message. This package
was built and validated against **Geant4 11.3.2 as released**. Geant4 was not
patched. GATE was: the patch in `integration/` adds the physics-list hook and
the fixes that let GATE 9.4.1 build with GCC 15 and CMake 4 and run under
Geant4 11.3.2 without visualisation (`integration/PATCH_NOTES.md`).

## Install

Five commands. `$GATE_SRC` is a GATE 9.4.1 source tree, `$PKG` is this package.

```bash
# 1. put the physics into the GATE tree and apply the patch
$PKG/integration/install_into_gate.sh $GATE_SRC

# 2. configure GATE against your Geant4 and ROOT
cmake -B build -S $GATE_SRC \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/gate-pli \
      -DGeant4_DIR=$GEANT4_INSTALL/lib/cmake/Geant4 \
      -DROOT_DIR=$ROOTSYS/cmake

# 3. build and install
cmake --build build -j $(nproc)
cmake --install build

# 4. put it on PATH
export PATH=$HOME/gate-pli/bin:$PATH

# 5. tell the physics where the materials database is
export PLI_MATERIALS_JSON=$PKG/materials/pli_materials.json
```

`install_into_gate.sh` copies five `.cc` and five `.hh` files into
`source/physics/`, which GATE compiles automatically because it globs that
directory, and applies `integration/gate-9.4.1-pli.patch`. The patch is 11
files: one hook that lets a macro select the physics list by name, and ten
portability fixes for GCC 15, CMake 4 and Geant4 11.3.2.

The production configuration (the GATE options of `ENVIRONMENT_MANIFEST.json`)
is in `integration/variants/build_gate_inflight.sh`. Measured for v1.1.0 on
the 16-core workstation, from a clean GATE tree at the base commit: configure,
build with 12 jobs and install in 95 s, build tree 86 MB, installed GATE 12 MB.

A GATE built from this release does not have the hash of the production
binary (e37be8b4), because the binary embeds the install paths of its machine
and one unused class that is not shipped. `ENVIRONMENT_MANIFEST.json` shows
the equivalence at the level of the object files instead.

## Quickstart

Three commands to a number.

```bash
export PLI_MATERIALS_JSON=$PKG/materials/pli_materials.json
cd $PKG/examples/quickstart
./run_quickstart.sh
```

It simulates 200,000 positrons stopping in a small water tube whose o-Ps
lifetime is set to 1.9215 ns, then recovers that lifetime from the delays the
physics emitted. Expect:

```
  tau_oPs = 1.9072 +/- 0.0107 ns
expected            : 1.9215 ns (configured value)
difference          : -0.0143 ns   tolerance +/-0.0500
RESULT: PASS
```

**Measured: 3 to 4 s wall**, negligible disk. There is deliberately no scanner in
this example. The check is on the physics, not on a detector readout, so the
quickstart needs none of the inputs that are not included.

## Tests

```bash
$PKG/tests/run_tests.sh
```

Nine entries, each printing PASS/FAIL and one line of what it verifies. Read
them as **seven checks and two recorded results**: `three_gamma_branching` and
`null_transport_vs_source` are marked `[recorded result]` because they validate
numbers measured on the reference run rather than executing package code on
yours. Two of the seven need things built first: `ore_powell_spectrum` needs
`custom_physics/build/test_orepowell`, and `fixed_seed_regression` needs the
quickstart to have been run. Both SKIP rather than fail if absent.

Three checks (`tao_eldrup_numeric`, `o2_quenching_conversion`,
`pedestal_analytic_identity`) import `tao_eldrup.py` and `pedestal.py` from the
analysis part. The suite looks for them in `analysis/pli_analysis/`,
`analysis/` and `analysis/legacy_v1.0.1/`, or in the folder named by
`PLI_ANALYSIS_PY`. `o2_quenching_conversion` requires the oxygen constant of
that model, of `tests/reference_values.json` and of the shipped database to be
the same (2.5704e-5 ns^-1 mmHg^-1, database version 8).

**Measured: 1 s** once the quickstart has run and `custom_physics` is built.

## Running the paper's simulations

`simulation/README.md` explains the variables, the runner and the run table.
In short:

```bash
export PLI_OUT=/scratch/pli-out                      # output root
export PLI_SCANNER_MACRO=/path/to/quadra_detector.mac # from the authors
export PLI_STACK=$HOME/pli-stack/install             # Geant4, ROOT, GATE (docs/BOOTSTRAP.md)
bash $PKG/simulation/runs/c0_eight_tube.sh
```

The production runs used 12 to 30 parallel shards on a 16-core, 78 GB
workstation, and the runner's memory caps and thread counts are set for it.

## Third-party inputs not included

| input | used by | why not included |
|---|---|---|
| Siemens Biograph Vision Quadra geometry macro (`quadra_detector.mac`) | every scanner simulation | derived from J. Peter's Musiré simulation package, licence for redistribution not cleared |
| GATE material database `gate-materials.db` | NEMA IQ, VIP-Man and stock-versus-custom runs | Musiré copy of a GATE material file, licence not cleared |
| VIP-Man voxel phantom (`vipman-atlas-*.mhd/.raw`) and its material and activity tables (`.dat`) | VIP-Man and stock-versus-custom runs | VIP-Man (Xu et al. 2000) redistribution terms unclear, the tables are Musiré inputs |

These inputs are available from the authors (dhritman@dhar.dev). Every macro
names them through a GATE alias (`{PLI_SCANNER_MACRO}`,
`{PLI_GATE_MATERIALS_DB}`, `{PLI_VIPMAN_DIR}`) read from the environment, and
`ENVIRONMENT_MANIFEST.json` records the sha256 of the files used, so copies
can be checked. The quickstart and the tests need none of them.

Two analysis-side scripts restate scanner parameters as numbers: the singles
generator `simulation/lu176/lu176_bg_generator.py` (crystal layout, for the
176Lu background) and `analysis/results/2026-08-23_acceptance/sigma_geo.py`
(ring radii, for the a priori acceptance). These are the scanner dimensions the
paper reports (Sec. 2.1 and Supplementary Note S20), not the GATE macro. The
scanner sensitivity image used by the VIP-Man reconstruction and the VIP-Man
tissue-label map at 4 mm are derived from the excluded inputs and are not
included either.

## How to add a material

Add an object to the `materials` array in `materials/pli_materials.json`. The
name must match the GATE material name in your macro exactly, or the run
aborts.

```json
{
  "name": "MyTissue",
  "aliases": [],
  "R_nm": 0.2800,
  "pO2_mmHg": 40.0,
  "I3": 0.27,
  "I3_convention": "ops_intensity",
  "tau_TE_ns": 1.9215,
  "tau_source": "Author Year, journal reference"
}
```

| field | required | meaning |
|---|---|---|
| `name` | yes | must equal the GATE material name |
| `tau_TE_ns` | yes | **effective** o-Ps lifetime `tau_eff`, ns. Must be > 0. Despite the name this is NOT the pure Tao-Eldrup lifetime, see below |
| `I3` | yes | positronium intensity. Read `I3_convention` before setting it |
| `I3_convention` | **yes** | `ops_intensity` or `total_formation` |
| `R_nm`, `pO2_mmHg` | no | recorded for provenance, not used at run time |
| `aliases` | no | alternative GATE material names mapping to this entry |

**`tau_TE_ns` holds `tau_eff`, not `tau_TE`.** The name is historical. It is
kept so that the shipped databases and the binaries that produced the paper's
numbers stay in correspondence. The stored value is the *effective* lifetime,
which includes the three-photon and oxygen-quenching terms that the Tao-Eldrup
model itself omits:

    tau_eff = 1 / ( lambda_po(R) + kappa * pO2 + lambda_3g )
    tau_TE  = 1 / lambda_po(R)                     <- pure Tao-Eldrup

For `water_shibuya2020_n2` the difference is 26.4 ps:

    tau_TE_ns   (stored)   1.92150 ns
    pure tau_TE (derived)  1.94786 ns

So a check of this column against a published Tao-Eldrup curve will disagree by
19 to 83 ps across the database's range of R, and that disagreement is expected.
To recover the pure value, invert the rate equation and subtract the terms TE
does not contain, for example with `lambda_po` of the analysis part's
`tao_eldrup.py`: `tau_TE = 1 / lambda_po(entry["R_nm"])`, 1.94786 ns here.

**`I3_convention` has no default and omitting it is a fatal error.** PALS
papers usually quote the o-Ps intensity directly (`ops_intensity`), and some
quote the total formation probability (`total_formation`). The two differ by
4/3, so guessing would introduce a 33 % error in the o-Ps yield silently. The
loader refuses rather than choose.

Likewise, **a material absent from the database is a fatal error**, not a
default. If you genuinely want the old permissive behaviour, set
`PLI_ALLOW_MATERIAL_FALLBACK=1` explicitly.

## Runtime switches

All are environment variables, read once at first use. The defaults are the
production behaviour. `ENVIRONMENT_MANIFEST.json` lists the file that reads
each one.

| variable | default | effect |
|---|---|---|
| `PLI_MATERIALS_JSON` | search: `materials/`, `../materials/`, `../../materials/`, `../../../materials/` | path of the materials database |
| `PLI_DIRECT_TAU_PS` | `0` | direct-annihilation lifetime in ps. `0` = instantaneous |
| `PLI_ALLOW_MATERIAL_FALLBACK` | `0` | non-zero = warn instead of aborting on an unknown material (tau 2.0 ns, formation probability 0.35) |
| `PLI_DISABLE_3GAMMA` | unset | `1` = two-photon-only o-Ps decay |
| `PLI_TAU_3GAMMA_NS` | `142.05` | vacuum three-photon lifetime, ns |
| `PLI_LIFETIME_MODE` | `transport` | `transport` / `source` / `fixed` |
| `PLI_FIXED_TAU_NS` | `2.0` | lifetime in ns when the mode is `fixed` |
| `PLI_LIFETIME_LOG` | unset | `1` = one `PLI_LIFETIME:` record per annihilation on stdout. The simulation runner sets it to `1` unless told otherwise |
| `PLI_LOG_XMAT` | unset | `1` = one `PLI_R2:` record per annihilation with emission and annihilation material |
| `PLI_INFLIGHT_STOCK` | unset | in-flight variant only (`integration/variants/`): `1` = keep Geant4's annihilation for positrons that annihilate in flight |

## Known limitations

Each of these is expanded in `docs/PHYSICS.md` and `docs/VALIDATION.md`.

* **Direct annihilation is instantaneous by default.** With a realistic finite
  value the recovered single-exponential lifetime changes substantially. This
  is the most consequential limitation here.
* **No positronium transport.** Ps annihilates where it forms.
* **Two-photon emission is collinear.** There is no intrinsic acollinearity.
* **Free volumes are spherical**, and for tissues `R` is back-derived from the
  published lifetime, so those materials cannot test the model.
* **No temperature dependence** in the lifetime model.
* **A recovered-lifetime offset of 9 to 16 ps** exists in the full analysis
  chain and is not explained. The physics sampler itself is exact to
  sub-picosecond.
* **A 2-sigma trend of −11.7 ± 6.0 ps/ns** in the calibration closure across
  lifetime is not explained. It is treated as a systematic uncertainty rather
  than statistical noise. The fitted floor is 3.9 ps, carrying 23 % of the
  residual variance. Releases up to v1.0.0 quoted −19.11 ± 6.27 ps/ns from a
  superseded analysis pipeline (`CHANGELOG.md`).
* **The vacuum o-Ps lifetime has two values in this package.**
  `custom_physics/src/PLIOrePowell.cc` uses 142.05 ns (Vallery et al., PRL 90,
  203402) for the three-photon branching. The Tao-Eldrup model of the analysis
  part, `custom_physics/include/PLIMaterialLoader.hh` and the `lambda_3g_ns`
  constant in the databases use 142.0 ns. The effect is small and quantified:
  0.0092 ps on `tau_eff` for water, and 3.5e-4 relative on `P(3 gamma)`. The
  values are deliberately **not** reconciled, because the paper's numbers come
  from binaries carrying both, and changing either would break that
  correspondence for an effect four orders of magnitude below the experimental
  noise floor. See `custom_physics/include/PLIConstants.hh`.
* **The lifetime model is standard Tao-Eldrup and is valid only for
  R < 1 nm.** Every database entry sits inside that range (largest R = 0.393 nm).
  Mesoporous materials require the extended model and are outside this
  implementation's scope.
* **Not every run can be rebuilt bit for bit.** Runs made before 2026-08-14
  14:54 used builds that were not retained, and the archived build 3522d567
  has no source snapshot (`simulation/README.md` lists which runs). The
  materials database 2a26b68e of the closure reruns and the earlier ones are
  not retained either. Two shipped variant databases
  (`materials/variants/lookup_baseline.json`,
  `materials/variants/r2a_equal_tau_materials.json`) predate the required
  `I3_convention` field, so this release cannot load them. They are kept as
  the record of what those runs were given.
* **The scanner geometry is not included**, so no scanner simulation runs
  without the inputs from the authors.

## A pitfall when reproducing rate-dependent results

If you use this package to study anything that depends on the **rate** of
decays (accidental coincidences, prompt competition, pile-up, randoms
subtraction), read this first. It cost us three runs.

GATE offers two ways to say how much to simulate:

```
# primaries-driven — DO NOT USE for rate-dependent work
/gate/application/setTotalNumberOfPrimaries 250000000

# time-driven — use this
/gate/source/mysource/setActivity 1.0e6 becquerel
/gate/application/setTimeStop 2.2 s
```

With `setTotalNumberOfPrimaries` GATE spaces the decays **deterministically**:
it divides the run duration by the primary count and emits one decay per slot.
The inter-decay time then has essentially zero variance instead of being
exponentially distributed. In our Sim-1 configuration this produced a hard
70.6 ns floor between consecutive decays, and the coefficient of variation of
the inter-decay interval collapsed from the Poisson value of 1.0 to about 0.

Every rate-dependent quantity is then wrong, and wrong in a way that looks
plausible. With no clustering there are almost no accidental coincidences, so
an artefact that scales with prompt competition simply vanishes. We measured a
null result twice before finding the cause.

**Diagnostic.** Compute the coefficient of variation of the inter-decay
intervals. Poisson emission gives CV = 1. Ours read 1.005 once the run was
time-driven, and effectively 0 before.

Note that matching the total counts between the two modes is not enough. The
counts can agree exactly while the time structure differs completely, which is
what makes this failure quiet. The eight-tube runs of the paper (C0, C2 and the
in-flight arm) are time-driven. The calibration closure (Sim-1) and the seed
study ran in primaries mode (`simulation/README.md`).

## Citation

See `CITATION.cff`. The paper is submitted: D. Dhar and J. Peter,
"Material-dependent positronium lifetimes in GATE and prompt-selection bias",
Physics in Medicine and Biology (submitted). The paper reference will be
added on publication.

## Licence

Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ).

The package is free software under the GNU Lesser General Public License
version 3 only (SPDX `LGPL-3.0-only`). `LICENSE` holds the LGPL-3.0 text, and
`COPYING` the GPL-3.0 text that the LGPL incorporates. The C++ sources, tests and
top-level scripts carry the SPDX identifier and the copyright line. The analysis,
figure and registration files are shipped as they ran and are covered by the same
licence through this file.

* The bundled nlohmann/json header (`custom_physics/third_party/nlohmann/json.hpp`)
  is MIT-licensed, copyright Niels Lohmann. Its notice is retained in the file
  and must be kept.
* `integration/gate-9.4.1-pli.patch` modifies GATE, which is distributed under
  LGPL-3.0, and is provided under GATE's licence.
* The third-party inputs listed above are not part of the release and are not
  covered by this licence.

## Contact

Questions about the package or the physics: dhritman@dhar.dev
