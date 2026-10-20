# Changelog

## v1.1.0 — 2026-10-20

The code release that accompanies the paper "Material-dependent positronium
lifetimes in GATE and prompt-selection bias" (D. Dhar and J. Peter, submitted
to Physics in Medicine and Biology). The physics code is unchanged.

**Licence.** The package is released under LGPL-3.0-only (`LICENSE`, with the
GPL-3.0 text it incorporates in `COPYING`). v1.0.x granted no licence. The
project's own physics, integration, materials, simulation, example and test
sources begin with an SPDX identifier and the copyright line.

**Physics sources.** Same text as v1.0.1 apart from the two header lines
(`tail -n +3 <file>` gives the v1.0.1 file). `ENVIRONMENT_MANIFEST.json` now
documents that this text is code-equivalent to the production binary
e37be8b4: built from this release, `PLIMaterialLoader.cc` and
`PLIPhysicsList.cc` give object files byte-identical to the production ones,
and the other three differ only in the lengths of reworded message strings.

**Materials database version 8 is the current database.**
`materials/pli_materials.json` and `materials/pli_materials_vipman.json` are the
files used by the paper's version-8 runs (sha256 23f607ea and 1e59f1cf as used),
and `materials/previous/` holds version 7 (f6d5a53c and 09e75eee as used). The
release copies of these and of the variant databases differ from the files used
only in one reworded phrase of a `comment` field that the loader does not read.
`ENVIRONMENT_MANIFEST.json` and `materials/PROVENANCE.md` give both hashes. At run time
the two versions differ only in the five `TumorTissue_pO2_*` entries (at most
0.11 ps), and version 8 changes the oxygen constant to 2.5704e-5 ns^-1 mmHg^-1
and adds I3 provenance classes. The v1.0.x copies (302f2645, a44c09b6), version 7
with reworded provenance text, are removed: version 8 carries every field they
carried and corrects attributions they still gave. `materials/variants/` holds
the lookup-arm, null-control and matched-yield databases of three paper runs.
`materials/make_provenance.py` writes the file hashes and the version-8 changes
into `materials/PROVENANCE.md`.

**New: `simulation/`.** The GATE macros, the parallel runner and the retained run
scripts of the simulations in the paper (`simulation/RUNS.json` marks the runs whose
macro or launch script was not retained), with a run table (binary, macro, script,
batches and shards, seeds, switches, database) in `simulation/README.md` and
`simulation/RUNS.json`. Paths are variables, and every edit to the production
files is listed in `simulation/EDITS.md`. The scanner geometry, the GATE
material database and the VIP-Man data it was used with are not included
(licence not cleared) and are available from the authors.

**New: `integration/variants/`.** The in-flight annihilation variant used for
one sensitivity arm, as a patch against the shipped `PLIPhysicsList.cc`, with
its build script and its physics-level test.

**New: `analysis/`, `figures/`, `registration/`.** The analysis code, the
figure scripts with their source data, and the pre-run registration records of
the paper. Each has its own README.

**Manifest rewritten.** `ENVIRONMENT_MANIFEST.json` records the production
build, the two other binaries and which runs used them, the GATE patch (GATE
was patched for build compatibility, which v1.0.x stated in two contradictory
ways), the shipped, v1.0.1 and as-compiled source hashes with the object-file
evidence, the database hashes, and all nine runtime switches with their
defaults (v1.0.x listed five).

**Tests.** `tests/test_pli.py` finds the Python models in the layouts the
analysis part has used, and its oxygen check now requires the analysis model,
`tests/reference_values.json` and the shipped database to agree on the oxygen
constant (version 8, 2.5704e-5).

**Documentation.** README rewritten for the release contents. `docs/PHYSICS.md`
gives the version-8 oxygen constant and the default of `PLI_FIXED_TAU_NS`
(2.0 ns, previously shown as none).

## v1.0.1 — 2026-08-27

Documentation only. No code, test, macro or database file changed.

Three files differed from v1.0.0: `README.md`, `docs/VALIDATION.md` and
`CITATION.cff`. The calibration-closure trend had been quoted from a superseded
pre-correction analysis pipeline. Every number in that section was replaced by
the post-correction rerun that the paper quotes:

    was   -19.11 +/- 6.27 ps/ns  (3 sigma),   -0.745 +/- 0.273 %/ns
          sign change near 2.01 ns, mean signed residual -2.54 +/- 2.28 ps
    now   -11.70 +/-  5.97 ps/ns (2.0 sigma), -0.448 +/- 0.259 %/ns
          sign change near 1.90 ns, mean signed residual -2.85 +/- 2.01 ps
          plus a fitted systematic floor of 3.9 ps, 23 % of the residual variance

The downgrade from a three-sigma to a two-sigma trend came from our own
pipeline correction. `CITATION.cff` gained `version: "1.0.1"`.

Tarball `pli-gate-v1.0.1.tar.gz`, sha256
`ee2d6da398cc85647f84f7371b2a948bcf775a4571c0b3b085dd4abc92a0f3df`, 43 files.

## v1.0.0 — 2026-08

First release, 43 files: the positronium physics for GATE 9.4.1
(`custom_physics/`), the install script and the 11-file GATE patch
(`integration/`), materials database version 7 with reworded provenance text
(`materials/`), a quickstart that recovers a configured lifetime in seconds
(`examples/quickstart/`), a nine-entry test suite (`tests/`), the physics,
validation and bootstrap documents (`docs/`) and an environment manifest. No
licence was granted.
