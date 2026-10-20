# What this package models

A short statement of the physics, the equations as implemented, and, as
important, the things that are deliberately not modelled.

## The problem

A positron entering matter either annihilates directly or first forms
positronium. Positronium comes in two spin states. Para-positronium (p-Ps)
decays to two photons with a lifetime of 125 ps. Ortho-positronium (o-Ps)
would live 142 ns in vacuum, but in matter it is quenched by pick-off
annihilation against surrounding electrons and lives a few nanoseconds. That
o-Ps lifetime depends on the size of the free volume the positronium sits in
and on the oxygen dissolved there, which is why it carries information about
tissue microstructure.

Positronium lifetime imaging measures that lifetime per location. This package
gives GATE the physics needed to simulate it.

## What is modelled

### 1. Formation

Each annihilation is assigned a channel by a single random draw against the
material's positronium formation probability:

    direct annihilation   probability 1 - P(form)
    p-Ps                  probability P(form) / 4
    o-Ps                  probability P(form) * 3 / 4

the 1:3 ratio being the spin-state multiplicity.

**The `I3_convention` field decides what the database number means, and it is
required.** PALS papers usually quote `I3` as the o-Ps intensity directly. Some
sources instead quote the total formation probability. The two differ by 4/3,
which is a 33 % error in the o-Ps yield if confused. So:

| `I3_convention` | meaning | loader does |
|---|---|---|
| `ops_intensity` | the source reports P(o-Ps) | multiplies by 4/3, so P(o-Ps) = I3 exactly |
| `total_formation` | the source reports P(form any Ps) | uses the value unconverted |

An entry that omits the field is a **hard error**, not a default. A default here
would silently pick one of two answers that differ by a third.

### 2. The o-Ps lifetime: Tao-Eldrup with oxygen quenching

The lifetime is a rate sum:

    1 / tau_eff(R, pO2) = lambda_po(R) + kappa * pO2 + lambda_3g

* `lambda_po(R)`, the pick-off rate, is the Tao-Eldrup model for a spherical
  free volume of radius R:

      lambda_po(R) = lambda_S * [ 1 - R/R_inf + sin(2*pi*R/R_inf) / (2*pi) ]
      R_inf        = R + delta_R

  with `lambda_S` = 2.0 ns^-1 the spin-averaged self-annihilation rate and
  `delta_R` = 0.1656 nm the empirical electron-layer thickness
  (Tao 1972, Eldrup 1981).

* `kappa` = 2.5704e-5 ns^-1 mmHg^-1 is the oxygen quenching constant, from
  Shibuya 2020 converted to a partial-pressure basis with the O2 solubility of
  plasma at 37 °C (1.26 µmol L^-1 mmHg^-1, Christmas and Bassingthwaighte
  2017). This is the value of materials database version 8. Up to version 7
  the conversion gave 2.45e-5, and the solubility of water (1.40) would give
  2.856e-5. The constant reaches the simulation only through the stored
  lifetimes, so between versions 7 and 8 only the five `TumorTissue_pO2_*`
  entries changed (by at most 0.11 ps).

* `lambda_3g` = 1/142 ns^-1 is the vacuum three-photon rate.

The arithmetic is in `tao_eldrup.py` of the analysis part of the release
(see `analysis/README.md`), which is the single source of truth. **The C++ does not re-derive it**: the database carries the already
computed `tau_TE_ns` and the loader reads it. That means a change to the model
requires regenerating the database, not recompiling.

Note that `tau_TE_ns` stores `tau_eff` (the full rate above, including the
three-photon and oxygen terms) and not the pure Tao-Eldrup lifetime
`1/lambda_po(R)`, which the TE model alone would give. The field name is
historical. See the schema table in `README.md` for the difference (26.4 ps on
water) and for how to recover the pure value.

### 3. Decay: two-photon and three-photon channels

o-Ps decays by two competing routes: pick-off, giving two photons, and
intrinsic self-annihilation, giving three. Because the routes compete, the
decay *time* is exponential with mean `tau_eff` regardless of channel, while
the *channel* is chosen with

    P(3 gamma) = lambda_3g * tau_eff = tau_eff / 142.05 ns

For `tau_eff` = 2.0 ns that is 1.41 %.

Three-photon energies follow the Ore-Powell spectrum. With
`x_i = E_i / m_e c^2` and `sum x_i = 2`, the Dalitz density is

    f(x1,x2,x3) = ((1-x1)/(x2 x3))^2 + ((1-x2)/(x1 x3))^2 + ((1-x3)/(x1 x2))^2

bounded by exactly 2 on the physical region, so plain rejection sampling with
envelope 2 is exact. `custom_physics/tests/test_orepowell.cc` checks energy and momentum
conservation and that the envelope is never violated.

## What is NOT modelled

Read this section before drawing conclusions from a simulation.

* **Direct annihilation is prompt by default.** Direct annihilations are
  emitted at exactly zero delay. Real direct annihilation has a finite
  lifetime of a few hundred picoseconds. Set `PLI_DIRECT_TAU_PS` to sample it
  from an exponential instead. This matters more than it looks: with a
  realistic value, a fraction of direct annihilations survives a 1 ns analysis
  censor and a single-exponential lifetime readout is pulled down
  substantially. See `docs/VALIDATION.md`.
* **No positronium transport.** Positronium annihilates where it is formed. Ps
  diffusion before annihilation is not simulated.
* **Two-photon emission is collinear.** No momentum smearing from the
  annihilating pair, so simulated two-photon events have no intrinsic
  acollinearity. A detector-level acollinearity model must be added separately
  if it matters for the study.
* **No temperature dependence.** Tao-Eldrup parameters are used at their
  tabulated values. The database records a measurement temperature where the
  source gives one, but the model does not use it.
* **Free volumes are spherical.** The Tao-Eldrup model assumes a spherical
  well. Real pores are not spherical, and R should be read as an effective
  radius.
* **For tissues, R is back-derived from the published lifetime.** Those
  materials therefore cannot be used to test the model: it would be predicting
  its own input.

## Runtime switches

All are environment variables, read once at first use, and all default to the
production behaviour.

| variable | default | effect |
|---|---|---|
| `PLI_MATERIALS_JSON` | search path | absolute path to the materials database |
| `PLI_DIRECT_TAU_PS` | `0` | direct-annihilation lifetime in ps, `0` means prompt |
| `PLI_ALLOW_MATERIAL_FALLBACK` | `0` | a non-zero integer lets a material absent from the database warn and use tau 2.0 ns and formation probability 0.35 instead of aborting |
| `PLI_DISABLE_3GAMMA` | unset | `1` restores two-photon-only o-Ps decay |
| `PLI_TAU_3GAMMA_NS` | `142.05` | override the vacuum three-photon lifetime |
| `PLI_LIFETIME_MODE` | `transport` | `transport` uses the material where the positron annihilates, `source` uses where it was created, and `fixed` uses one global lifetime |
| `PLI_FIXED_TAU_NS` | `2.0` | the lifetime in ns used when `PLI_LIFETIME_MODE=fixed` |
| `PLI_LIFETIME_LOG` | unset | `1` writes one `PLI_LIFETIME:` record per annihilation to stdout (`simulation/runner/run_gate_parallel.sh` sets it to `1` unless told otherwise) |
| `PLI_LOG_XMAT` | unset | `1` writes one `PLI_R2:` record per annihilation carrying both the emission and the annihilation material |

The in-flight variant of `integration/variants/` reads one more switch,
`PLI_INFLIGHT_STOCK`. With `1` it keeps Geant4's own annihilation process for
positrons that annihilate in flight. Unset, it behaves like the shipped list.
`ENVIRONMENT_MANIFEST.json` lists every switch with the file that reads it.

### The per-annihilation record

With `PLI_LIFETIME_LOG=1`:

    PLI_LIFETIME:<tau_ns>:<material>:<x_mm>:<y_mm>:<z_mm>:<global_time_ns>

`tau_ns` is the delay assigned to that annihilation: 0 for direct, a draw from
125 ps for p-Ps, a draw from the material's `tau_eff` for o-Ps. At end of run
each process prints a channel tally:

    PLI_CHANNEL_TALLY:<process>:direct2g=..:pPs2g=..:oPs2g=..:oPs3g=..:envelope_violations=..

`envelope_violations` must be zero. A non-zero value means the rejection
sampler's bound was exceeded and the spectrum is not trustworthy.

## Relationship to GATE's own positronium classes

GATE 9.4.1 ships `GateOrthoPositronium`, `GateParaPositronium` and related
classes. This package does not use or modify them. `PLIPhysicsList` removes the
standard `G4eplusAnnihilation` process and installs its own PostStep and AtRest
processes, which together cover every positron annihilation. Selecting
`PLI_Physics` therefore replaces the annihilation physics wholesale rather than
extending GATE's.
