// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// PLIOrePowell.hh: shared final-state physics for the PLI annihilation
//                  processes (three-photon o-Ps decay + channel selection)
//******************************************************************************
//
// Both PLI annihilation processes (G4PositroniumAnnihilation, PostStep, and
// G4PsAtRestProcess, AtRest) must apply IDENTICAL physics so the oracle path
// stays single-source. Everything that defines the final state therefore
// lives here and is called from both.
//
// Three-photon channel 
// -----------------------------------------------------------
// o-Ps decays by two competing routes:
//   pick-off (2 gamma, rate lambda_po, material dependent) and
//   intrinsic self-annihilation (3 gamma, rate lambda_3g, a vacuum constant).
// The observed effective lifetime satisfies 1/tau_eff = lambda_po + lambda_3g,
// so by the standard competing-risks result:
//   * the decay TIME is exponential with mean tau_eff, independent of channel,
//   * the CHANNEL is chosen with P(3 gamma) = lambda_3g * tau_eff
//                                           = tau_eff / tau_3g.
// tau_3g = 142.05 ns  [Vallery, Zitzewitz & Gidley, PRL 90, 203402 (2003)].
// For tau_eff = 2.0 ns this gives P(3 gamma) = 1.41 %.
//
// Energies follow the Ore-Powell spectrum [Ore & Powell, Phys. Rev. 75, 1696
// (1949)]. With x_i = E_i / m_e c^2 and sum_i x_i = 2, the Dalitz density is
//
//   f(x1,x2,x3) = ((1-x1)/(x2 x3))^2 + ((1-x2)/(x1 x3))^2 + ((1-x3)/(x1 x2))^2
//
// NOTE on a form that appears in some write-ups: the density is sometimes
// given as x_i(1-x_i)/(x_j x_k)^2. That form does NOT reproduce the analytic
// Ore-Powell single-photon spectrum (deviations up to 593 %, and it diverges
// at the Dalitz corner). The form above reproduces the analytic spectrum to
// numerical precision and is bounded by 2, and tests/test_orepowell.cc checks
// this.
//
// f is bounded by exactly 2 on the physical region (the bound is attained on
// the edges x_i = 1), so plain 2D rejection with envelope 2 is exact.
//
// Runtime switches (environment variables, read once):
//   PLI_DISABLE_3GAMMA=1     restore the two-gamma-only behaviour
//   PLI_TAU_3GAMMA_NS=<ns>   override tau_3g (sensitivity studies)
//   PLI_LIFETIME_MODE=transport|source|fixed   (alternative lifetime-source models)
//   PLI_FIXED_TAU_NS=<ns>    global tau for PLI_LIFETIME_MODE=fixed
//   PLI_LOG_XMAT=1           emit PLI_R2 per-annihilation records
//
//******************************************************************************

#ifndef PLIOrePowell_hh
#define PLIOrePowell_hh 1

#include "globals.hh"
#include "G4ThreeVector.hh"

/// Final-state channel selected for one positron annihilation.
enum class PLIChannel {
    kDirect    = 0,  ///< no Ps formed        -> 2 gamma, no delay
    kParaPs    = 1,  ///< p-Ps                -> 2 gamma, tau = 125 ps
    kOrthoPs2g = 2,  ///< o-Ps pick-off       -> 2 gamma, tau = tau_eff
    kOrthoPs3g = 3   ///< o-Ps self-annihil.  -> 3 gamma, tau = tau_eff
};

/// Which material supplies tau_eff (alternative lifetime-source models).
enum class PLILifetimeMode {
    kTransport = 0,  ///< material at the ANNIHILATION point (default, physical)
    kSource    = 1,  ///< material at the positron's CREATION vertex
    kFixed     = 2   ///< one global tau for the whole geometry
};

namespace PLIOrePowell {

/// Vacuum three-photon lifetime of o-Ps (Geant4 units).
G4double VacuumTau3Gamma();

/// P(3 gamma | o-Ps) = tau_eff / tau_3g, clamped to [0,1] with a one-shot
/// warning if tau_eff exceeds tau_3g (unphysical configuration).
G4double ThreeGammaProbability(G4double tauEff);

/// False when PLI_DISABLE_3GAMMA=1 (two-gamma-only regression reference).
G4bool ThreeGammaEnabled();

/// PLI_LIFETIME_MODE selector, from PLI_LIFETIME_MODE.
PLILifetimeMode LifetimeMode();

/// Global tau for PLI_LIFETIME_MODE=fixed (Geant4 units).
G4double FixedTau();

/// True when PLI_LOG_XMAT=1 (emit PLI_R2 per-annihilation records).
G4bool LogCrossMaterial();

/// is the unknown-material fallback permitted?
///
/// A material absent from the JSON silently received tau = 2.0 ns and the
/// scalar I3, which is a quiet default where a loud error belongs: the run
/// completes, the numbers look ordinary, and nothing in the output says a
/// compartment was never configured. It is now FATAL by default. Set
/// PLI_ALLOW_MATERIAL_FALLBACK=1 to restore the old permissive behaviour
/// (warn once, tally at destruction) for deliberate probes.
G4bool AllowMaterialFallback();

/// Ore-Powell Dalitz density, 0 outside the physical region.
G4double DalitzDensity(G4double x1, G4double x2, G4double x3);

/// Sample x_i = E_i / m_e c^2 with sum_i x_i = 2 (2D rejection, envelope 2).
void SampleFractions(G4double& x1, G4double& x2, G4double& x3);

/// Complete three-photon final state: energies (Geant4 units) and unit
/// directions satisfying sum_i E_i d_i = 0, with isotropic plane orientation.
void SampleThreeGamma(G4double energy[3], G4ThreeVector dir[3]);

/// Count of rejection-envelope violations (must stay 0, exposed for tests).
long EnvelopeViolations();

/// Human-readable channel name, for tallies and logs.
const char* ChannelName(PLIChannel c);

}  // namespace PLIOrePowell

#endif  // PLIOrePowell_hh
