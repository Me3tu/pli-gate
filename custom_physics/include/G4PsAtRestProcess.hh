// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// G4PsAtRestProcess.hh
//******************************************************************************
//
// AtRest-only process for positronium annihilation in PLI simulations.
//
// PURPOSE
// -------
// G4PositroniumAnnihilation (the existing PostStep process) intercepts
// near-thermalized positrons via PostStepDoIt when KE < 1 eV.  However,
// positrons produced by β+ decay and tracked with a 1 mm production cut
// are terminated by the CSDA range mechanism at ~300–400 keV, long before
// reaching KE < 1 eV, and sent directly to the AtRest queue.
//
// This class is an AtRest-only companion process.  It fires when a positron
// enters the AtRest stage (i.e. it was stopped by the range cut before the
// PostStep path could intercept it), applies the same positronium lifetime
// physics, and generates the two 511 keV gammas.
//
// REGISTRATION
// ------------
// As PLIPhysicsList.cc really does it:
//
//   G4PositroniumAnnihilation* psPost = new G4PositroniumAnnihilation();
//   G4PsAtRestProcess*         psRest = new G4PsAtRestProcess();
//   psRest->SetOrthoPsFraction(0.75);
//
//   // One loader call fills the per-material {tau_eff, I3} tables of BOTH
//   // processes, so they cannot drift apart.
//   PLIMaterialLoader::loadInto(psPost, psRest);
//
//   pManager->AddProcess(psPost, -1, -1,  1);   // PostStep only
//   pManager->AddProcess(psRest,  0, -1, -1);   // AtRest only
//
// Only one will fire for each positron track:
//   • PostStep path: β+ steps down to KE < 1 eV → psPost fires
//   • AtRest  path: β+ CSDA-range < cut at ~300 keV → psRest fires
//
// PHYSICS
// -------
// Identical to G4PositroniumAnnihilation::PostStepDoIt:
//   1. Sample Ps formation (prob 35%)
//   2. If Ps forms: sample ortho-Ps fraction (75%), sample exponential decay
//      from material-dependent lifetime (1.75–2.55 ns)
//   3. If direct annihilation: delay = 0
//   4. Generate two back-to-back 511 keV gammas with delay
//   5. Log PLI_LIFETIME oracle line (if PLI_LIFETIME_LOG=1)
//
// Project: Positronium Lifetime Imaging (PLI) - Siemens Quadra
//
//******************************************************************************

#ifndef G4PsAtRestProcess_hh
#define G4PsAtRestProcess_hh 1

#include "G4VRestProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "PLIOrePowell.hh"
#include <map>
#include <string>

class G4PsAtRestProcess : public G4VRestProcess
{
public:
    explicit G4PsAtRestProcess(
        const G4String& processName = "PsAnnihilAtRest");
    virtual ~G4PsAtRestProcess();

    // --- G4VProcess interface -------------------------------------------------

    /// Only applicable to positrons
    virtual G4bool IsApplicable(
        const G4ParticleDefinition& particle) override;

    // --- G4VRestProcess interface ---------------------------------------------

    /// Pure virtual in G4VRestProcess — returns 0 to fire immediately
    virtual G4double GetMeanLifeTime(
        const G4Track& track,
        G4ForceCondition* condition) override;

    /// Returns 0 — always fire immediately when positron enters AtRest stage
    virtual G4double AtRestGetPhysicalInteractionLength(
        const G4Track& track,
        G4ForceCondition* condition) override;

    /// Ps lifetime sampling + gamma generation (same physics as PostStep path)
    virtual G4VParticleChange* AtRestDoIt(
        const G4Track& track,
        const G4Step&  step) override;

    // --- Configuration -------------------------------------------------------

    /// Set material-specific o-Ps lifetime [Geant4 time units = ns by default]
    void SetMaterialLifetime(const G4String& materialName, G4double lifetime);

    /// Copy the entire lifetime map from another source
    void CopyLifetimeMap(const std::map<G4String, G4double>& lifetimes);

    /// Per-material Ps formation probability (I3). Loaded from the JSON
    /// parameter file by PLIMaterialLoader at physics-list init.
    void SetPsFormationProbability(const G4String& materialName, G4double prob);

    /// Scalar fallback for materials NOT in the per-material map.
    /// Warned about once per run.
    void SetPsFormationProbability(G4double prob)  { fPsFormationProb   = prob; }
    void SetOrthoPsFraction(G4double frac)         { fOrthoPosFraction  = frac; }
    void SetVerboseLevel(G4int level)              { fVerboseLevel      = level; }

private:
    // Physics parameters
    std::map<G4String, G4double> fMaterialLifetimeMap;
    std::map<G4String, G4double> fPsFormationMap;   // material -> I3 (JSON-loaded)
    G4double fPsFormationProb;    // Scalar fallback for materials not in fPsFormationMap
    G4double fOrthoPosFraction;   // Ortho-Ps fraction from spin statistics (0.75)
    G4double fParaPsLifetime;     // Para-Ps intrinsic lifetime (0.125 ns)
    G4int    fVerboseLevel;
    bool     fUnknownMaterialWarned;  // one-shot warning for unmapped materials
    mutable bool fLifetimeFallbackWarned = false;  // one-shot G4Exception for
                                                   // unknown-material lifetime
    mutable long fLifetimeFallbackCount = 0;       // run-end fallback tally

    // Per-channel tallies, reported at destruction. Used by the
    // branching test to check P(3 gamma) = tau_eff / tau_3g.
    mutable long fChannelCount[4] = {0, 0, 0, 0};

    // Helpers
    G4double GetLifetime(const G4Material* material) const;
    G4double SamplePsLifetime(const G4Material* material, PLIChannel& channel);
    void     GenerateGammas(G4VParticleChange* change,
                            const G4Track& track,
                            G4double delay,
                            PLIChannel channel);
    /// Material supplying tau_eff / I3 under the active PLI_LIFETIME_MODE
    /// (PLI_LIFETIME_MODE). Default 'transport' returns the annihilation material.
    const G4Material* LifetimeMaterial(const G4Track& track) const;
};

#endif // G4PsAtRestProcess_hh
