// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// G4PositroniumAnnihilation.hh
//******************************************************************************
//
// Custom discrete process for positron annihilation with positronium formation
// and material-dependent pick-off decay for PLI (Positronium Lifetime Imaging)
//
// Project: PLI Quadra Simulation
//
//******************************************************************************

#ifndef G4PositroniumAnnihilation_hh
#define G4PositroniumAnnihilation_hh 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "PLIOrePowell.hh"
#include <map>
#include <cstdlib>
#include <string>
#include <algorithm>

/**
 * Custom discrete process for positron annihilation with positronium formation
 *
 * Physics Model:
 * 1. Positron thermalizes (handled by standard EM physics)
 * 2. Forms positronium with probability ~35% in tissue
 * 3. Ps state: 75% ortho-Ps (o-Ps), 25% para-Ps (p-Ps)
 * 4. Decay lifetime:
 *    - Para-Ps: 125 ps (intrinsic)
 *    - Ortho-Ps: Material-dependent pick-off (1.8-3 ns in tissue)
 * 5. Annihilation: 2×511 keV back-to-back gammas with time delay
 */
class G4PositroniumAnnihilation : public G4VDiscreteProcess
{
public:
    // Constructor/Destructor
    explicit G4PositroniumAnnihilation(const G4String& processName = "PositroniumAnnihil");
    virtual ~G4PositroniumAnnihilation();

    // Required overrides for G4VDiscreteProcess
    virtual G4double GetMeanFreePath(const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) override;

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                           const G4Step& step) override;

    virtual G4bool IsApplicable(const G4ParticleDefinition& particle) override;

    // Configuration methods
    void SetMaterialLifetime(const G4String& materialName, G4double lifetime);

    // Per-material Ps formation probability (I3). Loaded from the JSON
    // material parameter file by PLIMaterialLoader at physics-list init.
    void SetPsFormationProbability(const G4String& materialName, G4double prob);

    // Scalar fallback used for materials NOT present in fPsFormationMap
    // (e.g. structural geometry materials that were never registered in
    // the JSON). Warned about once per run via the oracle log.
    void SetPsFormationProbability(G4double prob) { fPsFormationProb = prob; }

    void SetOrthoPsFraction(G4double frac) { fOrthoPosFraction = frac; }
    void SetVerboseLevel(G4int level) {
        // Don't override if PLI_LIFETIME_LOG env var is set
        const char* pliLog = std::getenv("PLI_LIFETIME_LOG");
        if (pliLog && std::string(pliLog) == "1") {
            fVerboseLevel = std::max(level, 1);  // At least level 1 for logging
        } else {
            fVerboseLevel = level;
        }
    }

    // Static lifetime-data accessors. These exist as a hook for an external
    // timing actor that is not shipped in this package; in-package the only
    // consumer is PLISteppingAction.
    static G4bool GetLastLifetimeData(G4double& lifetime,
                                      G4String& materialName,
                                      G4ThreeVector& position,
                                      G4int& trackID);
    static void ResetLifetimeData();

private:
    // Physics parameters
    std::map<G4String, G4double> fMaterialLifetimeMap;  // Material-specific o-Ps lifetimes
    std::map<G4String, G4double> fPsFormationMap;       // Material-specific I3 (TE+quenching JSON)
    G4double fPsFormationProb;                          // Scalar fallback for unknown materials
    G4double fOrthoPosFraction;                         // o-Ps fraction (0.75 from spin statistics)
    G4double fParaPsLifetime;                           // p-Ps intrinsic lifetime (125 ps)
    G4int fVerboseLevel;
    bool fUnknownMaterialWarned;                        // one-shot warning for unmapped materials
    mutable bool fLifetimeFallbackWarned = false;       // one-shot unknown-material lifetime warning
    mutable long fLifetimeFallbackCount  = 0;           // run-end fallback tally

    // Per-channel tallies, reported at destruction (see tests/test_orepowell.cc).
    mutable long fChannelCount[4] = {0, 0, 0, 0};

    // Helper methods
    G4double GetLifetime(const G4Material* material) const;
    G4double SamplePsLifetime(const G4Material* material, PLIChannel& channel);
    void GenerateGammas(G4VParticleChange* particleChange,
                       const G4Track& track,
                       G4double delay,
                       PLIChannel channel);
    /// Material supplying tau_eff / I3 under the active PLI_LIFETIME_MODE.
    const G4Material* LifetimeMaterial(const G4Track& track) const;
};

#endif // G4PositroniumAnnihilation_hh
