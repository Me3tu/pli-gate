// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// G4PsAtRestProcess.cc — AtRest positronium annihilation process
// ============================================================================
//
// Companion to G4PositroniumAnnihilation. Handles positrons that are
// stopped by Geant4's CSDA range mechanism (residual range < production cut)
// before reaching the 1 eV thermalization threshold of the PostStep process.
//
// In practice with a 1 mm positron production cut in water, approximately
// 88% of Sc-44 β+ positrons are stopped via CSDA and processed here,
// while ~12% thermalize gradually to < 1 eV and are handled by the
// PostStep companion. Both processes implement identical Ps physics.
//
// Registration: PLIPhysicsList registers this with ordAtRest=0, ensuring
//               it fires when the positron enters the AtRest state.
//
// Physics parameters: Set by PLIPhysicsList at runtime via
//                     SetMaterialLifetime() / CopyLifetimeMap().
//                     See G4PositroniumAnnihilation.cc header for
//                     parameter values and literature references.
//
// ============================================================================

#include <cmath>
#include "G4PsAtRestProcess.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4LogicalVolume.hh"
#include "PLIOrePowell.hh"
#include <cstdlib>

// ============================================================================
// Constructor
// ============================================================================

G4PsAtRestProcess::G4PsAtRestProcess(const G4String& processName)
    : G4VRestProcess(processName, fElectromagnetic),
      fPsFormationProb(0.35),
      fOrthoPosFraction(0.75),
      fParaPsLifetime(0.125 * ns),
      fVerboseLevel(0),
      fUnknownMaterialWarned(false)
{
    const char* pliLog = std::getenv("PLI_LIFETIME_LOG");
    if (pliLog && std::string(pliLog) == "1") {
        fVerboseLevel = 1;
    }

    // Fallback only — operational lifetimes set by PLIPhysicsList
    // (see G4PositroniumAnnihilation.cc header for rationale)

    SetProcessSubType(fAnnihilation);   // SubType = 5

    if (fVerboseLevel > 0) {
        G4cout << "G4PsAtRestProcess: AtRest positronium process constructed" << G4endl;
    }
}

// ============================================================================
// Destructor
// ============================================================================

G4PsAtRestProcess::~G4PsAtRestProcess()
{
    if (fLifetimeFallbackCount > 0) {
        G4cout << "[PLI] G4PsAtRestProcess: unknown-material lifetime "
               << "fallback (2.0 ns) was used " << fLifetimeFallbackCount
               << " times this run." << G4endl;
    }
    const long total = fChannelCount[0] + fChannelCount[1]
                     + fChannelCount[2] + fChannelCount[3];
    if (total > 0) {
        // Machine-parseable channel tally: the branching test reads
        // oPs3g / (oPs2g + oPs3g) from this line.
        G4cout << "PLI_CHANNEL_TALLY:G4PsAtRestProcess"
               << ":direct2g=" << fChannelCount[0]
               << ":pPs2g="    << fChannelCount[1]
               << ":oPs2g="    << fChannelCount[2]
               << ":oPs3g="    << fChannelCount[3]
               << ":envelope_violations=" << PLIOrePowell::EnvelopeViolations()
               << G4endl;
    }
}

// ============================================================================
// IsApplicable
// ============================================================================

G4bool G4PsAtRestProcess::IsApplicable(const G4ParticleDefinition& particle)
{
    return (&particle == G4Positron::Positron());
}

// ============================================================================
// GetMeanLifeTime / AtRestGetPhysicalInteractionLength
//   Both return 0 → process fires immediately when positron stops.
// ============================================================================

G4double G4PsAtRestProcess::GetMeanLifeTime(const G4Track&,
                                            G4ForceCondition* condition)
{
    *condition = NotForced;
    return 0.0;
}

G4double G4PsAtRestProcess::AtRestGetPhysicalInteractionLength(
    const G4Track&, G4ForceCondition* condition)
{
    *condition = NotForced;
    return 0.0;
}

// ============================================================================
// AtRestDoIt — identical physics to PostStepDoIt
// ============================================================================

G4VParticleChange* G4PsAtRestProcess::AtRestDoIt(const G4Track& track,
                                                 const G4Step&)
{
    aParticleChange.Initialize(track);

    const G4Material* material    = track.GetMaterial();      // annihilation site
    const G4Material* tauMaterial = LifetimeMaterial(track);   // PLI_LIFETIME_MODE selector
    PLIChannel channel = PLIChannel::kDirect;
    G4double delay = SamplePsLifetime(tauMaterial, channel);
    ++fChannelCount[static_cast<int>(channel)];

    // Oracle log — format UNCHANGED (six fields after the tag) so every
    // existing parser keeps working. The annihilation material is reported,
    // which is the truth quantity the analysis needs.
    if (fVerboseLevel >= 1) {
        G4ThreeVector pos = track.GetPosition();
        G4cout << "PLI_LIFETIME:"
               << delay / ns            << ":"
               << material->GetName()   << ":"
               << pos.x() / mm         << ":"
               << pos.y() / mm         << ":"
               << pos.z() / mm         << ":"
               << track.GetGlobalTime() / ns
               << G4endl;
    }

    // Opt-in record carrying BOTH the emission and the annihilation
    // material, which is what the emission->annihilation transition matrix
    // and the per-arm delay spectra are built from.
    if (PLIOrePowell::LogCrossMaterial()) {
        const G4LogicalVolume* lv = track.GetLogicalVolumeAtVertex();
        const G4String srcName = (lv && lv->GetMaterial())
                                     ? lv->GetMaterial()->GetName()
                                     : G4String("unknown");
        G4cout << "PLI_R2:" << delay / ns << ":" << srcName << ":"
               << material->GetName() << ":"
               << PLIOrePowell::ChannelName(channel) << G4endl;
    }

    GenerateGammas(&aParticleChange, track, delay, channel);

    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.ProposeLocalEnergyDeposit(0.0);

    return &aParticleChange;
}

// ============================================================================
// Physics helpers — identical to G4PositroniumAnnihilation
// ============================================================================

G4double G4PsAtRestProcess::GetLifetime(const G4Material* material) const
{
    // PLI_LIFETIME_MODE=fixed: one global tau for the whole geometry, emulating a
    // single-value lookup table. Materials are otherwise untouched.
    if (PLIOrePowell::LifetimeMode() == PLILifetimeMode::kFixed) {
        return PLIOrePowell::FixedTau();
    }
    auto it = fMaterialLifetimeMap.find(material->GetName());
    if (it != fMaterialLifetimeMap.end()) {
        return it->second;
    }
    if (!PLIOrePowell::AllowMaterialFallback()) {
        G4ExceptionDescription msg;
        msg << "Material '" << material->GetName() << "' has no configured "
            << "o-Ps lifetime in the PLI materials JSON. Earlier versions of this code "
            << "silently used a 2.0 ns fallback, so a compartment that was "
            << "never configured produced ordinary-looking numbers and nothing "
            << "in the output said so. Add the material to the JSON, or set "
            << "PLI_ALLOW_MATERIAL_FALLBACK=1 if the fallback is what you "
            << "actually want.";
        G4Exception("G4PsAtRestProcess::GetLifetime", "pli_missing_material_lifetime", FatalException, msg);
    }
    ++fLifetimeFallbackCount;
    if (!fLifetimeFallbackWarned) {
        G4ExceptionDescription msg;
        msg << "Material '" << material->GetName()
            << "' has no configured o-Ps lifetime; using the 2.0 ns "
            << "fallback. This warning prints once; a total fallback "
            << "count is reported at process destruction.";
        G4Exception("G4PsAtRestProcess::GetLifetime", "pli_fallback_lifetime",
                    JustWarning, msg);
        fLifetimeFallbackWarned = true;
    }
    return 2.0 * ns;  // Unknown material fallback (value unchanged)
}

const G4Material* G4PsAtRestProcess::LifetimeMaterial(const G4Track& track) const
{
    // PLI_LIFETIME_MODE. 'transport' (default) uses the material where the positron
    // actually annihilates; 'source' emulates a source-attached model by
    // using the material of the volume where the positron was created.
    if (PLIOrePowell::LifetimeMode() == PLILifetimeMode::kSource) {
        const G4LogicalVolume* lv = track.GetLogicalVolumeAtVertex();
        if (lv && lv->GetMaterial()) return lv->GetMaterial();
    }
    return track.GetMaterial();
}

// ---------------------------------------------------------------------------
// switchable early-time model for DIRECT annihilation.
// Default 0 reproduces the previous behaviour BIT-EXACTLY: the extra uniform
// draw is short-circuited, so the random stream is untouched when the switch is
// off. Set PLI_DIRECT_TAU_PS=388 to emit direct annihilation from an
// exponential with tau = 388 ps (mercolli2025iodine Table 3 fit constant)
// instead of at exactly zero delay.
// ---------------------------------------------------------------------------
namespace {
G4double PLIDirectTau()
{
    static G4double cached = -1.0;
    if (cached < 0.0) {
        cached = 0.0;
        if (const char* s = std::getenv("PLI_DIRECT_TAU_PS")) {
            const double v = std::atof(s);
            if (v > 0.0) cached = v * CLHEP::picosecond;
        }
        G4cout << "[PLI] direct-annihilation tau = " << cached / CLHEP::picosecond
               << " ps (PLI_DIRECT_TAU_PS; 0 = instantaneous, the default)"
               << G4endl;
    }
    return cached;
}
}  // namespace


G4double G4PsAtRestProcess::SamplePsLifetime(const G4Material* material,
                                             PLIChannel& channel)
{
    // Per-material I3 lookup with one-shot fallback warning. See
    // G4PositroniumAnnihilation::SamplePsLifetime for the rationale —
    // both processes must apply identical physics so the oracle path
    // remains single-source.
    G4double prob = fPsFormationProb;
    auto it = fPsFormationMap.find(material->GetName());
    if (it != fPsFormationMap.end()) {
        prob = it->second;
    } else {
        if (!PLIOrePowell::AllowMaterialFallback()) {
            G4ExceptionDescription msg;
            msg << "Material '" << material->GetName() << "' is not in the PLI "
                << "JSON formation-probability map. Earlier versions of this code "
                << "silently used the scalar I3 fallback. Add the material to "
                << "the JSON, or set PLI_ALLOW_MATERIAL_FALLBACK=1.";
            G4Exception("G4PsAtRestProcess::SamplePsLifetime", "pli_missing_material_i3", FatalException, msg);
        }
        if (!fUnknownMaterialWarned) {
            G4cout << "[PLI] WARNING (G4PsAtRestProcess): material '"
                   << material->GetName()
                   << "' not in JSON formation-probability map; using scalar "
                   << "fallback I3 = " << fPsFormationProb
                   << " (this warning prints once per run)." << G4endl;
            fUnknownMaterialWarned = true;
        }
    }

    if (G4UniformRand() > prob) {
        channel = PLIChannel::kDirect;
        {
            const G4double tauD = PLIDirectTau();
            // Short-circuit when off: consumes no random number, so the
            // pre-V2 stream is reproduced bit-exactly.
            if (tauD <= 0.0) return 0.0;
            return -tauD * std::log(G4UniformRand());
        }
    }
    if (G4UniformRand() < fOrthoPosFraction) {
        const G4double tau   = GetLifetime(material);
        const G4double delay = -tau * std::log(G4UniformRand());
        // Competing risks (see PLIOrePowell.hh): the decay TIME is
        // exponential with the total rate 1/tau_eff whatever the channel;
        // the CHANNEL is then drawn independently. Short-circuit ordering
        // guarantees that PLI_DISABLE_3GAMMA=1 consumes no extra random
        // number, so the two-gamma-only stream is reproduced bit-exactly.
        channel = PLIChannel::kOrthoPs2g;
        if (PLIOrePowell::ThreeGammaEnabled() &&
            G4UniformRand() < PLIOrePowell::ThreeGammaProbability(tau)) {
            channel = PLIChannel::kOrthoPs3g;
        }
        return delay;
    } else {
        channel = PLIChannel::kParaPs;
        return -fParaPsLifetime * std::log(G4UniformRand());
    }
}

void G4PsAtRestProcess::GenerateGammas(G4VParticleChange* change,
                                       const G4Track& track,
                                       G4double delay,
                                       PLIChannel channel)
{
    G4ThreeVector position   = track.GetPosition();
    G4double      globalTime = track.GetGlobalTime() + delay;

    if (channel == PLIChannel::kOrthoPs3g) {
        G4double      e[3];
        G4ThreeVector d[3];
        PLIOrePowell::SampleThreeGamma(e, d);
        change->SetNumberOfSecondaries(3);
        for (G4int i = 0; i < 3; ++i) {
            G4DynamicParticle* g =
                new G4DynamicParticle(G4Gamma::Gamma(), d[i], e[i]);
            G4Track* t = new G4Track(g, globalTime, position);
            t->SetTouchableHandle(track.GetTouchableHandle());
            change->AddSecondary(t);
        }
        if (fVerboseLevel >= 1) {
            G4cout << "PLI_3GAMMA:" << e[0] / keV << ":" << e[1] / keV << ":"
                   << e[2] / keV << ":" << track.GetMaterial()->GetName()
                   << G4endl;
        }
        return;
    }

    G4ThreeVector direction = G4RandomDirection();
    G4double energy = electron_mass_c2;

    G4DynamicParticle* g1 = new G4DynamicParticle(G4Gamma::Gamma(), direction, energy);
    G4Track* t1 = new G4Track(g1, globalTime, position);
    t1->SetTouchableHandle(track.GetTouchableHandle());
    change->AddSecondary(t1);

    G4DynamicParticle* g2 = new G4DynamicParticle(G4Gamma::Gamma(), -direction, energy);
    G4Track* t2 = new G4Track(g2, globalTime, position);
    t2->SetTouchableHandle(track.GetTouchableHandle());
    change->AddSecondary(t2);
}

// ============================================================================
// Configuration — called by PLIPhysicsList
// ============================================================================

void G4PsAtRestProcess::SetMaterialLifetime(const G4String& materialName,
                                            G4double lifetime)
{
    fMaterialLifetimeMap[materialName] = lifetime;
    if (fVerboseLevel > 0) {
        G4cout << "G4PsAtRestProcess: " << materialName
               << " → τ_oPs = " << lifetime / ns << " ns" << G4endl;
    }
}

void G4PsAtRestProcess::CopyLifetimeMap(
    const std::map<G4String, G4double>& lifetimes)
{
    fMaterialLifetimeMap = lifetimes;
}

// ============================================================================
// SetPsFormationProbability — per-material I3 setter
//   Mirrors G4PositroniumAnnihilation::SetPsFormationProbability. Called by
//   PLIMaterialLoader::loadInto for every material+alias from the JSON.
// ============================================================================

void G4PsAtRestProcess::SetPsFormationProbability(const G4String& materialName,
                                                  G4double prob)
{
    fPsFormationMap[materialName] = prob;
    if (fVerboseLevel > 0) {
        G4cout << "G4PsAtRestProcess: " << materialName
               << " → I3 = " << prob << G4endl;
    }
}

// End of G4PsAtRestProcess.cc
