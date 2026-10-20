// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// G4PositroniumAnnihilation.cc: PostStep positronium annihilation process
// ============================================================================
//
// Custom Geant4 discrete process for positronium lifetime imaging (PLI).
// Intercepts thermalized positrons (KE < 1 eV) and applies:
//   1. Ps formation probability (default 35% in tissue)
//   2. o-Ps / p-Ps branching (3:1 spin statistics)
//   3. Material-dependent o-Ps pick-off lifetime (exponential sampling)
//   4. Two back-to-back 511 keV gamma generation with proper time delay
//
// Trigger: G4VDiscreteProcess::PostStepDoIt fires when GetMeanFreePath
//          returns 0 (i.e., when positron KE drops below 1 eV).
//
// Companion: G4PsAtRestProcess handles positrons stopped by the CSDA
//            range cut before reaching 1 eV. Both processes implement
//            identical Ps physics. See PLIPhysicsList for registration.
//
// Oracle: When environment variable PLI_LIFETIME_LOG=1, each annihilation
//         event is logged to stdout in machine-parseable format:
//           PLI_LIFETIME:delay_ns:material:posX_mm:posY_mm:posZ_mm:thermalTime_ns
//         where delay_ns is the sampled Ps lifetime and thermalTime_ns is the
//         positron's global time at thermalization (before adding delay).
//
// Physics parameters and literature references:
//   Ps formation:  35% in condensed matter [Moskal et al., Phys Med Biol 2019]
//   o-Ps fraction:  75% (spin statistics: triplet:singlet = 3:1)
//   p-Ps lifetime:  0.125 ns (vacuum value, negligible medium correction)
//   o-Ps lifetimes: Material-dependent, set by PLIPhysicsList at runtime.
//                   Constructor defaults are fallbacks only.
//
// Process type: fElectromagnetic (type 2), fAnnihilation (subtype 5)
//   Matches standard G4eplusAnnihilation signature so that GATE's
//   digitizer and actor systems recognize the resulting gammas.
//
// ============================================================================

#include <cmath>
#include "G4PositroniumAnnihilation.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4LogicalVolume.hh"
#include "G4Exception.hh"
#include "PLIOrePowell.hh"
#include <mutex>
#include <cstdlib>
#include <string>

// Thread-safe static storage for last Ps decay information
// Read by GetLastLifetimeData(), a hook for an external timing actor that
// is NOT part of this package. In-package the only consumer is
// PLISteppingAction. Safe for single-threaded GATE.
static std::mutex g_lifetimeMutex;
static G4double g_lastPsLifetime = 0.0;
static G4String g_lastMaterialName = "";
static G4ThreeVector g_lastPosition;
static G4int g_lastTrackID = -1;
static G4bool g_lifetimeRecorded = false;

// ============================================================================
// Constructor
// ============================================================================

G4PositroniumAnnihilation::G4PositroniumAnnihilation(const G4String& processName)
    : G4VDiscreteProcess(processName),
      fPsFormationProb(0.35),
      fOrthoPosFraction(0.75),
      fParaPsLifetime(0.125 * ns),
      fVerboseLevel(0),
      fUnknownMaterialWarned(false)
{
    // Oracle logging via environment variable
    const char* pliLog = std::getenv("PLI_LIFETIME_LOG");
    if (pliLog && std::string(pliLog) == "1") {
        fVerboseLevel = 1;
        G4cout << "G4PositroniumAnnihilation: PLI lifetime logging ENABLED" << G4endl;
    }

    // Fallback material lifetimes, used ONLY for materials not explicitly
    // configured by PLIPhysicsList. In production, PLIPhysicsList::
    // AddCustomPositronAnnihilation() calls SetMaterialLifetime() for all
    // relevant materials, overriding these defaults.
    //
    // The default of 2.0 ns (generic soft tissue) is intentionally
    // conservative. Any material not in the map will use this value and
    // will be identifiable in the oracle output by its material name.
    fMaterialLifetimeMap["default_fallback"] = 2.0 * ns;

    // Match standard annihilation process signature for GATE compatibility
    SetProcessType(fElectromagnetic);    // ProcessType = 2
    SetProcessSubType(fAnnihilation);    // SubType = 5
}

// ============================================================================
// Destructor
// ============================================================================

G4PositroniumAnnihilation::~G4PositroniumAnnihilation()
{
    if (fLifetimeFallbackCount > 0) {
        G4cout << "[PLI] G4PositroniumAnnihilation: unknown-material lifetime "
               << "fallback (2.0 ns) was used " << fLifetimeFallbackCount
               << " times this run." << G4endl;
    }
    const long total = fChannelCount[0] + fChannelCount[1]
                     + fChannelCount[2] + fChannelCount[3];
    if (total > 0) {
        G4cout << "PLI_CHANNEL_TALLY:G4PositroniumAnnihilation"
               << ":direct2g=" << fChannelCount[0]
               << ":pPs2g="    << fChannelCount[1]
               << ":oPs2g="    << fChannelCount[2]
               << ":oPs3g="    << fChannelCount[3]
               << ":envelope_violations=" << PLIOrePowell::EnvelopeViolations()
               << G4endl;
    }
}

// ============================================================================
// IsApplicable: only positrons
// ============================================================================

G4bool G4PositroniumAnnihilation::IsApplicable(const G4ParticleDefinition& particle)
{
    return (&particle == G4Positron::Positron());
}

// ============================================================================
// GetMeanFreePath: triggers at thermalization (KE < 1 eV)
// ============================================================================

G4double G4PositroniumAnnihilation::GetMeanFreePath(const G4Track& track,
                                                    G4double,
                                                    G4ForceCondition* condition)
{
    *condition = NotForced;
    if (track.GetKineticEnergy() < 1.0 * eV) {
        return 0.0;   // Immediate interaction: positron is thermalized
    }
    return DBL_MAX;    // Transparent: let other processes handle transport
}

// ============================================================================
// GetLifetime: material-dependent o-Ps pick-off lifetime lookup
// ============================================================================

G4double G4PositroniumAnnihilation::GetLifetime(const G4Material* material) const
{
    // PLI_LIFETIME_MODE=fixed: one global tau for the whole geometry (see the AtRest
    // twin: both processes must apply identical physics).
    if (PLIOrePowell::LifetimeMode() == PLILifetimeMode::kFixed) {
        return PLIOrePowell::FixedTau();
    }
    G4String matName = material->GetName();
    auto it = fMaterialLifetimeMap.find(matName);
    if (it != fMaterialLifetimeMap.end()) {
        return it->second;
    }
    // Unknown material: conservative default, logged once (mirrors the
    // AtRest twin so a missing table entry can never pass silently).
    if (!PLIOrePowell::AllowMaterialFallback()) {
        G4ExceptionDescription msg;
        msg << "Material '" << material->GetName() << "' has no configured "
            << "o-Ps lifetime in the PLI materials JSON. Earlier versions of this code "
            << "silently used a 2.0 ns fallback, so a compartment that was "
            << "never configured produced ordinary-looking numbers and nothing "
            << "in the output said so. Add the material to the JSON, or set "
            << "PLI_ALLOW_MATERIAL_FALLBACK=1 if the fallback is what you "
            << "actually want.";
        G4Exception("G4PositroniumAnnihilation::GetLifetime", "pli_missing_material_lifetime", FatalException, msg);
    }
    ++fLifetimeFallbackCount;
    if (!fLifetimeFallbackWarned) {
        G4ExceptionDescription msg;
        msg << "Material '" << matName << "' has no configured o-Ps lifetime; "
            << "using the 2.0 ns fallback. This warning prints once; a total "
            << "fallback count is reported at process destruction.";
        G4Exception("G4PositroniumAnnihilation::GetLifetime",
                    "pli_fallback_lifetime", JustWarning, msg);
        fLifetimeFallbackWarned = true;
    }
    return 2.0 * ns;
}

const G4Material* G4PositroniumAnnihilation::LifetimeMaterial(
    const G4Track& track) const
{
    if (PLIOrePowell::LifetimeMode() == PLILifetimeMode::kSource) {
        const G4LogicalVolume* lv = track.GetLogicalVolumeAtVertex();
        if (lv && lv->GetMaterial()) return lv->GetMaterial();
    }
    return track.GetMaterial();
}

// ============================================================================
// SamplePsLifetime: Monte Carlo sampling of the Ps decay pathway
//
// Branching (per material, I3 from JSON via PLIMaterialLoader):
//   P(direct annihilation) = 1 - I3
//   P(o-Ps)                = I3 × fOrthoPosFraction (≈ 0.75 by spin stats)
//                            → delay ~ Exp(1/τ_material)
//   P(p-Ps)                = I3 × (1 - fOrthoPosFraction)
//                            → delay ~ Exp(1/0.125 ns)
//
// Materials not present in fPsFormationMap fall back to fPsFormationProb
// (the scalar default, 0.35) and emit a one-shot warning so the operator
// can spot a missing JSON entry.
//
// Sampling: inverse CDF method, -τ × ln(U), U ~ Uniform(0,1)
// ============================================================================

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


G4double G4PositroniumAnnihilation::SamplePsLifetime(const G4Material* material,
                                                     PLIChannel& channel)
{
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
            G4Exception("G4PositroniumAnnihilation::SamplePsLifetime", "pli_missing_material_i3", FatalException, msg);
        }
        if (!fUnknownMaterialWarned) {
            G4cout << "[PLI] WARNING (G4PositroniumAnnihilation): material '"
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
        // Ortho-Ps: material-dependent effective lifetime. The decay TIME is
        // exponential with 1/tau_eff whatever the channel. The CHANNEL is
        // drawn independently with P(3 gamma) = tau_eff / tau_3g (competing
        // risks, see PLIOrePowell.hh). Short-circuit ordering keeps the
        // two-gamma-only random stream bit-exact under PLI_DISABLE_3GAMMA=1.
        const G4double tau   = GetLifetime(material);
        const G4double delay = -tau * std::log(G4UniformRand());
        channel = PLIChannel::kOrthoPs2g;
        if (PLIOrePowell::ThreeGammaEnabled() &&
            G4UniformRand() < PLIOrePowell::ThreeGammaProbability(tau)) {
            channel = PLIChannel::kOrthoPs3g;
        }
        return delay;
    } else {
        // Para-Ps: intrinsic 125 ps lifetime
        channel = PLIChannel::kParaPs;
        return -fParaPsLifetime * std::log(G4UniformRand());
    }
}

// ============================================================================
// GenerateGammas: two back-to-back 511 keV gammas at annihilation position
// ============================================================================

void G4PositroniumAnnihilation::GenerateGammas(G4VParticleChange* particleChange,
                                               const G4Track& track,
                                               G4double delay,
                                               PLIChannel channel)
{
    G4ThreeVector position = track.GetPosition();
    G4double globalTime = track.GetGlobalTime() + delay;

    if (channel == PLIChannel::kOrthoPs3g) {
        G4double      e[3];
        G4ThreeVector d[3];
        PLIOrePowell::SampleThreeGamma(e, d);
        particleChange->SetNumberOfSecondaries(3);
        for (G4int i = 0; i < 3; ++i) {
            G4DynamicParticle* g =
                new G4DynamicParticle(G4Gamma::Gamma(), d[i], e[i]);
            G4Track* t = new G4Track(g, globalTime, position);
            t->SetTouchableHandle(track.GetTouchableHandle());
            particleChange->AddSecondary(t);
        }
        if (fVerboseLevel >= 1) {
            G4cout << "PLI_3GAMMA:" << e[0] / keV << ":" << e[1] / keV << ":"
                   << e[2] / keV << ":" << track.GetMaterial()->GetName()
                   << G4endl;
        }
        return;
    }

    G4ThreeVector direction = G4RandomDirection();
    G4double energy = electron_mass_c2;  // 511 keV

    G4DynamicParticle* gamma1 = new G4DynamicParticle(
        G4Gamma::Gamma(), direction, energy);
    G4Track* track1 = new G4Track(gamma1, globalTime, position);
    track1->SetTouchableHandle(track.GetTouchableHandle());
    particleChange->AddSecondary(track1);

    G4DynamicParticle* gamma2 = new G4DynamicParticle(
        G4Gamma::Gamma(), -direction, energy);
    G4Track* track2 = new G4Track(gamma2, globalTime, position);
    track2->SetTouchableHandle(track.GetTouchableHandle());
    particleChange->AddSecondary(track2);
}

// ============================================================================
// PostStepDoIt: main entry point when positron thermalizes
// ============================================================================

G4VParticleChange* G4PositroniumAnnihilation::PostStepDoIt(const G4Track& track,
                                                           const G4Step&)
{
    aParticleChange.Initialize(track);

    const G4Material* material    = track.GetMaterial();     // annihilation site
    const G4Material* tauMaterial = LifetimeMaterial(track);  // PLI_LIFETIME_MODE selector
    PLIChannel channel = PLIChannel::kDirect;
    G4double delay = SamplePsLifetime(tauMaterial, channel);
    ++fChannelCount[static_cast<int>(channel)];

    // Record for the external-actor hook (see GetLastLifetimeData)
    {
        std::lock_guard<std::mutex> lock(g_lifetimeMutex);
        g_lastPsLifetime = delay;
        g_lastMaterialName = material->GetName();
        g_lastPosition = track.GetPosition();
        g_lastTrackID = track.GetTrackID();
        g_lifetimeRecorded = true;
    }

    // Oracle log: machine-parseable lifetime data
    if (fVerboseLevel >= 1) {
        G4ThreeVector pos = track.GetPosition();
        G4cout << "PLI_LIFETIME:"
               << delay / ns << ":"
               << material->GetName() << ":"
               << pos.x() / mm << ":"
               << pos.y() / mm << ":"
               << pos.z() / mm << ":"
               << track.GetGlobalTime() / ns
               << G4endl;
    }

    // Opt-in record carrying BOTH emission and annihilation material.
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
// Static accessors for the external-actor hook
// ============================================================================

G4bool G4PositroniumAnnihilation::GetLastLifetimeData(G4double& lifetime,
                                                       G4String& materialName,
                                                       G4ThreeVector& position,
                                                       G4int& trackID)
{
    std::lock_guard<std::mutex> lock(g_lifetimeMutex);
    if (!g_lifetimeRecorded) return false;
    lifetime = g_lastPsLifetime;
    materialName = g_lastMaterialName;
    position = g_lastPosition;
    trackID = g_lastTrackID;
    g_lifetimeRecorded = false;
    return true;
}

void G4PositroniumAnnihilation::ResetLifetimeData()
{
    std::lock_guard<std::mutex> lock(g_lifetimeMutex);
    g_lifetimeRecorded = false;
}

// ============================================================================
// SetMaterialLifetime: called by PLIPhysicsList to configure lifetimes
// ============================================================================

void G4PositroniumAnnihilation::SetMaterialLifetime(const G4String& materialName,
                                                    G4double lifetime)
{
    fMaterialLifetimeMap[materialName] = lifetime;
    if (fVerboseLevel > 0) {
        G4cout << "G4PositroniumAnnihilation: " << materialName
               << " → τ_oPs = " << lifetime / ns << " ns" << G4endl;
    }
}

// ============================================================================
// SetPsFormationProbability: per-material I3 setter
//   Called by PLIMaterialLoader::loadInto for every material+alias from
//   the JSON parameter file at physics-list init.
// ============================================================================

void G4PositroniumAnnihilation::SetPsFormationProbability(
    const G4String& materialName, G4double prob)
{
    fPsFormationMap[materialName] = prob;
    if (fVerboseLevel > 0) {
        G4cout << "G4PositroniumAnnihilation: " << materialName
               << " → I3 = " << prob << G4endl;
    }
}

// End of G4PositroniumAnnihilation.cc
