// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// PLIPhysicsList.cc — Custom physics list for Positronium Lifetime Imaging
// ============================================================================
//
// Registers all physics needed for PLI simulation with GATE:
//   1. G4EmStandardPhysics_option4 — accurate EM for medical physics
//   2. G4DecayPhysics — standard particle decays
//   3. G4RadioactiveDecayPhysics — β+ decay of PET isotopes (Sc-44 etc.)
//   4. G4PositroniumAnnihilation — custom PostStep Ps annihilation (KE < 1 eV)
//   5. G4PsAtRestProcess — custom AtRest Ps annihilation (CSDA-stopped)
//
// The standard G4eplusAnnihilation is removed and replaced by processes
// 4 and 5, which together cover 100% of positron annihilations.
//
// PER-MATERIAL PHYSICS PARAMETERS
//   Material lifetimes (tau_eff from Tao-Eldrup + O2 quenching) and
//   formation probabilities (I3 from PALS literature) are loaded at
//   physics-list init by PLIMaterialLoader from
//
//       materials/pli_materials.json
//
//   The JSON is the single source of truth for the parameter set. To revise
//   a material, edit the JSON; the C++ side picks up the new values at the
//   next run and no code change in this file is required. See
//   docs/PHYSICS.md for the required fields and materials/PROVENANCE.md for
//   where each published value comes from.
//
// GATE integration:
//   Loaded via: /gate/physics/addPhysicsList PLI_Physics
//   Requires GateRunManager mapping: "PLI_Physics" → PLIPhysicsList()
//
// Build: Compile as part of a shared library loaded by GATE at runtime,
//        or integrate directly into GATE source tree and recompile.
//        See CMakeLists.txt for standalone build instructions.
//
// ============================================================================

#include "PLIPhysicsList.hh"
#include "G4PositroniumAnnihilation.hh"
#include "G4PsAtRestProcess.hh"
#include "PLIMaterialLoader.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NuclideTable.hh"

#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4Positron.hh"
#include "G4eplusAnnihilation.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include <stdexcept>

// ============================================================================
// Constructor
// ============================================================================

PLIPhysicsList::PLIPhysicsList()
    : G4VModularPhysicsList(),
      fVerboseLevel(0)
{
    defaultCutValue = 0.1 * mm;
    SetVerboseLevel(fVerboseLevel);

    G4cout << "[PLI] Custom physics list for positronium lifetime imaging" << G4endl;

    RegisterPhysics(new G4EmStandardPhysics_option4(fVerboseLevel));
    RegisterPhysics(new G4DecayPhysics(fVerboseLevel));

    // Lower G4NuclideTable mean-life threshold so Lu-176
    // (T_½ = 3.78×10¹⁰ y, mean life 1.19×10¹⁸ s, just above Geant4's
    // default 1×10¹⁸ s threshold) is registered as unstable and
    // tracked through G4RadioactiveDecay. Without this, Lu-176 ions
    // are treated as stable and never decay; /particle/select Lu176
    // fails.
    //
    // Threshold value: 1.0e10 s -- the value every validation run in
    // docs/VALIDATION.md was produced with.
    // A mean-life threshold excludes nuclides shorter-lived than it, so a very
    // small threshold (e.g. the 1*ns rdecay01 default) enables essentially the
    // whole nuclide table and the per-process memory cost becomes prohibitive.
    // Lu-176's mean life is ~1.7e18 s, far above 1.0e10 s, so it is registered
    // here; that is the isotope this threshold exists to keep available.
    G4NuclideTable::GetInstance()->SetMeanLifeThreshold(1.0e10 * second);

    RegisterPhysics(new G4RadioactiveDecayPhysics(fVerboseLevel));
}

PLIPhysicsList::~PLIPhysicsList() {}

// ============================================================================
// ConstructParticle / ConstructProcess
// ============================================================================

void PLIPhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();
}

void PLIPhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();
    AddCustomPositronAnnihilation();
}

// ============================================================================
// AddCustomPositronAnnihilation
//   Removes standard G4eplusAnnihilation and registers both PostStep and
//   AtRest custom processes with the material lifetime table above.
// ============================================================================

void PLIPhysicsList::AddCustomPositronAnnihilation()
{
    G4cout << "[PLI] Registering custom positronium annihilation processes" << G4endl;

    G4Positron*       positron = G4Positron::Positron();
    G4ProcessManager* pManager = positron->GetProcessManager();

    if (!pManager) {
        G4Exception("PLIPhysicsList::AddCustomPositronAnnihilation",
                    "pli0001", FatalException,
                    "Positron has no process manager!");
        return;
    }

    // ── Remove standard G4eplusAnnihilation ──────────────────────────────
    G4ProcessVector* processList = pManager->GetProcessList();
    for (G4int i = 0; i < (G4int)processList->size(); ++i) {
        if ((*processList)[i]->GetProcessName() == "annihil") {
            G4cout << "[PLI] Removing standard G4eplusAnnihilation" << G4endl;
            pManager->RemoveProcess((*processList)[i]);
            break;
        }
    }

    // ── PostStep process (fires at KE < 1 eV) ───────────────────────────
    G4PositroniumAnnihilation* psPost = new G4PositroniumAnnihilation();
    psPost->SetOrthoPsFraction(0.75);   // spin statistics, universal

    // ── AtRest process (fires for CSDA-stopped positrons) ────────────────
    G4PsAtRestProcess* psRest = new G4PsAtRestProcess();
    psRest->SetOrthoPsFraction(0.75);

    // Respect PLI_LIFETIME_LOG environment variable for both processes
    {
        const char* pliLog = std::getenv("PLI_LIFETIME_LOG");
        G4int oracleVerbose = (pliLog && std::string(pliLog) == "1") ? 1 : fVerboseLevel;
        psPost->SetVerboseLevel(oracleVerbose);
        psRest->SetVerboseLevel(oracleVerbose);
    }

    // ── Load per-material {tau_eff, I3} from JSON into BOTH processes ────
    // tau_eff_ns is precomputed in the JSON from the Tao-Eldrup + O2
    // quenching model; docs/PHYSICS.md gives the equations and
    // materials/PROVENANCE.md the source of every published value.
    try {
        PLIMaterialLoader::loadInto(psPost, psRest);
    } catch (const std::exception& ex) {
        G4Exception("PLIPhysicsList::AddCustomPositronAnnihilation",
                    "pli0002", FatalException, ex.what());
        return;
    }

    // ── Register: PostStep (ordPostStep=1) + AtRest (ordAtRest=0) ────────
    pManager->AddProcess(psPost, -1, -1,  1);   // PostStep only
    pManager->AddProcess(psRest,  0, -1, -1);   // AtRest only

    G4cout << "[PLI] Registered: PostStep (KE<1eV) + AtRest (CSDA-stopped)" << G4endl;

    // ── Verify ───────────────────────────────────────────────────────────
    G4ProcessVector* finalList = pManager->GetProcessList();
    G4cout << "[PLI] Final positron processes (" << finalList->size() << "):" << G4endl;
    for (size_t i = 0; i < finalList->size(); ++i) {
        G4cout << "[PLI]   " << i << ": " << (*finalList)[i]->GetProcessName() << G4endl;
    }
}

// ============================================================================
// ConfigureMaterialLifetimes — no-op (configuration done in constructor)
// ============================================================================

void PLIPhysicsList::ConfigureMaterialLifetimes() {}

// ============================================================================
// SetCuts — production thresholds
// ============================================================================

void PLIPhysicsList::SetCuts()
{
    // Default 0.1 mm — may be overridden per-region in the simulation macro.
    // Note: the production macro typically sets 1 mm for the world region,
    // which causes ~88% of positrons to be CSDA-stopped (→ AtRest path).
    // This is intentional and does not affect lifetime accuracy.
    SetCutValue(0.1 * mm, "gamma");
    SetCutValue(0.1 * mm, "e-");
    SetCutValue(0.1 * mm, "e+");
    SetCutValue(0.1 * mm, "proton");
}

// ============================================================================
// SetVerboseLevel
// ============================================================================

void PLIPhysicsList::SetVerboseLevel(G4int level)
{
    fVerboseLevel = level;
    G4VModularPhysicsList::SetVerboseLevel(level);
}

// End of PLIPhysicsList.cc
