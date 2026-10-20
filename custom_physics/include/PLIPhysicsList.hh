// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// PLIPhysicsList.hh
//******************************************************************************
//
// Custom Geant4 physics list for Positronium Lifetime Imaging (PLI)
//
// Purpose:
//   Provides a complete physics list that includes:
//     1. Standard electromagnetic physics (Option4 - high accuracy)
//     2. Decay physics (for radioactive sources)
//     3. Radioactive decay physics (for beta+ emitters like 44Sc, 68Ga)
//     4. Custom positronium pick-off annihilation model
//
// Usage in GATE:
//   /gate/physics/addPhysicsList PLI_Physics
//
// Key Features:
//   - Full positron transport (range ~3 mm for 44Sc, ~2 mm for 68Ga)
//   - Material-dependent positronium lifetime (1.8-3.0 ns)
//   - Ortho-Ps vs para-Ps discrimination
//   - Time-delayed gamma emission for PLI reconstruction
//
// Project: Positronium Lifetime Imaging (PLI) - Siemens Quadra
//
//******************************************************************************

#ifndef PLIPhysicsList_hh
#define PLIPhysicsList_hh 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//******************************************************************************
// PLIPhysicsList Class Declaration
//******************************************************************************

class PLIPhysicsList : public G4VModularPhysicsList
{
public:
    // Constructor/Destructor
    explicit PLIPhysicsList();
    virtual ~PLIPhysicsList();

    // Override from base class
    virtual void SetCuts() override;

    // Configuration methods
    void SetVerboseLevel(G4int level);

    // Material-specific lifetime configuration
    // Call this after construction to customize lifetimes
    void ConfigureMaterialLifetimes();

private:
    // Private methods for physics construction
    void ConstructParticle() override;
    void ConstructProcess() override;

    // Helper to add custom positron annihilation
    void AddCustomPositronAnnihilation();

    G4int fVerboseLevel;
};

#endif // PLIPhysicsList_hh
