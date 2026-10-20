// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// PLISteppingAction.hh
//******************************************************************************
//
// Geant4 stepping action to record Ps lifetime data for PLI analysis
//
// This class intercepts positronium annihilation events and records:
// - Ps lifetime (delay before gamma emission)
// - Material name where annihilation occurred
// - Position of annihilation
// - Event ID for correlation with GATE output
//
// Project: PLI Quadra Simulation
//
//******************************************************************************

#ifndef PLISteppingAction_hh
#define PLISteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <string>
#include <fstream>
#include <mutex>

// Data structure for Ps lifetime record
struct PLILifetimeRecord {
    int eventID;
    double lifetime_ns;       // Ps lifetime in nanoseconds
    double posX, posY, posZ;  // Position in mm
    std::string materialName;
    double globalTime_ns;     // Global simulation time in ns
};

class PLISteppingAction : public G4UserSteppingAction
{
public:
    PLISteppingAction();
    virtual ~PLISteppingAction();

    virtual void UserSteppingAction(const G4Step* step) override;

    // Data access
    const std::vector<PLILifetimeRecord>& GetRecords() const { return fRecords; }

    // File output
    void WriteToFile(const std::string& filename);
    void WriteToROOT(const std::string& filename);

    // Control
    void Clear() { fRecords.clear(); }
    void SetVerbose(bool v) { fVerbose = v; }

    // Singleton access (for global registration)
    static PLISteppingAction* GetInstance();
    static void SetInstance(PLISteppingAction* instance);

private:
    std::vector<PLILifetimeRecord> fRecords;
    bool fVerbose;
    std::mutex fMutex;

    static PLISteppingAction* fgInstance;
};

#endif // PLISteppingAction_hh
