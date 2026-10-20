// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// PLISteppingAction.cc
//******************************************************************************

#include "PLISteppingAction.hh"
#include "G4PositroniumAnnihilation.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include <iostream>
#include <iomanip>

// Singleton instance
PLISteppingAction* PLISteppingAction::fgInstance = nullptr;

//******************************************************************************
// Constructor
//******************************************************************************

PLISteppingAction::PLISteppingAction()
    : G4UserSteppingAction(),
      fVerbose(false)
{
    fRecords.reserve(10000);  // Pre-allocate for efficiency
}

//******************************************************************************
// Destructor
//******************************************************************************

PLISteppingAction::~PLISteppingAction()
{
    if (fgInstance == this) {
        fgInstance = nullptr;
    }
}

//******************************************************************************
// Singleton access
//******************************************************************************

PLISteppingAction* PLISteppingAction::GetInstance()
{
    return fgInstance;
}

void PLISteppingAction::SetInstance(PLISteppingAction* instance)
{
    fgInstance = instance;
}

//******************************************************************************
// UserSteppingAction
//******************************************************************************

void PLISteppingAction::UserSteppingAction(const G4Step* step)
{
    // Check if this step produced gammas from positronium annihilation
    const G4Track* track = step->GetTrack();

    // Look for newly created gammas from PositroniumAnnihil process
    const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();

    if (secondaries && secondaries->size() >= 2) {
        // Check if these are annihilation gammas
        bool hasGammaPair = false;
        int gammaCount = 0;

        for (const G4Track* sec : *secondaries) {
            if (sec->GetParticleDefinition() == G4Gamma::Gamma()) {
                // Check if energy is ~511 keV
                G4double energy = sec->GetKineticEnergy();
                if (energy > 500*keV && energy < 520*keV) {
                    gammaCount++;
                }
            }
        }

        if (gammaCount >= 2) {
            // This is likely a positronium annihilation event
            // Try to get lifetime data from G4PositroniumAnnihilation
            G4double lifetime;
            G4String materialName;
            G4ThreeVector position;
            G4int eventID;

            if (G4PositroniumAnnihilation::GetLastLifetimeData(lifetime, materialName,
                                                               position, eventID)) {
                // Get current event ID from run manager
                const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();
                int gateEventID = event ? event->GetEventID() : -1;

                // Create record
                PLILifetimeRecord record;
                record.eventID = gateEventID;
                record.lifetime_ns = lifetime / ns;
                record.posX = position.x() / mm;
                record.posY = position.y() / mm;
                record.posZ = position.z() / mm;
                record.materialName = materialName;
                record.globalTime_ns = step->GetPostStepPoint()->GetGlobalTime() / ns;

                // Thread-safe storage
                {
                    std::lock_guard<std::mutex> lock(fMutex);
                    fRecords.push_back(record);
                }

                if (fVerbose) {
                    G4cout << "PLI: Event " << record.eventID
                           << " Lifetime=" << std::fixed << std::setprecision(4)
                           << record.lifetime_ns << " ns"
                           << " Material=" << record.materialName
                           << " Pos=(" << record.posX << "," << record.posY << "," << record.posZ << ") mm"
                           << G4endl;
                }
            }
        }
    }
}

//******************************************************************************
// WriteToFile - CSV output
//******************************************************************************

void PLISteppingAction::WriteToFile(const std::string& filename)
{
    std::ofstream file(filename);

    if (!file.is_open()) {
        G4cerr << "PLISteppingAction: Cannot open file " << filename << G4endl;
        return;
    }

    // Write header
    file << "eventID,lifetime_ns,posX_mm,posY_mm,posZ_mm,material,globalTime_ns\n";

    // Write data
    for (const auto& rec : fRecords) {
        file << rec.eventID << ","
             << std::fixed << std::setprecision(6) << rec.lifetime_ns << ","
             << std::setprecision(3) << rec.posX << ","
             << rec.posY << ","
             << rec.posZ << ","
             << rec.materialName << ","
             << std::setprecision(6) << rec.globalTime_ns << "\n";
    }

    file.close();

    G4cout << "PLISteppingAction: Wrote " << fRecords.size()
           << " Ps lifetime records to " << filename << G4endl;
}

//******************************************************************************
// WriteToROOT - ROOT TTree output
//******************************************************************************

void PLISteppingAction::WriteToROOT(const std::string& filename)
{
    // This requires ROOT headers - implement if ROOT is available
    // For now, just write CSV and use Python to convert to ROOT

    // Alternative: write binary format that can be read by Python/ROOT
    std::string csvFile = filename;
    size_t dotPos = csvFile.rfind('.');
    if (dotPos != std::string::npos) {
        csvFile = csvFile.substr(0, dotPos) + ".csv";
    } else {
        csvFile += ".csv";
    }

    WriteToFile(csvFile);

    G4cout << "PLISteppingAction: ROOT output requested, wrote CSV to " << csvFile << G4endl;
    G4cout << "  See docs/PHYSICS.md for the record format" << G4endl;
}

//******************************************************************************
// End of PLISteppingAction.cc
//******************************************************************************
