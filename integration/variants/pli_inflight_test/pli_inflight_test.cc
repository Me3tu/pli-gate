// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// pli_inflight_test.cc  (2026-09-30, feedback round 3: in-flight annihilation)
// Derived from pli_boundary_test.cc (referee round 2). Changes: physics mode
// 'hybrid' (shipped list + stock annihil kept for IN-FLIGHT annihilation via
// PLI_INFLIGHT_STOCK=1, see PLIPhysicsList_inflight.cc), source 'i124'
// (two main beta+ branches), and per-photon energies in the output.
// Original header follows.
//
//
// Standalone Geant4 application that links the SHIPPED PLI process sources
// (the same .cc files, sha256-identical to the copies compiled into GATE) and
// audits, per annihilation:
//   * which process fired (PostStep "PositroniumAnnihil" / AtRest
//     "PsAnnihilAtRest", stock "annihil" in the reference mode)
//   * the positron kinetic energy at the annihilation step and the step length
//   * the material the process used (pre-step point of the annihilation step,
//     i.e. exactly what G4Track::GetMaterial() returns inside the DoIt)
//   * the TRUE material at the vertex, located by an independent G4Navigator
//   * the signed distance of the vertex from a planar interface at x = 0
//   * the emitted photons: number, sum of energies, |sum of momenta|,
//     emission time relative to the positron time, emission position
//
// Geometry: world half-space x < 0 = material A, daughter box x in [0, 200] mm
// = material B, one exactly planar interface at x = 0.
//
// Usage:
//   pli_boundary_test <physics pli|stock> <source> <N> <cut_mm> <seed> <out.csv>
//   source = sc44                 44Sc beta+ spectrum, x0 ~ U(-1, 1) mm, isotropic
//          | ctrl:<d_mm>:<side>:<KE_eV>
//                                  every positron starts at x0 = side*d_mm with
//                                  kinetic energy KE_eV (0 -> AtRest path,
//                                  0.5 -> PostStep path), i.e. annihilation is
//                                  forced at a controlled distance from the plane
// ============================================================================
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserSteppingAction.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4UImanager.hh"
#include "PLIPhysicsList.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace {
G4String gMatA = "water_shibuya2020_n2";   // x < 0
G4String gMatB = "adipose_moskal2023";     // x > 0
std::FILE* gOut = nullptr;
}

// ---------------------------------------------------------------- geometry
class Detector : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    // Materials carry the exact names of the PLI database entries, and the
    // lookup is by name, so these are the entries the processes will use.
    G4Material* a = nist->BuildMaterialWithNewDensity(gMatA, "G4_WATER", 1.00 * g / cm3);
    G4Material* b = nist->BuildMaterialWithNewDensity(gMatB, "G4_ADIPOSE_TISSUE_ICRP", 0.95 * g / cm3);
    auto* sW = new G4Box("world", 200 * mm, 200 * mm, 200 * mm);
    auto* lW = new G4LogicalVolume(sW, a, "world");
    auto* pW = new G4PVPlacement(nullptr, {}, lW, "world", nullptr, false, 0, true);
    auto* sB = new G4Box("half_B", 100 * mm, 200 * mm, 200 * mm);
    auto* lB = new G4LogicalVolume(sB, b, "half_B");
    new G4PVPlacement(nullptr, G4ThreeVector(100 * mm, 0, 0), lB, "half_B", lW, false, 0, true);
    return pW;
  }
};

// ---------------------------------------------------------------- stock physics (reference)
class StockPhysics : public G4VModularPhysicsList {
public:
  StockPhysics() {
    RegisterPhysics(new G4EmStandardPhysics_option4(0));
    RegisterPhysics(new G4DecayPhysics(0));
  }
};

// ---------------------------------------------------------------- primaries
class Primary : public G4VUserPrimaryGeneratorAction {
public:
  explicit Primary(const std::string& src) {
    fGun = new G4ParticleGun(1);
    fGun->SetParticleDefinition(G4Positron::Positron());
    if (src.rfind("ctrl:", 0) == 0) {
      std::stringstream ss(src.substr(5));
      std::string t; std::vector<double> v;
      while (std::getline(ss, t, ':')) v.push_back(std::atof(t.c_str()));
      fCtrl = true; fD = v.at(0) * mm; fSide = v.at(1); fKE = v.at(2) * eV;
    }
    fI124 = (src == "i124");
  }
  ~Primary() override { delete fGun; }

  // 44Sc beta+ kinetic energy: allowed shape with a non-relativistic Fermi
  // factor for a positron leaving 44Ca (Z = 20), endpoint 1.474 MeV.
  G4double SampleSc44() {
    const G4double me = electron_mass_c2, E0 = 1.474 * MeV, Z = 20., al = fine_structure_const;
    auto f = [&](G4double T) {
      const G4double W = T / me + 1., W0 = E0 / me + 1., p = std::sqrt(W * W - 1.);
      if (p <= 0) return 0.;
      const G4double eta = -Z * al * W / p, x = 2. * pi * eta;
      const G4double F = x / (1. - std::exp(-x));
      return p * W * (W0 - W) * (W0 - W) * F;
    };
    if (fFmax < 0) { for (int i = 1; i < 2000; ++i) fFmax = std::max(fFmax, f(E0 * i / 2000.)); fFmax *= 1.05; }
    while (true) { const G4double T = E0 * G4UniformRand(); if (G4UniformRand() * fFmax < f(T)) return T; }
  }

  // 124I beta+: branches E0 = 2.1376 MeV (10.7 %) and 1.5349 MeV (11.7 %),
  // daughter 124Te (Z = 52). Same allowed shape and Fermi approximation.
  G4double SampleI124() {
    const G4double me = electron_mass_c2, Z = 52., al = fine_structure_const;
    const G4double E0 = (G4UniformRand() < 10.7 / (10.7 + 11.7)) ? 2.1376 * MeV : 1.5349 * MeV;
    auto f = [&](G4double T) {
      const G4double W = T / me + 1., W0 = E0 / me + 1., p = std::sqrt(W * W - 1.);
      if (p <= 0) return 0.;
      const G4double eta = -Z * al * W / p, x = 2. * pi * eta;
      const G4double F = x / (1. - std::exp(-x));
      return p * W * (W0 - W) * (W0 - W) * F;
    };
    G4double fmax = -1; for (int i = 1; i < 2000; ++i) fmax = std::max(fmax, f(E0 * i / 2000.)); fmax *= 1.05;
    while (true) { const G4double T = E0 * G4UniformRand(); if (G4UniformRand() * fmax < f(T)) return T; }
  }

  void GeneratePrimaries(G4Event* ev) override {
    if (fCtrl) {
      fGun->SetParticlePosition(G4ThreeVector(fSide * fD, 40 * mm * (G4UniformRand() - .5), 40 * mm * (G4UniformRand() - .5)));
      fGun->SetParticleEnergy(fKE);
    } else {
      fGun->SetParticlePosition(G4ThreeVector(2 * mm * (G4UniformRand() - .5), 40 * mm * (G4UniformRand() - .5), 40 * mm * (G4UniformRand() - .5)));
      fGun->SetParticleEnergy(fI124 ? SampleI124() : SampleSc44());
    }
    fGun->SetParticleMomentumDirection(G4RandomDirection());
    fGun->GeneratePrimaryVertex(ev);
  }
private:
  G4ParticleGun* fGun;
  bool fCtrl = false, fI124 = false; G4double fD = 0, fSide = 1, fKE = 0, fFmax = -1;
};

// ---------------------------------------------------------------- per-annihilation record
class Stepping : public G4UserSteppingAction {
public:
  void UserSteppingAction(const G4Step* step) override {
    const G4Track* trk = step->GetTrack();
    if (trk->GetDefinition() != G4Positron::Positron()) return;
    const G4StepPoint* pre = step->GetPreStepPoint();
    const G4StepPoint* post = step->GetPostStepPoint();
    const G4VProcess* proc = post->GetProcessDefinedStep();
    if (!proc) return;
    const G4String& pn = proc->GetProcessName();
    if (pn != "PositroniumAnnihil" && pn != "PsAnnihilAtRest" && pn != "annihil") return;

    if (!fNav) {
      fNav = new G4Navigator();
      fNav->SetWorldVolume(G4TransportationManager::GetTransportationManager()
                               ->GetNavigatorForTracking()->GetWorldVolume());
    }
    const G4ThreeVector pos = post->GetPosition();
    const G4VPhysicalVolume* pv = fNav->LocateGlobalPointAndSetup(pos, nullptr, false, true);
    const G4String truthMat = pv ? pv->GetLogicalVolume()->GetMaterial()->GetName() : G4String("none");
    const G4String usedMat = pre->GetMaterial()->GetName();   // == track.GetMaterial() inside DoIt
    const G4bool atRest = (post->GetStepStatus() == fAtRestDoItProc);

    G4double sumE = 0, tmin = 1e300, tmax = -1e300, dmax = 0;
    G4ThreeVector sumP;
    G4int ng = 0;
    std::vector<G4double> eg;
    if (const auto* sec = step->GetSecondaryInCurrentStep()) {
      for (const G4Track* s : *sec) {
        if (s->GetDefinition() != G4Gamma::Gamma()) continue;
        ++ng; sumE += s->GetKineticEnergy(); sumP += s->GetMomentum(); eg.push_back(s->GetKineticEnergy());
        tmin = std::min(tmin, s->GetGlobalTime()); tmax = std::max(tmax, s->GetGlobalTime());
        dmax = std::max(dmax, (s->GetPosition() - pos).mag());
      }
    }
    const G4double tPos = post->GetGlobalTime();
    std::sort(eg.begin(), eg.end(), std::greater<G4double>());
    while (eg.size() < 3) eg.push_back(0.);
    std::fprintf(gOut, "%s,%d,%s,%s,%.15g,%.15g,%.15g,%.9g,%.9g,%.9g,%d,%.12g,%.6g,%.9g,%.9g,%.3g,%.6g,%.6g,%.6g\n",
                 pn.c_str(), atRest ? 1 : 0, usedMat.c_str(), truthMat.c_str(),
                 pos.x() / mm, pos.y() / mm, pos.z() / mm,
                 pre->GetKineticEnergy() / eV, post->GetKineticEnergy() / eV, step->GetStepLength() / mm,
                 ng, sumE / keV, sumP.mag() / keV,
                 ng ? (tmin - tPos) / ns : -1.0, ng ? (tmax - tmin) / ns : -1.0, dmax / mm,
                 eg[0] / keV, eg[1] / keV, eg[2] / keV);
  }
private:
  G4Navigator* fNav = nullptr;
};

class Actions : public G4VUserActionInitialization {
public:
  explicit Actions(const std::string& src) : fSrc(src) {}
  void Build() const override {
    SetUserAction(new Primary(fSrc));
    SetUserAction(new Stepping());
  }
private:
  std::string fSrc;
};

int main(int argc, char** argv) {
  if (argc < 7) { std::fprintf(stderr, "usage: %s pli|hybrid|stock sc44|i124|ctrl:d:s:KE N cut_mm seed out.csv\n", argv[0]); return 2; }
  const std::string phys = argv[1], src = argv[2];
  const long N = std::atol(argv[3]);
  const double cut = std::atof(argv[4]);
  const long seed = std::atol(argv[5]);
  gOut = std::fopen(argv[6], "w");
  std::fprintf(gOut, "process,at_rest,used_mat,truth_mat,x_mm,y_mm,z_mm,preKE_eV,postKE_eV,steplen_mm,"
                     "n_gamma,sumE_keV,sumP_keV,delay_ns,dt_spread_ns,emit_offset_mm,e1_keV,e2_keV,e3_keV\n");
  CLHEP::HepRandom::setTheSeed(seed);

  auto* rm = new G4RunManager();
  rm->SetUserInitialization(new Detector());
  if (phys == "hybrid") setenv("PLI_INFLIGHT_STOCK", "1", 1);
  else unsetenv("PLI_INFLIGHT_STOCK");
  if (phys == "stock") rm->SetUserInitialization(new StockPhysics());
  else rm->SetUserInitialization(new PLIPhysicsList());
  rm->SetUserInitialization(new Actions(src));
  auto* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/run/setCut " + std::to_string(cut) + " mm");
  rm->Initialize();
  // Print the positron process table (exactly one annihilation process per mode)
  auto* pm = G4Positron::Positron()->GetProcessManager();
  auto* pl = pm->GetProcessList();
  for (size_t i = 0; i < pl->size(); ++i) {
    G4VProcess* p = (*pl)[i];
    G4cout << "POSITRON_PROCESS:" << p->GetProcessName()
           << ":atRest=" << pm->GetAtRestIndex(p)
           << ":alongStep=" << pm->GetAlongStepIndex(p)
           << ":postStep=" << pm->GetPostStepIndex(p) << G4endl;
  }
  rm->BeamOn(N);
  std::fclose(gOut);
  delete rm;
  return 0;
}
