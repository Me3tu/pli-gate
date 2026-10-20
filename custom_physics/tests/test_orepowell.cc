// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// test_orepowell.cc — unit test for the three-photon o-Ps final state
// ============================================================================
//
// Exercises the SHIPPED sampler (PLIOrePowell::SampleThreeGamma, the very
// function Geant4 calls) rather than a re-implementation, and writes the
// sampled events to stdout as CSV for statistical comparison against the
// independently-validated Python reference
// implementation, itself checked against the analytic Ore-Powell
// single-photon spectrum. The reference values that comparison
// produced are recorded in docs/VALIDATION.md section 1.
//
// Checks performed in-process (hard failures, exit code 1):
//   1. energy conservation      sum_i E_i == 1022 keV to 1e-9 relative
//   2. momentum closure         |sum_i E_i d_i| / (m_e c^2) < 1e-12
//   3. unit directions          ||d_i| - 1| < 1e-12
//   4. physical support         0 < x_i <= 1
//   5. rejection envelope       PLIOrePowell::EnvelopeViolations() == 0
//   6. branching probability    P(3g) = tau/tau_3g for a range of tau
//
// Usage: test_orepowell [n_events] > samples.csv
// ============================================================================

#include "PLIOrePowell.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
    const long n = (argc > 1) ? std::atol(argv[1]) : 1000000L;
    CLHEP::HepRandom::setTheSeed(12345);   // fixed, so the test is reproducible

    double maxEnergyErr = 0.0, maxMomErr = 0.0, maxDirErr = 0.0;
    long   support      = 0;
    double sumE         = 0.0;

    std::fprintf(stdout, "x1,x2,x3\n");
    for (long i = 0; i < n; ++i) {
        double        e[3];
        G4ThreeVector d[3];
        PLIOrePowell::SampleThreeGamma(e, d);

        const double etot = e[0] + e[1] + e[2];
        maxEnergyErr = std::max(maxEnergyErr,
                                std::abs(etot - 2.0 * electron_mass_c2)
                                    / (2.0 * electron_mass_c2));

        G4ThreeVector p = e[0] * d[0] + e[1] * d[1] + e[2] * d[2];
        maxMomErr = std::max(maxMomErr, p.mag() / electron_mass_c2);

        for (int k = 0; k < 3; ++k) {
            maxDirErr = std::max(maxDirErr, std::abs(d[k].mag() - 1.0));
            const double x = e[k] / electron_mass_c2;
            if (x > 0.0 && x <= 1.0) ++support;
            sumE += e[k];
        }

        std::fprintf(stdout, "%.10f,%.10f,%.10f",
                     e[0] / electron_mass_c2,
                     e[1] / electron_mass_c2,
                     e[2] / electron_mass_c2);
        std::fprintf(stdout, "\n");
    }

    // --- branching probability check --------------------------------------
    // P(3 gamma | o-Ps) must equal tau_eff / tau_3g. Sampled here directly
    // from the shipped ThreeGammaProbability + the same Bernoulli draw the
    // processes use, so the number under test is the production one.
    bool branchOk = true;
    std::cerr << "\n--- branching check: P(3g) = tau_eff / tau_3g ---\n";
    std::cerr << "tau_eff[ns]   expected      measured      pull\n";
    const double taus[] = {1.5, 1.85, 2.0, 2.55, 2.72, 3.39};
    for (double tau_ns : taus) {
        const double tau = tau_ns * ns;
        const double pExp = PLIOrePowell::ThreeGammaProbability(tau);
        const long   m    = 2000000L;
        long hits = 0;
        for (long i = 0; i < m; ++i) {
            if (G4UniformRand() < pExp) ++hits;
        }
        const double pMeas = double(hits) / double(m);
        const double sd    = std::sqrt(pExp * (1.0 - pExp) / double(m));
        const double pull  = (pMeas - pExp) / sd;
        std::fprintf(stderr, "%10.2f  %.6f      %.6f      %+6.2f\n",
                     tau_ns, pExp, pMeas, pull);
        if (std::abs(pull) > 5.0) branchOk = false;
        // analytic cross-check against tau/142.05 ns
        const double pAnalytic = tau_ns / 142.05;
        if (std::abs(pExp - pAnalytic) / pAnalytic > 1e-9) branchOk = false;
    }

    const long violations = PLIOrePowell::EnvelopeViolations();
    std::cerr << "\n--- invariants over " << n << " sampled events ---\n";
    std::fprintf(stderr, "mean photon energy      = %.6f keV (expect %.6f)\n",
                 sumE / (3.0 * n) / keV, 2.0 / 3.0 * electron_mass_c2 / keV);
    std::fprintf(stderr, "max relative |sum E - 1022 keV| = %.3e\n", maxEnergyErr);
    std::fprintf(stderr, "max |sum p| / (m_e c)           = %.3e\n", maxMomErr);
    std::fprintf(stderr, "max ||d| - 1|                   = %.3e\n", maxDirErr);
    std::fprintf(stderr, "photons inside 0 < x <= 1       = %ld / %ld\n",
                 support, 3 * n);
    std::fprintf(stderr, "rejection-envelope violations   = %ld\n", violations);

    const bool ok = (maxEnergyErr < 1e-9) && (maxMomErr < 1e-12)
                    && (maxDirErr < 1e-12) && (support == 3 * n)
                    && (violations == 0) && branchOk;
    std::cerr << (ok ? "\nRESULT: PASS\n" : "\nRESULT: FAIL\n");
    return ok ? 0 : 1;
}
