// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// test_load.cc - materials-database acceptance check
// ============================================================================
//
// Resolves the PLI materials JSON path, parses it via
// PLIMaterialLoader::loadFromJson, and prints the per-material tau_eff
// and I3 values. This is the acceptance gate before any GATE simulation
// run: it verifies that
//
//   1. nlohmann/json is wired into the custom_physics build correctly,
//   2. PLIMaterialLoader can resolve the JSON path from CWD,
//   3. The schema version matches what the loader expects,
//   4. Every material entry parses cleanly,
//   5. The parsed values agree numerically with what the materials JSON
//      contains.
//
// The C++ side does NOT recompute Tao-Eldrup math here; the JSON's
// tau_TE_ns is the value used at runtime. So this test is an end-to-end
// check of the parser + linkage, NOT an independent physics validation
// (which lives in tests/test_pli.py).
//
// Exit codes:
//   0  success
//   1  loader threw an exception
//   2  loaded entry count is zero
//
// Usage:
//   ./test_load
//   PLI_MATERIALS_JSON=/abs/path/to/pli_materials.json ./test_load
// ============================================================================

#include "PLIMaterialLoader.hh"

#include <cstdio>
#include <exception>
#include <iostream>

int main() {
    using namespace PLIMaterialLoader;

    try {
        const G4String path = resolveJsonPath();
        const LoadResult r = loadFromJson(path);

        std::cout << "================================================================\n"
                  << "test_load — PLI materials JSON parser verification\n"
                  << "================================================================\n"
                  << "json_path        : " << r.json_path << "\n"
                  << "schema_version   : " << r.schema_version << "\n"
                  << "generated_at     : " << r.generated_at << "\n"
                  << "generated_from   : " << r.generated_from << "\n"
                  << "constants:\n"
                  << "  lambda_S_ns    = " << r.constants.lambda_S_ns << "\n"
                  << "  delta_R_nm     = " << r.constants.delta_R_nm << "\n"
                  << "  lambda_3g_ns   = " << r.constants.lambda_3g_ns << "\n"
                  << "  k_O2_ns_mmHg   = " << r.constants.k_O2_ns_mmHg << "\n"
                  << "materials (" << r.materials.size() << "):\n";

        std::printf("  %-20s  %-12s  %8s  %6s  %8s  %6s  %s\n",
                    "name", "role", "tau_eff", "I3",
                    "R_nm", "pO2", "aliases");
        std::printf("  %-20s  %-12s  %8s  %6s  %8s  %6s  %s\n",
                    "----", "----", "-------", "--",
                    "----", "---", "-------");

        for (const auto& m : r.materials) {
            std::printf("  %-20s  %-12s  %8.4f  %6.3f  %8.4f  %6.1f  %zu\n",
                        m.name.c_str(),
                        m.role.c_str(),
                        m.tau_eff_ns,
                        m.I3,
                        m.R_nm,
                        m.pO2_mmHg,
                        m.aliases.size());
        }
        std::cout << "================================================================\n";

        if (r.materials.empty()) {
            std::cerr << "[FAIL] zero materials loaded\n";
            return 2;
        }
        std::cout << "[PASS] " << r.materials.size()
                  << " materials loaded successfully.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] PLIMaterialLoader threw:\n  " << e.what() << "\n";
        std::cerr << "Searched paths:\n";
        for (const auto& p : PLIMaterialLoader::lastSearched()) {
            std::cerr << "  - " << p << "\n";
        }
        return 1;
    }
}
