// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// PLIMaterialLoader.cc — JSON-driven per-material PLI parameter loader
// ============================================================================
//
// See PLIMaterialLoader.hh for the design rationale and the contract.
//
// This file is the only place in custom_physics/ that includes
// <nlohmann/json.hpp>. JSON parsing is contained here; the rest of the
// library only sees plain G4 types.
//
// ============================================================================

#include "PLIMaterialLoader.hh"
#include "G4PositroniumAnnihilation.hh"
#include "G4PsAtRestProcess.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace PLIMaterialLoader {

namespace {

constexpr int kSupportedSchemaVersion = 1;

// Candidates searched (in order) when PLI_MATERIALS_JSON is not set.
// Resolved relative to the current working directory of the running
// process — keep this in sync with the layout under
// ~/research/pli-quadra-simulation/.
const std::vector<G4String> kCandidates = {
    "materials/pli_materials.json",
    "../materials/pli_materials.json",
    "../../materials/pli_materials.json",
    "../../../materials/pli_materials.json",
};

std::vector<G4String> g_lastSearched;

bool fileExists(const G4String& path) {
    if (path.empty()) return false;
    std::ifstream f(path);
    return f.good();
}

} // namespace

// ----------------------------------------------------------------------------
// Path resolution
// ----------------------------------------------------------------------------

const std::vector<G4String>& lastSearched() {
    return g_lastSearched;
}

G4String resolveJsonPath() {
    g_lastSearched.clear();

    if (const char* env = std::getenv("PLI_MATERIALS_JSON")) {
        G4String envPath(env);
        g_lastSearched.push_back(G4String("[env PLI_MATERIALS_JSON] ") + envPath);
        if (fileExists(envPath)) {
            return envPath;
        }
    }

    for (const auto& cand : kCandidates) {
        g_lastSearched.push_back(cand);
        if (fileExists(cand)) {
            return cand;
        }
    }
    return G4String("");
}

// ----------------------------------------------------------------------------
// JSON parsing
// ----------------------------------------------------------------------------

LoadResult loadFromJson(const G4String& jsonPath) {
    using nlohmann::json;

    if (jsonPath.empty()) {
        std::string msg =
            "PLIMaterialLoader: no materials JSON found. Searched:\n";
        for (const auto& p : g_lastSearched) {
            msg += "  - " + std::string(p) + "\n";
        }
        msg += "Set the PLI_MATERIALS_JSON environment variable to an "
               "absolute path, or run from a directory where the "
               "materials/ folder is reachable.";
        throw std::runtime_error(msg);
    }

    std::ifstream in(jsonPath);
    if (!in) {
        throw std::runtime_error(
            std::string("PLIMaterialLoader: cannot open ") + std::string(jsonPath));
    }

    json doc;
    try {
        in >> doc;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("PLIMaterialLoader: JSON parse error in ")
            + std::string(jsonPath) + ": " + e.what());
    }

    LoadResult out;
    out.json_path      = jsonPath;
    out.schema_version = doc.value("_schema_version", 0);
    if (out.schema_version != kSupportedSchemaVersion) {
        throw std::runtime_error(
            std::string("PLIMaterialLoader: unsupported _schema_version ")
            + std::to_string(out.schema_version)
            + " in " + std::string(jsonPath)
            + "; expected " + std::to_string(kSupportedSchemaVersion)
            + ". Bump kSupportedSchemaVersion in PLIMaterialLoader.cc "
              "if a new schema is intentional.");
    }
    out.generated_at   = doc.value("_generated_at",   std::string("(none)"));
    out.generated_from = doc.value("_generated_from", std::string("(none)"));

    // constants block
    if (!doc.contains("constants")) {
        throw std::runtime_error(
            "PLIMaterialLoader: missing 'constants' block in " + std::string(jsonPath));
    }
    const auto& cj = doc.at("constants");
    out.constants.lambda_S_ns  = cj.at("lambda_S_ns").get<double>();
    out.constants.delta_R_nm   = cj.at("delta_R_nm").get<double>();
    out.constants.lambda_3g_ns = cj.at("lambda_3g_ns").get<double>();
    out.constants.k_O2_ns_mmHg = cj.at("k_O2_ns_mmHg").get<double>();

    // materials array
    if (!doc.contains("materials") || !doc["materials"].is_array()) {
        throw std::runtime_error(
            "PLIMaterialLoader: missing or non-array 'materials' in "
            + std::string(jsonPath));
    }

    for (const auto& m : doc.at("materials")) {
        MaterialEntry e;
        e.name = m.at("name").get<std::string>();

        if (m.contains("aliases") && m["aliases"].is_array()) {
            for (const auto& a : m["aliases"]) {
                if (a.is_string()) {
                    e.aliases.push_back(a.get<std::string>());
                }
            }
        }

        e.role           = m.value("role",          std::string(""));
        e.R_nm           = m.value("R_nm",          0.0);
        e.pO2_mmHg       = m.value("pO2_mmHg",      0.0);
        e.I3_raw         = m.value("I3",            0.0);
        e.tau_eff_ns     = m.value("tau_TE_ns",     0.0);
        e.tau_pub_ns     = m.value("tau_pub_ns",    0.0);
        e.tau_source     = m.value("tau_source",    std::string(""));
        e.I3_source      = m.value("I3_source",     std::string(""));
        e.I3_flag        = m.value("I3_flag",       std::string(""));
        // I3_convention is REQUIRED. It used to default to "total_formation"
        // when absent, which silently gave 17 of the then-30 database entries a
        // convention nobody had asserted -- and at least three of those were
        // demonstrably wrong, understating their o-Ps yield by 25 % (see the
        // convention audit). An undeclared convention is now a
        // hard error, matching the existing behaviour for an unrecognised one.
        if (!m.contains("I3_convention")) {
            throw std::runtime_error(
                std::string("PLIMaterialLoader: material '") + std::string(e.name)
                + "' does not declare I3_convention. This field is REQUIRED: the "
                  "o-Ps yield differs by a factor 4/3 between the two conventions, "
                  "so it must never be inferred from a default. Add "
                  "\"I3_convention\": \"ops_intensity\" (the source reports the o-Ps "
                  "intensity I3 directly, as PALS 3-component decompositions do) or "
                  "\"total_formation\" (the source reports P(form any positronium)).");
        }
        e.I3_convention  = m.value("I3_convention", std::string("total_formation"));

        // Apply convention conversion. Our G4PositroniumAnnihilation /
        // G4PsAtRestProcess use I3 as P(form any positronium). The 3:1
        // ortho:para spin split is applied AFTER the formation draw, so
        //   P(o-Ps) = I3_internal * 0.75
        //   P(p-Ps) = I3_internal * 0.25
        //   P(direct) = 1 - I3_internal
        // PALS literature commonly reports the o-Ps intensity directly:
        //   I3_PALS = P(o-Ps)
        // Convention: ops_intensity means P(o-Ps) = I3 exactly. The 3:1
        // ortho:para split is applied after the formation draw, so the
        // stored I3 is scaled by 4/3 at load to give total Ps formation.
        //   "ops_intensity -> P(o-Ps) = I3 directly, 3:1 split independent of I3"
        // we set I3_internal = I3_raw * 4/3 for ops_intensity entries, which
        // makes P(o-Ps) = (4/3 I3) * 0.75 = I3 exactly. P(p-Ps) becomes I3/3.
        if (e.I3_convention == "total_formation") {
            e.I3 = e.I3_raw;
        } else if (e.I3_convention == "ops_intensity") {
            e.I3 = e.I3_raw * (4.0 / 3.0);
            if (e.I3 > 1.0) {
                throw std::runtime_error(
                    std::string("PLIMaterialLoader: material '") + std::string(e.name)
                    + "' has I3 (ops_intensity) = " + std::to_string(e.I3_raw)
                    + " which exceeds the 0.75 ceiling. After 4/3 conversion the "
                      "implied total formation probability would be > 1, which is "
                      "unphysical. Either lower the JSON value or recheck the "
                      "convention assignment.");
            }
        } else {
            throw std::runtime_error(
                std::string("PLIMaterialLoader: material '") + std::string(e.name)
                + "' has unrecognised I3_convention = '"
                + std::string(e.I3_convention)
                + "'. Valid values: 'total_formation', 'ops_intensity'.");
        }

        if (e.tau_eff_ns <= 0.0) {
            throw std::runtime_error(
                std::string("PLIMaterialLoader: material '") + std::string(e.name)
                + "' has non-positive tau_TE_ns in JSON; aborting.");
        }
        if (e.I3 < 0.0 || e.I3 > 1.0) {
            throw std::runtime_error(
                std::string("PLIMaterialLoader: material '") + std::string(e.name)
                + "' has out-of-range I3 = " + std::to_string(e.I3));
        }

        out.materials.push_back(std::move(e));
    }

    return out;
}

// ----------------------------------------------------------------------------
// Process population
// ----------------------------------------------------------------------------

LoadResult loadInto(G4PositroniumAnnihilation* psPost,
                    G4PsAtRestProcess*         psRest)
{
    G4String path = resolveJsonPath();
    LoadResult result = loadFromJson(path);

    G4cout << "[PLI] PLIMaterialLoader loaded " << result.materials.size()
           << " materials from " << result.json_path << G4endl;
    G4cout << "[PLI]   schema v" << result.schema_version
           << ", generated_at=" << result.generated_at
           << ", generated_from=" << result.generated_from << G4endl;
    G4cout << "[PLI]   constants: lambda_S=" << result.constants.lambda_S_ns
           << " ns^-1, delta_R=" << result.constants.delta_R_nm
           << " nm, lambda_3g=" << result.constants.lambda_3g_ns
           << " ns^-1, k_O2=" << result.constants.k_O2_ns_mmHg
           << " ns^-1 mmHg^-1" << G4endl;

    auto applyOne = [&](const G4String& key,
                        G4double tau_ns_value, G4double i3_value)
    {
        const G4double tau_g4 = tau_ns_value * ns;
        if (psPost) {
            psPost->SetMaterialLifetime(key, tau_g4);
            psPost->SetPsFormationProbability(key, i3_value);
        }
        if (psRest) {
            psRest->SetMaterialLifetime(key, tau_g4);
            psRest->SetPsFormationProbability(key, i3_value);
        }
    };

    for (const auto& m : result.materials) {
        applyOne(m.name, m.tau_eff_ns, m.I3);
        for (const auto& alias : m.aliases) {
            applyOne(alias, m.tau_eff_ns, m.I3);
        }
        G4cout << "[PLI]   " << m.name
               << " [" << m.role << "]"
               << "  tau_eff=" << m.tau_eff_ns << " ns,"
               << "  I3=" << m.I3
               << " (raw " << m.I3_raw << " " << m.I3_convention << ")"
               << "  (R=" << m.R_nm << " nm,"
               << " pO2=" << m.pO2_mmHg << " mmHg)";
        if (!m.aliases.empty()) {
            G4cout << "  +" << m.aliases.size() << " alias(es)";
        }
        G4cout << G4endl;
    }

    return result;
}

} // namespace PLIMaterialLoader
