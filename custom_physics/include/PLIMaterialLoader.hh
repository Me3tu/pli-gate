// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
// ============================================================================
// PLIMaterialLoader.hh: JSON-driven per-material PLI parameter loader
// ============================================================================
//
// Reads materials/pli_materials.json (the file produced by
// the materials JSON) and populates G4PositroniumAnnihilation +
// G4PsAtRestProcess with per-material o-Ps lifetime and I3 (formation
// probability) values.
//
// Design notes
// ------------
// 1. Python is the single source of truth for the Tao-Eldrup + O2
//    quenching math. The JSON carries the *already-computed* tau_TE_ns
//    per material. This C++ loader does NOT re-derive R or tau. That
//    eliminates floating-point divergence between the Python validation
//    and the simulation runtime.
//
// 2. Each material entry in the JSON may carry an "aliases" array.
//    The loader populates the process maps for the canonical name AND
//    every alias, so downstream G4Material::GetName() lookups will hit
//    regardless of which naming convention the GATE macro used.
//
// 3. Path resolution order:
//      a) environment variable PLI_MATERIALS_JSON (absolute path)
//      b) candidate search relative to CWD (materials/, ../materials/, ...)
//    If nothing is found, loadFromJson throws std::runtime_error with a
//    diagnostic listing every path that was tried.
//
// 4. Schema version is enforced. The loader refuses to load a JSON whose
//    _schema_version differs from kSupportedSchemaVersion (currently 1).
//    Bump the constant in PLIMaterialLoader.cc when the schema changes.
//
// ============================================================================

#ifndef PLIMaterialLoader_hh
#define PLIMaterialLoader_hh 1

#include "globals.hh"
#include <vector>

class G4PositroniumAnnihilation;
class G4PsAtRestProcess;

namespace PLIMaterialLoader {

/// Physics constants block (mirrors the JSON `constants` object).
/// Stored for diagnostic printing only. The loader does not perform
/// any TE math itself.
struct PLIConstants {
    G4double lambda_S_ns   = 2.0;
    G4double delta_R_nm    = 0.1656;
    G4double lambda_3g_ns  = 1.0 / 142.0;
    G4double k_O2_ns_mmHg  = 2.45e-5;
};

/// One material entry parsed from the JSON.
struct MaterialEntry {
    G4String              name;          // canonical name (matches G4Material if possible)
    std::vector<G4String> aliases;       // additional GATE/G4 material names that map here
    G4String              role;          // "validation" | "calibration"
    G4double              R_nm        = 0.0;   // diagnostic, not used at runtime
    G4double              pO2_mmHg    = 0.0;
    G4double              I3          = 0.0;       // value passed to G4PositroniumAnnihilation
                                                   // (post-conversion if convention != total_formation)
    G4double              I3_raw      = 0.0;       // value as read from JSON, before convention conversion
    G4String              I3_convention;           // "total_formation" (default) | "ops_intensity"
    G4double              tau_eff_ns  = 0.0;       // = tau_TE_ns from JSON, the value used at runtime
    G4double              tau_pub_ns  = 0.0;       // diagnostic
    G4String              tau_source;
    G4String              I3_source;
    G4String              I3_flag;     // "measured" | "estimated"
};

/// Aggregate result of loading the parameter file.
struct LoadResult {
    G4String                    json_path;
    int                         schema_version = 0;
    G4String                    generated_at;
    G4String                    generated_from;
    PLIConstants                constants;
    std::vector<MaterialEntry>  materials;
};

/// Locate the JSON file. See header notes for resolution order.
/// Returns an empty G4String if no candidate exists. The caller should
/// then inspect lastSearched() and emit a diagnostic.
G4String resolveJsonPath();

/// List of paths that the most recent resolveJsonPath() call attempted.
/// Reset on every call to resolveJsonPath().
const std::vector<G4String>& lastSearched();

/// Parse the JSON at the given path. Throws std::runtime_error on:
///   - empty/missing path
///   - file open failure
///   - JSON parse error
///   - unsupported _schema_version
///   - missing required fields
LoadResult loadFromJson(const G4String& jsonPath);

/// Convenience: resolve the JSON path, parse it, and populate both
/// processes' lifetime + Ps-formation maps. Iterates over each
/// material's canonical name AND every alias. Either pointer may be
/// nullptr (the loader skips it). Returns the parsed LoadResult so the
/// caller can print a summary.
///
/// This is the entry point that PLIPhysicsList::AddCustomPositronAnnihilation
/// calls at physics-list initialisation.
LoadResult loadInto(G4PositroniumAnnihilation* psPost,
                    G4PsAtRestProcess*         psRest);

} // namespace PLIMaterialLoader

#endif // PLIMaterialLoader_hh
