// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
#ifndef PLICONSTANTS_HH
#define PLICONSTANTS_HH
// ---------------------------------------------------------------------------
// Single-constant header. FOR FUTURE WORK -- deliberately NOT yet wired in.
//
// The vacuum ortho-positronium lifetime currently has TWO values in this
// package:
//
//   PLIOrePowell.cc:59        142.05 ns   three-photon branching
//   PLIMaterialLoader.hh:52   142.0  ns   lifetime model
//   analysis/tao_eldrup.py    142.0  ns   lifetime model
//   pli_materials*.json       1/142.0     lambda_3g_ns constant
//
// The accepted value is 142.05 ns:
//   Vallery, Zitzewitz & Gidley, Phys. Rev. Lett. 90, 203402 (2003).
//
// Quantified impact of the discrepancy, so nobody has to re-derive it:
//   tau_eff (water, pO2 = 0)   1.921496643 ns at 142.0
//                              1.921505795 ns at 142.05
//                              difference 0.0092 ps
//   P(3 gamma) = tau/tau_3g    relative difference 3.5e-4
//
// For scale: the demonstrated per-voxel lifetime uncertainty on the modelled
// scanner is about 530 ps, so the discrepancy is roughly 1 part in 58,000 of
// it.
//
// WHY IT IS NOT SIMPLY FIXED. Every validation number in docs/VALIDATION.md
// was produced by ONE binary that carried both values. Changing either one
// regenerates the shipped material databases and breaks the claim that the
// numbers and the binary correspond. That provenance is worth more than
// 0.0092 ps.
//
// The intended migration, when a rebuild is happening anyway for other
// reasons: include this header in PLIOrePowell.cc, PLIMaterialLoader.hh and
// the database generator, regenerate the databases, and re-run VALIDATION.md
// end to end. Do not do it piecemeal -- a half-migration would leave the
// package with two values AND a header claiming there is one.
// ---------------------------------------------------------------------------
namespace pli {

/// Vacuum o-Ps three-photon lifetime, ns. Vallery et al., PRL 90, 203402.
constexpr double kTau3GammaNs = 142.05;

/// Vacuum o-Ps three-photon rate, ns^-1.
constexpr double kLambda3GammaNs = 1.0 / kTau3GammaNs;

}  // namespace pli
#endif  // PLICONSTANTS_HH
