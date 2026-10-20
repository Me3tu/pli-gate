// SPDX-License-Identifier: LGPL-3.0-only
// Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
//******************************************************************************
// PLIOrePowell.cc: see PLIOrePowell.hh for the physics and the references.
//******************************************************************************

#include "PLIOrePowell.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Exception.hh"
#include "Randomize.hh"
#include "G4ios.hh"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <string>

namespace {

std::atomic<long> gEnvelopeViolations{0};

G4double EnvDouble(const char* name, G4double fallback)
{
    const char* s = std::getenv(name);
    if (!s) return fallback;
    try {
        return std::stod(std::string(s));
    } catch (...) {
        G4cout << "[PLI] WARNING: could not parse " << name << "='" << s
               << "'; using " << fallback << G4endl;
        return fallback;
    }
}


// unknown-material fallback gate. Read once, like every other
// PLI runtime switch. Absent or "0" -> fatal, "1" -> the old permissive path.
G4bool AllowFallbackCached()
{
    static int cached = -1;
    if (cached < 0) {
        const char* v = std::getenv("PLI_ALLOW_MATERIAL_FALLBACK");
        cached = (v && std::atoi(v) != 0) ? 1 : 0;
        if (cached) {
            G4cout << "[PLI] PLI_ALLOW_MATERIAL_FALLBACK=1: a material absent "
                   << "from the JSON will WARN and use tau = 2.0 ns and the "
                   << "scalar I3 instead of aborting the run." << G4endl;
        }
    }
    return cached != 0;
}

}  // namespace

G4double PLIOrePowell::VacuumTau3Gamma()
{
    // Vallery, Zitzewitz & Gidley, Phys. Rev. Lett. 90, 203402 (2003).
    static const G4double tau = EnvDouble("PLI_TAU_3GAMMA_NS", 142.05) * ns;
    return tau;
}

G4double PLIOrePowell::ThreeGammaProbability(G4double tauEff)
{
    if (tauEff <= 0.0) return 0.0;
    const G4double p = tauEff / VacuumTau3Gamma();
    if (p > 1.0) {
        static G4bool warned = false;
        if (!warned) {
            G4ExceptionDescription msg;
            msg << "tau_eff = " << tauEff / ns << " ns exceeds the vacuum "
                << "three-photon lifetime tau_3g = " << VacuumTau3Gamma() / ns
                << " ns, which is unphysical (o-Ps cannot live longer than "
                << "its vacuum self-annihilation limit). Clamping "
                << "P(3 gamma) to 1. Check the material lifetime table.";
            G4Exception("PLIOrePowell::ThreeGammaProbability",
                        "pli_tau_exceeds_vacuum", JustWarning, msg);
            warned = true;
        }
        return 1.0;
    }
    return p;
}

G4bool PLIOrePowell::ThreeGammaEnabled()
{
    static const G4bool enabled = []() {
        const char* s = std::getenv("PLI_DISABLE_3GAMMA");
        const G4bool off = (s && std::string(s) == "1");
        if (off) {
            G4cout << "[PLI] three-photon o-Ps channel DISABLED "
                   << "(PLI_DISABLE_3GAMMA=1) — two-gamma-only reference behaviour"
                   << G4endl;
        }
        return !off;
    }();
    return enabled;
}

PLILifetimeMode PLIOrePowell::LifetimeMode()
{
    static const PLILifetimeMode mode = []() {
        const char* s = std::getenv("PLI_LIFETIME_MODE");
        const std::string v = s ? std::string(s) : std::string("transport");
        PLILifetimeMode m = PLILifetimeMode::kTransport;
        if (v == "source") m = PLILifetimeMode::kSource;
        else if (v == "fixed") m = PLILifetimeMode::kFixed;
        else if (v != "transport" && s) {
            G4cout << "[PLI] WARNING: unknown PLI_LIFETIME_MODE='" << v
                   << "'; using 'transport'" << G4endl;
        }
        if (m != PLILifetimeMode::kTransport) {
            G4cout << "[PLI] lifetime mode = " << v
                   << " (alternative lifetime-source model; NOT the default physics)"
                   << G4endl;
        }
        return m;
    }();
    return mode;
}

G4double PLIOrePowell::FixedTau()
{
    static const G4double tau = EnvDouble("PLI_FIXED_TAU_NS", 2.0) * ns;
    return tau;
}

G4bool PLIOrePowell::LogCrossMaterial()
{
    static const G4bool on = []() {
        const char* s = std::getenv("PLI_LOG_XMAT");
        return (s && std::string(s) == "1");
    }();
    return on;
}

G4bool PLIOrePowell::AllowMaterialFallback()
{
    return AllowFallbackCached();
}

G4double PLIOrePowell::DalitzDensity(G4double x1, G4double x2, G4double x3)
{
    if (x1 <= 0.0 || x2 <= 0.0 || x3 <= 0.0 ||
        x1 > 1.0  || x2 > 1.0  || x3 > 1.0) {
        return 0.0;   // outside the physical (triangle-inequality) region
    }
    const G4double a = (1.0 - x1) / (x2 * x3);
    const G4double b = (1.0 - x2) / (x1 * x3);
    const G4double c = (1.0 - x3) / (x1 * x2);
    return a * a + b * b + c * c;
}

void PLIOrePowell::SampleFractions(G4double& x1, G4double& x2, G4double& x3)
{
    // sup(DalitzDensity) on the physical region is exactly 2. Proof: with
    // x1+x2+x3 = 2 the Ore-Powell density is
    //     f = sum_i ((1-x_i)/(x_j x_k))^2 ... normalised so f -> 2 as any
    // x_i -> 1, where the corresponding term vanishes and the other two each
    // tend to 1. f is continuous on the closed physical triangle and attains
    // its supremum on the boundary x_i = 1, so sup f = 2 exactly. Rejection
    // with this envelope is therefore exact, not approximate. Acceptance ~33 %.
    // tests/test_orepowell.cc asserts zero envelope violations over 1e6 draws.
    constexpr G4double kEnvelope = 2.0;

    for (G4int trial = 0; trial < 1000; ++trial) {
        G4double u = G4UniformRand();
        G4double v = G4UniformRand();
        // Fold the unit square onto the Dalitz triangle {u+v >= 1}, uniformly.
        if (u + v < 1.0) { u = 1.0 - u; v = 1.0 - v; }
        const G4double w = 2.0 - u - v;

        const G4double f = DalitzDensity(u, v, w);
        if (f <= 0.0) continue;
        if (f > kEnvelope) gEnvelopeViolations.fetch_add(1);

        if (G4UniformRand() * kEnvelope < f) {
            x1 = u; x2 = v; x3 = w;
            return;
        }
    }
    // Unreachable in practice: P(1000 consecutive rejections) ~ 1e-176.
    x1 = x2 = x3 = 2.0 / 3.0;
}

void PLIOrePowell::SampleThreeGamma(G4double energy[3], G4ThreeVector dir[3])
{
    G4double x1, x2, x3;
    SampleFractions(x1, x2, x3);

    energy[0] = x1 * electron_mass_c2;
    energy[1] = x2 * electron_mass_c2;
    energy[2] = x3 * electron_mass_c2;

    // Photon 1 isotropic. The decay plane then gets a uniform azimuth about it.
    dir[0] = G4RandomDirection();
    const G4ThreeVector ref = (std::abs(dir[0].x()) < 0.9)
                                  ? G4ThreeVector(1.0, 0.0, 0.0)
                                  : G4ThreeVector(0.0, 1.0, 0.0);
    const G4ThreeVector e1 = dir[0].cross(ref).unit();
    const G4ThreeVector e2 = dir[0].cross(e1).unit();
    const G4double phi = twopi * G4UniformRand();
    const G4ThreeVector perp = std::cos(phi) * e1 + std::sin(phi) * e2;

    // Momentum closure fixes the opening angle between photons 1 and 2:
    //   x3^2 = |x1 d1 + x2 d2|^2 = x1^2 + x2^2 + 2 x1 x2 cos(theta_12)
    G4double cos12 = (x3 * x3 - x1 * x1 - x2 * x2) / (2.0 * x1 * x2);
    if (cos12 >  1.0) cos12 =  1.0;
    if (cos12 < -1.0) cos12 = -1.0;
    const G4double sin12 = std::sqrt(std::max(0.0, 1.0 - cos12 * cos12));

    dir[1] = (cos12 * dir[0] + sin12 * perp).unit();
    dir[2] = -(x1 * dir[0] + x2 * dir[1]) / x3;
    dir[2] = dir[2].unit();
}

long PLIOrePowell::EnvelopeViolations()
{
    return gEnvelopeViolations.load();
}

const char* PLIOrePowell::ChannelName(PLIChannel c)
{
    switch (c) {
        case PLIChannel::kDirect:    return "direct2g";
        case PLIChannel::kParaPs:    return "pPs2g";
        case PLIChannel::kOrthoPs2g: return "oPs2g";
        case PLIChannel::kOrthoPs3g: return "oPs3g";
    }
    return "unknown";
}
