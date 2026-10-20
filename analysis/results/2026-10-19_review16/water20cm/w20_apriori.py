#!/usr/bin/env python3
"""A priori acceptance readout of the 20 cm water cylinder (Table 5, Sec. 3.1), written for release v1.1.0.

The published w20_apriori.json had no script on disk (inventory flag F2). This script reproduces it.
(1) sigma_geo, the a priori path-length spread of the box variable, with the method and the frozen detector model of
    2026-08-23_acceptance/sigma_geo.py (hit_distance, C) and 2026-08-26_A2_cond2/sigma_geo_sim1.py: source uniform in a
    water cylinder, isotropic pair axis and independent isotropic prompt, all three photons inside the axial
    acceptance, 2e6 decays. One generator, seed 20261019, used in turn for the 20 mm cylinder (R 10 mm, half length
    150 mm), the 20 cm cylinder (R 100 mm, half length 150 mm) and the geometry of sigma_geo.py (R 25 mm, 150 mm).
    sigma_e = hypot(sigma_geo, sigma_t / sqrt 2) with sigma_t = 90.9 ps (64.3 ps pair-timing term).
(2) The windowed-moment readout of VOI_15x25 inverted with the acceptance of width sigma_e (20 cm), next to the
    sharp-edge production readout, with the builder, inversion and shard bootstrap (4 x 250 resamples, seed 11) of
    w20_extract.py (its definitions are executed unchanged from that file).
    Without arguments: the all-prompts per-shard histograms are those of the production all-prompts subset
    (|s| <= 5 ns, all truth classes) stored in 2026-10-20_review18/w20_edge_hists.npz, identical to the w20_extract
    histograms in the (1, 5) ns fit window. The closest-prompt histograms were not kept, so the closest-prompt point
    value is obtained from the windowed mean recovered from the stored sharp readout (w20_extract.json) and its
    bootstrap uncertainty is not computed (null).
    With the eight batch files as arguments (raw singles, see analysis/README.md), both rules are rebuilt from the
    singles exactly as w20_extract.py does.
Usage: w20_apriori.py [water20_batch_200.h5 ... water20_batch_207.h5]
Writes w20_apriori.json next to this script."""
import json, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(RES / "2026-08-23_acceptance"))
from sigma_geo import hit_distance, C  # noqa: E402  (the frozen detector model, and the module prints its own table)

SIG_T_PAIR = 0.0909 / np.sqrt(2.0)
TAU_CFG, TAU_STEINBERGER = 1.84, 1.815


def sigma_geo(rng, R, hz, n=2_000_000):
    r = R * np.sqrt(rng.random(n)); th = rng.uniform(0, 2 * np.pi, n)
    p = np.c_[r * np.cos(th), r * np.sin(th), rng.uniform(-hz, hz, n)]

    def iso(m):
        c = rng.uniform(-1, 1, m); s = np.sqrt(1 - c * c)
        ph = rng.uniform(0, 2 * np.pi, m)
        return np.c_[s * np.cos(ph), s * np.sin(ph), c]
    u = iso(n)
    d1, ok1 = hit_distance(p, u)
    d2, ok2 = hit_distance(p, -u)
    v = iso(n)
    dp, okp = hit_distance(p, v)
    ok = ok1 & ok2 & okp
    return float(((dp[ok] - 0.5 * (d1[ok] + d2[ok])) / C).std(ddof=1))


files = sys.argv[1:]
for fn in files:
    if not Path(fn).exists():
        raise SystemExit(f"{fn}: raw data not included in the release, see analysis/README.md")

out = {}
rng = np.random.default_rng(20261019)
for key, R, hz in (("cyl20mm", 10.0, 150.0), ("cyl20cm", 100.0, 150.0), ("sigma_geo_py_water", 25.0, 150.0)):
    s = sigma_geo(rng, R, hz)
    out[key] = dict(sigma_geo_ps=1000 * s, sigma_e_ps=1000 * float(np.hypot(s, SIG_T_PAIR)))
    print(f"{key:20s} sigma_geo {1000 * s:6.1f} ps  sigma_e {out[key]['sigma_e_ps']:6.2f} ps")
SE = out["cyl20cm"]["sigma_e_ps"] / 1000

W = {"__file__": str(HERE / "w20_extract.py"), "__name__": "w20_extract_defs"}
exec(compile((HERE / "w20_extract.py").read_text().split("files = sys.argv[1:]")[0], "w20_extract_defs", "exec"), W)
J = json.loads((HERE / "w20_extract.json").read_text())
IV = 1 + [v["name"] for v in W["VOIS"]].index("VOI_15x25")
curves = {"sharp": W["curve"](), "apriori": W["curve"](SE)}
for rule in ("all", "closest"):
    HS = None
    if files:
        HS = W["run"](files, rule)[0]
    elif rule == "all":
        Z = np.load(RES / "2026-10-20_review18/w20_edge_hists.npz")
        E0 = int(np.flatnonzero(np.isclose(Z["edges"], W["X"].TAU_BINS[0]))[0]); NB = len(W["X"].TAU_BINS) - 1
        assert np.allclose(Z["edges"][E0:E0 + NB + 1], W["X"].TAU_BINS)
        HS = Z["H"][:, :, :, 1].sum(2)[..., E0:E0 + NB].transpose(0, 2, 1, 3).astype(float)   # [shard, region, arm, bin]
    sharp_pub = J[rule]["vois"]["VOI_15x25"]["sharp"][0]
    if HS is not None:
        T0 = W["invert"](HS.sum(0), curves)
        assert abs(T0["sharp"][IV] - sharp_pub) < 1e-9, (rule, T0["sharp"][IV], sharp_pub)
        rng = np.random.default_rng(11); S = len(HS); flat = HS.reshape(S, -1); TB = []
        for _ in range(4):
            w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(250)]).astype(float)
            TB.append(W["invert"]((w @ flat).reshape((250,) + HS.shape[1:]), curves))
        sharp, ap = float(T0["sharp"][IV]), float(T0["apriori"][IV])
        err = float(np.nanstd(np.concatenate([tb["apriori"][:, IV] for tb in TB]), ddof=1))
    else:
        mean = np.interp(sharp_pub, W["tg"], curves["sharp"])      # windowed mean on (1, 5) ns
        sharp, ap, err = float(sharp_pub), float(np.interp(mean, curves["apriori"], W["tg"])), None
    out[rule] = dict(sharp=sharp, apriori=ap, apriori_err=err, dev_cfg=100 * (ap - TAU_CFG) / TAU_CFG,
                     dev_meas=100 * (ap - TAU_STEINBERGER) / TAU_STEINBERGER)
    print(f"{rule:8s} VOI_15x25 sharp {sharp:.4f}  a priori ({1000 * SE:.2f} ps) {ap:.4f}"
          + (f" +- {err:.4f}" if err is not None else " (uncertainty needs the raw singles)")
          + f"  vs configured {out[rule]['dev_cfg']:+.2f} %, vs Steinberger {out[rule]['dev_meas']:+.2f} %")
json.dump(out, open(HERE / "w20_apriori.json", "w"), indent=1)
print("wrote", HERE / "w20_apriori.json")
