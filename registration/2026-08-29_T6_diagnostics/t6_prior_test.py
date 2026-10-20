#!/usr/bin/env python3
"""Mechanism test: is T6's -6.26 % the Dirichlet prior pulling o-Ps intensity down?

Registered prediction BEFORE running (reviewer 84):
  P-T6-1  Removing the prior shrinks |dev| on T6 by more than half.
  P-T6-2  The shrink is LARGEST on T6 of all eight tubes -- i.e. the effect is
          specific to it, not a uniform shift.
  P-T6-3  Fitted I_oPs rises above its prior-constrained value on T6.
If P-T6-1 fails, the prior is NOT the mechanism and the finding stays open.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
REPO = Path("/home/dhritman/research/pli-quadra-simulation")
RES = REPO / "reconstruction/results"
sys.path.insert(0, str(RES / "2026-08-26_A2_cond2"))
import fit3_sim1 as F
E = F.E
z = np.load(RES / "2026-08-26_B2_paired/hists_sim1_388.npz", allow_pickle=True)
names = [str(s) for s in z["voi_names"]]; cfg = z["tau_cfg"]

def nll_noprior(th, counts):
    if not (0.3 < th[0] < 8.0 and 0.01 < th[2] < 1.0 and abs(th[1]) < 2.0): return 1e12
    w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
    p = np.zeros(len(F.EDGES) - 1)
    for Ik, tk in zip(I, (E.TAU_PPS, E.TAU_DIR, th[0])):
        p += Ik * np.diff(E.emg_cdf(F.EDGES, tk, th[1], th[2]))
    p = np.clip(p * F.ACCB, 1e-15, None); p /= p.sum()
    mu = counts.sum() * p
    return float(np.sum(mu - counts * np.log(mu)))

def fit(counts, obj):
    best = None
    for s in F.STARTS:
        st = np.array([s, 0.0, .0909, np.log(.15/.23), np.log(.62/.23)])
        r = minimize(obj, st, args=(counts,), method="Nelder-Mead",
                     options=dict(maxiter=30000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    th = best.x; w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
    return float(th[0]), I

print(f"{'tube':<28} {'cfg':>6} {'dev PRIOR':>10} {'dev NOPRIOR':>12} {'shrink':>8} "
      f"{'I_oPs pri':>10} {'I_oPs no':>9}")
out=[]
for ki, name in enumerate(names):
    if not name.startswith("T"): continue
    c = F.spectrum(z["all_sig"][:, ki, :].sum(0), z["all_rnd"][:, ki, :].sum(0))
    t_p, I_p = fit(c, F._nll)
    t_n, I_n = fit(c, nll_noprior)
    k = float(cfg[ki])
    dp = 100*(t_p-k)/k; dn = 100*(t_n-k)/k
    shrink = (abs(dp)-abs(dn))/abs(dp)*100 if dp else float('nan')
    out.append(dict(tube=name, dev_prior=dp, dev_noprior=dn, shrink_pct=shrink,
                    I_ops_prior=float(I_p[2]), I_ops_noprior=float(I_n[2])))
    print(f"{name:<28} {k:6.3f} {dp:10.2f} {dn:12.2f} {shrink:7.1f}% "
          f"{I_p[2]:10.3f} {I_n[2]:9.3f}")

t6 = next(r for r in out if "brain_healthy" in r["tube"])
others = [r for r in out if r is not t6]
print(f"\n  P-T6-1  |dev| shrinks by >50 % on T6 : "
      f"{t6['shrink_pct']:.1f} %  -> {'PASS' if t6['shrink_pct']>50 else 'FAIL'}")
print(f"  P-T6-2  shrink is largest on T6      : T6 {t6['shrink_pct']:.1f} % vs "
      f"max other {max(r['shrink_pct'] for r in others):.1f} %  -> "
      f"{'PASS' if t6['shrink_pct']>max(r['shrink_pct'] for r in others) else 'FAIL'}")
print(f"  P-T6-3  I_oPs rises without the prior: {t6['I_ops_prior']:.3f} -> "
      f"{t6['I_ops_noprior']:.3f}  -> "
      f"{'PASS' if t6['I_ops_noprior']>t6['I_ops_prior'] else 'FAIL'}")
(Path(__file__).parent/"t6_prior_test.json").write_text(json.dumps(out, indent=2))
