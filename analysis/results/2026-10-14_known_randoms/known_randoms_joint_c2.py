#!/usr/bin/env python3
"""Review round 14 (referee 2): the known-content test under the realistic timing model with the estimator used there,
the count-level joint likelihood of the three-component model (comp_joint_prod.fit_joint, fit3_sim1 shape, fixed
lifetimes, prior and acceptance unchanged), on the labeled C2 replay (known_randoms_c2_hists.npz).
Per tube, rule and rung (k = 0.15, 1.00, 3.00):
  tau_sub = joint fit of (signal arm, randoms arm)              -- the clinical readout
  tau_tc  = joint fit of (G + PA_corr, no delayed arm)          -- what an exact removal of the uncorrelated content leaves
  tau_G   = joint fit of (G, no delayed arm)                    -- genuine triples only
  tau_kn  = joint fit of (signal arm, AP + PA_unc as the delayed arm) -- the same objective with the known uncorrelated
            content in place of the randoms arm (an independent realization of what the randoms arm estimates)
  dR = 100 (tau_sub - tau_tc) / tau_cfg   (randoms-correction error under the component fit, includes the profiling
       behaviour of the objective, since tau_tc has no delayed arm)
  dK = 100 (tau_sub - tau_kn) / tau_cfg   (randoms arm against the known uncorrelated content, same objective)
  bG = 100 (tau_G - tau_cfg) / tau_cfg    (estimator bias on genuine triples)
Shard bootstrap B = 100, one resample shared by all fits of a tube and, in this version, by the three rungs of that tube
(seed 20261014 + 100 rule + voi), so that differences between rungs are paired; warm starts from the point fit.
(The first run used per-rung seeds; its output is kept as known_randoms_joint_c2_unpaired.json.)
Writes known_randoms_joint_c2.json."""
import os, sys, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
NPZ = str(RES / '2026-10-14_known_randoms/known_randoms_c2_hists.npz')
OUT = RES / '2026-10-14_known_randoms/known_randoms_joint_c2.json'
KSEL = [0.15, 1.00, 3.00]
B = 100


def job(args):
    ik, k, ir, rule, vi, name, tau, npz, b = args
    sys.path.insert(0, str(RES / '2026-09-30_referee2')); sys.path.insert(0, str(RES / '2026-08-26_A2_cond2'))
    import comp_joint_prod as CJ
    F = CJ.F
    Z = np.load(npz)
    H = Z['hist'][:, ik, ir, vi].astype(float)            # [shard, cls, bin]
    CL = list(Z['classes']); iS, iR, iG, iPc, iAP, iPu = (CL.index(c) for c in ('sig', 'rnd', 'G', 'PA_corr', 'AP', 'PA_unc'))
    zero = np.zeros(int(F.M.sum()))

    def fits(h, x0=None):
        nt, nd = h[iS], h[iR]
        ts, xs = CJ.fit_joint(nt[F.M], nd[F.M], x0=None if x0 is None else x0[0])
        tt, xt = CJ.fit_joint((h[iG] + h[iPc])[F.M], zero, x0=None if x0 is None else x0[1])
        tg, xg = CJ.fit_joint(h[iG][F.M], zero, x0=None if x0 is None else x0[2])
        tk, xk = CJ.fit_joint(nt[F.M], (h[iAP] + h[iPu])[F.M], x0=None if x0 is None else x0[3])
        return (ts, tt, tg, tk), (xs, xt, xg, xk)
    (ts, tt, tg, tk), x0 = fits(H.sum(0))
    rng = np.random.default_rng(20261014 + 100 * ir + vi)
    S = H.shape[0]; bd = np.full((b, 4), np.nan)
    for i in range(b):
        w = np.bincount(rng.integers(0, S, S), minlength=S).astype(float)
        (a, c, d, e), _ = fits(np.tensordot(w, H, axes=(0, 0)), x0)
        bd[i] = [100 * (a - c) / tau, 100 * (d - tau) / tau, 100 * (a - tau) / tau, 100 * (a - e) / tau]
    sd = lambda v: float(0.5 * np.subtract(*np.nanpercentile(v, [84, 16])))
    return dict(k=k, rule=rule, tube=name, tau_cfg=tau, tau_sub=ts, tau_tc=tt, tau_G=tg,
                dR=100 * (ts - tt) / tau, dR_sd=sd(bd[:, 0]), bG=100 * (tg - tau) / tau, bG_sd=sd(bd[:, 1]),
                bsub=100 * (ts - tau) / tau, bsub_sd=sd(bd[:, 2]), dK=100 * (ts - tk) / tau, dK_sd=sd(bd[:, 3]), tau_kn=tk,
                boot_dR=bd[:, 0].tolist(), boot_bG=bd[:, 1].tolist(), boot_dK=bd[:, 3].tolist(),
                n_nan=int(np.isnan(bd).any(1).sum()))


if __name__ == '__main__':
    Z = np.load(NPZ); KS = list(Z['ks']); RULES = list(Z['rules']); names = list(Z['voi_names']); tc = Z['tau_cfg']
    jobs = [(KS.index(k), k, ir, r, vi, str(names[vi]), float(tc[vi]), NPZ, B)
            for k in KSEL for ir, r in enumerate(RULES) for vi in range(len(names)) if str(names[vi]).startswith('T')]
    t0 = time.time(); rows = []
    with ProcessPoolExecutor(max_workers=24) as ex:
        for i, r in enumerate(ex.map(job, jobs)):
            rows.append(r)
            print(f"{r['k']:.2f} {r['rule']:7} {r['tube'][:22]:22} dR {r['dR']:+6.2f}+-{r['dR_sd']:.2f}  bG {r['bG']:+6.2f}+-{r['bG_sd']:.2f}  "
                  f"bsub {r['bsub']:+6.2f}  dK {r['dK']:+6.2f}+-{r['dK_sd']:.2f}  nan {r['n_nan']}  {time.time() - t0:.0f} s", flush=True)
    summ = {}
    for k in KSEL:
        for r_ in RULES:
            sel = [r for r in rows if r['k'] == k and r['rule'] == r_]
            dRb = np.nanmean(np.array([r['boot_dR'] for r in sel]), axis=0); bGb = np.nanmean(np.array([r['boot_bG'] for r in sel]), axis=0)
            dKb = np.nanmean(np.array([r['boot_dK'] for r in sel]), axis=0)
            sd = lambda v: float(0.5 * np.subtract(*np.nanpercentile(v, [84, 16])))
            summ[f'{k:.2f}_{r_}'] = dict(dR_mean=float(np.mean([r['dR'] for r in sel])), dR_sd=sd(dRb),
                                        bG_mean=float(np.mean([r['bG'] for r in sel])), bG_sd=sd(bGb),
                                        bsub_mean=float(np.mean([r['bsub'] for r in sel])),
                                        dR_range=[float(min(r['dR'] for r in sel)), float(max(r['dR'] for r in sel))],
                                        max_abs_z=float(max(abs(r['dR']) / r['dR_sd'] for r in sel)),
                                        dK_mean=float(np.mean([r['dK'] for r in sel])), dK_sd=sd(dKb),
                                        dK_range=[float(min(r['dK'] for r in sel)), float(max(r['dK'] for r in sel))],
                                        dK_max_abs_z=float(max(abs(r['dK']) / r['dK_sd'] for r in sel)))
            print(k, r_, {a: (round(b_, 3) if isinstance(b_, float) else b_) for a, b_ in summ[f'{k:.2f}_{r_}'].items()}, flush=True)
    # paired change of the all-prompts and closest-prompt eight-tube means from the lowest to the highest rung
    for r_ in RULES:
        lo = [r for r in rows if r['k'] == KSEL[0] and r['rule'] == r_]; hi = [r for r in rows if r['k'] == KSEL[-1] and r['rule'] == r_]
        lo.sort(key=lambda r: r['tube']); hi.sort(key=lambda r: r['tube'])
        for q in ('dR', 'dK'):
            dl = np.nanmean(np.array([r[f'boot_{q}'] for r in lo]), axis=0); dh = np.nanmean(np.array([r[f'boot_{q}'] for r in hi]), axis=0)
            d0 = float(np.mean([r[q] for r in hi]) - np.mean([r[q] for r in lo])); dd = dh - dl
            summ[f'rise_{q}_{r_}'] = dict(change=d0, sd=float(0.5 * np.subtract(*np.nanpercentile(dd, [84, 16]))))
            print('rise', q, r_, round(d0, 2), '+-', round(summ[f'rise_{q}_{r_}']['sd'], 2), flush=True)
    for r in rows: r.pop('boot_dR'); r.pop('boot_bG'); r.pop('boot_dK')
    json.dump(dict(doc=__doc__, B=B, rows=rows, summary=summ), open(OUT, 'w'), indent=1)
    print(f'done {time.time() - t0:.0f} s')
