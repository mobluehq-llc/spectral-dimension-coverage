#!/usr/bin/env python3
"""Noise-free check of E7 (Amendment 8): zero-width box at the exact q must return a band containing the exact d(t).
(a) atoms = exact distinct eigenvalues + the 400/799/1597 log grid (feasible by construction): tests the LP/Charnes-Cooper implementation.
(b) atoms = log grid only (400/799/1597), box widened by relative 1e-6 only: reports feasibility and containment of the grid-only arm."""
import sys, numpy as np, pandas as pd, time
import p22lib as L, p22amend8 as A8
import run_p22 as R

def exact_atoms(g):
    lam = g.lam; lam = lam[lam > 1e-14]; s = np.unique(np.round(-np.log(lam), 13)); return s

def band_with(ts, lo, hi, N, idx, s_atoms):
    S = A8.lp_setup(ts, lo, hi, s_atoms, N)
    res = [A8.lp_one(S, j) for j in idx]
    return np.array([r[0] for r in res]), np.array([r[1] for r in res]), max(r[2] for r in res)

cells = [("ring", 256, (256,)), ("torus2", 24, (24, 24)), ("torus3", 16, (16,) * 3), ("aniso2_w16", 128, (128, 16)), ("aniso3_w8", 32, (32, 8, 8)),
         ("aniso2_w16", 256, (256, 16)), ("aniso3_w8", 96, (96, 8, 8))]
rows = []
for fam, l, dims in cells:
    g = L.Graph("torus", 0.5, dims=dims); N = g.N
    big = L.make_grid(24, 34, 3_000_000); q, d = g.exact(big)
    cap = big[np.nonzero(N * q >= 4.5)[0][-1]]; ts = L.make_grid(24, 34, int(cap)); q, d = g.exact(ts)
    J = L.window_len(N * q * 10 ** 3, q, N, 6.0)         # exact window by N q >= 6 only
    idx = list(range(J))
    for tol in (1e-9,):
        lo = q * (1 - tol); hi = q * (1 + tol)
        se = exact_atoms(g); t0 = time.time()
        sa = np.concatenate([[0.0], se, np.exp(np.linspace(np.log(0.02 / ts[-1]), 0, 400)), [np.inf]])
        a = band_with(ts, lo, hi, N, idx, sa); ta = time.time() - t0
        rows.append(dict(cell=f"{fam}_{l}", mode="exact-eig+grid", J=J, n_atoms=len(sa), status=a[2], contains=bool(((d[:J] >= a[0] - 1e-6) & (d[:J] <= a[1] + 1e-6)).all()) if a[2] == 0 else None,
                         max_viol=float(np.maximum(a[0] - d[:J], d[:J] - a[1]).max()) if a[2] == 0 else None, mean_width=float((a[1] - a[0]).mean()) if a[2] == 0 else None, seconds=round(ta, 1)))
        for n in A8.NATOMS:
            for tolg in (1e-9, 1e-6):
                lo = q * (1 - tolg); hi = q * (1 + tolg)
                b = band_with(ts, lo, hi, N, idx, A8.atoms_grid(ts[-1], n))
                rows.append(dict(cell=f"{fam}_{l}", mode=f"grid{n}_tol{tolg:g}", J=J, n_atoms=n + 2, status=b[2], contains=bool(((d[:J] >= b[0] - 1e-6) & (d[:J] <= b[1] + 1e-6)).all()) if b[2] == 0 else None,
                                 max_viol=float(np.maximum(b[0] - d[:J], d[:J] - b[1]).max()) if b[2] == 0 else None, mean_width=float((b[1] - b[0]).mean()) if b[2] == 0 else None, seconds=0))
    print(rows[-7:], flush=True)
df = pd.DataFrame(rows); df.to_csv(R.RES + "/amend8_noisefree.csv", index=False); print(df.to_string())
