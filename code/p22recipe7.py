#!/usr/bin/env python3
"""P22 Amendment 7 predictions (P1 frozen J-only recipe; P2 noise-scaled). Uses only exact curves + amend6_recipe_params.csv.
Usage: p22recipe7.py  -> results/amend7_predictions.csv, results/amend7_cells_exact.csv"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results"); sys.path.insert(0, HERE)
import p22lib as L
G3 = [32, 48, 64, 80, 96, 112, 128, 160, 192, 256]; G2 = [128, 160, 192, 224, 256, 320, 384, 512, 768]
FAM = {"aniso3_w8": (G3, lambda l: (l, 8, 8)), "aniso2_w16": (G2, lambda l: (l, 16))}
EST = ["E0", "E1q", "E4c", "E6"]; CF = dict(E0=1.0, E1q=0.5, E4c=0.5, E6=0.5)   # count factor of the window-selection data

def exactJ(q, ts, N, wpn, cf, thr=L.THR, cmin=L.CMIN, t0=24):
    J = 0
    for t, qq in zip(ts, q):
        if t < t0: continue
        if N * qq >= thr and wpn * cf * N * qq >= cmin: J += 1
        else: break
    return J

def minL_idx(det, grid):
    for i in range(len(grid)):
        if all(det[l] for l in grid[i:]): return i
    return len(grid)

def main():
    par = pd.read_csv(RES + "/amend6_recipe_params.csv").set_index("est")
    rows, prows = [], []
    for fam, (grid, dims) in FAM.items():
        for l in grid:
            g = L.Graph("torus", 0.5, dims=dims(l)); ts = L.make_grid(24, 34, 3_000_000)
            q, d = g.exact(ts); N = g.N
            run = [i for i in range(len(ts)) if ts[i] >= 24][:]
            k = [i for i in run]; jj = []
            for i in k:
                if N * q[i] >= L.THR: jj.append(i)
                else: break
            S = float((d[jj[0]] - d[jj]).max())
            for wpn in (8, 32, 128):
                for e in EST:
                    rows.append(dict(fam=fam, L=l, N=N, wpn=wpn, est=e, J=exactJ(q, ts, N, wpn, CF[e]), S=S))
    C = pd.DataFrame(rows); C.to_csv(RES + "/amend7_cells_exact.csv", index=False)
    for (fam, wpn, e), g in C.groupby(["fam", "wpn", "est"]):
        p = par.loc[e]; grid = FAM[fam][0]
        for model in ("P1", "P2"):
            det = {}
            for _, r in g.iterrows():
                hw = p.a * (r.J / 7.0) ** (-p.delta) * ((wpn / 32.0) ** -0.5 if model == "P2" else 1.0)
                det[int(r.L)] = bool(r.J >= p.Jmin and (r.S - 0.25) / hw >= p["T"])
            i = minL_idx(det, grid)
            prows.append(dict(fam=fam, wpn=wpn, est=e, model=model, pred_idx=i, pred_minL=grid[i] if i < len(grid) else f">{grid[-1]}",
                              detail=";".join(f"{l}:J{int(g[g.L==l].J.iloc[0])}{'Y' if v else 'n'}" for l, v in det.items())))
    P = pd.DataFrame(prows); P.to_csv(RES + "/amend7_predictions.csv", index=False)
    print(P.drop(columns="detail").pivot_table(index=["fam", "est"], columns=["model", "wpn"], values="pred_minL", aggfunc="first").to_string())

if __name__ == "__main__": main()
