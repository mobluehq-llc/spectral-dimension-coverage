"""EXPLORATORY (not pre-registered; added after E3 returned the pinned grid edge 3.6 on every causet): E3 with the d search range widened to 0.5..14."""
import pickle, glob, os, numpy as np, pandas as pd
import p22lib as L, causet_lib as cl
rows = [pickle.load(open(f, "rb")) for f in glob.glob("results/causet/r_*.pkl") if os.path.getsize(f) > 0]
by = {}
for r in rows: by.setdefault((r["dim"], r["N"]), []).append(r)
idx = {int(t): i for i, t in enumerate(cl.GRID)}
Ns = [500, 1000, 2000, 4000, 8000]; out = []
for dim in (1, 3):
    for con in ("k0", "k3", "k15", "rgg0", "rgg3", "rgg15"):
        for gname, grid in (("G2", cl.GRID2), ("G4", cl.GRID4)):
            ci = [idx[int(t)] for t in grid]
            for i in range(2, 5):
                trip = [Ns[i], Ns[i - 1], Ns[i - 2]]
                fam = []
                for m in trip:
                    rs = sorted(by[(dim, m)], key=lambda x: x["rep"]); fam.append((np.array([r[f"q_{con}"][ci] for r in rs]), 32.0 * m / len(rs), grid, m))
                rng = np.random.default_rng(np.random.SeedSequence([22, 779, dim, trip[0], len(con), len(gname)]))
                r = L.est_E3(fam, rng, None, B=80, thr=1.05, t0=float(grid[0]), dmax=14.0)
                out.append(dict(dim=dim, con=con, grid=gname, N=trip[0], **{k: r.get(k) for k in ("abstain", "d_early", "d_late", "decl_neg", "decl_pos", "hw")}))
pd.DataFrame(out).to_csv("results/causet_E3_wide_exploratory.csv", index=False); print(pd.DataFrame(out).query("grid=='G4'").round(2).to_string())
