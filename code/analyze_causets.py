"""Analysis for PLAN.md AMENDMENT 2 (causal sets). Reads results/causet/*.pkl, writes results/causet_*.csv/json and figures."""
import glob, pickle, json, sys, os
import numpy as np, pandas as pd
import p22lib as L, p22amend as A
import causet_lib as cl

CAP = {}
def _wrap(orig):
    def sc(dhat, Dstar, dtrue, tag, dpop=None):
        r = orig(dhat, Dstar, dtrue, tag, dpop)
        se = np.maximum(Dstar.std(0, ddof=1), 1e-9)
        c = np.quantile((np.abs(Dstar - dhat) / se).max(1), 0.95)
        CAP[tag] = (np.asarray(dhat).copy(), c * se)
        return r
    return sc
_orig = L.score
L.score = _wrap(_orig); A.score = L.score
DUM = dict(d_exact=lambda t: np.zeros(len(np.atleast_1d(t))), q_exact=lambda t: np.ones(len(np.atleast_1d(t))))
CONS = ["k0", "k3", "k15", "rgg0", "rgg3", "rgg15"]
GRIDS = {"G2": (cl.GRID2, 2), "G4": (cl.GRID4, 4)}
THRS = (6.0, 1.3)
B = 200
DIR = "results/causet_smoke" if "--smoke" in sys.argv else "results/causet"
SUF = "_smoke" if "--smoke" in sys.argv else ""


def load():
    rows = [pickle.load(open(f, "rb")) for f in glob.glob(DIR + "/r_*.pkl") if os.path.getsize(f) > 0]  # size-0 files are skip-markers (3+1 N=8000 reps 9-11 cut for time)
    return rows


def cols(grid):
    idx = {int(t): i for i, t in enumerate(cl.GRID)}
    return [idx[int(t)] for t in grid]


def row_of(est, r, dim, con, N, gname, thr, nunits, ts):
    out = dict(dim=dim, con=con, N=N, grid=gname, thr=thr, est=est, R=nunits)
    out["abstain"] = int(r.get("abstain", 1))
    if out["abstain"]:
        out["why"] = r.get("why", ""); return out
    for k in ("J", "tlast", "decl_neg", "decl_pos", "decl_zero", "dmin", "dmax", "hw", "h", "d_early", "d_late", "tsplit"):
        if k in r: out[k] = r[k]
    return out


def analyse(rows):
    out, curves = [], []
    by = {}
    for r in rows:
        by.setdefault((r["dim"], r["N"]), []).append(r)
    Ns = sorted({k[1] for k in by})
    for dim in (1, 3):
        for con in CONS:
            for gname, (grid, step) in GRIDS.items():
                ci = cols(grid)
                for thr in THRS:
                    for N in Ns:
                        if (dim, N) not in by: continue
                        rs = sorted(by[(dim, N)], key=lambda x: x["rep"])
                        U = np.array([r[f"q_{con}"][ci] for r in rs]) # units x K
                        n = len(rs); mw = 32.0 * N / n
                        rng = np.random.default_rng(np.random.SeedSequence([22, 777, dim, N, CONS.index(con), step, int(thr * 10)]))
                        def rec(est, res, tag=None):
                            out.append(row_of(est, res, dim, con, N, gname, thr, n, grid))
                            if tag and tag in CAP and not res.get("abstain"):
                                d, hw = CAP[tag]
                                J = len(d)
                                tt = grid[:J] if est == "E1q" else (np.sqrt(grid[:J] * grid[step:J + step].astype(float)) if est == "E2" else None)
                                curves.append(dict(dim=dim, con=con, N=N, grid=gname, thr=thr, est=est, t=None if tt is None else [float(x) for x in tt],
                                                   d=[float(x) for x in d], hw=[float(x) for x in hw]))
                        CAP.clear()
                        res = A.est_E1b([U], mw, grid, N, rng, DUM, deg=2, under=1.0, thr=thr, tag="E1q"); rec("E1q", res, "E1q")
                        CAP.clear()
                        res = L.est_E0([U], mw, grid, N, rng, DUM, B=B, thr=thr)
                        rec("E0", res[1], "E0")
                        if "E0" in CAP and not res[1].get("abstain"):
                            curves.append(dict(dim=dim, con=con, N=N, grid=gname, thr=thr, est="E0seg", t=None, d=[float(x) for x in CAP["E0"][0]], hw=[float(x) for x in CAP["E0"][1]]))
                        CAP.clear()
                        res = L.est_E2([U], grid, N, rng, DUM, B=400, thr=thr, step=step)
                        rec("E2", res[0], "E2")
                    if thr != 6.0: continue
                    # E3: triplets of consecutive sizes, labelled by the largest
                    for i in range(2, len(Ns)):
                        trip = [Ns[i], Ns[i - 1], Ns[i - 2]]
                        if any((dim, m) not in by for m in trip): continue
                        fam = []
                        for m in trip:
                            rs = sorted(by[(dim, m)], key=lambda x: x["rep"])
                            Um = np.array([r[f"q_{con}"][ci] for r in rs]); fam.append((Um, 32.0 * m / len(rs), grid, m))
                        rng = np.random.default_rng(np.random.SeedSequence([22, 778, dim, trip[0], CONS.index(con), step, int(thr * 10)]))
                        res = L.est_E3(fam, rng, None, B=80, thr=1.05 if thr == 6.0 else 1.05, t0=2.0)
                        o = row_of("E3", res, dim, con, trip[0], gname, thr, len(rs), grid); o["triplet"] = str(trip); out.append(o)
    return pd.DataFrame(out), curves


if __name__ == "__main__":
    rows = load()
    print("loaded", len(rows))
    df, curves = analyse(rows)
    df.to_csv(f"results/causet_estimates{SUF}.csv", index=False)
    json.dump(curves, open(f"results/causet_curves{SUF}.json", "w"))
    meta = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("q_")} for r in rows])
    meta.to_csv(f"results/causet_graph_meta{SUF}.csv", index=False)
    print(df.groupby(["est", "con"]).abstain.mean())
