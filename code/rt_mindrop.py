#!/usr/bin/env python3
"""SA3 pick 2 (minimum compatible drop by bisection on delta). Pre-registration: PREREG_REDTEAM_FINAL.md, amendment 11-24 PM (P6).
Imports rt_nodrop (frozen, not edited) and sets its module DELTA per call.
Usage: rt_mindrop.py exact | data [--ndata N] [--ncpu K] | controls [--ndata N]"""
import os, sys, time
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np, pandas as pd
import multiprocessing as mp
import rt_nodrop as RT
import p22lib as L, p22amend8 as A8
import run_p22 as R
sys.argv = _argv
args = sys.argv[1:]
TOL = 0.005


class Tester:
    def __init__(self, ts, lo, hi, N, Jw):
        self.ts, self.lo, self.hi, self.N, self.Jw = ts, lo, hi, N, Jw
        self.S = A8.atoms_grid and RT.A8.lp_setup(ts, lo, hi, A8.atoms_grid(float(ts[-1]), RT.NATOMS), N)
        self.S2 = None
        self.d0min, self.d0max, self.st = A8.lp_one(self.S, 0)
        self.nlp = 0; self.log = []

    def test(self, delta):
        """outcome at this delta with the 1597-atom gate on DROP (same logic as rt_nodrop.global_nodrop)"""
        RT.DELTA = delta
        res = RT.find_nodrop(self.S, self.lo, self.hi, self.N, self.Jw, self.d0min, self.d0max); self.nlp += res["nlp"]
        out = res["outcome"]
        if out == "DROP":
            if self.S2 is None:
                self.S2 = RT.A8.lp_setup(self.ts, self.lo, self.hi, A8.atoms_grid(float(self.ts[-1]), 2 * (RT.NATOMS - 1) + 1), self.N)
                self.b2 = A8.lp_one(self.S2, 0)
            a0, a1, st2 = self.b2
            if st2 != 0: out = "INCONCLUSIVE"
            else:
                r2 = RT.find_nodrop(self.S2, self.lo, self.hi, self.N, self.Jw, a0, a1); self.nlp += r2["nlp"]
                if r2["outcome"] != "DROP": out = r2["outcome"] if r2["outcome"] != "WITNESS" else "WITNESS"
        self.log.append((round(delta, 5), out))
        return out

    def bisect(self):
        if self.st != 0: return dict(D_lo=np.nan, D_hi=np.nan, status="BOXINFEASIBLE", steps=0, n_incon=0)
        top = float(np.ceil(self.d0max / TOL) * TOL)
        o = self.test(0.0)
        if o != "DROP":
            return dict(D_lo=0.0, D_hi=0.0, status="NODROP0_" + o, steps=1, n_incon=int(o == "INCONCLUSIVE"))
        lo, hi = 0.0, top
        while hi - lo > TOL:
            m = 0.5 * (lo + hi); o = self.test(m)
            if o == "DROP": lo = m
            else: hi = m
        ninc = sum(1 for d, o in self.log if o == "INCONCLUSIVE")
        return dict(D_lo=lo, D_hi=hi, status="OK", steps=len(self.log), n_incon=ninc)


def data_box(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    p = U.mean(0); Wt = U.shape[0] * m; cnt = p * Wt
    Jw = L.window_len(cnt, p, N, thr)
    if Jw < 2: return None
    qhat, lo, hi, Ps, binfo = A8.make_box(U, m, rngs[7])
    return ts, lo, hi, N, Jw


CELLS = None
def work_data(task):
    cid, ds = task; c = CELLS[cid]; t0 = time.time()
    base = dict(cell=cid, fam=c["fam"], L=c["L"], ds=ds, mode="data", exact_drop=float(c["max_drop"]))
    b = data_box(c, ds)
    if b is None: return dict(base, status="NOWINDOW")
    ts, lo, hi, N, Jw = b
    T = Tester(ts, lo, hi, N, Jw); r = T.bisect()
    return dict(base, Jw=Jw, d0max=T.d0max, nlp=T.nlp, seconds=time.time() - t0, trace=";".join(f"{d}:{o[0]}" for d, o in T.log), **r)


def work_exact(task):
    cid, eps = task; c = CELLS[cid]; t0 = time.time()
    ts, N = c["ts_w"], c["N"]; Jw = int(c["exact_J"]); q = c["graph"].exact(ts)[0]
    T = Tester(ts, q * (1 - eps), q * (1 + eps), N, Jw); r = T.bisect()
    return dict(cell=cid, fam=c["fam"], L=c["L"], ds=-1, mode="exact", eps=eps, exact_drop=float(c["max_drop"]), Jw=Jw, d0max=T.d0max, nlp=T.nlp,
                seconds=time.time() - t0, trace=";".join(f"{d}:{o[0]}" for d, o in T.log), **r)


def main():
    global CELLS
    mode = args[0]; t0 = time.time()
    want = RT.CONTROLS if mode == "controls" else (RT.TARGETS + RT.CONTROLS if mode == "exact" else RT.TARGETS)
    sp_all = [s for s in RT.RI.core_specs() if s["kind"] in ("control", "crossover")]
    ids = {s["id"] for s in sp_all if (s["fam"], s["L"]) in want}
    sp = RT.specs(ids); CELLS = {s["id"]: s for s in sp}
    print("prep %.0fs" % (time.time() - t0), [(s["id"], s["fam"], s["L"], round(float(s["max_drop"]), 4)) for s in sp], flush=True)
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 5; assert ncpu <= 5
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (30 if mode == "controls" else 100)
    if mode == "exact":
        tasks = [(s["id"], e) for s in sp for e in (1e-3, 1e-2)]; fn = work_exact; out = "mindrop_exact.csv"
    else:
        tasks = [(s["id"], ds) for ds in range(nd) for s in sp]; fn = work_data; out = "mindrop_controls.csv" if mode == "controls" else "mindrop_data.csv"
    out = os.path.join(R.RES, out); rows = []
    with mp.get_context("fork").Pool(min(ncpu, len(tasks))) as pool:
        for k, r in enumerate(pool.imap_unordered(fn, tasks)):
            rows.append(r)
            if k % 5 == 0: print(f"{k+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True); pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False)
    if mode == "exact":
        for r in rows: print({k: r[k] for k in ("fam", "L", "eps", "exact_drop", "D_lo", "D_hi", "status", "steps", "seconds")}, r["trace"])
    print("done", len(rows), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
