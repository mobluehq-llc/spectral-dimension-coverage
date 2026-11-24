#!/usr/bin/env python3
"""EXPLORATORY follow-up of exploratory-log row 29 (Amendment 12, wpn=128, E0g misses on cells 9212 / 9213 = torus3 L=28 / L=32, controls).
NOT part of the pre-registered programme; no confirmatory claim. Written after the Amendment 12 results were seen.
Uses the frozen run_amend12.run_dataset UNCHANGED (same seeds SeedSequence([22,12000+id,ds])) so ds 0..199 must reproduce the stored rows,
and ds 200..799 are new replicates (3x the original 200). Extra: a mirror of est_E0's global arm with the window length FIXED (J = 5..10 and the
data-driven J) on the same walker data, to separate window selection from bias and SE miscalibration.
usage: row29_explore.py OUTDIR NCPU DS_LO DS_HI"""
import os, sys, time, pickle, json
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k] = "1"
OUT, NCPU, LO, HI = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np, multiprocessing as mp
sys.argv = ["run_amend12.py", "b128"]
import run_amend12 as RA, run_p22 as R, run_p22i as RI, p22lib as L
RA.verify()
RI.GROUP = "b128"; SP = [s for s in RI.select() if s["id"] in (9212, 9213)]
rr = np.random.default_rng(2200)
for s in SP: R.prep(s, rr)
CELLS = {s["id"]: s for s in SP}; RA.CELLS = CELLS
FIXJ = (5, 6, 7, 8, 9, 10)

def fixedJ_rows(c, ds):
    """same walker data as run_dataset (rngs[1]); E0g-style global weighted fit with bootstrap SE, for data-driven J and fixed J."""
    ss = np.random.SeedSequence([22, 12000 + c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    p = L.grand_mean([U]); W = U.shape[0] * m; cnt = p * W; Jd = L.window_len(cnt, p, N, thr)
    rng = np.random.default_rng(np.random.SeedSequence([22, 29000 + c["id"], ds]))   # own stream: extra variants must not disturb nothing stored
    Ps_all = L.boot_means([U], 200, rng); floor = 0.5 / W
    rows = []
    for tag, J in [("dataJ", Jd)] + [(f"J{j}", j) for j in FIXJ]:
        if J < 3 or J > len(ts): continue
        x = np.log(ts[:J].astype(float)); w = cnt[:J]; xb = (w * x).sum() / w.sum(); row = w * (x - xb) / (w * (x - xb) ** 2).sum()
        y = np.log(np.maximum(p[:J], floor)); Ys = np.log(np.maximum(Ps_all[:, :J], floor))
        dh = float(-2 * (row @ y)); Ds = -2 * Ys @ row; se = float(Ds.std(0, ddof=1))
        tref = np.exp(x.mean()); dt = float(c["d_exact"](np.array([tref]))[0]); qe = np.log(c["q_exact"](ts[:J])); dpop = float(-2 * (row @ qe))
        rows.append(dict(cell=c["id"], ds=ds, variant=tag, J=J, Jdata=Jd, dhat=dh, se_boot=se, d_true=dt, d_pop=dpop,
                         n_cnt_min=float(cnt[:J].min()), Np_min=float((N * p[:J]).min())))
    return rows

def work(task):
    cid, lo, hi = task
    rows, ex = [], []
    for ds in range(lo, hi):
        rows += RA.run_dataset(CELLS[cid], ds); ex += fixedJ_rows(CELLS[cid], ds)
    return cid, lo, hi, rows, ex

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True); CH = 5
    tasks = [(cid, lo, min(lo + CH, HI)) for lo in range(LO, HI, CH) for cid in CELLS]
    todo = [t for t in tasks if not os.path.exists(f"{OUT}/t_{t[0]}_{t[1]}_{t[2]}.pkl")]
    print(time.strftime("%F %T"), f"{len(tasks)} tasks, {len(todo)} todo", flush=True); t0 = time.time()
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, r in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f = f"{OUT}/t_{r[0]}_{r[1]}_{r[2]}.pkl"; pickle.dump(r, open(f + ".tmp", "wb")); os.replace(f + ".tmp", f)
            if k % 5 == 0: print(time.strftime("%F %T"), f"{k+1}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("DONE", flush=True)
