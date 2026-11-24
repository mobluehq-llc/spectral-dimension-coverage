#!/usr/bin/env python3
"""P22 amendment runner: E1u/E1q/E1qu/E3s on the 14 control + 13 crossover cells, same seeds/walker blocks as run_p22.py.
Usage: .venv/bin/python run_p22b.py [--smoke] [--ndata N]. Output results/raw_amend.csv (+ results/amend_repro_check.csv in smoke)."""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A
import run_p22 as R

SMOKE = "--smoke" in sys.argv
NCPU = 9
CELLS = None


def run_dataset(c, ds, check=False):
    ss = np.random.SeedSequence([22, c["id"], ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(6)]       # identical stream layout to run_p22.py
    rb = [np.random.default_rng(s) for s in np.random.SeedSequence([23, c["id"], ds]).spawn(4)]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts_w, N = c["ts_w"], c["N"]
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts_w, L.NBLOCK)
    U = [C / m]
    out = [A.est_E1b(U, m, ts_w, N, rb[0], truth, 1, 0.5, 200, c["thr"], "E1u"),
           A.est_E1b(U, m, ts_w, N, rb[1], truth, 2, 1.0, 200, c["thr"], "E1q"),
           A.est_E1b(U, m, ts_w, N, rb[2], truth, 2, 0.7, 200, c["thr"], "E1qu"),
           A.est_E3s(U, m, ts_w, N, rb[3], truth)]
    if check:
        out.append(L.est_E1(U, m, ts_w, N, rngs[2], truth, True, 200, c["thr"], "E1"))
    for r in out:
        r.update(cell=c["id"], ds=ds)
    return out


def worker(task):
    cid, lo, hi, check = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi):
        rows += run_dataset(c, ds, check)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    t0 = time.time()
    specs = R.cell_specs()          # ids identical to the first run
    nd = 8 if SMOKE else 200
    if "--ndata" in sys.argv: nd = int(sys.argv[sys.argv.index("--ndata") + 1])
    specs = [s for s in specs if s["kind"] in ("control", "crossover") and s["note"] == ""]
    rr = np.random.default_rng(2200)
    for s in specs: R.prep(s, rr)
    CELLS = {s["id"]: s for s in specs}
    print(f"prep {time.time()-t0:.0f}s, {len(specs)} cells", flush=True)
    chunk = 4 if SMOKE else 10
    tasks = [(s["id"], lo, min(lo + chunk, nd), SMOKE) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows_all)
    out = "raw_amend_smoke.csv" if SMOKE else "raw_amend.csv"
    df.to_csv(os.path.join(R.RES, out), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, "cell_cpu_seconds_amend%s.csv" % ("_smoke" if SMOKE else "")))
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
