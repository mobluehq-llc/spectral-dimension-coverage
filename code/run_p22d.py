#!/usr/bin/env python3
"""P22 Amendment 4 runner: E5qu/E5q/E6 on the 27 control+crossover cells. Usage: run_p22d.py [--smoke] [--ndata N]"""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22shape as S
import run_p22 as R
SMOKE = "--smoke" in sys.argv
NCPU = 9
CELLS = None

def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(6)]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * c["N"], c["ts_w"], L.NBLOCK)
    out = S.est_E5E6([C / m], m, c["ts_w"], c["N"], rngs[2], rngs[3], truth, c["thr"])
    for r in out: r.update(cell=c["id"], ds=ds)
    return out

def worker(task):
    cid, lo, hi = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0

def main():
    global CELLS
    t0 = time.time(); specs = R.cell_specs()
    nd = 8 if SMOKE else 200
    if "--ndata" in sys.argv: nd = int(sys.argv[sys.argv.index("--ndata") + 1])
    specs = [s for s in specs if s["kind"] in ("control", "crossover") and s["note"] == ""]
    rr = np.random.default_rng(2200)
    for s in specs: R.prep(s, rr)
    CELLS = {s["id"]: s for s in specs}
    chunk = 4 if SMOKE else 10
    tasks = [(s["id"], lo, min(lo + chunk, nd)) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows_all)
    df.to_csv(os.path.join(R.RES, "raw_amend4_smoke.csv" if SMOKE else "raw_amend4.csv"), index=False)
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__": main()
