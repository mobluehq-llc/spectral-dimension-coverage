#!/usr/bin/env python3
"""P22 Amendment 3 runner: E4o/E4c/E4r on the 27 control+crossover cells. Usage: .venv/bin/python run_p22c.py [--smoke] [--ndata N] [--cells fam:L,...]"""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22honest as H
import run_p22 as R

SMOKE = "--smoke" in sys.argv
NCPU = 9
CELLS = None


def run_dataset(c, ds, check=False):
    ss = np.random.SeedSequence([22, c["id"], ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(6)]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts_w, N = c["ts_w"], c["N"]
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts_w, L.NBLOCK)
    U = [C / m]
    out = H.est_E4(U, m, ts_w, N, rngs[2], truth, c["Mo_cum"], c["thr"])
    if check:
        out.append(L.est_E1(U, m, ts_w, N, np.random.default_rng(ss.spawn(6)[2]) if False else np.random.default_rng(np.random.SeedSequence([22, c["id"], ds]).spawn(6)[2]), truth, True, 200, c["thr"], "E1"))
    for r in out:
        r.update(cell=c["id"], ds=ds)
    return out


def worker(task):
    cid, lo, hi, check = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi):
        rows += run_dataset(c, ds, check)
    return cid, rows, time.time() - t0, dict(H.STATS)


def main():
    global CELLS
    t0 = time.time()
    specs = R.cell_specs()
    nd = 8 if SMOKE else 200
    if "--ndata" in sys.argv: nd = int(sys.argv[sys.argv.index("--ndata") + 1])
    specs = [s for s in specs if s["kind"] in ("control", "crossover") and s["note"] == ""]
    if "--cells" in sys.argv:
        want = set(sys.argv[sys.argv.index("--cells") + 1].split(","))
        specs = [s for s in specs if f'{s["fam"]}:{s["L"]}' in want]
    rr = np.random.default_rng(2200)
    allspecs = [s for s in R.cell_specs() if s["kind"] in ("control", "crossover") and s["note"] == ""]
    # prep consumes rr only for rrg cells (not in this set); same rr seed as run_p22b
    for s in specs:
        R.prep(s, rr)
        s["Mo_cum"] = H.exact_g3_cummax(s["graph"], s["ts_w"])
    CELLS = {s["id"]: s for s in specs}
    print(f"prep {time.time()-t0:.0f}s, {len(specs)} cells", flush=True)
    chunk = 4 if SMOKE else 10
    tasks = [(s["id"], lo, min(lo + chunk, nd), SMOKE) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}; st = {"fail": 0, "solves": 0}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt, s2) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows_all)
    out = "raw_amend3_smoke.csv" if SMOKE else "raw_amend3.csv"
    df.to_csv(os.path.join(R.RES, out), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, "cell_cpu_seconds_amend3%s.csv" % ("_smoke" if SMOKE else "")))
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
