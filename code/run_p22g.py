#!/usr/bin/env python3
"""P22 Amendment 7 runner (walker budget). Usage: run_p22g.py [--smoke] [--ndata N] [--ncpu K] [--wpn 8,128,32]
Cells: wpn 8 (ids 500+), 128 (600+): 8 controls + 19 crossover; wpn 32 fill-in (700+): 7 crossover. Output results/raw_amend7[_smoke].csv"""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A, p22amend6 as A6
import run_p22 as R

args = sys.argv[1:]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 16
assert NCPU <= 16
WPNS = [int(x) for x in args[args.index("--wpn") + 1].split(",")] if "--wpn" in args else [8, 128, 32]
CTRL = [("ring", 256, (256,)), ("ring", 1024, (1024,)), ("torus2", 24, (24, 24)), ("torus2", 48, (48, 48)), ("torus2", 96, (96, 96)),
        ("torus3", 16, (16,) * 3), ("torus3", 24, (24,) * 3), ("torus3", 32, (32,) * 3)]
G3 = [32, 48, 64, 80, 96, 112, 128, 160, 192, 256]; G2 = [128, 160, 192, 224, 256, 320, 384, 512, 768]
FILL3 = [80, 112, 160, 256]; FILL2 = [160, 224, 320]
CELLS = None

def T(d): return L.Graph("torus", 0.5, dims=tuple(d))

def specs():
    S = []
    def add(fam, l, dims, kind, wpn, base):
        S.append(dict(fam=fam, L=l, build=lambda d=dims: T(d), sizes=None, plateau=None, kind=kind, grid=(24, 34), lazy=0.5, note="", wpn=wpn, id=base + len([s for s in S if s["wpn"] == wpn])))
    for wpn, base in ((8, 500), (128, 600)):
        if wpn not in WPNS: continue
        for fam, l, d in CTRL: add(fam, l, d, "control", wpn, base)
        for l in G3: add("aniso3_w8", l, (l, 8, 8), "crossover", wpn, base)
        for l in G2: add("aniso2_w16", l, (l, 16), "crossover", wpn, base)
    if 32 in WPNS:
        for l in FILL3: add("aniso3_w8", l, (l, 8, 8), "crossover", 32, 700)
        for l in FILL2: add("aniso2_w16", l, (l, 16), "crossover", 32, 700)
    return S

def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(9)]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts = c["ts_w"]; N = c["N"]
    C, m = c["graph"].sample_blocks(rngs[1], c["wpn"] * N, ts, L.NBLOCK); U = [C / m]
    out = list(L.est_E0(U, m, ts, N, rngs[5], truth, 200, c["thr"]))
    out.append(A.est_E1b(U, m, ts, N, rngs[6], truth, 2, 1.0, 200, c["thr"], "E1q"))
    out += A6.est_all(U, m, ts, N, rngs[2], rngs[3], rngs[4], truth, c["thr"], ("E4c", "E4r"))
    for r in out:
        r.update(cell=c["id"], ds=ds, wpn=c["wpn"])
        for k in ("e6lo", "e6hi", "e4lo", "e4hi"): r.pop(k, None)
    return out

def worker(task):
    cid, lo, hi = task; c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0

def main():
    global CELLS
    t0 = time.time(); sp = specs()
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (4 if SMOKE else 200)
    rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    CELLS = {s["id"]: s for s in sp}
    sfx = "_smoke" if SMOKE else ""
    keys = ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(R.RES, f"cells_meta_amend7{sfx}.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(sp)} cells, ncpu {NCPU}", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = sorted([(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[0]]["N"] * CELLS[t[0]]["wpn"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
            if k % 100 == 0: pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend7{sfx}.partial.csv"), index=False)
    pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend7{sfx}.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}_w{CELLS[k]["wpn"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, f"cell_cpu_seconds_amend7{sfx}.csv"))
    print(f"done {len(rows_all)} rows in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__": main()
