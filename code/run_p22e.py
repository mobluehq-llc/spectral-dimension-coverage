#!/usr/bin/env python3
"""P22 Amendment 5 runner. Usage: run_p22e.py [--smoke] [--ndata N] [--cells fam:L,...]. Output results/raw_amend5.csv, cells_meta_amend5.csv."""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A, p22honest as H, p22amend5 as A5
import run_p22 as R

SMOKE = "--smoke" in sys.argv
NCPU = 16  # as run on host-T (local copy used 9; the only change)
CELLS = None


def new_specs():
    S = []
    def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))
    def add(fam, Lval, build, sizes, kind):
        S.append(dict(fam=fam, L=Lval, build=build, sizes=sizes, plateau=None, kind=kind, grid=(24, 34), lazy=0.5, note="", new=True))
    for l in (128, 192, 256, 320):
        add("twostep_8_32", l, lambda l=l: T(l, 8, 32), [lambda l=l: T(l, 8, 32), lambda l=l: T(int(.75 * l), 8, 32), lambda l=l: T(l // 2, 8, 32)], "twostep")
    for lb, lf in ((8, 24), (12, 16), (16, 16)):
        add(f"brush_f{lf}", lb, lambda lb=lb, lf=lf: A5.brush(lb, lf),
            [lambda lb=lb, lf=lf: A5.brush(lb, lf), lambda lb=lb, lf=lf: A5.brush(int(.75 * lb), lf), lambda lb=lb, lf=lf: A5.brush(lb // 2, lf)], "brush")
    for n in (5, 6, 7):
        add("gasket", n, lambda n=n: L.gasket(n), None, "gasket")
    for i, s in enumerate(S): s["id"] = 100 + i
    return S


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds])
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch[:6]]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts_w, N, thr = c["ts_w"], c["N"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts_w, L.NBLOCK)
    U = [C / m]
    out = [A5.est_E6(U, m, ts_w, N, rngs[2], rngs[3], truth, thr)]
    out += A5.est_M(U, m, ts_w, N, np.random.default_rng(ch[7]), truth, thr, plateau=c["plateau"])
    if c.get("new"):
        out.append(A.est_E1b(U, m, ts_w, N, np.random.default_rng(ch[8]), truth, 2, 1.0, 200, thr, "E1q"))
        e0 = L.est_E0(U, m, ts_w, N, np.random.default_rng(ch[6]), truth, 200, thr); out.append(e0[1])
        e4 = H.est_E4(U, m, ts_w, N, np.random.default_rng(ch[4]), truth, c["Mo_cum"], thr); out.append(e4[1])
    if c.get("size_graphs") is not None:
        fam = [(U[0], m, ts_w, N)]
        for sg, cap in zip(c["size_graphs"][1:], c["cap_sizes"][1:]):
            tsi = L.make_grid(*c["grid"], int(cap))
            Ci, mi = sg.sample_blocks(rngs[5], R.WPN * sg.N, tsi, L.NBLOCK)
            fam.append((Ci / mi, mi, tsi, sg.N))
        out.append(A5.est_CL(fam, rngs[5], truth, c["plateau"]))
    for r in out: r.update(cell=c["id"], ds=ds)
    return out


def worker(task):
    cid, lo, hi = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    t0 = time.time()
    old = [s for s in R.cell_specs() if s["kind"] in ("control", "crossover") and s["note"] == ""]
    specs = old + new_specs()
    if "--cells" in sys.argv:
        want = set(sys.argv[sys.argv.index("--cells") + 1].split(","))
        specs = [s for s in specs if f'{s["fam"]}:{s["L"]}' in want]
    nd = 8 if SMOKE else 200
    if "--ndata" in sys.argv: nd = int(sys.argv[sys.argv.index("--ndata") + 1])
    rr = np.random.default_rng(2200)
    for s in specs:
        R.prep(s, rr)
        if s.get("new"): s["Mo_cum"] = H.exact_g3_cummax(s["graph"], s["ts_w"])
    CELLS = {s["id"]: s for s in specs}
    keys = ("id", "fam", "L", "N", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "cap_e2")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in specs]).to_csv(os.path.join(R.RES, "cells_meta_amend5%s.csv" % ("_smoke" if SMOKE else "")), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(specs)} cells", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = [(s["id"], lo, min(lo + chunk, nd)) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows_all)
    df.to_csv(os.path.join(R.RES, "raw_amend5%s.csv" % ("_smoke" if SMOKE else "")), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, "cell_cpu_seconds_amend5%s.csv" % ("_smoke" if SMOKE else "")))
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__": main()
