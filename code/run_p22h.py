#!/usr/bin/env python3
"""P22 Amendment 8 runner (E7 moment-cone LP band, E8 secant bracket, A1' alarm). Usage:
  run_p22h.py GROUP [--smoke] [--ndata N] [--ncpu K] [--cells id,id,...] [--verify]
GROUP: core (wpn 32: 14 controls, 13 old crossovers + 3 big2d + 7 fill-in, 10 misspecified), w8 / w128 (Amendment-7 cells, ungated).
Walker data: SeedSequence([22, id, ds]).spawn(9)[1] -> sample_blocks, exactly as run_p22e/f/g. E6 (for A1') on controls + misspecified cells
(+ all cells with --verify) with streams [2],[3] as in run_p22f; E7/E8 bootstrap uses child [7]. Output results/raw_amend8_<GROUP>.csv"""
import os, sys, time
_argv = sys.argv; sys.argv = [_argv[0], "--wpn", "8,128,32"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend6 as A6, p22amend8 as A8
import run_p22 as R, run_p22e as RE
import run_p22g as RG
sys.argv = _argv

args = sys.argv[1:]
GROUP = args[0]
SMOKE = "--smoke" in args; VERIFY = "--verify" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 16
assert NCPU <= 16
CELLS = None
E6_KINDS = ("control", "twostep", "brush", "gasket")


def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))


def big_specs():
    return [dict(fam="aniso2_w16", L=l, build=lambda l=l: T(l, 16), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note="", id=400 + i, wpn=32)
            for i, l in enumerate((384, 512, 768))]


def select():
    if GROUP == "core":
        old = [s for s in R.cell_specs() if s["kind"] in ("control", "crossover") and s["note"] == ""]
        for s in old: s["wpn"] = 32
        new = RE.new_specs()
        for s in new: s["wpn"] = 32
        fill = [s for s in RG.specs() if s["wpn"] == 32]
        return old + new + big_specs() + fill
    w = int(GROUP[1:])
    return [s for s in RG.specs() if s["wpn"] == w]


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds])
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts, N, wpn = c["ts_w"], c["N"], c["wpn"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    out = []; J6 = 0; e6lo = e6hi = None
    if c["kind"] in E6_KINDS or VERIFY:
        rows = A6.est_all([U], m, ts, N, rngs[2], rngs[3], rngs[4], truth, c["thr"], ())
        r6 = rows[0]
        if r6.get("abstain", 1) == 0:
            e6lo = np.array([float(v) for v in r6["e6lo"].split(";")]); e6hi = np.array([float(v) for v in r6["e6hi"].split(";")]); J6 = len(e6lo)
        out.append(r6)
    out += A8.run_e7_e8(U, m, ts, N, rngs[7], truth, c["thr"], J6, e6lo, e6hi)
    for r in out: r.update(cell=c["id"], ds=ds, wpn=wpn)
    return out


def worker(task):
    cid, lo, hi = task; c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    t0 = time.time(); sp = select()
    if "--cells" in args:
        want = {int(x) for x in args[args.index("--cells") + 1].split(",")}; sp = [s for s in sp if s["id"] in want]
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (8 if SMOKE else 200)
    rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    CELLS = {s["id"]: s for s in sp}
    sfx = "%s%s" % (GROUP, "_smoke" if SMOKE else "")
    keys = ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(R.RES, f"cells_meta_amend8_{sfx}.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(sp)} cells, group {GROUP}, ncpu {NCPU}", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = sorted([(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[0]]["N"] * CELLS[t[0]]["wpn"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
            if k % 100 == 0: pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend8_{sfx}.partial.csv"), index=False)
    pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend8_{sfx}.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}_w{CELLS[k]["wpn"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, f"cell_cpu_seconds_amend8_{sfx}.csv"))
    print(f"done {len(rows_all)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__": main()
