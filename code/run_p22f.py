#!/usr/bin/env python3
"""P22 Amendment 6 runner. Usage: run_p22f.py GROUP [--smoke] [--ndata N] [--cells fam:L,...] [--ncpu K]
GROUP: orig (37 Amendment-5 cells, seeds [22,id,ds]); fresh (replication cells, seeds [22,id,1000+ds]); sep (twostep C256xC8xCb, ids 300-303); big2d (aniso2_w16 L=384,512,768, ids 400-402; E4c,E4r,E6).
Output results/raw_amend6_<GROUP>.csv"""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend6 as A6, p22honest as H
import run_p22 as R, run_p22e as RE

args = sys.argv[1:]
GROUP = args[0]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 16
assert NCPU <= 16
CELLS = None
SEEDOFF = 1000 if GROUP == "fresh" else 0


def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))


def sep_specs():
    S = []
    for i, b in enumerate((12, 16, 20, 24)):
        S.append(dict(fam=f"twostep_8_{b}", L=256, build=lambda b=b: T(256, 8, b), sizes=None, plateau=None, kind="twostep", grid=(24, 34), lazy=0.5, note="", id=300 + i))
    return S


def big_specs():
    S = []
    for i, l in enumerate((384, 512, 768)):
        S.append(dict(fam="aniso2_w16", L=l, build=lambda l=l: T(l, 16), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note="", id=400 + i))
    return S


def all_specs():
    old = [s for s in R.cell_specs() if s["kind"] in ("control", "crossover") and s["note"] == ""]
    return old + RE.new_specs()


def select():
    if GROUP == "sep": return sep_specs()
    if GROUP == "big2d": return big_specs()
    sp = all_specs()
    if GROUP == "fresh":
        keep = lambda s: (s["fam"] == "aniso2_w16" and s["L"] >= 128) or (s["fam"] == "aniso3_w8" and s["L"] >= 48) or s["fam"] == "twostep_8_32" or s["fam"].startswith("brush")
        sp = [s for s in sp if keep(s)]
    return sp


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], SEEDOFF + ds])
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch[:6]]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * c["N"], c["ts_w"], L.NBLOCK)
    tags = ("E4c", "E4r") if GROUP == "big2d" else ("E4c",)
    out = A6.est_all([C / m], m, c["ts_w"], c["N"], rngs[2], rngs[3], rngs[4], truth, c["thr"], tags)
    for r in out: r.update(cell=c["id"], ds=SEEDOFF + ds)
    return out


def worker(task):
    cid, lo, hi = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    t0 = time.time(); specs = select()
    if "--cells" in args:
        want = set(args[args.index("--cells") + 1].split(","))
        specs = [s for s in specs if f'{s["fam"]}:{s["L"]}' in want]
    nd = 8 if SMOKE else 200
    if "--ndata" in args: nd = int(args[args.index("--ndata") + 1])
    rr = np.random.default_rng(2200)
    for s in specs: R.prep(s, rr)
    CELLS = {s["id"]: s for s in specs}
    keys = ("id", "fam", "L", "N", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end")
    sfx = "%s%s" % (GROUP, "_smoke" if SMOKE else "")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in specs]).to_csv(os.path.join(R.RES, f"cells_meta_amend6_{sfx}.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(specs)} cells, group {GROUP}, ncpu {NCPU}", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = [(s["id"], lo, min(lo + chunk, nd)) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows_all)
    df.to_csv(os.path.join(R.RES, f"raw_amend6_{sfx}.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, f"cell_cpu_seconds_amend6_{sfx}.csv"))
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__": main()
