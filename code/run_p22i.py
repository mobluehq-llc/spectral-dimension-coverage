#!/usr/bin/env python3
"""P22 Amendment 9 runner. Usage: run_p22i.py GROUP [--smoke] [--ndata N] [--ncpu K] [--cells id,id,...]
GROUP:
  a      9a: E7j joint change test on the Amendment-8 core cells (controls + crossovers), ORIGINAL data seeds SeedSequence([22, id, ds]) (same data as Amendment 8)
  amis   9a on the ten misspecified cells (ungated)
  b8 / b32 / b128   9b: frozen published arms (E0, E1q, E4c, E4r, E6, E7, E8) at that walker budget, fresh seeds SeedSequence([22, id, 2000 + ds])
  c      9c: second design (grid ratio 2^(1/4), count >= 100, N q >= 10, J >= 8), arms E0, E1q, E4c, E6, E7, E8 at wpn 32, fresh seeds, ids 9300+
Output results/raw_amend9_<GROUP>.csv, cells_meta_amend9_<GROUP>.csv"""
import os, sys, time
_argv = sys.argv; sys.argv = [_argv[0], "--wpn", "32"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A, p22amend6 as A6, p22amend8 as A8, p22amend9 as A9
import run_p22 as R, run_p22e as RE
import run_p22g as RG
sys.argv = _argv

args = sys.argv[1:]
GROUP = args[0]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 16
assert NCPU <= 16
CELLS = None
G3 = RG.G3; G2 = RG.G2
SEED_OFF = 0 if GROUP in ("a", "amis") else 2000
ARMS_C = ("E0", "E1q", "E4c", "E6", "E7", "E8")
IDBASE = {"b32": 9000, "b8": 9100, "b128": 9200, "c": 9300}


def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))


def big_specs():
    return [dict(fam="aniso2_w16", L=l, build=lambda l=l: T(l, 16), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note="", id=400 + i, wpn=32)
            for i, l in enumerate((384, 512, 768))]


def core_specs():
    old = [s for s in R.cell_specs() if s["kind"] in ("control", "crossover") and s["note"] == ""]
    for s in old: s["wpn"] = 32
    new = RE.new_specs()
    for s in new: s["wpn"] = 32
    RG.WPNS[:] = [32]
    fill = [s for s in RG.specs() if s["wpn"] == 32]
    return old + new + big_specs() + fill


def controls():
    c = [s for s in R.cell_specs() if s["kind"] == "control" and s["note"] == ""]
    assert len(c) == 14
    return c


def crossovers_real():
    S = []
    for l in G3:
        if l >= 48: S.append(dict(fam="aniso3_w8", L=l, build=lambda l=l: T(l, 8, 8), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note=""))
    for l in G2: S.append(dict(fam="aniso2_w16", L=l, build=lambda l=l: T(l, 16), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note=""))
    return S


def select():
    if GROUP == "a":
        return [s for s in core_specs() if s["kind"] in ("control", "crossover")]
    if GROUP == "amis":
        return [s for s in core_specs() if s["kind"] in ("twostep", "brush", "gasket")]
    if GROUP in ("b8", "b32", "b128"):
        w = int(GROUP[1:]); S = controls() + crossovers_real()
        for i, s in enumerate(S): s["wpn"] = w; s["id"] = IDBASE[GROUP] + i
        return S
    if GROUP == "c":
        S = controls()
        for l in G3: S.append(dict(fam="aniso3_w8", L=l, build=lambda l=l: T(l, 8, 8), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note=""))
        for l in G2: S.append(dict(fam="aniso2_w16", L=l, build=lambda l=l: T(l, 16), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note=""))
        for i, s in enumerate(S): s["wpn"] = 32; s["id"] = IDBASE["c"] + i
        return S
    raise SystemExit("bad group")


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], SEED_OFF + ds])
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    if GROUP in ("a", "amis"):
        out = [A9.run_e7j(U, m, ts, N, rngs[7], thr)]
    else:
        arms = ARMS_C if GROUP == "c" else ("E0", "E1q", "E4c", "E4r", "E6", "E7", "E8")
        out = []
        if "E0" in arms:
            out += [r for r in L.est_E0([U], m, ts, N, rngs[5], truth, 200, thr) if r["est"] == "E0"]
        out.append(A.est_E1b([U], m, ts, N, rngs[6], truth, 2, 1.0, 200, thr, "E1q"))
        tags = ("E4c",) if GROUP == "c" else ("E4c", "E4r")
        rows = A6.est_all([U], m, ts, N, rngs[2], rngs[3], rngs[4], truth, thr, tags)
        J6 = 0; e6lo = e6hi = None
        r6 = rows[0]
        if r6.get("abstain", 1) == 0:
            e6lo = np.array([float(v) for v in r6["e6lo"].split(";")]); e6hi = np.array([float(v) for v in r6["e6hi"].split(";")]); J6 = len(e6lo)
        out += rows
        e78 = A8.run_e7_e8(U, m, ts, N, rngs[7], truth, thr, J6, e6lo, e6hi)
        for r in e78:
            r.pop("e7lo", None); r.pop("e7hi", None)
        out += e78
    for r in out:
        r.update(cell=c["id"], ds=ds, wpn=wpn)
        for k in ("e6lo", "e6hi", "e4lo", "e4hi"): r.pop(k, None)
    return out


def worker(task):
    cid, lo, hi = task; c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    if GROUP == "c":
        A9.patch_design2(100, 8)
    t0 = time.time(); sp = select()
    if "--cells" in args:
        want = {int(x) for x in args[args.index("--cells") + 1].split(",")}; sp = [s for s in sp if s["id"] in want]
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (8 if SMOKE else 200)
    rr = np.random.default_rng(2200)
    for s in sp:
        R.prep(s, rr)
        if GROUP == "c": s["thr"] = 10.0
    CELLS = {s["id"]: s for s in sp}
    sfx = "%s%s" % (GROUP, "_smoke" if SMOKE else "")
    keys = ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "thr")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(R.RES, f"cells_meta_amend9_{sfx}.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(sp)} cells, group {GROUP}, ncpu {NCPU}", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = sorted([(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[0]]["N"] * CELLS[t[0]]["wpn"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
            if k % 100 == 0: pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend9_{sfx}.partial.csv"), index=False)
    pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend9_{sfx}.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}_w{CELLS[k]["wpn"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, f"cell_cpu_seconds_amend9_{sfx}.csv"))
    print(f"done {len(rows_all)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__": main()
