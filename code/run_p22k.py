#!/usr/bin/env python3
"""P22 Amendment 11 runner: generalisation holdout on graph families never used in Amendments 1-10.
Usage:
  run_p22k.py exact                      -> results/amend11_cells_exact.json (exact curves only, NO simulation); prints its sha256
  run_p22k.py run [--smoke] [--ndata N] [--ncpu K] [--out NAME]
Cells (lazy walk 1/2, grid (24,34), same window rule/threshold as the main benchmark):
  crossover 2D->1D: C_L x C_12, L in 96,128,192,256,384  (fam aniso2_w12)
  crossover 3D->1D: C_L x C_6 x C_6, L in 48,64,96,128,192 (fam aniso3_w6)
  controls: ring C_L L in 384,768; torus2 L in 40,80; torus3 L in 18,26
wpn in 8,32,128; 200 datasets; seeds SeedSequence([22, 11000 + id, ds]); id = 100*wpn_index + cell_index (wpn_index 0,1,2 = 8,32,128).
Arms: E1q (est_E1b), E4c, E4r (A6.est_all), E7, E8 (A8.run_e7_e8); E6 row kept (feeds the ungated A1' alarm only)."""
import os, sys, time, json, hashlib
_argv = sys.argv; sys.argv = [_argv[0], "--wpn", "32"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A, p22amend6 as A6, p22amend8 as A8
import run_p22 as R
sys.argv = _argv

args = sys.argv[1:]
MODE = args[0]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 10
assert NCPU <= 10
WPNS = (8, 32, 128)
CELLS = None


def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))


def base_specs():
    S = []
    def add(fam, l, build, kind):
        S.append(dict(fam=fam, L=l, build=build, sizes=None, plateau=None, kind=kind, grid=(24, 34), lazy=0.5, note=""))
    for l in (96, 128, 192, 256, 384): add("aniso2_w12", l, lambda l=l: T(l, 12), "crossover")
    for l in (48, 64, 96, 128, 192): add("aniso3_w6", l, lambda l=l: T(l, 6, 6), "crossover")
    for l in (384, 768): add("ring", l, lambda l=l: T(l), "control")
    for l in (40, 80): add("torus2", l, lambda l=l: T(l, l), "control")
    for l in (18, 26): add("torus3", l, lambda l=l: T(l, l, l), "control")
    for i, s in enumerate(S): s["cidx"] = i
    return S


def specs():
    out = []
    for wi, w in enumerate(WPNS):
        for s in base_specs():
            s["wpn"] = w; s["id"] = 100 * wi + s["cidx"]; out.append(s)
    return out


def exact_table(S):
    rr = np.random.default_rng(2200); rows = []
    for s in S:
        R.prep(s, rr)
        rows.append(dict(cidx=s["cidx"], fam=s["fam"], L=s["L"], N=int(s["N"]), kind=s["kind"], exact_window_end=int(s["exact_window_end"]), exact_J=int(s["exact_J"]),
                         max_drop=None if not np.isfinite(s["max_drop"]) else round(float(s["max_drop"]), 6), d_start=None if not np.isfinite(s["d_start"]) else round(float(s["d_start"]), 6),
                         d_end=None if not np.isfinite(s["d_end"]) else round(float(s["d_end"]), 6), realised=bool(s["real"])))
    return rows


def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, 11000 + c["id"], ds])
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    out = [A.est_E1b([U], m, ts, N, rngs[6], truth, 2, 1.0, 200, thr, "E1q")]
    rows = A6.est_all([U], m, ts, N, rngs[2], rngs[3], rngs[4], truth, thr, ("E4c", "E4r"))
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


def main_exact():
    t0 = time.time(); rows = exact_table(base_specs())
    p = os.path.join(R.RES, "amend11_cells_exact.json")
    json.dump(rows, open(p, "w"), indent=1, sort_keys=True)
    print(f"exact table in {time.time()-t0:.0f}s; sha256 {sha(p)}")
    for r in rows: print(r)


def main_run():
    global CELLS
    t0 = time.time(); sp = specs()
    pj = os.path.join(R.RES, "amend11_cells_exact.json"); assert os.path.exists(pj), "run `exact` first"
    print(f"cells_exact sha256 at run start: {sha(pj)}", flush=True)
    rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    ex = {(r["cidx"]): r for r in json.load(open(pj))}
    for s in sp:
        e = ex[s["cidx"]]
        assert bool(s["real"]) == e["realised"] and s["exact_window_end"] == e["exact_window_end"], f"MISMATCH exact table cell {s['cidx']}"
    CELLS = {s["id"]: s for s in sp}
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (4 if SMOKE else 200)
    sfx = args[args.index("--out") + 1] if "--out" in args else ("smoke" if SMOKE else "full")
    keys = ("id", "cidx", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "thr", "real")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(R.RES, f"cells_meta_amend11_{sfx}.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(sp)} cells, nd {nd}, ncpu {NCPU}", flush=True)
    chunk = 2 if SMOKE else 5
    tasks = sorted([(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[0]]["N"] * CELLS[t[0]]["wpn"] * len(CELLS[t[0]]["ts_w"]))
    rows_all = []; ct = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            rows_all += rows; ct[cid] = ct.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks elapsed {time.time()-t0:.0f}s", flush=True)
            if k % 100 == 0: pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend11_{sfx}.partial.csv"), index=False)
    pd.DataFrame(rows_all).to_csv(os.path.join(R.RES, f"raw_amend11_{sfx}.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}_w{CELLS[k]["wpn"]}': v for k, v in ct.items()}).to_csv(os.path.join(R.RES, f"cell_cpu_seconds_amend11_{sfx}.csv"))
    print(f"done {len(rows_all)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main_exact() if MODE == "exact" else main_run()
