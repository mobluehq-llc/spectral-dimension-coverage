#!/usr/bin/env python3
"""P22 runner. Usage: .venv/bin/python run_p22.py [--smoke] [--only fam,fam_L,...] [--ndata N]. Outputs to results/."""
import os, sys, time
import numpy as np
import pandas as pd
import multiprocessing as mp
import p22lib as L

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)
SMOKE = "--smoke" in sys.argv
NDATA = 200
ONLY = None
for i, a in enumerate(sys.argv):
    if a == "--ndata": NDATA = int(sys.argv[i + 1])
    if a == "--only": ONLY = sys.argv[i + 1].split(",")
WPN, PROBES, NW_CAP = 32, 32, 4096
NCPU = 9


def torus(dims, lazy=0.5): return L.Graph("torus", lazy, dims=tuple(dims))


def cell_specs():
    S = []
    def add(fam, Lval, build, sizes, plateau=None, kind="control", grid=(24, 34), lazy=0.5, note=""):
        S.append(dict(fam=fam, L=Lval, build=build, sizes=sizes, plateau=plateau, kind=kind, grid=grid, lazy=lazy, note=note))
    def fam3(mk, l):  # three sizes of the long axis (isotropic families scale all axes)
        return [lambda: mk(l), lambda: mk(int(.75 * l)), lambda: mk(l // 2)]
    for l in (128, 256, 512, 1024):
        add("ring", l, lambda l=l: torus((l,)), fam3(lambda x: torus((x,)), l), 1.0)
    for l in (24, 32, 48, 64, 96):
        add("torus2", l, lambda l=l: torus((l, l)), fam3(lambda x: torus((x, x)), l), 2.0)
    for l in (16, 20, 24, 28, 32):
        add("torus3", l, lambda l=l: torus((l,) * 3), fam3(lambda x: torus((x,) * 3), l), 3.0)
    for l in (32, 48, 64, 96, 128, 192, 256):
        add("aniso2_w16", l, lambda l=l: torus((l, 16)), fam3(lambda x: torus((x, 16)), l), None, "crossover")
    for l in (32, 48, 64, 96, 128, 192):
        add("aniso3_w8", l, lambda l=l: torus((l, 8, 8)), fam3(lambda x: torus((x, 8, 8)), l), None, "crossover")
    for n in (4, 5, 6, 7):
        add("gasket", n, lambda n=n: L.gasket(n), None, None, "descriptive")
    for N in (512, 2048):
        add("rrg3", N, None, None, None, "descriptive", grid=(2, 3))
    for fam, l, dims, pl in (("ring", 256, (256,), 1.0), ("torus2", 48, (48, 48), 2.0)):
        for lazy, grid, tag in ((0.0, (24, 34), "nonlazy_evengrid"), (0.0, (25, 35), "nonlazy_oddevengrid"), (0.5, (25, 35), "lazy_oddevengrid")):
            add(fam, l, lambda dims=dims, lazy=lazy: torus(dims, lazy), None, pl, "ablation", grid=grid, lazy=lazy, note=tag)
    for i, s in enumerate(S): s["id"] = i
    return S


def prep(c, rng_ref):
    """exact truth objects + grids for a cell (parent process)."""
    if c["fam"] == "rrg3":
        N = c["L"]
        refs = [L.rrg(N, 3, rng_ref) for _ in range(24 if N <= 512 else 12)]
        lam = np.concatenate([g.lam for g in refs])
        g0 = L.Graph("step", 0.5, nbr=refs[0].nbr, deg=refs[0].deg, lam=lam); g0.N = N
        c["graph"] = g0; c["thr"] = 1.3
    else:
        c["graph"] = c["build"](); c["thr"] = L.THR
    g = c["graph"]; c["N"] = g.N
    a0, b0 = c["grid"]
    big = L.make_grid(a0, b0, 3_000_000 if g.kind == "torus" else NW_CAP)
    if c["kind"] == "ablation" and c["lazy"] == 0.0:
        # truth for the non-lazy walk is replaced by the plateau (no analytic continuation in t)
        c["q_exact"] = lambda t, p=c["plateau"]: np.asarray(t, float) ** (-p / 2)
        c["d_exact"] = lambda t, p=c["plateau"]: np.full(len(np.atleast_1d(t)), p)
        c["cap_e2"] = int({"ring": 4000, "torus2": 1200}[c["fam"]]); c["ts_w"] = L.make_grid(a0, b0, c["cap_e2"]); c["ts_e2"] = c["ts_w"]
        c.update(exact_window_end=0, exact_J=0, max_drop=np.nan, max_rise=np.nan, d_start=np.nan, d_end=np.nan, max_abs_delta=np.nan, real=False)
        return c
    q, d = g.exact(big)
    Nq = g.N * q
    cap_idx = np.nonzero(Nq >= (c["thr"] - 1.5 if c["fam"] != "rrg3" else 1.05))[0]
    cap_e2 = big[cap_idx[-1]] if len(cap_idx) else big[0]
    cap_e2 = min(cap_e2, NW_CAP * 4 if g.kind == "torus" else NW_CAP)
    cap_w = cap_e2
    if c["sizes"] is not None and c["kind"] != "ablation":
        c["size_graphs"] = [f() for f in c["sizes"]]
        capsw = []
        for sg in c["size_graphs"]:
            qs, _ = sg.exact(big)
            idx = np.nonzero(sg.N * qs >= 1.08)[0]
            capsw.append(big[idx[-1]] if len(idx) else big[0])
        c["cap_sizes"] = capsw
        cap_w = max(cap_w, capsw[0])
    c["cap_e2"] = int(cap_e2)
    c["ts_w"] = L.make_grid(a0, b0, int(cap_w)); c["ts_e2"] = L.make_grid(a0, b0, int(cap_e2))
    c["q_exact"] = lambda t, g=g: g.exact(t)[0]
    c["d_exact"] = lambda t, g=g: g.exact(t)[1]
    run = []
    for i in range(len(big)):
        if big[i] < a0: continue
        if Nq[i] >= c["thr"]: run.append(i)
        else: break
    run = np.array(run, dtype=int)
    if len(run) >= 2:
        c["exact_window_end"] = int(big[run[-1]]); c["exact_J"] = len(run)
        c["max_drop"] = float((d[run[0]] - d[run]).max()); c["max_rise"] = float((d[run] - d[run[0]]).max())
        c["d_start"] = float(d[run[0]]); c["d_end"] = float(d[run[-1]])
        c["max_abs_delta"] = float(np.abs(d[run] - d[run[0]]).max())
    else:
        c.update(exact_window_end=0, exact_J=len(run), max_drop=np.nan, max_rise=np.nan, d_start=np.nan, d_end=np.nan, max_abs_delta=np.nan)
    c["real"] = bool(c["kind"] == "crossover" and c["max_drop"] >= 0.6)
    return c


CELLS = None


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(6)]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts_w, ts_e2 = c["ts_w"], c["ts_e2"]
    if c["fam"] == "rrg3":
        G = 4; graphs = [L.rrg(c["N"], 3, rngs[0]) for _ in range(G)]
    else:
        G = 1; graphs = [c["graph"]]
    N = c["N"]; nb = L.NBLOCK // G
    U = []; mw = None
    for g in graphs:
        C, m = g.sample_blocks(rngs[1], WPN * N // G, ts_w, nb); U.append(C / m); mw = m
    out = [L.est_E1(U, mw, ts_w, N, rngs[2], truth, True, 200, c["thr"], "E1"),
           L.est_E1(U, mw, ts_w, N, rngs[2], truth, False, 200, c["thr"], "E1x")]
    out += L.est_E0(U, mw, ts_w, N, rngs[3], truth, 200, c["thr"])
    X = [g.sample_probes(rngs[4], PROBES // G, ts_e2) for g in graphs]
    out += L.est_E2(X, ts_e2, N, rngs[3], truth, 400, c["thr"])
    if c.get("size_graphs") is not None:
        fam = [(U[0], mw, ts_w, N)]
        for sg, cap in zip(c["size_graphs"][1:], c["cap_sizes"][1:]):
            tsi = L.make_grid(*c["grid"], int(cap))
            Ci, mi = sg.sample_blocks(rngs[5], WPN * sg.N, tsi, L.NBLOCK)
            fam.append((Ci / mi, mi, tsi, sg.N))
        out.append(L.est_E3(fam, rngs[5], c["plateau"]))
    for r in out:
        r.update(cell=c["id"], ds=ds)
    return out


def worker(task):
    cid, lo, hi = task
    c = CELLS[cid]; t0 = time.time(); rows = []
    for ds in range(lo, hi):
        rows += run_dataset(c, ds)
    return cid, rows, time.time() - t0


def main():
    global CELLS
    t0 = time.time()
    specs = cell_specs()
    if ONLY: specs = [s for s in specs if s["fam"] in ONLY or f'{s["fam"]}_{s["L"]}' in ONLY]
    rr = np.random.default_rng(2200)
    for s in specs: prep(s, rr)
    CELLS = {s["id"]: s for s in specs}
    keys = ("id", "fam", "L", "N", "kind", "lazy", "note", "real", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "cap_e2")
    pd.DataFrame([{k: s[k] for k in keys} | {"grid": str(s["grid"])} for s in specs]).to_csv(os.path.join(RES, "cells_meta.csv"), index=False)
    print(f"prep {time.time()-t0:.0f}s, {len(specs)} cells", flush=True)
    nd = 8 if SMOKE else NDATA
    chunk = 4 if SMOKE else 10
    tasks = [(s["id"], lo, min(lo + chunk, nd)) for s in specs for lo in range(0, nd, chunk)]
    tasks.sort(key=lambda t: -CELLS[t[0]]["N"] * (len(CELLS[t[0]]["ts_w"]) + CELLS[t[0]]["cap_e2"] / 20))
    allrows = []; cell_time = {}
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, rows, dt) in enumerate(pool.imap_unordered(worker, tasks)):
            allrows += rows; cell_time[cid] = cell_time.get(cid, 0) + dt
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} tasks  elapsed {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(allrows)
    df.to_csv(os.path.join(RES, "raw_smoke.csv" if SMOKE else "raw.csv"), index=False)
    pd.Series({f'{CELLS[k]["fam"]}_{CELLS[k]["L"]}_{CELLS[k]["note"]}': v for k, v in cell_time.items()}).to_csv(os.path.join(RES, "cell_cpu_seconds_smoke.csv" if SMOKE else "cell_cpu_seconds.csv"))
    print(f"done {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
