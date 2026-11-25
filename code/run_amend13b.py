#!/usr/bin/env python3
"""P22 Amendment 13(b) runner: random-graph arm (Watts-Strogatz small-world, supercritical 2D site-percolation clusters).
Usage:
  run_amend13b.py ens   [--smoke] [--ncpu K]   ensemble truth: R_ENS independent realisations per cell, eigen-decomposition, exact q_r(t), dq_r(t) on the grid; writes
                                               results/amend13b_cells_exact.json (cell table + ensemble-mean curves). Done BEFORE any dataset is generated; its sha256 is logged.
  run_amend13b.py data  [--smoke] [--ncpu K] [--ndata N] [--cells id,id]  datasets: one fresh graph realisation per (cell, ds) shared across the three budgets, walkers fresh per budget.
Arms: E1q (frozen p22amend.est_E1b), E4c (V0, frozen construction), E4c_V2 (independent split), E7 (frozen p22amend8.run_e7_e8; E8 row dropped).
Every arm is scored against BOTH the single-realisation exact d_s(t) (columns without suffix) and the ensemble-mean curve (columns with suffix _E).
Resumable (one pickle per task, atomic), manifests with sha256, run under nice, <= 6 workers."""
import os, sys, time, json, hashlib, pickle, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
ARGV = list(sys.argv); STEP = ARGV[1]
SMOKE = "--smoke" in ARGV
NCPU = int(ARGV[ARGV.index("--ncpu") + 1]) if "--ncpu" in ARGV else 6
assert NCPU <= 10
import numpy as np, pandas as pd
import multiprocessing as mp
sys.argv = [ARGV[0], "--wpn", "32"]
import p22lib as L, p22amend as A, p22amend8 as A8, p22amend13 as A13
import run_p22 as R
sys.argv = ARGV
RES = R.RES
SFX = "_smoke" if SMOKE else ""
R_ENS = 6 if SMOKE else 300
WPNS = (8, 32, 128)
THR = L.THR
NCAP = 4096
BIG = L.make_grid(24, 34, NCAP)
FROZEN = {"p22lib.py": "a4517f3e056246266a6eaae2db0063b0d5fb1e2fb903946b40eb4248b3d43b91",
          "p22amend.py": "f496087a9dad23ee92ca76f28724cc829f50b50c2011c36db86ac94d9b6806d8",
          "p22amend8.py": "6b0964bffc9e69c8ec73d1573e2798378486f3ca638d8e5ad5ceec4f0dc8ebf3",
          "p22honest.py": "871ed8b3e2a5049f23fa8c84c1a8f50e8a5556363f1c3e3f35096925ff192216",
          "run_p22.py": "4a1937dff73cc1cb6f52bb67d5f5e7ec8f0b0bcd9087053f05feba375bc9848d"}

# cell list (fixed in PLAN.md Amendment 13b): id, family, parameters
CELLS_DEF = [
    dict(id=0, fam="ws", name="ws_N1000_p0.05", par=dict(N=1000, k=4, p=0.05)),
    dict(id=1, fam="ws", name="ws_N2000_p0.01", par=dict(N=2000, k=4, p=0.01)),
    dict(id=2, fam="ws", name="ws_N2000_p0.03", par=dict(N=2000, k=4, p=0.03)),
    dict(id=3, fam="ws", name="ws_N4000_p0.05", par=dict(N=4000, k=4, p=0.05)),
    dict(id=4, fam="perc", name="perc_L40_p0.65", par=dict(Lside=40, p=0.65)),
    dict(id=5, fam="perc", name="perc_L64_p0.65", par=dict(Lside=64, p=0.65)),
    dict(id=6, fam="perc", name="perc_L40_p0.75", par=dict(Lside=40, p=0.75)),
    dict(id=7, fam="perc", name="perc_L64_p0.75", par=dict(Lside=64, p=0.75)),
    dict(id=8, fam="ws", name="ws_N4000_p0.02", par=dict(N=4000, k=4, p=0.02)),
    dict(id=9, fam="perc", name="perc_L64_p0.70", par=dict(Lside=64, p=0.70)),
]
NMAX = 4000
E7_ON = True
CT = None


def log(*a): print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def verify():
    here = os.path.dirname(os.path.abspath(__file__)); bad = [(f, h, sha(os.path.join(here, f))) for f, h in FROZEN.items() if sha(os.path.join(here, f)) != h]
    if bad: log("MISMATCH frozen files", bad); sys.exit(2)
    log("frozen estimator files byte-identical to the Amendment 8-12 runs")


def graph_realisation(cd, ss_child):
    """one realisation; rejects (and redraws, counted) realisations with N > NMAX (cannot occur for the registered parameters at the sizes used, asserted)."""
    rng = np.random.default_rng(ss_child)
    for tries in range(20):
        g = A13.make_random_graph(cd["fam"], cd["par"], rng)
        if g is not None and g.N <= NMAX: return g, tries
    raise RuntimeError("no realisation with N<=NMAX")


# ------------------------------------------------------------------------------ ensemble step
def ens_one(task):
    cid, r = task
    cd = CT_DEF[cid]
    ss = np.random.SeedSequence([22, 17000 + cid, r])
    g, tries = graph_realisation(cd, ss.spawn(1)[0])
    lam = g.lam; q, dq = A13.exact_curves(lam, g.N, BIG)
    return dict(cid=cid, r=r, N=int(g.N), q=q, dq=dq, tries=tries, lam_min=float(lam.min()), lam_max=float(lam.max()))


CT_DEF = {c["id"]: c for c in CELLS_DEF}


def step_ens():
    verify()
    defs = [c for c in CELLS_DEF if (not SMOKE or c["id"] in (0, 4))]
    d = os.path.join(RES, f"ens_amend13b{SFX}"); os.makedirs(d, exist_ok=True)
    tasks = [(c["id"], r) for c in defs for r in range(R_ENS)]
    todo = [t for t in tasks if not os.path.exists(os.path.join(d, f"{t[0]}_{t[1]}.pkl"))]
    log(f"[ens] {len(tasks)} realisations, {len(tasks)-len(todo)} done, ncpu {NCPU}")
    t0 = time.time()
    if todo:
        with mp.get_context("fork").Pool(NCPU) as pool:
            for k, out in enumerate(pool.imap_unordered(ens_one, sorted(todo, key=lambda t: -CT_DEF[t[0]]["par"].get("N", CT_DEF[t[0]]["par"].get("Lside", 0) ** 2)), chunksize=1)):
                p = os.path.join(d, f"{out['cid']}_{out['r']}.pkl"); pickle.dump(out, open(p + ".tmp", "wb")); os.replace(p + ".tmp", p)
                if k % 25 == 0: log(f"[ens] {k+1}/{len(todo)} {time.time()-t0:.0f}s")
    rows = []
    for c in defs:
        outs = [pickle.load(open(os.path.join(d, f"{c['id']}_{r}.pkl"), "rb")) for r in range(R_ENS)]
        N = np.array([o["N"] for o in outs]); Q = np.stack([o["q"] for o in outs]); DQ = np.stack([o["dq"] for o in outs])
        qbar = Q.mean(0); dqbar = DQ.mean(0); dbar = -2 * BIG * dqbar / qbar
        dr = -2 * BIG[None, :] * DQ / Q                        # per-realisation d_s(t)
        Nq = (N[:, None] * Q).mean(0)                           # mean over realisations of N_r q_r(t)
        cap_idx = np.nonzero(Nq >= THR - 1.5)[0]; cap = int(BIG[cap_idx[-1]]) if len(cap_idx) else int(BIG[0])
        ts_w = L.make_grid(24, 34, min(cap, NCAP)); K = len(ts_w)
        run = []
        for i in range(K):
            if Nq[i] >= THR: run.append(i)
            else: break
        run = np.array(run, dtype=int)
        if len(run) >= 2:
            ewe = int(BIG[run[-1]]); ej = len(run); md = float(dbar[run[0]] - dbar[run].min()); mr = float(dbar[run].max() - dbar[run[0]]); mad = float(np.abs(dbar[run] - dbar[run[0]]).max())
            sdd = float(dr[:, run].std(0, ddof=1).mean()); ds_, de_ = float(dbar[run[0]]), float(dbar[run[-1]])
        else:
            ewe, ej, md, mr, mad, sdd, ds_, de_ = 0, len(run), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        rows.append(dict(id=c["id"], name=c["name"], fam=c["fam"], par=c["par"], R_ens=R_ENS, N_mean=float(N.mean()), N_min=int(N.min()), N_max=int(N.max()), N_sd=float(N.std(ddof=1)),
                         ts_w=[int(t) for t in ts_w], exact_window_end=ewe, exact_J=ej, max_drop=md, max_rise=mr, max_abs_delta=mad, d_start=ds_, d_end=de_,
                         sd_d_between_graphs_mean_over_window=sdd, crossover_like=bool(np.isfinite(mad) and mad >= 0.6), flat_like=bool(np.isfinite(mad) and mad < 0.25),
                         dbar=[float(v) for v in dbar[:K]], qbar=[float(v) for v in qbar[:K]], tries_max=int(max(o["tries"] for o in outs))))
    p = os.path.join(RES, f"amend13b_cells_exact{SFX}.json"); json.dump(rows, open(p, "w"), indent=1, sort_keys=True)
    log("ensemble table written, sha256", sha(p))
    for r in rows:
        log("  ", {k: r[k] for k in ("id", "name", "N_mean", "N_min", "N_max", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "sd_d_between_graphs_mean_over_window", "crossover_like", "flat_like")})


# ------------------------------------------------------------------------------ data step
def run_dataset(c, ds):
    cd = CT_DEF[c["id"]]
    ss = np.random.SeedSequence([22, 16000 + c["id"], ds])
    g, tries = graph_realisation(cd, ss.spawn(1)[0])
    t0 = time.time(); lam = g.lam; teig = time.time() - t0
    ts = np.array(c["ts_w"], dtype=np.int64); q, dq = A13.exact_curves(lam, g.N, ts)
    dS = -2 * ts * dq / q; dE = np.array(c["dbar"])
    fS = lambda t, dS=dS: dS[:len(np.atleast_1d(t))]
    fE = lambda t, dE=dE: dE[:len(np.atleast_1d(t))]
    qf = lambda t, q=q: q[:len(np.atleast_1d(t))]
    truthS = dict(q_exact=qf, d_exact=fS); truthE = dict(q_exact=lambda t, qb=np.array(c["qbar"]): qb[:len(np.atleast_1d(t))], d_exact=fE)
    # exact excursion within the realisation's own exact window (diagnostic for the false-change scoring)
    Nq = g.N * q; run = []
    for i in range(len(ts)):
        if Nq[i] >= THR: run.append(i)
        else: break
    ex = dict(N_r=int(g.N), tries=tries, eig_s=teig, exact_J_r=len(run),
              exc_S=float(np.abs(dS[run] - dS[run[0]]).max()) if len(run) >= 2 else np.nan, exc_E=float(np.abs(dE[run] - dE[run[0]]).max()) if len(run) >= 2 else np.nan,
              gap_SE=float(np.abs(dS[run] - dE[run]).max()) if len(run) >= 1 else np.nan)
    rows = []
    for wpn in WPNS:
        ssw = np.random.SeedSequence([22, 16000 + c["id"], ds, wpn]); rr = [np.random.default_rng(s) for s in ssw.spawn(8)]
        Cc, m = g.sample_blocks(rr[0], wpn * g.N, ts, L.NBLOCK); U = Cc / m
        out = []
        r1 = A.est_E1b([U], m, ts, g.N, np.random.default_rng(np.random.SeedSequence([22, 16000 + c["id"], ds, wpn, 1])), truthS, 2, 1.0, 200, THR, "E1q")
        r1e = A.est_E1b([U], m, ts, g.N, np.random.default_rng(np.random.SeedSequence([22, 16000 + c["id"], ds, wpn, 1])), truthE, 2, 1.0, 200, THR, "E1q")
        if r1.get("abstain", 1) == 0 and r1e.get("abstain", 1) == 0:
            for k in ("cov_simult", "cov_point", "delta_cov", "bias", "mae"):
                if k in r1e: r1[k + "_E"] = r1e[k]
        out.append(r1)
        out += A13.est_E4v2([U], m, ts, g.N, rr[3], rr[4], [("S", fS), ("E", fE)], None, THR, tags=("E4c",), variants=("V0", "V2"), keep_bands=False)
        if E7_ON:
            try:
                e78 = A8.run_e7_e8(U, m, ts, g.N, rr[5], truthS, THR)
                r7 = e78[0]; lo7 = r7.pop("e7lo", ""); hi7 = r7.pop("e7hi", "")
                if r7.get("nowin", 0) == 0 and r7.get("infeasible", 0) == 0 and lo7:
                    Jw = r7["J"]; lo = np.array([float(v) for v in lo7.split(";")])[:Jw]; hi = np.array([float(v) for v in hi7.split(";")])[:Jw]
                    cov = (dE[:Jw] >= lo - 1e-12) & (dE[:Jw] <= hi + 1e-12); r7["cov_simult_E"] = int(cov.all()); r7["cov_point_E"] = float(cov.mean())
                elif r7.get("nowin", 0) == 0:
                    r7["cov_simult_E"] = 0; r7["cov_point_E"] = 0.0
                out.append(r7)
            except Exception as e:
                out.append(dict(est="E7", abstain=1, why="error:" + repr(e)[:80]))
        for r in out: r.update(cell=c["id"], ds=ds, wpn=wpn, **ex)
        rows += out
    return rows


def worker(task):
    key, (cid, lo, hi) = task
    try:
        rows = []
        for ds in range(lo, hi): rows += run_dataset(CT[cid], ds)
        return key, rows, None
    except Exception:
        return key, None, traceback.format_exc()[-800:]


def step_data():
    global CT
    verify()
    p = os.path.join(RES, f"amend13b_cells_exact{SFX}.json"); assert os.path.exists(p), "run ens first"
    log("amend13b exact/ensemble table sha256 at data start:", sha(p))
    tab = json.load(open(p)); CT = {r["id"]: r for r in tab}
    want = [int(v) for v in ARGV[ARGV.index("--cells") + 1].split(",")] if "--cells" in ARGV else sorted(CT)
    nd = int(ARGV[ARGV.index("--ndata") + 1]) if "--ndata" in ARGV else (3 if SMOKE else 300)
    name = "amend13b" + SFX + ("_c" + "".join(map(str, want)) if "--cells" in ARGV else "")
    d = os.path.join(RES, f"items_{name}"); os.makedirs(d, exist_ok=True)
    chunk = 1 if SMOKE else 2
    tasks = sorted([(f"{c}_{lo}_{min(lo+chunk, nd)}", (c, lo, min(lo + chunk, nd))) for c in want for lo in range(0, nd, chunk)], key=lambda t: -CT[t[1][0]]["N_mean"])
    todo = [t for t in tasks if not os.path.exists(os.path.join(d, t[0] + ".pkl"))]
    log(f"[{name}] {len(want)} cells, {len(tasks)} tasks, {len(tasks)-len(todo)} done, ncpu {NCPU}")
    errors = []; t0 = time.time()
    if todo:
        with mp.get_context("fork").Pool(NCPU) as pool:
            for k, (key, rows, err) in enumerate(pool.imap_unordered(worker, todo, chunksize=1)):
                if err: errors.append(key); log("ERROR", key, err); continue
                tmp = os.path.join(d, key + ".tmp"); pickle.dump(rows, open(tmp, "wb")); os.replace(tmp, os.path.join(d, key + ".pkl"))
                if k % 10 == 0: log(f"[{name}] {k+1}/{len(todo)} new, {time.time()-t0:.0f}s")
    rows = []
    for key, _ in tasks:
        pp = os.path.join(d, key + ".pkl")
        if os.path.exists(pp): rows += pickle.load(open(pp, "rb"))
    csv = os.path.join(RES, f"raw_{name}.csv"); pd.DataFrame(rows).to_csv(csv, index=False)
    items = sorted(f for f in os.listdir(d) if f.endswith(".pkl")); here = os.path.dirname(os.path.abspath(__file__))
    man = dict(name=name, tasks_expected=len(tasks), items_present=len(items), errors=errors, complete=(len(items) == len(tasks) and not errors), n_rows=len(rows),
               csv_sha256=sha(csv), items_sha256=hashlib.sha256("".join(f"{f} {sha(os.path.join(d, f))}\n" for f in items).encode()).hexdigest(), cells_exact_sha256=sha(p),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"), frozen=FROZEN, p22amend13_sha256=sha(os.path.join(here, "p22amend13.py")), runner_sha256=sha(os.path.abspath(__file__)), numpy=np.__version__, e7_on=E7_ON)
    json.dump(man, open(os.path.join(RES, f"manifest_{name}.json"), "w"), indent=1, sort_keys=True)
    log(f"[{name}] manifest complete={man['complete']} rows {len(rows)} csv {man['csv_sha256'][:16]}")


if __name__ == "__main__":
    {"ens": step_ens, "data": step_data}[STEP]()
