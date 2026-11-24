#!/usr/bin/env python3
"""P22 Amendment 10 runner. Usage:
  run_p22j.py validate                     : 10b validation on 2 tiny cases (finer atom grids stay inside the certified outer band) + exact-q containment
  run_p22j.py select                       : write results/amend10b_selection.json (sha256 logged) BEFORE data generation
  run_p22j.py b [--ncpu K] [--smoke]       : 10b certified E7/E7j on the selected datasets (needs selection.json)
  run_p22j.py exactq                       : 10b exact-q versions of the 7 checked cells
  run_p22j.py windows                      : 10c fixed windows from exact curves -> results/amend10c_windows.json (sha256 logged)
  run_p22j.py c [--ncpu K] [--smoke] [--cells id,..]  : 10c fixed-window budget comparison (needs windows.json)
Walker data: SeedSequence([22, id, ds]).spawn(9)[1] and box from child [7] exactly as run_p22h/i (10b); 10c uses fresh seeds [22, id, 3000+ds]."""
import os, sys, time, json, hashlib
_argv = sys.argv; sys.argv = [_argv[0], "--wpn", "32"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A, p22amend6 as A6, p22amend8 as A8, p22amend9 as A9, p22amend10 as A10
import run_p22 as R, run_p22e as RE
import run_p22g as RG
sys.argv = _argv
args = sys.argv[1:]
MODE = args[0] if args else ""
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 12
assert NCPU <= 16
RES = R.RES
CELLS = {}


def T(*d): return L.Graph("torus", 0.5, dims=tuple(d))
def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()


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


def get_cell(cid):
    if cid not in CELLS:
        for s in core_specs():
            if s["id"] == cid:
                R.prep(s, np.random.default_rng(2200)); CELLS[cid] = s; break
    return CELLS[cid]


def box_for(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    ts, N, wpn = c["ts_w"], c["N"], c["wpn"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    nb, K = U.shape; Wt = nb * m; p = U.mean(0); cnt = p * Wt
    Jw = L.window_len(cnt, p, N, c["thr"])
    qhat, lo, hi, Ps, binfo = A8.make_box(U, m, rngs[7])
    return Jw, lo, hi


# --------------------------------------------------------------- validation
def validate():
    A10_out = []
    for cid, ds in ((6, 0), (1, 0)):          # torus2 L=32 (N=1024), ring 256: tiny
        c = get_cell(cid); ts, N = c["ts_w"], c["N"]
        Jw, lo, hi = box_for(c, ds); idx = list(range(Jw)); dtrue = c["d_exact"](ts[:Jw])
        print(f"cell {cid} {c['fam']} L={c['L']} N={N} Jw={Jw} K={len(ts)} tmax={ts[-1]}", flush=True)
        t0 = time.time(); bB = A10.band_B(ts, lo, hi, N, idx); bA = A10.band_A(ts, lo, hi, N, idx); print(f"  certified bands in {time.time()-t0:.1f}s, status B {bB[2]} A {bA[2]}", flush=True)
        for n in (400, 1597, 6385, 12769):
            lo_f, hi_f, st = A8.lp_band(ts, lo, hi, N, idx, float(ts[-1]), n)
            if st != 0: print(f"  atoms {n}: LP status {st}"); continue
            inB = bool(((bB[0] <= lo_f + 1e-6) & (bB[1] >= hi_f - 1e-6)).all()); inA = bool(((bA[0] <= lo_f + 1e-6) & (bA[1] >= hi_f - 1e-6)).all())
            print(f"  fine grid {n:6d} atoms: band inside outer B: {inB}  inside outer A: {inA}  max inward violation B {max((bB[0]-lo_f).max(), (hi_f-bB[1]).max()):+.2e}", flush=True)
            A10_out.append((cid, n, inB, inA))
        print("  B band", np.round(bB[0], 3), np.round(bB[1], 3)); print("  A band", np.round(bA[0], 3), np.round(bA[1], 3)); print("  exact d", np.round(dtrue, 3))
        # random continuous measures (atoms anywhere) that satisfy the box, via a random-support LP at 3 random supports
        rng = np.random.default_rng(10)
        for rep in range(3):
            s_atoms = np.concatenate([[0.0], np.sort(np.exp(rng.uniform(np.log(0.02 / ts[-1]), 0, 700))), [np.inf]])
            S = A8.lp_setup(ts, lo, hi, s_atoms, N); res = [A8.lp_one(S, j) for j in idx]
            if max(r[2] for r in res) != 0: print(f"  random-support {rep}: infeasible"); continue
            lo_r = np.array([r[0] for r in res]); hi_r = np.array([r[1] for r in res])
            print(f"  random off-grid support {rep}: inside B {bool(((bB[0] <= lo_r + 1e-6) & (bB[1] >= hi_r - 1e-6)).all())} inside A {bool(((bA[0] <= lo_r + 1e-6) & (bA[1] >= hi_r - 1e-6)).all())}", flush=True)
    ok = all(x[2] for x in A10_out)
    print("VALIDATION (finer / random grids stay inside outer B on both cases):", "PASS" if ok else "FAIL")


# --------------------------------------------------------------- selection
def select():
    d8 = pd.read_csv(os.path.join(RES, "raw_amend8_core.csv"), low_memory=False); e7 = d8[d8.est == "E7"]
    sens = e7[e7.gate_fail.fillna(0) == 1][["cell", "ds"]]
    a9 = pd.read_csv(os.path.join(RES, "raw_amend9_a.csv"), low_memory=False)
    wh = a9[a9.gate_flag.fillna(0) > 0][["cell", "ds"]]
    stable = e7[(e7.gate_fail.fillna(0) == 0) & (e7.infeasible.fillna(0) == 0) & (e7.nowin.fillna(0) == 0)][["cell", "ds"]].sort_values(["cell", "ds"]).reset_index(drop=True)
    rs = np.random.default_rng(2210); pick = np.sort(rs.choice(len(stable), 300, replace=False)); st = stable.iloc[pick]
    minl = [(c, ds) for c in (702, 26, 703) for ds in range(200)]
    sel = dict(sens=[[int(a), int(b)] for a, b in sens.values], withheld=[[int(a), int(b)] for a, b in wh.values], stable300=[[int(a), int(b)] for a, b in st.values],
               minL3d=[[a, b] for a, b in minl], seed_stable=2210, n_stable_pool=int(len(stable)),
               note="sens: E7 rows of raw_amend8_core with gate_fail==1; withheld: raw_amend9_a gate_flag>0; stable300: default_rng(2210).choice over feasible gate-passing E7 rows sorted by (cell,ds); minL3d: all 200 datasets of the 3D cells L=160,192,256 (ids 702,26,703) added so that the 3D min L = 160 bar can be assessed")
    p = os.path.join(RES, "amend10b_selection.json"); json.dump(sel, open(p, "w"), indent=1)
    print("selection written", p, "sha256", sha(p), {k: len(v) for k, v in sel.items() if isinstance(v, list)})


# --------------------------------------------------------------- 10b run
def run_one(task):
    cid, ds, grp = task; c = get_cell(cid); t0 = time.time()
    ts, N = c["ts_w"], c["N"]
    Jw, lo, hi = box_for(c, ds)
    row = dict(cell=cid, ds=ds, group=grp, fam=c["fam"], L=c["L"], kind=c["kind"], J=Jw)
    if Jw < 2: row["nowin"] = 1; return row
    dtrue = c["d_exact"](ts[:Jw])
    try:
        row.update(A10.run_cert(ts, lo, hi, N, Jw, dtrue))
    except Exception as ex:
        row["error"] = repr(ex)[:200]
    row["seconds"] = time.time() - t0
    return row


def run_b():
    sel = json.load(open(os.path.join(RES, "amend10b_selection.json")))
    print("selection sha256", sha(os.path.join(RES, "amend10b_selection.json")), flush=True)
    tasks = {}
    for g in ("sens", "withheld", "stable300", "minL3d"):
        for cid, ds in sel[g]: tasks.setdefault((cid, ds), []).append(g)
    tl = [(cid, ds, "+".join(g)) for (cid, ds), g in tasks.items()]
    if SMOKE: tl = tl[:6] + [t for t in tl if t[0] in (702,)][:2]
    tl.sort(key=lambda t: -t[0] % 7)               # mix cells
    for cid in sorted({t[0] for t in tl}): get_cell(cid)   # prep in parent (forked workers inherit)
    print(f"{len(tl)} datasets, ncpu {NCPU}", flush=True)
    t0 = time.time(); rows = []; sfx = "_smoke" if SMOKE else ""
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, r in enumerate(pool.imap_unordered(run_one, tl, chunksize=1)):
            rows.append(r)
            if k % 10 == 0: print(f"{k+1}/{len(tl)} elapsed {time.time()-t0:.0f}s", flush=True)
            if k % 25 == 0: pd.DataFrame(rows).to_csv(os.path.join(RES, f"raw_amend10b{sfx}.partial.csv"), index=False)
    pd.DataFrame(rows).to_csv(os.path.join(RES, f"raw_amend10b{sfx}.csv"), index=False)
    print(f"done {len(rows)} rows in {time.time()-t0:.0f}s", flush=True)


# --------------------------------------------------------------- exact-q cells (the 7 of noisefree_e7)
EXQ = [("ring", 256, (256,)), ("torus2", 24, (24, 24)), ("torus3", 16, (16,) * 3), ("aniso2_w16", 128, (128, 16)), ("aniso3_w8", 32, (32, 8, 8)),
       ("aniso2_w16", 256, (256, 16)), ("aniso3_w8", 96, (96, 8, 8))]


def exq_one(spec):
    fam, l, dims = spec; g = L.Graph("torus", 0.5, dims=dims); N = g.N
    big = L.make_grid(24, 34, 3_000_000); q, d = g.exact(big)
    cap = big[np.nonzero(N * q >= 4.5)[0][-1]]; ts = L.make_grid(24, 34, int(cap)); q, d = g.exact(ts)
    Jw = L.window_len(N * q * 10 ** 3, q, N, 6.0)
    out = dict(cell=f"{fam}_{l}", J=Jw, N=N)
    for tol in (1e-9,):
        lo = q * (1 - tol); hi = q * (1 + tol)
        r = A10.run_cert(ts, lo, hi, N, Jw, d[:Jw], do_j=True, do_A=True)
        out.update(r); out["tol"] = tol
        out["true_drop"] = float((d[:Jw] - d[0]).min()); out["true_rise"] = float((d[:Jw] - d[0]).max())
    return out


def run_exq():
    t0 = time.time()
    with mp.get_context("fork").Pool(min(NCPU, 7)) as pool: rows = pool.map(exq_one, EXQ)
    pd.DataFrame(rows).to_csv(os.path.join(RES, "raw_amend10b_exactq.csv"), index=False)
    print(pd.DataFrame(rows)[["cell", "J", "o_status", "B_status", "o_cov", "B_cov", "B_width", "o_width", "A_status", "A_cov", "A_nopen", "B_contains_orig", "true_drop", "Bj_neg", "Bj_amb", "Bj_certnon"]].to_string())
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    {"validate": validate, "select": select, "b": run_b, "exactq": run_exq}.get(MODE, lambda: print("10c modes below"))()
