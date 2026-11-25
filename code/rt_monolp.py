#!/usr/bin/env python3
"""Red-team R1 (PREREG_REDTEAM_GEMINI.md): MonoLP (E7 / E7j) with exact structural constraints q(1)=1/2 (C1) and q(2) exact (C2).
Frozen files are imported, never edited; constraints are injected by wrapping A8.lp_setup / A8.atoms_grid at run time.
Usage: rt_monolp.py [--smoke] [--ndata N] [--ncpu K] [--cells id,id] [--only2d]
Data: group a of run_p22i (SeedSequence([22, id, ds]), walkers child [1], box bootstrap child [7]); the box is built ONCE per dataset and shared by all variants."""
import os, sys, time, json
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend8 as A8, p22amend9 as A9
import run_p22 as R, run_p22i as RI
sys.argv = _argv
args = sys.argv[1:]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 5
assert NCPU <= 5
EPS = 1e-9
VAR = dict(name="B0", ext=False, cons=())          # cons: tuple of (m, rhs): sum w e^{-m s} = rhs
_orig_atoms, _orig_setup = A8.atoms_grid, A8.lp_setup
EXT = np.exp(np.linspace(0.0, np.log(30.0), 61))[1:]    # 60 atoms on (1, 30]


def atoms_ext(tmax, n):
    g = _orig_atoms(tmax, n)                 # [0, grid..., inf]
    return np.concatenate([g[:-1], EXT, [np.inf]]) if VAR["ext"] else g


def setup_c(ts, lo, hi, s_atoms, N):
    S = _orig_setup(ts, lo, hi, s_atoms, N)
    if not VAR["cons"]:
        return S
    n = S["n"]; rows = []
    sf = np.where(np.isinf(s_atoms), np.inf, s_atoms)
    for m_, rhs in VAR["cons"]:
        c = np.exp(-m_ * sf)                 # inf -> 0
        a = np.zeros(n + 1); a[:n] = c; a[n] = -rhs
        rows.append(a); rows.append(-a)
    S["Aub"] = np.vstack([S["Aub"], np.array(rows)]); S["bub"] = np.concatenate([S["bub"], np.full(len(rows), EPS)])
    return S


A8.atoms_grid = atoms_ext; A8.lp_setup = setup_c


def secant_bounds(ts, lo, hi, Jw):
    """certified (from box only, log-convexity of q in t) lower bound on d(t0) and upper bound on d(t_j); best certified drop over window."""
    t = ts.astype(float); f_lo = np.log(np.maximum(lo, 1e-300)); f_hi = np.log(hi)
    j0 = 0
    Lb = max(-2 * t[0] * (f_hi[k] - f_lo[0]) / (t[k] - t[0]) for k in range(1, len(t)))
    Ub = np.array([min(-2 * t[j] * (f_lo[j] - f_hi[a]) / (t[j] - t[a]) for a in range(j)) if j > 0 else np.inf for j in range(Jw)])
    drop = Lb - Ub[1:]
    return float(Lb), float(drop.max()) if len(drop) else np.nan


def variants(c):
    g = c["graph"]; q2 = float(g.exact(np.array([2.0]))[0][0])
    v = [dict(name="B0", ext=False, cons=()),
         dict(name="C1", ext=True, cons=((1.0, 0.5),)),
         dict(name="C2", ext=True, cons=((1.0, 0.5), (2.0, q2)))]
    if c["fam"] == "aniso2_w16":
        v.insert(1, dict(name="B0x", ext=True, cons=()))
    return v


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    nb, K = U.shape; Wt = nb * m; p = U.mean(0); cnt = p * Wt
    Jw = L.window_len(cnt, p, N, thr)
    base = dict(cell=c["id"], fam=c["fam"], L=c["L"], kind=c["kind"], ds=ds, Jw=Jw)
    if Jw < 2:
        return [dict(base, variant="B0", nowin=1)]
    qhat, lo, hi, Ps, binfo = A8.make_box(U, m, rngs[7])
    dtrue = c["d_exact"](ts[:Jw])
    sb = secant_bounds(ts, lo, hi, Jw)
    out = []
    for v in variants(c):
        VAR.update(v)
        r = dict(base, variant=v["name"], nowin=0, c_box=binfo["c_box"], n_cp=binfo["n_cp"], sec_lb0=sb[0], sec_drop=sb[1])
        t0 = time.time()
        dmin, dmax, st, gi = A8.e7_gated(ts, lo, hi, N, list(range(Jw)))
        r.update(e7_status=st, natoms=gi["natoms"], gate_fail=gi["gate_fail"], t_e7=time.time() - t0)
        if st != 0:
            r.update(infeasible=1); out.append(r); continue
        lo7, hi7 = dmin[:Jw], dmax[:Jw]
        cov = (dtrue >= lo7 - 1e-12) & (dtrue <= hi7 + 1e-12)
        dl_hi = hi7[1:] - lo7[0]; dl_lo = lo7[1:] - hi7[0]
        r.update(infeasible=0, cov_simult=int(cov.all()), width=float((hi7 - lo7).mean()), width_last=float(hi7[-1] - lo7[-1]), width_first=float(hi7[0] - lo7[0]),
                 e7_neg=int((dl_hi < -L.DELTA).any()), e7_pos=int((dl_lo > L.DELTA).any()), best_endpoint_drop=float(-dl_hi.min()))
        t1 = time.time()
        j = A9.e7j_decide(ts, lo, hi, N, Jw, gi["natoms"])
        r.update(j_infeasible=j["j_infeasible"], j_drop=j["drop"], j_rise=j["rise"], j_gate_flag=j["gate_flag"], j_nbad=j["nbad"], j_amb=j["n_amb"],
                 j_minU_drop=j["minU_drop"], j_maxL_drop=j["maxL_drop"], t_joint=time.time() - t1)
        out.append(r)
    VAR.update(name="B0", ext=False, cons=())
    return out


CELLS = None


def worker(task):
    cid, lo, hi = task; c = CELLS[cid]; rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows


def main():
    global CELLS
    t0 = time.time(); sp = [s for s in RI.core_specs() if s["kind"] in ("control", "crossover")]
    if "--only2d" in args: sp = [s for s in sp if s["fam"] == "aniso2_w16"]
    if "--cells" in args:
        want = {int(x) for x in args[args.index("--cells") + 1].split(",")}; sp = [s for s in sp if s["id"] in want]
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (4 if SMOKE else 200)
    rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    CELLS = {s["id"]: s for s in sp}
    print(f"prep {time.time()-t0:.0f}s {len(sp)} cells", flush=True)
    sfx = "_smoke" if SMOKE else ""
    chunk = 2 if SMOKE else 5
    tasks = sorted([(s["id"], lo, min(lo + chunk, nd)) for s in sp for lo in range(0, nd, chunk)], key=lambda t: (CELLS[t[0]]["kind"] == "control", -CELLS[t[0]]["N"]))
    rows = []
    out = os.path.join(R.RES, f"rt_monolp{sfx}{args[args.index('--tag') + 1] if '--tag' in args else ''}.csv")
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, r) in enumerate(pool.imap_unordered(worker, tasks)):
            rows += r
            if k % 25 == 0:
                print(f"{k+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True); pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"done {len(rows)} rows {time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    main()
