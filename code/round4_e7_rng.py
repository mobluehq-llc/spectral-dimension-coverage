#!/usr/bin/env python3
"""ROUND 4 EXPLORATORY check 1 (spec: ROUND4_CHECKS_SPEC.md section 1): E7 bootstrap RNG determinism and the c_box drift. usage: round4_e7_rng.py {single|scan} [--out DIR]
single: dataset 0 of the REPRO_CHECK 9b cell. scan: float-vs-integer count>=30 disagreement over 9b b32 cells (sampling only)."""
import os, sys, json, re
ARGV = list(sys.argv); MODE = ARGV[1]
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = ARGV[ARGV.index("--out") + 1] if "--out" in ARGV else os.path.join(HERE, "round4_checks"); os.makedirs(OUT, exist_ok=True)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(v, "1")
sys.argv = [ARGV[0], "verify"]
import numpy as np, pandas as pd
import run_amend1011 as M
sys.argv = ARGV
R, RI, L = M.R, M.RI, M.L
import p22amend8 as A8

def setup():
    RI.GROUP = "b32"; RI.SEED_OFF = 2000
    sp = RI.select(); s = [s for s in sp if s["kind"] == "control" and s["fam"] == "ring" and s["L"] == 128][0]
    R.prep(s, np.random.default_rng(2200)); RI.CELLS = {s["id"]: s}; return s

def draw(s, ds):
    ss = np.random.SeedSequence([22, s["id"], RI.SEED_OFF + ds]); ch = ss.spawn(9)
    rngs = [np.random.default_rng(x) for x in ch]
    ts, N, wpn = s["ts_w"], s["N"], s["wpn"]
    C, m = s["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK)
    return rngs, C, m, ts, N

def logbox(U, m, rng, label):
    """replicates make_box's bootstrap with the same calls, returns the logged quantities (uses make_box itself for c)."""
    st = rng.bit_generator.state["state"]["state"]
    probe = np.random.Generator(np.random.PCG64()); probe.bit_generator.state = rng.bit_generator.state
    n = U.shape[0]; mult = probe.multinomial(n, np.full(n, 1.0 / n), size=5)       # first five draws only (same stream position as make_box's first five rows)
    qhat, lo, hi, Ps, info = A8.make_box(U, m, rng)
    big = qhat * (U.shape[0] * m) >= L.CMIN
    se = np.maximum(Ps[:, big].std(0, ddof=1), 1e-300); T = (np.abs(Ps[:, big] - qhat[big]) / se).max(1)
    return dict(label=label, state_before=str(st), first5_mult_head=[m_[:12].tolist() for m_ in mult], first5_mult_sha=[__import__("hashlib").sha256(m_.tobytes()).hexdigest()[:12] for m_ in mult],
                first5_T=T[:5].tolist(), c_box=info["c_box"], n_boot=info["n_boot"], n_cp=info["n_cp"], T_sha=__import__("hashlib").sha256(T.tobytes()).hexdigest()[:16], hw_head=(hi - lo)[:6].tolist())

def single():
    s = setup(); res = {"cell": s["id"], "numpy": np.__version__, "CMIN": L.CMIN}
    rngs, C, m, ts, N = draw(s, 0); U = C / m
    # (b) two full run_dataset in one process
    rows_a = RI.run_dataset(s, 0); rows_b = RI.run_dataset(s, 0)
    e7 = lambda rows: {k: r.get(k) for r in rows if r["est"] == "E7" for k in ("c_box", "n_boot", "n_cp", "hw", "width_mean", "a1p_maxgap", "gate1")}
    res["run_dataset_twice_same_process_E7"] = [e7(rows_a), e7(rows_b)]; res["twice_identical"] = e7(rows_a) == e7(rows_b)
    # (c) isolated make_box with fresh rngs[7] (no other estimator first), +/- global np.random seeding
    outs = {}
    for lab, seed in (("isolated", None), ("global_seed0", 0), ("global_seed1", 1)):
        rngs, C, m, ts, N = draw(s, 0)
        if seed is not None: np.random.seed(seed)
        outs[lab] = logbox(C / m, m, rngs[7], lab)
    # E7 as inside the full sequence: run the other estimators first (their rngs are separate children) then E7 with rngs[7]
    rngs, C, m, ts, N = draw(s, 0); U = C / m
    truth = dict(q_exact=s["q_exact"], d_exact=s["d_exact"])
    _ = L.est_E0([U], m, ts, N, rngs[5], truth, 200, s["thr"]); _ = M.A.est_E1b([U], m, ts, N, rngs[6], truth, 2, 1.0, 200, s["thr"], "E1q")
    outs["after_other_estimators"] = logbox(U, m, rngs[7], "after_other_estimators")
    res["make_box_logs"] = outs
    keys = ["first5_mult_sha", "first5_T", "c_box", "n_boot", "T_sha"]
    res["all_make_box_runs_identical"] = all(all(outs[k][q] == outs["isolated"][q] for q in keys) for k in outs)
    # (d) near-tie
    rngs, C, m, ts, N = draw(s, 0); U = C / m; nb = U.shape[0]; Wt = nb * m
    qhat = U.mean(0); cntf = qhat * Wt; cnti = C.sum(0)
    tie = [dict(k=int(k), t=float(ts[k]), cnt_float=repr(float(cntf[k])), cnt_int=float(cnti[k]), diff=float(cntf[k] - 30)) for k in range(len(ts)) if abs(cntf[k] - 30) < 1e-4 or abs(cnti[k] - 30) < 1e-9]
    res["points_near_30"] = tie
    res["n_big_float_rule"] = int((cntf >= L.CMIN).sum()); res["n_big_int_rule"] = int((cnti >= L.CMIN).sum())
    res["m"] = float(m); res["nb"] = int(nb); res["cnt_float_all"] = [repr(float(v)) for v in cntf]; res["cnt_int_all"] = [float(v) for v in cnti]
    # c_box under both rules, same Ps draws (rngs[7] fresh)
    rngs, C, m, ts, N = draw(s, 0); U = C / m
    Ps = L.boot_means([U], A8.BOOT, rngs[7])
    def cbox(big):
        se = np.maximum(Ps[:, big].std(0, ddof=1), 1e-300); T = (np.abs(Ps[:, big] - qhat[big]) / se).max(1); return float(np.quantile(T, 1 - A8.A_BOX / 4)), int(big.sum())
    res["c_box_float_rule"] = cbox(cntf >= L.CMIN); res["c_box_int_rule"] = cbox(cnti >= L.CMIN)
    # also: different summation orders of the same mean
    alt = {}
    for lab, f in (("np.mean_axis0", lambda: U.mean(0)), ("math.fsum", lambda: np.array([__import__("math").fsum(U[:, k]) / nb for k in range(U.shape[1])])), ("sum_C_over_m_over_nb", lambda: C.sum(0) / m / nb), ("float32", lambda: U.astype(np.float32).mean(0).astype(np.float64))):
        q = f(); alt[lab] = int(((q * Wt) >= L.CMIN).sum())
    res["n_big_by_mean_method"] = alt
    json.dump(res, open(os.path.join(OUT, "e7_rng_single.json"), "w"), indent=1, default=str); print(json.dumps({k: v for k, v in res.items() if not k.startswith("cnt_")}, indent=1, default=str)[:6000])

def scan():
    """sampling only: save raw C per dataset for the 9b groups under the running numpy (compare two envs afterwards with --compare)."""
    import multiprocessing as mp
    global _SP
    groups = ARGV[ARGV.index("--groups") + 1].split(",") if "--groups" in ARGV else ["b8", "b32", "b128"]
    ncpu = int(ARGV[ARGV.index("--ncpu") + 1]) if "--ncpu" in ARGV else 4; nd = int(ARGV[ARGV.index("--nd") + 1]) if "--nd" in ARGV else 200
    store = {}
    for g in groups:
        RI.GROUP = g; RI.SEED_OFF = 2000
        sp = RI.select(); rr = np.random.default_rng(2200)
        for s in sp: R.prep(s, rr)
        _SP = {(g, s["id"]): s for s in sp}
        tasks = [(g, s["id"], lo, min(lo + 20, nd)) for s in sp for lo in range(0, nd, 20)]
        with mp.get_context("fork").Pool(ncpu) as pool:
            for part in pool.imap_unordered(_scan_task, tasks):
                store.update(part)
        print(g, "cells", len(sp), "datasets so far", len(store), flush=True)
    np.savez_compressed(os.path.join(OUT, f"e7_walks_np{np.__version__}.npz"), **store)

def _scan_task(t):
    g, cid, lo, hi = t; s = _SP[(g, cid)]; out = {}
    for ds in range(lo, hi):
        rngs, C, m, ts, N = draw(s, ds); out[f"{g}_{cid}_{ds}"] = C.astype(np.int16)
    return out

def compare():
    a = np.load(ARGV[ARGV.index("--a") + 1]); b = np.load(ARGV[ARGV.index("--b") + 1]); rows = []
    for k in a.files:
        x, y = a[k], b[k]
        diff = (x != y).any(0); first = int(np.argmax(diff)) if diff.any() else -1
        g, cid, ds = k.split("_")[0], int(k.split("_")[1]), int(k.split("_")[2])
        rows.append(dict(group=g, cell=cid, ds=ds, differs=int(diff.any()), first_diff_t=first, n_diff_t=int(diff.sum()), K=x.shape[1]))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(OUT, "e7_walk_diff_by_dataset.csv"), index=False)
    for g, t in d.groupby("group"):
        print(g, "datasets", len(t), "differ", int(t.differs.sum()), "cells with any diff", int(t.groupby("cell").differs.max().sum()), "of", t.cell.nunique())
        print("  first differing time-index distribution (from end: K-1-idx):", (t[t.differs == 1].K - 1 - t[t.differs == 1].first_diff_t).value_counts().sort_index().to_dict())

if __name__ == "__main__":
    {"single": single, "scan": scan, "compare": compare}[MODE]()
