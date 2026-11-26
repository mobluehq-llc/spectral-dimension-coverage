#!/usr/bin/env python3
"""BD (PLAN.md): exact Clopper-Pearson box, near-null controls, power map. Frozen p22 files are imported, never edited.
Usage: rt_cpbox.py MODE [--ncpu K] [--ndata N]    MODE in harness | b-freeze | a | b | c | d | e | f
Run from a directory holding the p22_spectral_dimension_flow/*.py files and results/ (deep ~/hp/BD).
Rows: one JSON line per dataset in results/bd_<mode>.jsonl (resumable: done keys are skipped)."""
import os, sys, time, json, hashlib
_argv = sys.argv; sys.argv = [_argv[0], "a"]
import numpy as np
import multiprocessing as mp
from scipy.stats import beta
import p22lib as L, p22amend8 as A8, rt_nodrop as RT
import run_p22 as R, run_p22i as RI
sys.argv = _argv
args = sys.argv[1:]
MODE = args[0]
DELTA = L.DELTA
ALPHA = 0.05
RES = R.RES
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 10
NDATA = int(args[args.index("--ndata") + 1]) if "--ndata" in args else 100
A_CELLS = [("aniso2_w16", 768), ("aniso3_w8", 160)]
A_CTRL = [("ring", 1024), ("torus2", 96), ("torus3", 24)]
A_WPN = (8, 32, 128)
C_WPN = (8, 32, 128, 512)
TARGET_DROPS = (0.15, 0.20, 0.23, 0.25, 0.27, 0.30, 0.35)
B_WPN = (32, 128)
HARNESS_N = 2000
CELLS = {}


# ---------------------------------------------------------------- the box
def cp_box(U, m, alpha=ALPHA):
    """Two-sided Clopper-Pearson at level alpha/K at each of the K grid times of the cell (K = full grid). Returns qhat, lo, hi."""
    K = U.shape[1]; Wt = U.shape[0] * m
    k = np.rint(U.mean(0) * Wt); a = alpha / K
    lo = np.where(k > 0, beta.ppf(a / 2, np.maximum(k, 1), Wt - k + 1), 0.0)
    hi = np.where(k < Wt, beta.ppf(1 - a / 2, k + 1, np.maximum(Wt - k, 1)), 1.0)
    return k / Wt, lo, hi


# ---------------------------------------------------------------- data (same seeds as paper A where stored)
def gen(c, wpn, ds):
    if wpn == 32 and c.get("stored", True): ss = np.random.SeedSequence([22, c["id"], ds])
    else: ss = np.random.SeedSequence([22, c["id"], 9000 + wpn, ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(9)]
    ts, N, thr = c["ts_w"], c["N"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    p = U.mean(0); Wt = U.shape[0] * m
    Jw = L.window_len(p * Wt, p, N, thr)
    return rngs, U, m, Jw


def secant(ts, lo, hi, Jw, dtrue):
    """convex-secant procedure on a box: log q is convex, so d_j >= 2 t_j (log lo_j - log hi_p)/(t_p - t_j) (p > j) and d_j <= 2 t_j (log hi_p - log lo_j)/(t_j - t_p) (p < j)."""
    K = len(ts); t = ts.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        llo = np.log(lo); lhi = np.log(hi)
    dl = np.zeros(Jw); du = np.full(Jw, np.inf)
    for j in range(Jw):
        for p in (j + 1, j + 2):
            if p < K and np.isfinite(llo[j]): dl[j] = max(dl[j], 2 * t[j] * (llo[j] - lhi[p]) / (t[p] - t[j]))
        for p in (j - 1, j - 2):
            if p >= 0 and np.isfinite(llo[j]): du[j] = min(du[j], 2 * t[j] * (lhi[p] - llo[j]) / (t[j] - t[p]))
    neg = int(((du[1:] - dl[0]) < -DELTA).any()) if Jw > 1 else 0
    cov = int(((dtrue[:Jw] >= dl - 1e-12) & (dtrue[:Jw] <= du + 1e-12)).all())
    return dict(s_neg=neg, s_cov=cov)


def key_of(t): return "%s|%s|%s|%s|%s" % (t["mode"], t["cid"], t["wpn"], t["ds"], t.get("box", "cp"))


def work(t):
    c = CELLS[t["cid"]]; wpn = t["wpn"]; ds = t["ds"]; t0 = time.time()
    base = dict(key=key_of(t), mode=t["mode"], cell=c["id"], fam=c["fam"], L=c["L"], N=c["N"], wpn=wpn, ds=ds, true_drop=float(c.get("max_drop", np.nan)))
    try:
        if t["mode"] == "harness":
            return harness_one(c, ds, base)
        if t["mode"] == "exact":
            ts, N = c["ts_w"], c["N"]; Jw = int(c["exact_J"]); q = c["graph"].exact(ts)[0]; eps = t["eps"]
            r = dict(base, eps=eps, Jw=Jw); r.update(RT.global_nodrop(ts, q * (1 - eps), q * (1 + eps), N, Jw)); r["wall"] = time.time() - t0; return r
        rngs, U, m, Jw = gen(c, wpn, ds)
        ts, N = c["ts_w"], c["N"]
        r = dict(base, Jw=Jw, K=len(ts), box=t.get("box", "cp"))
        if Jw < 2: r["g_outcome"] = "NOWINDOW"; r["wall"] = time.time() - t0; return r
        qn, lon, hin, Ps, binfo = A8.make_box(U, m, rngs[7])
        qh, locp, hicp = cp_box(U, m)
        qex = c["graph"].exact(ts)[0]
        r["box_cov"] = int(((locp <= qex) & (qex <= hicp)).all())
        r["hw_ratio"] = float(np.mean((hicp - locp)[:Jw] / np.maximum((hin - lon)[:Jw], 1e-300)))
        if t["mode"] == "f":
            CUR_Q3[0] = c["graph"].exact(np.array([1.0, 2.0, 3.0]))[0]
            r.update(RT.global_nodrop(ts, locp, hicp, N, Jw))
        elif t["mode"] == "a":
            r.update(RT.global_nodrop(ts, locp, hicp, N, Jw))
            t1 = time.time(); jr, gi = RT.joint(ts, locp, hicp, N, Jw); r.update(jr); r["j_seconds"] = time.time() - t1
            r.update(secant(ts, locp, hicp, Jw, c["d_exact"](ts)))
        elif t["mode"] in ("b", "c"):
            lo, hi = (locp, hicp) if t.get("box", "cp") == "cp" else (lon, hin)
            r.update(RT.global_nodrop(ts, lo, hi, N, Jw))
        r["wall"] = time.time() - t0
        return r
    except Exception as e:                       # never lose the pool; the error row is counted in the report
        return dict(base, g_outcome="ERROR", err=repr(e)[:300], wall=time.time() - t0)


def harness_one(c, ds, base):
    ss = np.random.SeedSequence([22, 9999, ds]); rng = np.random.default_rng(ss)
    ts, N = c["ts_w"], c["N"]
    C, m = c["graph"].sample_blocks(rng, 32 * N, ts, L.NBLOCK); U = C / m
    qh, lo, hi = cp_box(U, m); qex = c["graph"].exact(ts)[0]
    return dict(base, K=len(ts), cover_all=int(((lo <= qex) & (qex <= hi)).all()), n_miss=int(((lo > qex) | (qex > hi)).sum()))


# ---------------------------------------------------------------- cell construction
def torus(dims): return L.Graph("torus", 0.5, dims=tuple(dims))


def prod_spec(Lv, w, cid):
    return dict(fam="prod", L=Lv, w=w, build=lambda: torus((Lv, w)), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=0.5, note="", id=cid, wpn=32, stored=False)


def core_cells(pairs):
    sp_all = [s for s in RI.core_specs() if s["kind"] in ("control", "crossover")]
    ids = {s["id"] for s in sp_all if (s["fam"], s["L"]) in pairs}
    return RT.specs(ids)


def power_pairs():
    sp_all = [s for s in RI.core_specs() if s["kind"] in ("control", "crossover")]
    out = []
    for fam in ("aniso2_w16", "aniso3_w8"):
        Ls = sorted({s["L"] for s in sp_all if s["fam"] == fam}, reverse=True)[:4]
        out += [(fam, l) for l in Ls]
    return out


def freeze_b():
    rr = np.random.default_rng(2200); cand = []
    for Lv in range(24, 321):
        for w in range(2, 25):
            if Lv * w > 3000: continue
            c = R.prep(prod_spec(Lv, w, 8000), rr)
            if c.get("exact_J", 0) >= 2 and np.isfinite(c["max_drop"]):
                cand.append(dict(L=Lv, w=w, N=int(c["N"]), exact_J=int(c["exact_J"]), max_drop=float(c["max_drop"]), d_start=float(c["d_start"]), d_end=float(c["d_end"])))
    chosen = []; used = set()
    for tgt in TARGET_DROPS:
        best = sorted([x for x in cand if (x["L"], x["w"]) not in used], key=lambda x: (round(abs(x["max_drop"] - tgt), 3), x["N"]))[0]
        used.add((best["L"], best["w"])); chosen.append(dict(best, target=tgt, id=8100 + len(chosen)))
    json.dump(dict(targets=TARGET_DROPS, graphs=chosen, n_candidates=len(cand)), open(os.path.join(RES, "bd_b_graphs.json"), "w"), indent=1)
    for g in chosen: print(g, flush=True)
    print("n candidates", len(cand))


def load_b():
    g = json.load(open(os.path.join(RES, "bd_b_graphs.json")))["graphs"]; rr = np.random.default_rng(2200); out = []
    for x in g:
        c = R.prep(prod_spec(x["L"], x["w"], x["id"]), rr)
        assert abs(c["max_drop"] - x["max_drop"]) < 1e-9, "frozen graph list does not reproduce"
        c["target"] = x["target"]; out.append(c)
    return out


# ---------------------------------------------------------------- driver
def run_tasks(tasks, name):
    out = os.path.join(RES, "bd_%s.jsonl" % name)
    done = set()
    if os.path.exists(out):
        for ln in open(out):
            try: done.add(json.loads(ln)["key"])
            except Exception: pass
    todo = [t for t in tasks if key_of(t) not in done]
    cost = lambda t: -CELLS[t["cid"]]["N"] * t["wpn"] * len(CELLS[t["cid"]]["ts_w"])
    todo.sort(key=cost)
    print("%s: %d tasks, %d done, %d to run, ncpu %d" % (name, len(tasks), len(done), len(todo), NCPU), flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(NCPU) as pool, open(out, "a") as f:
        for k, r in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); f.flush()
            if k % 25 == 0: print("%d/%d %.0fs last %s %s" % (k + 1, len(todo), time.time() - t0, r.get("key"), r.get("g_outcome")), flush=True)
    print("done %s %.0fs" % (name, time.time() - t0), flush=True)


def gate():
    p = os.path.join(RES, "bd_harness_summary.json")
    if not os.path.exists(p) or not json.load(open(p)).get("pass"):
        raise SystemExit("HARNESS GATE (bar 1) not passed: %s missing or pass != true; refusing to run mode %s" % (p, MODE))


def main():
    global CELLS
    t0 = time.time(); os.makedirs(RES, exist_ok=True)
    if MODE == "b-freeze": freeze_b(); return
    if MODE == "harness":
        sp = core_cells([("ring", 1024)]); CELLS = {s["id"]: s for s in sp}; c = sp[0]
        run_tasks([dict(mode="harness", cid=c["id"], wpn=32, ds=i) for i in range(HARNESS_N)], "harness")
        rows = [json.loads(l) for l in open(os.path.join(RES, "bd_harness.jsonl"))]
        n = len(rows); k = sum(r["cover_all"] for r in rows); lowcp = float(beta.ppf(0.05, k, n - k + 1)) if k > 0 else 0.0
        s = dict(n=n, covered=k, rate=k / n, one_sided_cp_lower=lowcp, K=rows[0]["K"], pass_=bool(lowcp >= 0.94))
        s["pass"] = s.pop("pass_"); json.dump(s, open(os.path.join(RES, "bd_harness_summary.json"), "w"), indent=1); print(s); return
    gate()
    if MODE in ("a", "f"):
        sp = core_cells(A_CELLS + A_CTRL); CELLS = {s["id"]: s for s in sp}
        tasks = []
        for s in sp:
            wl = A_WPN if (s["fam"], s["L"]) in A_CELLS else (32,)
            if MODE == "f": wl = (32,)
            for w in wl:
                for ds in range(NDATA): tasks.append(dict(mode=MODE, cid=s["id"], wpn=w, ds=ds))
        if MODE == "f":
            tasks = [t for t in tasks if (CELLS[t["cid"]]["fam"], CELLS[t["cid"]]["L"]) in A_CELLS]
            orig = A8.lp_setup
            def lp_setup_exact3(ts, lo, hi, s_atoms, N, orig=orig):
                S = orig(ts, lo, hi, s_atoms, N)
                return add_low_order(S, s_atoms, N)
            A8.lp_setup = lp_setup_exact3
        run_tasks(tasks, MODE); return
    if MODE == "b":
        bc = load_b(); CELLS = {s["id"]: s for s in bc}
        tasks = [dict(mode="b", cid=s["id"], wpn=w, ds=ds, box=bx) for s in bc for w in B_WPN for ds in range(NDATA) for bx in ("cp", "nom")]
        run_tasks(tasks, "b"); return
    if MODE in ("c", "d"):
        sp = core_cells(power_pairs()); CELLS = {s["id"]: s for s in sp}
        print("POWER-MAP CELLS", [(s["fam"], s["L"], s["id"], s["N"], round(float(s["max_drop"]), 4)) for s in sp], flush=True)
        a_ids = {s["id"] for s in sp if (s["fam"], s["L"]) in A_CELLS}
        if MODE == "c":
            tasks = [dict(mode="c", cid=s["id"], wpn=w, ds=ds) for s in sp for w in C_WPN for ds in range(NDATA) if not (s["id"] in a_ids and w in A_WPN)]
        else:
            tasks = [dict(mode="exact", cid=s["id"], wpn=0, ds=-1, eps=e) for s in sp for e in (1e-4, 1e-3, 1e-2)]
        run_tasks(tasks, MODE); return
    if MODE == "e":
        import rt_certify as RC
        sp = core_cells(A_CELLS); CELLS = {s["id"]: s for s in sp}
        rows = [json.loads(l) for l in open(os.path.join(RES, "bd_a.jsonl"))]
        drops = sorted([r for r in rows if r.get("g_outcome") == "DROP"], key=lambda r: (r["cell"], r["wpn"], r["ds"]))
        rng = np.random.default_rng(2211); mpset = set(rng.choice(len(drops), min(10, len(drops)), replace=False).tolist()) if drops else set()
        tasks = [dict(mode="e", cid=r["cell"], wpn=r["wpn"], ds=r["ds"], mp=int(i in mpset), skip3d=0) for i, r in enumerate(drops)]
        run_tasks_cert(tasks, RC); return
    raise SystemExit("bad mode")


def add_low_order(S, s_atoms, N):
    """arm (f): add q(1), q(2), q(3) as (near-)equalities (exact values from degrees, edges, triangles) to the LP."""
    q3 = CUR_Q3[0]
    if q3 is None: return S
    tt = np.array([1.0, 2.0, 3.0])
    A = np.exp(-np.outer(tt, np.where(np.isinf(s_atoms), 0.0, s_atoms))); A[:, np.isinf(s_atoms)] = 0.0
    As3 = A / q3[:, None]
    lo, hi = q3 * (1 - 1e-9), q3 * (1 + 1e-9); n = S["n"]
    rows = np.vstack([np.hstack([-As3, (lo / q3)[:, None]]), np.hstack([As3, -(hi / q3)[:, None]])])
    S["Aub"] = np.vstack([S["Aub"], rows]); S["bub"] = np.concatenate([S["bub"], np.zeros(rows.shape[0])])
    return S
CUR_Q3 = [None]


def cert_work(t):
    import rt_certify as RC
    c = CELLS[t["cid"]]; rngs, U, m, Jw = gen(c, t["wpn"], t["ds"]); ts, N = c["ts_w"], c["N"]
    qh, lo, hi = cp_box(U, m); t0 = time.time()
    try:
        r = RC.certify_dataset(ts, lo, hi, N, Jw, verify_mp=bool(t["mp"]))
    except Exception as e:
        r = dict(verdict="ERROR", err=repr(e)[:300])
    return dict(key=key_of(t), mode="e", cell=t["cid"], wpn=t["wpn"], ds=t["ds"], Jw=Jw, wall=time.time() - t0, **r)


def run_tasks_cert(tasks, RC):
    out = os.path.join(RES, "bd_e.jsonl"); done = set()
    if os.path.exists(out):
        for ln in open(out):
            try: done.add(json.loads(ln)["key"])
            except Exception: pass
    todo = [t for t in tasks if key_of(t) not in done]
    print("e: %d tasks, %d to run" % (len(tasks), len(todo)), flush=True)
    with mp.get_context("fork").Pool(NCPU) as pool, open(out, "a") as f:
        for k, r in enumerate(pool.imap_unordered(cert_work, todo, chunksize=1)):
            f.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); f.flush()
            print(k + 1, len(todo), r.get("key"), r.get("verdict"), flush=True)


if __name__ == "__main__":
    main()
