#!/usr/bin/env python3
"""P22 Amendments 10b / 10c / 11: resumable orchestrator (no estimator, bar, cell or seed is defined here; see PLAN.md and IMPLEMENTATION_NOTES_A10_A11.md).
Usage: run_amend1011.py STEP [--smoke] [--ncpu K]
Steps (each idempotent and resumable: one pickle per task in results/items_<name>/, a manifest results/manifest_<name>.json at the end):
  verify   : estimator files byte-identical to 1ba9ca6 (hard-coded sha256) -> else exit 2
  repro    : re-run Amendment 9b (b32) on 3 cells x 10 datasets and compare every row with results/raw_amend9_b32.csv
  b10      : 10b certified E7/E7j on the pre-registered selection (run_p22j.run_one)
  exq10    : 10b exact-q versions of the 7 checked cells
  win10c   : 10c fixed windows from exact curves -> results/amend10c_windows.json (sha256 logged) BEFORE any data
  c10c     : 10c fixed-window budget comparison (+ original-window rows) + E4 bound diagnosis
  k11exact : Amendment 11 exact table (no simulation)
  k11      : Amendment 11 generalisation holdout (run ONCE)
  all      : verify, repro, then (b10, exq10, win10c, c10c, k11exact, k11) in that order.
--smoke uses tiny subsets and writes to *_smoke names."""
import os, sys, time, json, hashlib, pickle, traceback
STEP = sys.argv[1] if len(sys.argv) > 1 else "all"
ARGV = list(sys.argv)
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
SMOKE = "--smoke" in ARGV
NCPU = int(ARGV[ARGV.index("--ncpu") + 1]) if "--ncpu" in ARGV else 10
assert NCPU <= 10, "brief: 10 worker processes"
import numpy as np, pandas as pd
import multiprocessing as mp
sys.argv = [ARGV[0], "verify"]          # neutral argv for the module-level parsing of the legacy runners
import p22lib as L, p22amend as A, p22amend6 as A6, p22amend8 as A8, p22honest as H
import run_p22 as R, run_p22e as RE, run_p22g as RG
import run_p22i as RI, run_p22j as RJ, run_p22k as RK
sys.argv = ARGV
RES = R.RES
SFX = "_smoke" if SMOKE else ""
FROZEN = {"p22lib.py": "a4517f3e056246266a6eaae2db0063b0d5fb1e2fb903946b40eb4248b3d43b91",
          "p22amend.py": "f496087a9dad23ee92ca76f28724cc829f50b50c2011c36db86ac94d9b6806d8",
          "p22amend6.py": "1d56c21774ff4402078942ccc88234bb4e9ed0d5a363ac14b0778b48b00bad18",
          "p22amend8.py": "6b0964bffc9e69c8ec73d1573e2798378486f3ca638d8e5ad5ceec4f0dc8ebf3",
          "p22honest.py": "871ed8b3e2a5049f23fa8c84c1a8f50e8a5556363f1c3e3f35096925ff192216",
          "run_p22.py": "4a1937dff73cc1cb6f52bb67d5f5e7ec8f0b0bcd9087053f05feba375bc9848d"}


def log(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def verify():
    here = os.path.dirname(os.path.abspath(__file__)); bad = []
    for f, h in FROZEN.items():
        got = sha(os.path.join(here, f))
        if got != h: bad.append((f, h, got))
    if bad:
        log("MISMATCH frozen estimator files:", bad); sys.exit(2)
    log("frozen estimator files byte-identical to 1ba9ca6:", {f: h[:12] for f, h in FROZEN.items()})


# ----------------------------------------------------------------------------------------- resumable pool
_FN = None


def _wrap(task):
    key, t = task
    try:
        return key, _FN(t), None
    except Exception:
        return key, None, traceback.format_exc()[-600:]


def pool_run(name, tasks, fn, ncpu=None):
    """tasks: list of (key, payload). fn(payload) -> list of dict rows. Items are pickled to results/items_<name>/<key>.pkl (atomic), skipped when present.
    Returns (rows, n_tasks, n_done_before, errors)."""
    global _FN
    ncpu = ncpu or NCPU
    d = os.path.join(RES, f"items_{name}"); os.makedirs(d, exist_ok=True)
    todo = [t for t in tasks if not os.path.exists(os.path.join(d, t[0] + ".pkl"))]
    done_before = len(tasks) - len(todo); errors = []
    log(f"[{name}] {len(tasks)} tasks, {done_before} already done, {len(todo)} to run, ncpu {ncpu}")
    _FN = fn; t0 = time.time()
    if todo:
        with mp.get_context("fork").Pool(ncpu) as pool:
            for k, (key, rows, err) in enumerate(pool.imap_unordered(_wrap, todo, chunksize=1)):
                if err is not None:
                    errors.append((key, err)); log(f"[{name}] ERROR {key}: {err}")
                    continue
                tmp = os.path.join(d, key + ".tmp"); pickle.dump(rows, open(tmp, "wb")); os.replace(tmp, os.path.join(d, key + ".pkl"))
                if k % 10 == 0: log(f"[{name}] {k+1}/{len(todo)} new tasks, elapsed {time.time()-t0:.0f}s")
    rows = []
    for key, _ in tasks:
        p = os.path.join(d, key + ".pkl")
        if os.path.exists(p): rows += pickle.load(open(p, "rb"))
    return rows, len(tasks), done_before, errors


def finish(name, rows, ntasks, done_before, errors, extra=None):
    d = os.path.join(RES, f"items_{name}")
    csv = os.path.join(RES, f"raw_{name}.csv"); pd.DataFrame(rows).to_csv(csv, index=False)
    items = sorted(f for f in os.listdir(d) if f.endswith(".pkl"))
    isha = {f: sha(os.path.join(d, f)) for f in items}
    man = dict(name=name, smoke=SMOKE, tasks_expected=ntasks, items_present=len(items), tasks_resumed_from_disk=done_before, errors=[k for k, _ in errors],
               complete=(len(items) == ntasks and not errors), n_rows=len(rows), csv=os.path.basename(csv), csv_sha256=sha(csv),
               items_sha256_of_sha256s=hashlib.sha256("".join(f"{f} {isha[f]}\n" for f in items).encode()).hexdigest(),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"), frozen_files=FROZEN)
    if extra: man.update(extra)
    json.dump(man, open(os.path.join(RES, f"manifest_{name}.json"), "w"), indent=1, sort_keys=True)
    open(os.path.join(RES, f"items_{name}.sha256"), "w").write("".join(f"{isha[f]}  items_{name}/{f}\n" for f in items))
    log(f"[{name}] manifest: complete={man['complete']} items {len(items)}/{ntasks} rows {len(rows)} csv sha256 {man['csv_sha256']}")
    return man


def chunk_tasks(cells, nd, chunk=5):
    ts = []
    for s in cells:
        for lo in range(0, nd, chunk): ts.append((s["id"], lo, min(lo + chunk, nd)))
    return ts


# ----------------------------------------------------------------------------------------- repro (Amendment 9b)
def repro():
    RI.GROUP = "b32"; RI.SEED_OFF = 2000
    sp = RI.select()
    want = [s for s in sp if (s["kind"] == "control" and s["fam"] == "ring" and s["L"] == 128) or (s["fam"] == "aniso3_w8" and s["L"] == 48) or (s["fam"] == "aniso2_w16" and s["L"] == 128)]
    rr = np.random.default_rng(2200)
    for s in want: R.prep(s, rr)
    RI.CELLS = {s["id"]: s for s in want}
    nd = 10; rows = []
    for s in want:
        for ds in range(nd): rows += RI.run_dataset(s, ds)
    new = pd.DataFrame(rows)
    old = pd.read_csv(os.path.join(RES, "raw_amend9_b32.csv"), low_memory=False)
    old = old[old.cell.isin([s["id"] for s in want]) & (old.ds < nd)]
    # RI.run_dataset records ds as the data index (0..nd-1); the 2000 offset enters only the seed. old csv has the same convention.
    num = [c for c in ("abstain", "J", "hw", "cov_simult", "cov_point", "decl_neg", "decl_pos", "dmin", "dmax", "M") if c in new.columns and c in old.columns]
    m = new.merge(old, on=["cell", "ds", "est"], suffixes=("_n", "_o"))
    nmis = 0; worst = 0.0; ndisc = 0; DISC = ("abstain", "J", "cov_simult", "cov_point", "decl_neg", "decl_pos", "M")
    for c in num:
        a = m[c + "_n"].astype(float).values; b = m[c + "_o"].astype(float).values
        both = ~(np.isnan(a) & np.isnan(b)); dif = np.where(both, np.abs(np.nan_to_num(a, nan=1e9) - np.nan_to_num(b, nan=1e9)), 0)
        bad = dif > 1e-9; nmis += int(bad.sum())
        if c in DISC and c != "M": ndisc += int(bad.sum())
        if c == "M": ndisc += int((dif > 1e-6).sum())
        worst = max(worst, float(dif.max()) if len(dif) else 0.0)
        if bad.any(): log(f"  mismatch column {c}: {int(bad.sum())} rows, by est: {m.loc[bad, 'est'].value_counts().to_dict()}, max {float(dif.max()):.3e}")
    log(f"repro: {len(new)} new rows, {len(old)} old rows, {len(m)} matched, columns {num}, mismatching cells {nmis}, max abs diff {worst:.3e}")
    ok = len(m) == len(new) == len(old) and ndisc == 0          # discrete outcomes (abstain, J, coverage, declarations, M) must be identical; continuous solver outputs may differ at ties across platforms (logged)
    ex = m[m.est == "E1q"][["cell", "ds", "cov_simult_n", "cov_simult_o", "hw_n", "hw_o"]].head(3).to_string()
    log("repro sample E1q rows:\n" + ex)
    p = os.path.join(RES, f"repro_amend9_b32{SFX}.json")
    json.dump(dict(ok=bool(ok), n_new=len(new), n_old=len(old), n_matched=len(m), mismatching_cells=nmis, discrete_mismatches=ndisc, max_abs_diff=worst, cells=[s["id"] for s in want], nd=nd), open(p, "w"), indent=1)
    log("REPRO", "PASS" if ok else "FAIL (max diff above; inspect before launch)")
    if not ok: sys.exit(3)


# ----------------------------------------------------------------------------------------- 10b
def b10():
    selp = os.path.join(RES, "amend10b_selection.json"); sel = json.load(open(selp)); log("10b selection sha256", sha(selp))
    tasks = {}
    for g in ("sens", "withheld", "stable300", "minL3d"):
        for cid, ds in sel[g]: tasks.setdefault((cid, ds), []).append(g)
    tl = [(cid, ds, "+".join(g)) for (cid, ds), g in sorted(tasks.items())]
    if SMOKE:
        for cid in sorted({t[0] for t in tl}): RJ.get_cell(cid)
        tl = sorted(tl, key=lambda t: RJ.get_cell(t[0])["N"])[:2]
    else:
        for cid in sorted({t[0] for t in tl}): RJ.get_cell(cid)           # prep in the parent; forked workers inherit
    tl.sort(key=lambda t: -(RJ.get_cell(t[0])["N"]))
    pt = [(f"{c}_{d}", t) for t in tl for c, d in [(t[0], t[1])]]
    rows, n, db, er = pool_run("amend10b" + SFX, pt, lambda t: [RJ.run_one(t)])
    finish("amend10b" + SFX, rows, n, db, er, dict(selection_sha256=sha(selp)))


def exq10():
    rows_, n, db, er = pool_run("amend10b_exactq" + SFX, [(f"{f}_{l}", s) for s in (RJ.EXQ[:1] if SMOKE else RJ.EXQ) for f, l, _ in [s]], lambda t: [RJ.exq_one(t)], ncpu=min(NCPU, 7))
    finish("amend10b_exactq" + SFX, [{k: v for k, v in r.items()} for r in rows_], n, db, er)


# ----------------------------------------------------------------------------------------- 10c
def cells_10c():
    """the realised crossover cells of Amendment 9b at wpn 8/32/128 (ids IDBASE + index, as run_p22i), prepped."""
    out = {}
    rr = np.random.default_rng(2200)
    for w, g in ((8, "b8"), (32, "b32"), (128, "b128")):
        RI.GROUP = g
        sp = [s for s in RI.select() if s["kind"] == "crossover"]
        for s in sp: R.prep(s, rr)
        out[w] = sp
    return out


def win10c():
    cc = cells_10c(); wins = {}
    for w in (8, 32, 128):
        for s in cc[w]:
            key = f"{s['fam']}_{s['L']}"
            e = wins.setdefault(key, dict(fam=s["fam"], L=int(s["L"]), N=int(s["N"]), ids={}, realised=bool(s["real"]), exact_J_wpn32=None))
            e["ids"][str(w)] = int(s["id"])
            if w == 32: e["exact_J_wpn32"] = int(s["exact_J"]); e["exact_window_end"] = int(s["exact_window_end"]); e["max_drop"] = float(s["max_drop"]); e["realised"] = bool(s["real"])
    p = os.path.join(RES, "amend10c_windows.json")
    if SMOKE: p = os.path.join(RES, "amend10c_windows_smoke.json")
    if os.path.exists(p) and not SMOKE:
        log("10c windows already written, sha256", sha(p)); return
    json.dump(dict(rule="J fixed per cell = exact window (Nq >= thr from the exact curve) = wpn-32 exact_J, computed before any simulation", cells=wins), open(p, "w"), indent=1, sort_keys=True)
    log("10c windows written", p, "sha256", sha(p), {k: (v["exact_J_wpn32"], v["realised"]) for k, v in wins.items()})


_G3 = {}


def g3cm(c):
    k = (c["fam"], c["L"])
    if k not in _G3: _G3[k] = H.exact_g3_cummax(c["graph"], c["ts_w"])
    return _G3[k]


def _diag(rows, c):
    cm = g3cm(c)
    for r in rows:
        if r.get("est") in ("E4c", "E4r") and r.get("abstain", 1) == 0:
            g = float(cm[int(r["J"]) - 1]); r["g3_exact"] = g; r["m_viol"] = int(r["M"] < g)
    return rows


def _run_arms(U, m, ts, N, thr, truth, ss, tags):
    ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    out = []
    try:
        out.append(A.est_E1b([U], m, ts, N, rngs[6], truth, 2, 1.0, 200, thr, "E1q"))
    except Exception as ex:
        out.append(dict(est="E1q", abstain=1, why="error:" + repr(ex)[:80]))
    try:
        out += A6.est_all([U], m, ts, N, rngs[2], rngs[3], rngs[4], truth, thr, tags, do_e6=False)
    except Exception as ex:
        out += [dict(est=t, abstain=1, why="error:" + repr(ex)[:80]) for t in tags]
    return out


C10 = None


def run_10c(task):
    cid, lo, hi = task; c = C10[cid]; rows = []
    wl_A, wl_A6 = A.window_len, A6.window_len
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    tags = ("E4c", "E4r")
    for ds in range(lo, hi):
        ss = np.random.SeedSequence([22, c["id"], 3000 + ds])
        rngd = np.random.default_rng(ss.spawn(9)[1])
        Cc, m = c["graph"].sample_blocks(rngd, wpn * N, ts, L.NBLOCK); U = Cc / m
        arms = [("orig", None)] + ([("fixed", c["J_fixed"])] if c.get("J_fixed") else [])
        for mode, J in arms:
            if J is not None:
                A.window_len = A6.window_len = (lambda cnt, p, N_, thr_, cmin=L.CMIN, J=J: J)
            try:
                out = _run_arms(U, m, ts, N, thr, truth, ss, tags)
            finally:
                A.window_len, A6.window_len = wl_A, wl_A6
            _diag(out, c)
            nat = None
            for r in out:
                r.update(cell=c["id"], ds=ds, wpn=wpn, mode=mode, J_fixed=J if J is not None else -1)
                for k in ("e6lo", "e6hi", "e4lo", "e4hi"): r.pop(k, None)
            rows += out
    return rows


def c10c():
    global C10
    wp = os.path.join(RES, "amend10c_windows_smoke.json" if SMOKE else "amend10c_windows.json"); assert os.path.exists(wp), "run win10c first"
    log("10c windows.json sha256 (logged before data generation):", sha(wp))
    wins = json.load(open(wp))["cells"]
    cc = cells_10c(); C10 = {}; cells = []
    for w in (8, 32, 128):
        for s in cc[w]:
            e = wins[f"{s['fam']}_{s['L']}"]
            assert e["ids"][str(w)] == s["id"] and e["exact_J_wpn32"] is not None, "MISMATCH windows vs cells"
            if not e["realised"]: continue                     # non-realised cells are not part of the 'realised cells' set; listed in windows.json
            if SMOKE and not (s["fam"] == "aniso2_w16" and s["L"] == 128 and w in (8, 32)): continue
            s["J_fixed"] = int(e["exact_J_wpn32"]); C10[s["id"]] = s; cells.append(s)
    # separation-scan cells at ratios 2.5 (b=20, id 302) and 3.0 (b=24, id 303): wpn 32, original windows, E4 diagnosis only (no fixed window specified)
    rr = np.random.default_rng(2200)
    import run_p22f as RF
    for s in RF.sep_specs():
        if s["id"] in (302, 303) and not SMOKE:
            s["wpn"] = 32; R.prep(s, rr); s["J_fixed"] = None; C10[s["id"]] = s; cells.append(s)
    nd = 2 if SMOKE else 200
    pd.DataFrame([{k: s.get(k) for k in ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "J_fixed", "max_drop", "real", "thr")} for s in cells]).to_csv(os.path.join(RES, f"cells_meta_amend10c{SFX}.csv"), index=False)
    for s in cells: g3cm(s)                                   # exact g3 sup, computed in the parent, before data
    cost = lambda t: -(C10[t[0]]["N"] * C10[t[0]]["wpn"] * len(C10[t[0]]["ts_w"]))
    tl = sorted(chunk_tasks(cells, nd, 2 if SMOKE else 5), key=cost)
    rows, n, db, er = pool_run("amend10c" + SFX, [(f"{a}_{b}_{c_}", (a, b, c_)) for a, b, c_ in tl], run_10c)
    finish("amend10c" + SFX, rows, n, db, er, dict(windows_sha256=sha(wp)))


# ----------------------------------------------------------------------------------------- 11
def k11exact():
    p = os.path.join(R.RES, "amend11_cells_exact.json")
    if SMOKE: p_s = p
    if os.path.exists(p):
        log("amend11 exact table already present, sha256", sha(p)); return
    rows = RK.exact_table(RK.base_specs())
    json.dump(rows, open(p, "w"), indent=1, sort_keys=True)
    log("amend11 exact table written (BEFORE any simulation), sha256", sha(p))
    for r in rows: log("  ", r)


def k11():
    pj = os.path.join(RES, "amend11_cells_exact.json"); assert os.path.exists(pj), "run k11exact first"
    log("amend11 cells_exact sha256 at run start:", sha(pj))
    sp = RK.specs(); rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    ex = {r["cidx"]: r for r in json.load(open(pj))}
    for s in sp:
        e = ex[s["cidx"]]
        assert bool(s["real"]) == e["realised"] and s["exact_window_end"] == e["exact_window_end"], f"MISMATCH exact table cell {s['cidx']}"
    if SMOKE: sp = [s for s in sp if s["wpn"] == 8 and s["fam"] in ("ring", "torus2") and s["L"] in (384, 40)]
    RK.CELLS = {s["id"]: s for s in sp}
    nd = 2 if SMOKE else 200
    keys = ("id", "cidx", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "thr", "real")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(RES, f"cells_meta_amend11{SFX}.csv"), index=False)
    cost = lambda t: -(RK.CELLS[t[0]]["N"] * RK.CELLS[t[0]]["wpn"] * len(RK.CELLS[t[0]]["ts_w"]))
    tl = sorted(chunk_tasks(sp, nd, 2 if SMOKE else 5), key=cost)

    def fn(t):
        cid, lo, hi = t; c = RK.CELLS[cid]; rows = []
        for ds in range(lo, hi): rows += RK.run_dataset(c, ds)
        return rows
    rows, n, db, er = pool_run("amend11" + SFX, [(f"{a}_{b}_{c_}", (a, b, c_)) for a, b, c_ in tl], fn)
    finish("amend11" + SFX, rows, n, db, er, dict(cells_exact_sha256=sha(pj)))


STEPS = dict(verify=verify, repro=repro, b10=b10, exq10=exq10, win10c=win10c, c10c=c10c, k11exact=k11exact, k11=k11)

if __name__ == "__main__":
    verify()
    if STEP == "all":
        repro()
        for s in ("b10", "exq10", "win10c", "c10c", "k11exact", "k11"):
            log("=== STEP", s); STEPS[s]()
        log("ALL DONE")
    else:
        STEPS[STEP]()
