#!/usr/bin/env python3
"""P22 Amendment 13(a) runner: CONFIRMATORY pass of the independent-split bias-aware variants E4o_V2 and E4c_V2 (plus the frozen V0 arms E4o, E4c on the SAME data).
Usage: run_amend13a.py DESIGN WPN [--smoke] [--ncpu K] [--ndata N]
  DESIGN rep  = replication design: the Amendment 9b cell set (14 controls + crossover cells; ids 9000/9100/9200 + index, run_p22i.select)
  DESIGN hold = holdout design: the Amendment 11 cell set (run_p22k.specs, ids 100*wpn_index + cell_index)
  WPN 8 | 32 | 128 (walkers per node).  200 datasets per cell; fresh seeds SeedSequence([22, 15000 + id, ds]).
Resumable: one pickle per task chunk in results/items_amend13a_<design>_b<WPN>/ (atomic write), manifest results/manifest_amend13a_<design>_b<WPN>.json (counts, sha256s, frozen-file hashes),
raw csv results/raw_amend13a_<design>_b<WPN>.csv. At most 6 worker processes (brief: <= 10 cores per host, hosts are shared); run under nice."""
import os, sys, time, json, hashlib, pickle, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
ARGV = list(sys.argv); DESIGN = ARGV[1]; WPN = int(ARGV[2])
assert DESIGN in ("rep", "hold") and WPN in (8, 32, 128)
SMOKE = "--smoke" in ARGV
NCPU = int(ARGV[ARGV.index("--ncpu") + 1]) if "--ncpu" in ARGV else 6
assert NCPU <= 10
import numpy as np, pandas as pd
import multiprocessing as mp
sys.argv = [ARGV[0], "b32"]
import p22lib as L, p22honest as H, p22amend13 as A13
import run_p22 as R, run_p22i as RI, run_p22k as RK
sys.argv = ARGV
RES = R.RES
SFX = "_smoke" if SMOKE else ""
NAME = f"amend13a_{DESIGN}_b{WPN}{SFX}"
FROZEN = {"p22lib.py": "a4517f3e056246266a6eaae2db0063b0d5fb1e2fb903946b40eb4248b3d43b91",
          "p22honest.py": "871ed8b3e2a5049f23fa8c84c1a8f50e8a5556363f1c3e3f35096925ff192216",
          "run_p22.py": "4a1937dff73cc1cb6f52bb67d5f5e7ec8f0b0bcd9087053f05feba375bc9848d"}
CELLS = None


def log(*a): print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def verify():
    here = os.path.dirname(os.path.abspath(__file__)); bad = [(f, h, sha(os.path.join(here, f))) for f, h in FROZEN.items() if sha(os.path.join(here, f)) != h]
    if bad: log("MISMATCH frozen files", bad); sys.exit(2)
    log("frozen p22lib.py / p22honest.py / run_p22.py byte-identical to Amendment 11/12")


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, 15000 + c["id"], ds]); rngs = [np.random.default_rng(s) for s in ss.spawn(9)]
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    out = A13.est_E4v2([U], m, ts, N, rngs[2], rngs[3], [("S", c["d_exact"])], c["Mo_cum"], thr)
    for r in out: r.update(cell=c["id"], ds=ds, wpn=wpn, design=DESIGN)
    return out


def worker(task):
    key, (cid, lo, hi) = task
    try:
        rows = []
        for ds in range(lo, hi): rows += run_dataset(CELLS[cid], ds)
        return key, rows, None
    except Exception:
        return key, None, traceback.format_exc()[-800:]


def select():
    if DESIGN == "rep":
        RI.GROUP = f"b{WPN}"; return RI.select()
    return [s for s in RK.specs() if s["wpn"] == WPN]


def main():
    global CELLS
    verify()
    sp = select()
    if SMOKE: sp = [s for s in sp if (s["fam"] == "ring" and s["L"] in (128, 384)) or (s["fam"] in ("aniso2_w16", "aniso2_w12") and s["L"] in (128, 96))]
    rr = np.random.default_rng(2200); g3 = {}
    for s in sp:
        R.prep(s, rr); k = (s["fam"], s["L"])
        if k not in g3: g3[k] = H.exact_g3_cummax(s["graph"], s["ts_w"])
        s["Mo_cum"] = g3[k]
    CELLS = {s["id"]: s for s in sp}
    nd = int(ARGV[ARGV.index("--ndata") + 1]) if "--ndata" in ARGV else (4 if SMOKE else 200)
    d = os.path.join(RES, f"items_{NAME}"); os.makedirs(d, exist_ok=True)
    keys = ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "thr", "real")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(RES, f"cells_meta_{NAME}.csv"), index=False)
    chunk = 2 if SMOKE else 5
    tasks = sorted([(f"{s['id']}_{lo}_{min(lo+chunk, nd)}", (s["id"], lo, min(lo + chunk, nd))) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[1][0]]["N"] * CELLS[t[1][0]]["wpn"] * len(CELLS[t[1][0]]["ts_w"]))
    todo = [t for t in tasks if not os.path.exists(os.path.join(d, t[0] + ".pkl"))]
    log(f"[{NAME}] {len(sp)} cells, {len(tasks)} tasks, {len(tasks)-len(todo)} already done, ncpu {NCPU}")
    errors = []; t0 = time.time()
    if todo:
        with mp.get_context("fork").Pool(NCPU) as pool:
            for k, (key, rows, err) in enumerate(pool.imap_unordered(worker, todo, chunksize=1)):
                if err: errors.append(key); log("ERROR", key, err); continue
                tmp = os.path.join(d, key + ".tmp"); pickle.dump(rows, open(tmp, "wb")); os.replace(tmp, os.path.join(d, key + ".pkl"))
                if k % 20 == 0: log(f"[{NAME}] {k+1}/{len(todo)} new, {time.time()-t0:.0f}s")
    rows = []
    for key, _ in tasks:
        p = os.path.join(d, key + ".pkl")
        if os.path.exists(p): rows += pickle.load(open(p, "rb"))
    csv = os.path.join(RES, f"raw_{NAME}.csv"); pd.DataFrame(rows).to_csv(csv, index=False)
    items = sorted(f for f in os.listdir(d) if f.endswith(".pkl"))
    here = os.path.dirname(os.path.abspath(__file__))
    man = dict(name=NAME, tasks_expected=len(tasks), items_present=len(items), errors=errors, complete=(len(items) == len(tasks) and not errors), n_rows=len(rows),
               csv_sha256=sha(csv), items_sha256=hashlib.sha256("".join(f"{f} {sha(os.path.join(d, f))}\n" for f in items).encode()).hexdigest(),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"), frozen=FROZEN, p22amend13_sha256=sha(os.path.join(here, "p22amend13.py")), runner_sha256=sha(os.path.abspath(__file__)),
               numpy=np.__version__)
    json.dump(man, open(os.path.join(RES, f"manifest_{NAME}.json"), "w"), indent=1, sort_keys=True)
    log(f"[{NAME}] manifest complete={man['complete']} rows {len(rows)} csv {man['csv_sha256'][:16]}")


if __name__ == "__main__":
    main()
