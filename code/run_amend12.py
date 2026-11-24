#!/usr/bin/env python3
"""P22 Amendment 12 runner (Gaspar-Kukal estimator). Usage: run_amend12.py GROUP [--smoke] [--ncpu K] [--ndata N]   GROUP in b32 b8 b128 (walkers per node 32 / 8 / 128)
Cells = exactly the Amendment 9b cell set (14 controls + realised crossover cells; ids 9000/9100/9200 + index, run_p22i.select()).
Seeds SeedSequence([22, 12000 + id, ds]) (fresh). Arms: GKW GKL GKWm GKWb GK2W GK2L GK2Wb (p22amend12) + references E0g, E0 (p22lib, frozen) on the SAME walker data.
Resumable: one pickle per task chunk in results/items_amend12_<GROUP>/, atomic; manifest results/manifest_amend12_<GROUP>.json with counts and sha256; raw csv results/raw_amend12_<GROUP>.csv.
At most 6 worker processes (brief)."""
import os, sys, time, json, hashlib, pickle, traceback
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
ARGV = list(sys.argv); GROUP = ARGV[1]
SMOKE = "--smoke" in ARGV
NCPU = int(ARGV[ARGV.index("--ncpu") + 1]) if "--ncpu" in ARGV else 6
assert NCPU <= 6
import numpy as np, pandas as pd
import multiprocessing as mp
sys.argv = [ARGV[0], "b32"]
import p22lib as L, p22amend12 as G
import run_p22 as R, run_p22i as RI
sys.argv = ARGV
RES = R.RES
SFX = "_smoke" if SMOKE else ""
IDB = {"b32": 9000, "b8": 9100, "b128": 9200}
FROZEN = {"p22lib.py": "a4517f3e056246266a6eaae2db0063b0d5fb1e2fb903946b40eb4248b3d43b91",
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
    log("frozen p22lib.py / run_p22.py byte-identical to the Amendment 11 run")


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, 12000 + c["id"], ds]); ch = ss.spawn(9); rngs = [np.random.default_rng(s) for s in ch]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    C, m = c["graph"].sample_blocks(rngs[1], wpn * N, ts, L.NBLOCK); U = C / m
    nom = G.NOMINAL.get(c["fam"]) if c["kind"] == "control" else None
    out = G.est_gk([U], m, ts, N, rngs[5], truth, 200, thr, nom)
    out += L.est_E0([U], m, ts, N, rngs[6], truth, 200, thr)
    for r in out: r.update(cell=c["id"], ds=ds, wpn=wpn)
    return out


def worker(task):
    key, (cid, lo, hi) = task
    try:
        rows = []
        for ds in range(lo, hi): rows += run_dataset(CELLS[cid], ds)
        return key, rows, None
    except Exception:
        return key, None, traceback.format_exc()[-800:]


def main():
    global CELLS
    verify()
    RI.GROUP = GROUP; sp = RI.select()
    if SMOKE: sp = [s for s in sp if (s["fam"] == "ring" and s["L"] in (128, 256)) or (s["fam"] == "aniso2_w16" and s["L"] == 64) or (s["fam"] == "torus2" and s["L"] == 24)]
    rr = np.random.default_rng(2200)
    for s in sp: R.prep(s, rr)
    CELLS = {s["id"]: s for s in sp}
    nd = int(ARGV[ARGV.index("--ndata") + 1]) if "--ndata" in ARGV else (6 if SMOKE else 200)
    name = f"amend12_{GROUP}{SFX}"; d = os.path.join(RES, f"items_{name}"); os.makedirs(d, exist_ok=True)
    keys = ("id", "fam", "L", "N", "wpn", "kind", "exact_window_end", "exact_J", "max_drop", "max_rise", "d_start", "d_end", "thr", "real")
    pd.DataFrame([{k: s.get(k) for k in keys} for s in sp]).to_csv(os.path.join(RES, f"cells_meta_{name}.csv"), index=False)
    chunk = 3 if SMOKE else 10
    tasks = sorted([(f"{s['id']}_{lo}_{min(lo+chunk, nd)}", (s["id"], lo, min(lo + chunk, nd))) for s in sp for lo in range(0, nd, chunk)],
                   key=lambda t: -CELLS[t[1][0]]["N"] * CELLS[t[1][0]]["wpn"])
    todo = [t for t in tasks if not os.path.exists(os.path.join(d, t[0] + ".pkl"))]
    log(f"[{name}] {len(sp)} cells, {len(tasks)} tasks, {len(tasks)-len(todo)} already done, ncpu {NCPU}")
    errors = []; t0 = time.time()
    if todo:
        with mp.get_context("fork").Pool(NCPU) as pool:
            for k, (key, rows, err) in enumerate(pool.imap_unordered(worker, todo, chunksize=1)):
                if err: errors.append(key); log("ERROR", key, err); continue
                tmp = os.path.join(d, key + ".tmp"); pickle.dump(rows, open(tmp, "wb")); os.replace(tmp, os.path.join(d, key + ".pkl"))
                if k % 10 == 0: log(f"[{name}] {k+1}/{len(todo)} new, {time.time()-t0:.0f}s")
    rows = []
    for key, _ in tasks:
        p = os.path.join(d, key + ".pkl")
        if os.path.exists(p): rows += pickle.load(open(p, "rb"))
    csv = os.path.join(RES, f"raw_{name}.csv"); pd.DataFrame(rows).to_csv(csv, index=False)
    items = sorted(f for f in os.listdir(d) if f.endswith(".pkl"))
    man = dict(name=name, tasks_expected=len(tasks), items_present=len(items), errors=errors, complete=(len(items) == len(tasks) and not errors), n_rows=len(rows),
               csv_sha256=sha(csv), items_sha256=hashlib.sha256("".join(f"{f} {sha(os.path.join(d, f))}\n" for f in items).encode()).hexdigest(),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"), frozen=FROZEN, estimator_file_sha256=sha(os.path.join(os.path.dirname(os.path.abspath(__file__)), "p22amend12.py")))
    json.dump(man, open(os.path.join(RES, f"manifest_{name}.json"), "w"), indent=1, sort_keys=True)
    log(f"[{name}] manifest complete={man['complete']} rows {len(rows)} csv {man['csv_sha256'][:16]}")


if __name__ == "__main__":
    main()
