#!/usr/bin/env python3
"""Gemini round 2 (2026-11-26): full-chain Arb ball recheck (rt_certify_arbfull.certify_dataset2, the P9 check) of EVERY DROP
outcome of the exact Clopper-Pearson box run: the 423 of Table III (bd_a.jsonl) and all DROPs of Table S6 / cpmap (bd_a for the
two A cells at 8/32/128, bd_c for the rest). Frozen files (rt_cpbox.py, rt_certify*.py) are imported, never edited.
Datasets are regenerated from the frozen seeds (rt_cpbox.gen + cp_box); reproduction is asserted by comparing d0min/d0max of the
float LP with the stored mode-e rows where they exist (the 423).
Usage (dir holding the p22 *.py, rt_cpbox.py, results/{bd_a,bd_c,bd_e}.jsonl): rt_cpbox_arball.py [--ncpu K] [--limit N]
Output: results/bd_arball.jsonl (resumable)."""
import os, sys, json, time
_argv = sys.argv; sys.argv = [_argv[0], "e"]
import numpy as np
import multiprocessing as mp
import rt_cpbox as B
import rt_certify_arbfull as AF
sys.argv = _argv
args = sys.argv[1:]
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 8
LIMIT = int(args[args.index("--limit") + 1]) if "--limit" in args else None
RES = B.RES
CELLS = {}


def work(t):
    c = CELLS[t["cell"]]; t0 = time.time()
    try:
        rngs, U, m, Jw = B.gen(c, t["wpn"], t["ds"]); ts, N = c["ts_w"], c["N"]
        qh, lo, hi = B.cp_box(U, m)
        r = AF.certify_dataset2(ts, lo, hi, N, Jw)
    except Exception as e:
        r = dict(verdict="ERROR", err=repr(e)[:300])
    r.update(key=t["key"], src=t["src"], cell=t["cell"], fam=c["fam"], L=c["L"], wpn=t["wpn"], ds=t["ds"], wall=time.time() - t0)
    return r


def main():
    global CELLS
    sp = B.core_cells(B.power_pairs()); CELLS = {s["id"]: s for s in sp}; B.CELLS.update(CELLS)
    tasks = []
    for src in ("bd_a", "bd_c"):
        for ln in open(os.path.join(RES, src + ".jsonl")):
            r = json.loads(ln)
            if r.get("g_outcome") == "DROP":
                assert r["cell"] in CELLS, r["cell"]
                tasks.append(dict(key="%s|%s|%s|%s" % (src, r["cell"], r["wpn"], r["ds"]), src=src, cell=r["cell"], wpn=r["wpn"], ds=r["ds"]))
    tasks.sort(key=lambda t: t["key"])
    out = os.path.join(RES, "bd_arball.jsonl"); done = set()
    if os.path.exists(out):
        for ln in open(out):
            try: done.add(json.loads(ln)["key"])
            except Exception: pass
    todo = [t for t in tasks if t["key"] not in done]
    if LIMIT: todo = todo[:LIMIT]
    print("arball: %d DROP tasks (%d bd_a, %d bd_c), %d to run" % (len(tasks), sum(t["src"] == "bd_a" for t in tasks), sum(t["src"] == "bd_c" for t in tasks), len(todo)), flush=True)
    with mp.get_context("fork").Pool(NCPU) as pool, open(out, "a") as f:
        for k, r in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); f.flush()
            print(k + 1, len(todo), r["key"], r.get("verdict"), r.get("n_arbfull"), r.get("ncells"), round(r["wall"], 1), flush=True)


if __name__ == "__main__":
    main()
