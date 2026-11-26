#!/usr/bin/env python3
"""Gemini round 3 (after review): rigorous Clopper-Pearson endpoints for the full-chain Arb DROP recheck.
rt_cpbox.cp_box computes the CP endpoints in double precision (scipy beta.ppf). Here each stored endpoint is widened outward
(lo -> lo*(1-eps), hi -> hi*(1+eps), rounded outward, clamped to [0,1]) and the widening is PROVEN sufficient in Arb ball
arithmetic: with X ~ Bin(W,p), the CP lower endpoint p_L (k>0) satisfies P(X>=k | p_L) = a/2 and the upper endpoint p_U (k<W)
satisfies P(X<=k | p_U) = a/2, a = 0.05/K exactly, both tails monotone in p. So lo_w <= p_L if P(X>=k | lo_w) < a/2, and
hi_w >= p_U if P(X<=k | hi_w) < a/2. Each tail is bounded above rigorously: M terms summed exactly in Arb balls (log-gamma
binomial coefficient) plus a geometric bound on the remainder (the term ratio is monotone in j). eps starts at 2^-30 and is
multiplied by 2^6 until every endpoint verifies. The widened (true-CP-containing) box is then fed to the unchanged
rt_certify_arbfull.certify_dataset2 (P9). Frozen files are imported, never edited.
Usage: rt_cpbox_arball_rig.py [--ncpu K] [--limit N]   Output: results/bd_arball_rig.jsonl (resumable)."""
import os, sys, json, time
_argv = sys.argv; sys.argv = [_argv[0], "e"]
import numpy as np
import multiprocessing as mp
import rt_cpbox as B
import rt_certify_arbfull as AF
from flint import arb, ctx
sys.argv = _argv
args = sys.argv[1:]
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 8
LIMIT = int(args[args.index("--limit") + 1]) if "--limit" in args else None
RES = B.RES; CELLS = {}
ctx.prec = 128


def _term(W, j, p, lp, l1p):
    # C(W,j) p^j (1-p)^(W-j) as an Arb ball
    return (arb(W + 1).lgamma() - arb(j + 1).lgamma() - arb(W - j + 1).lgamma() + j * lp + (W - j) * l1p).exp()


def upper_tail_ub(W, k, p, M):
    """rigorous upper bound on P(X>=k), X~Bin(W,p), p an exact double (as arb)."""
    p = arb(p); lp = p.log(); l1p = (1 - p).log(); r = p / (1 - p)
    t = _term(W, k, p, lp, l1p); s = arb(0); j = k
    while j <= W and j < k + M:
        s += t; t = t * (W - j) / (j + 1) * r; j += 1
    if j > W: return s
    rho = arb(W - j) / (j + 1) * r          # ratio t_{i+1}/t_i for i >= j is <= rho
    if not (rho < 1): return arb(1)
    return s + t / (1 - rho)


def lower_tail_ub(W, k, p, M):
    """rigorous upper bound on P(X<=k)."""
    p = arb(p); lp = p.log(); l1p = (1 - p).log(); r = (1 - p) / p
    t = _term(W, k, p, lp, l1p); s = arb(0); j = k
    while j >= 0 and j > k - M:
        s += t; t = t * j / (W - j + 1) * r; j -= 1
    if j < 0: return s
    rho = arb(j) / (W - j) * r              # ratio t_{i-1}/t_i for i <= j is <= rho
    if not (rho < 1): return arb(1)
    return s + t / (1 - rho)


def rig_box(U, m, lo, hi, alpha=B.ALPHA):
    K = U.shape[1]; Wt = int(U.shape[0] * m)
    k = np.rint(U.mean(0) * Wt).astype(np.int64)
    half = arb(1) / 40 / K if abs(alpha - 0.05) < 1e-15 else None
    assert half is not None, alpha
    eps = 2.0 ** -30
    while eps < 1e-2:
        lo_w = np.where(lo > 0, np.nextafter(lo * (1 - eps), 0.0), 0.0); lo_w = np.maximum(lo_w, 0.0)
        hi_w = np.minimum(np.nextafter(hi * (1 + eps), 2.0), 1.0)
        ok = True
        for i in range(K):
            ki = int(k[i])
            if ki > 0:
                assert lo_w[i] <= lo[i]
                if not any(upper_tail_ub(Wt, ki, float(lo_w[i]), M) < half for M in (64, 4096)): ok = False; break
            if ki < Wt and hi_w[i] < 1.0:
                assert hi_w[i] >= hi[i]
                if not any(lower_tail_ub(Wt, ki, float(hi_w[i]), M) < half for M in (64, 4096)): ok = False; break
        if ok: return lo_w, hi_w, eps
        eps *= 64
    raise RuntimeError("CP endpoints not verified")


def work(t):
    c = CELLS[t["cell"]]; t0 = time.time()
    try:
        rngs, U, m, Jw = B.gen(c, t["wpn"], t["ds"]); ts, N = c["ts_w"], c["N"]
        qh, lo, hi = B.cp_box(U, m)
        lo_w, hi_w, eps = rig_box(U, m, lo, hi)
        r = AF.certify_dataset2(ts, lo_w, hi_w, N, Jw); r["cp_eps"] = eps
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
                tasks.append(dict(key="%s|%s|%s|%s" % (src, r["cell"], r["wpn"], r["ds"]), src=src, cell=r["cell"], wpn=r["wpn"], ds=r["ds"]))
    tasks.sort(key=lambda t: t["key"])
    out = os.path.join(RES, "bd_arball_rig.jsonl"); done = set()
    if os.path.exists(out):
        for ln in open(out):
            try: done.add(json.loads(ln)["key"])
            except Exception: pass
    todo = [t for t in tasks if t["key"] not in done]
    if LIMIT: todo = todo[:LIMIT]
    print("arball_rig: %d DROP tasks, %d to run" % (len(tasks), len(todo)), flush=True)
    with mp.get_context("fork").Pool(NCPU) as pool, open(out, "a") as f:
        for k, r in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            f.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); f.flush()
            print(k + 1, len(todo), r["key"], r.get("verdict"), r.get("n_arbfull"), r.get("ncells"), r.get("cp_eps"), round(r["wall"], 1), flush=True)


if __name__ == "__main__":
    main()
