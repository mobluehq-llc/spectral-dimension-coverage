#!/usr/bin/env python3
"""Red-team R2 (RBC, CCF style) and R3 (LocQuad bias decomposition) on the Amendment 10c fixed-window design. See PREREG_REDTEAM_GEMINI.md.
Frozen files imported, never edited. LocQuad is mirrored step for step from p22amend.est_E1b (identity with the frozen function is asserted per dataset
on cov_simult, h and hw), then RBC = local-cubic slope at LocQuad's half-A bandwidth with the SAME bootstrap resamples (stream [6] regenerated twice).
Usage: rt_rbc.py [--smoke] [--ndata N] [--ncpu K] [--wpn 8,32,128] [--cells L1,L2]"""
import os, sys, time, json
_argv = sys.argv; sys.argv = [_argv[0], "verify"]
import numpy as np, pandas as pd
import multiprocessing as mp
import p22lib as L, p22amend as A
import run_p22 as R, run_amend1011 as RM
sys.argv = _argv
args = sys.argv[1:]
SMOKE = "--smoke" in args
NCPU = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 5
assert NCPU <= 5


def fmt(v): return ";".join("%.5f" % x for x in v)


def core(U, m, ts, N, rngs_pair, truth, J, thr):
    """returns dict(LQ=..., RBC=..., diag) or None (abstain)."""
    UA = [U[0::2]]; UB = [U[1::2]]
    pA = L.grand_mean(UA); WA = UA[0].shape[0] * m; cntA = pA * WA
    if J < 6: return None
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    best = None
    for h in A.H_QUAD:
        val, _ = A.poly_rows(x, wA, h, 2); hd = np.diag(val)
        if np.any(np.abs(1 - hd) < 1e-8): continue
        cv = (wA * ((yA - val @ yA) / (1 - hd)) ** 2).sum()
        if best is None or cv < best[0]: best = (cv, h)
    if best is None: return None
    h = best[1]
    _, sl2 = A.poly_rows(x, wA, h, 2); _, sl3 = A.poly_rows(x, wA, h, 3); _, sl2u = A.poly_rows(x, wA, 0.7 * h, 2)
    WB = UB[0].shape[0] * m; floor = 0.5 / WB
    pB = np.maximum(L.grand_mean(UB), floor)[:J]
    lq = np.log(truth["q_exact"](ts[:J])); dtrue = truth["d_exact"](ts[:J])
    res = {}
    for tag, sl, rng in (("LQ", sl2, rngs_pair[0]), ("RBC", sl3, rngs_pair[1]), ("LQU", sl2u, rngs_pair[2])):
        dhat = -2 * sl @ np.log(pB)
        Ps = np.maximum(L.boot_means(UB, 200, rng)[:, :J], floor)
        Ds = -2 * np.log(Ps) @ sl.T
        dpop = -2 * sl @ lq
        r = A.score(dhat, Ds, dtrue, tag, dpop)
        se = np.maximum(Ds.std(0, ddof=1), 1e-9); b = dpop - dtrue
        c = float(np.quantile((np.abs(Ds - dhat) / se).max(1), 0.95))
        # recentred band: centre dhat - (dpop - d), half-width c*se; covers d iff |dhat - dpop| <= c se (identity, computed independently here)
        rec = int((np.abs(dtrue - (dhat - b)) <= c * se).all())
        r.update(h=(0.7 * h if tag == 'LQU' else h), J=J, worst_ratio=float((np.abs(b) / se).max()), worst_t=int((np.abs(b) / se).argmax()), crit=c, rec_cov=rec,
                 b_vec=fmt(b), dev_vec=fmt(dhat - dpop), se_vec=fmt(se), mean_bias=float(b.mean()))
        res[tag] = r
    return res


def run_dataset(c, ds):
    ss = np.random.SeedSequence([22, c["id"], 3000 + ds])
    rngd = np.random.default_rng(ss.spawn(9)[1])
    ts, N, wpn, thr = c["ts_w"], c["N"], c["wpn"], c["thr"]
    truth = dict(q_exact=c["q_exact"], d_exact=c["d_exact"])
    Cc, m = c["graph"].sample_blocks(rngd, wpn * N, ts, L.NBLOCK); U = Cc / m
    ss.spawn(9)                                           # consumed by mode "orig" in run_10c
    ch = ss.spawn(9)                                      # mode "fixed"
    J = c["J_fixed"]
    wl = A.window_len
    A.window_len = (lambda cnt, p, N_, thr_, cmin=L.CMIN, J=J: J)
    try:
        frozen = A.est_E1b([U], m, ts, N, np.random.default_rng(ch[6]), truth, 2, 1.0, 200, thr, "E1q")
    finally:
        A.window_len = wl
    out = []
    base = dict(cell=c["id"], fam=c["fam"], L=c["L"], wpn=wpn, ds=ds, J_fixed=J)
    res = core(U, m, ts, N, (np.random.default_rng(ch[6]), np.random.default_rng(ch[6]), np.random.default_rng(ch[6])), truth, J, thr)
    if res is None:
        return [dict(base, est="LQ", abstain=1, frozen_abstain=int(frozen.get("abstain", 0)))]
    lq = res["LQ"]
    ident = (frozen.get("abstain", 1) == 0 and int(frozen["cov_simult"]) == int(lq["cov_simult"]) and frozen["h"] == lq["h"] and abs(frozen["hw"] - lq["hw"]) < 1e-9)
    for tag in ("LQ", "RBC", "LQU"):
        r = res[tag]; r.update(base); r.update(abstain=0, ident_frozen=int(ident), frozen_cov=frozen.get("cov_simult"), frozen_hw=frozen.get("hw"), frozen_h=frozen.get("h"))
        out.append(r)
    return out


C10 = None


def worker(task):
    cid, lo, hi = task; c = C10[cid]; rows = []
    for ds in range(lo, hi): rows += run_dataset(c, ds)
    return cid, rows


def main():
    global C10
    t0 = time.time()
    wins = json.load(open(os.path.join("results_amend10_11", "amend10c_windows.json")))["cells"]
    cc = RM.cells_10c(); C10 = {}
    wl = [int(x) for x in args[args.index("--wpn") + 1].split(",")] if "--wpn" in args else [8, 32, 128]
    Ls = {(x.split(":")[0], int(x.split(":")[1])) for x in args[args.index("--cells") + 1].split(",")} if "--cells" in args else None
    for w in wl:
        for s in cc[w]:
            e = wins[f"{s['fam']}_{s['L']}"]
            assert e["ids"][str(w)] == s["id"]
            if not e["realised"]: continue
            if Ls is not None and (s["fam"], s["L"]) not in Ls: continue
            s["J_fixed"] = int(e["exact_J_wpn32"])
            if s["J_fixed"] < 6: continue                          # structurally unanswered (J<6)
            C10[s["id"]] = s
    nd = int(args[args.index("--ndata") + 1]) if "--ndata" in args else (3 if SMOKE else 200)
    print(f"prep {time.time()-t0:.0f}s {len(C10)} cell-budgets", flush=True)
    chunk = 5
    tasks = sorted([(i, lo, min(lo + chunk, nd)) for i in C10 for lo in range(0, nd, chunk)], key=lambda t: -C10[t[0]]["N"] * C10[t[0]]["wpn"])
    rows = []; out = os.path.join(R.RES, f"rt_rbc{'_smoke' if SMOKE else ''}.csv")
    with mp.get_context("fork").Pool(NCPU) as pool:
        for k, (cid, r) in enumerate(pool.imap_unordered(worker, tasks)):
            rows += r
            if k % 25 == 0:
                print(f"{k+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True); pd.DataFrame(rows).to_csv(out + ".partial", index=False)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"done {len(rows)} rows {time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    main()
