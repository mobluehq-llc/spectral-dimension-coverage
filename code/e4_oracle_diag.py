#!/usr/bin/env python3
"""EXPLORATORY E4 oracle-arm diagnosis (spec: E4_ORACLE_DIAGNOSIS_SPEC.md). Usage:
  e4_oracle_diag.py replay|fresh --cells 0,1,2,3 --n N [--ncpu K] --out DIR
Re-executes the frozen E4o branch of p22honest.est_E4 (same code calls) keeping the intermediate quantities.
Writes DIR/e4diag_<arm>_cell<id>.csv (one row per dataset) and DIR/e4diag_<arm>_perj_cell<id>.csv (one row per dataset x j)."""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p22lib as L, p22honest as H
import run_p22 as R

ARM = sys.argv[1]
args = sys.argv[2:]
CELLS = {}


def dataset(c, ds):
    sid = c["id"] if ARM == "replay" else 13000 + c["id"]
    ss = np.random.SeedSequence([22, sid, ds])
    rngs = [np.random.default_rng(s) for s in ss.spawn(6)]
    ts, N = c["ts_w"], c["N"]
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts, L.NBLOCK)
    Ublk = [C / m]
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = L.grand_mean(UA); WA = sum(u.shape[0] for u in UA) * m; cntA = pA * WA
    J = L.window_len(cntA, pA, N, c["thr"])
    base = dict(cell=c["id"], ds=ds, arm=ARM)
    if J < 6:
        return dict(base, abstain=1, J=J), []
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J])
    M = max(float(c["Mo_cum"][J - 1]), H.FLOORM)
    WB = sum(u.shape[0] for u in UB) * m; floor = 0.5 / WB
    yB = np.log(np.maximum(L.grand_mean(UB), floor)[:J])
    Ps = np.maximum(L.boot_means(UB, 200, rngs[2])[:, :J], floor); Ys = np.log(Ps)
    S = np.cov(Ys, rowvar=False)
    g = np.log(c["q_exact"](ts[:J])); beta_t = -c["d_exact"](ts[:J]) / 2.0
    Wm = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J)
    for j in range(J):
        r = H.honest_weights(x, x[j], S, M)
        if r is None: return dict(base, abstain=1, J=J, why="solver"), []
        Wm[j], sd[j], Bj[j] = r
    sd = np.maximum(sd, 1e-12)
    bh = Wm @ yB; Bs = Ys @ Wm.T
    crit = float(np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95))
    hwj = crit * sd + Bj
    err = bh - beta_t
    b = Wm @ g - beta_t                       # exact smoothing bias
    nz = Wm @ (yB - g)                        # stochastic part
    assert np.abs(err - (b + nz)).max() < 1e-9, "identity err = b + n failed"
    miss = bool((np.abs(err) > hwj).any())
    bound_viol = bool((np.abs(b) > Bj + 1e-12).any())
    tmax = float((np.abs(nz) / sd).max())
    row = dict(base, abstain=0, J=J, M=M, crit=crit, hw=float(2 * hwj.mean()), cov_simult=int(not miss), miss=int(miss), bound_viol=int(bound_viol),
               max_b_over_B=float((np.abs(b) / np.maximum(Bj, 1e-300)).max()), tmax_noise=tmax, noise_ok=int(tmax <= crit),
               miss_kind=("bound" if (miss and bound_viol) else ("calibration" if miss else "")),
               miss_also_noise_excess=int(miss and bound_viol and tmax > crit))
    pj = [dict(base, J=J, j=j, b=float(b[j]), B=float(Bj[j]), n=float(nz[j]), sd=float(sd[j]), err=float(err[j]), hwj=float(hwj[j])) for j in range(J)]
    return row, pj


def worker(task):
    cid, lo, hi = task
    c = CELLS[cid]; rows = []; pjs = []
    for ds in range(lo, hi):
        r, p = dataset(c, ds); rows.append(r); pjs += p
    return cid, rows, pjs


def main():
    global CELLS
    nd = int(args[args.index("--n") + 1]); out = args[args.index("--out") + 1]; os.makedirs(out, exist_ok=True)
    want = [int(v) for v in args[args.index("--cells") + 1].split(",")]
    ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 4
    assert ncpu <= 12
    rr = np.random.default_rng(2200)
    for s in R.cell_specs():
        if s["id"] in want:
            R.prep(s, rr); s["Mo_cum"] = H.exact_g3_cummax(s["graph"], s["ts_w"]); CELLS[s["id"]] = s
    for cid in want:
        c = CELLS[cid]; assert c["fam"] == "ring" and c["kind"] == "control", c["fam"]
    t0 = time.time()
    tasks = [(cid, lo, min(lo + 10, nd)) for cid in want for lo in range(0, nd, 10)]
    acc = {cid: ([], []) for cid in want}
    with mp.get_context("fork").Pool(ncpu) as pool:
        for k, (cid, rows, pjs) in enumerate(pool.imap_unordered(worker, tasks)):
            acc[cid][0].extend(rows); acc[cid][1].extend(pjs)
            if k % 20 == 0: print(f"{k+1}/{len(tasks)} {time.time()-t0:.0f}s", flush=True)
    for cid in want:
        pd.DataFrame(sorted(acc[cid][0], key=lambda r: r["ds"])).to_csv(os.path.join(out, f"e4diag_{ARM}_cell{cid}.csv"), index=False)
        pd.DataFrame(acc[cid][1]).to_csv(os.path.join(out, f"e4diag_{ARM}_perj_cell{cid}.csv"), index=False)
    print("done", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
