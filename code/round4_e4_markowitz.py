#!/usr/bin/env python3
"""ROUND 4 EXPLORATORY check 2 (spec: ROUND4_CHECKS_SPEC.md section 2): E4o on ring L=512 (cell id 2, wpn 32) with V0 baseline / V1 B=1000 / V2 independent-split covariance.
usage: round4_e4_markowitz.py {replay|fresh} --n N [--ncpu K<=6] --out DIR. Re-executes the frozen E4o branch of p22honest (same calls as e4_oracle_diag.py); only S and crit differ per variant."""
import os, sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(v, "1")
import p22lib as L, p22honest as H
import run_p22 as R
ARM = sys.argv[1]; args = sys.argv[2:]; CID = 2; CELL = {}

def variant(vn, c, J, x, ts, M, yB, g, beta_t, S, Ys):
    Wm = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J)
    for j in range(J):
        r = H.honest_weights(x, x[j], S, M)
        if r is None: return None
        Wm[j], sd[j], Bj[j] = r
    sd = np.maximum(sd, 1e-12); bh = Wm @ yB; Bs = Ys @ Wm.T
    crit = float(np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95)); hwj = crit * sd + Bj
    err = bh - beta_t; b = Wm @ g - beta_t; nz = Wm @ (yB - g)
    assert np.abs(err - (b + nz)).max() < 1e-9
    miss = bool((np.abs(err) > hwj).any()); bv = bool((np.abs(b) > Bj + 1e-12).any()); tmax = float((np.abs(nz) / sd).max())
    return dict(v=vn, cov=int(not miss), bound_viol=int(bv), crit=crit, noise_ok=int(tmax <= crit), tmax=tmax, hw=float(2 * hwj.mean())), [dict(v=vn, j=j, n=float(nz[j]), sd=float(sd[j]), b=float(b[j]), B=float(Bj[j])) for j in range(J)]

def dataset(ds):
    c = CELL; sid = c["id"] if ARM == "replay" else 13000 + c["id"]
    ss = np.random.SeedSequence([22, sid, ds]); rngs = [np.random.default_rng(s) for s in ss.spawn(6)]
    ts, N = c["ts_w"], c["N"]
    C, m = c["graph"].sample_blocks(rngs[1], R.WPN * N, ts, L.NBLOCK); Ublk = [C / m]
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = L.grand_mean(UA); WA = sum(u.shape[0] for u in UA) * m; J = L.window_len(pA * WA, pA, N, c["thr"])
    base = dict(ds=ds, arm=ARM)
    if J < 6: return [dict(base, v="all", abstain=1, J=J)], []
    x = np.log(ts[:J].astype(float)); M = max(float(c["Mo_cum"][J - 1]), H.FLOORM)
    WB = sum(u.shape[0] for u in UB) * m; floor = 0.5 / WB
    yB = np.log(np.maximum(L.grand_mean(UB), floor)[:J])
    g = np.log(c["q_exact"](ts[:J])); beta_t = -c["d_exact"](ts[:J]) / 2.0
    # V0 (frozen): 200 draws of half B, rngs[2]
    Ys0 = np.log(np.maximum(L.boot_means(UB, 200, rngs[2])[:, :J], floor))
    # V1: 1000 draws half B, own stream
    r1 = np.random.default_rng(np.random.SeedSequence([22, sid, ds, 1000]))
    Ys1 = np.log(np.maximum(L.boot_means(UB, 1000, r1)[:, :J], floor))
    # V2: S from half A (independent of yB), 200 draws; crit from V0's half-B draws
    r2 = np.random.default_rng(np.random.SeedSequence([22, sid, ds, 2000]))
    floorA = 0.5 / WA; YsA = np.log(np.maximum(L.boot_means(UA, 200, r2)[:, :J], floorA))
    rows = []; pjs = []
    for vn, S, Ys in (("V0_baseline", np.cov(Ys0, rowvar=False), Ys0), ("V1_B1000", np.cov(Ys1, rowvar=False), Ys1), ("V2_split_A", np.cov(YsA, rowvar=False), Ys0)):
        o = variant(vn, c, J, x, ts, M, yB, g, beta_t, S, Ys)
        if o is None: rows.append(dict(base, v=vn, abstain=1, J=J, why="solver")); continue
        r, p = o; rows.append(dict(base, abstain=0, J=J, **r)); pjs += [dict(q, ds=ds, J=J) for q in p]
    return rows, pjs

def worker(t):
    rows = []; pjs = []
    for ds in range(*t):
        r, p = dataset(ds); rows += r; pjs += p
    return rows, pjs

def cp(k, n):
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)

def main():
    global CELL
    nd = int(args[args.index("--n") + 1]); out = args[args.index("--out") + 1]; ncpu = int(args[args.index("--ncpu") + 1]) if "--ncpu" in args else 4; assert ncpu <= 6
    os.makedirs(out, exist_ok=True); rr = np.random.default_rng(2200)
    for s in R.cell_specs():
        if s["id"] in (0, 1, 2, 3):
            R.prep(s, rr)
            if s["id"] == CID: s["Mo_cum"] = H.exact_g3_cummax(s["graph"], s["ts_w"]); CELL = s
    assert CELL["fam"] == "ring" and CELL["L"] == 512 and CELL["kind"] == "control"
    t0 = time.time(); rows = []; pjs = []
    with mp.get_context("fork").Pool(ncpu) as pool:
        for r, p in pool.imap_unordered(worker, [(lo, min(lo + 10, nd)) for lo in range(0, nd, 10)]): rows += r; pjs += p
    d = pd.DataFrame(rows); pj = pd.DataFrame(pjs); d.to_csv(os.path.join(out, f"markowitz_{ARM}_rows.csv"), index=False); pj.to_csv(os.path.join(out, f"markowitz_{ARM}_perj.csv"), index=False)
    lines = [f"EXPLORATORY E4o Markowitz check, ring L=512 (cell {CID}), arm={ARM}, n={nd}, numpy {np.__version__}, {time.time()-t0:.0f}s"]
    nmin = 30 if nd <= 200 else 100
    for vn in ("V0_baseline", "V1_B1000", "V2_split_A"):
        a = d[(d.v == vn) & (d.abstain == 0)]; ab = nd - len(a); k = int(a["cov"].sum()); n = len(a); lo, hi = cp(k, n)
        p = pj[pj.v == vn]; res = []
        for J, gj in p.groupby("J"):
            if gj.ds.nunique() < nmin: continue
            for j, h in gj.groupby("j"): res.append(h.sd.median() / h.n.std(ddof=1))
        z = (p.n / p.sd).to_numpy()
        lines.append(f"{vn}: coverage {k}/{n} = {100*k/n:.1f}% (CP {100*lo:.1f}-{100*hi:.1f}); abstained {ab}; sd ratio (median over J-strata>={nmin}, j; strata used {len(res)}) {np.median(res) if res else float('nan'):.3f}; "
                     f"pooled 1/std(n/sd) {1/z.std(ddof=1):.3f}; stochastic-part coverage {100*a.noise_ok.mean():.1f}%; mean crit {a.crit.mean():.3f}; mean hw {a.hw.mean():.4f}; datasets with |b|>B: {int(a.bound_viol.sum())}")
    if ARM == "replay":
        r3 = pd.read_csv(os.path.join(HERE, "results", "raw_amend3.csv"), low_memory=False); r3 = r3[(r3.est == "E4o") & (r3.cell == CID)]
        v0 = d[(d.v == "V0_baseline") | (d.v == "all")].merge(r3, on="ds", suffixes=("", "_o"))
        both = v0[(v0.abstain == 0) & (v0.abstain_o == 0)]
        lines.append(f"replay check vs stored raw_amend3 E4o cell {CID}: rows {len(v0)}; abstain equal {(v0.abstain == v0.abstain_o).sum()}; answered {len(both)}; J equal {(both.J == both.J_o).sum()}; cov equal {(both['cov'] == both.cov_simult).sum()}; hw within 1e-6 {(abs(both.hw - both.hw_o) < 1e-6).sum()}; stored coverage {int(r3[r3.abstain == 0].cov_simult.sum())}/{int((r3.abstain == 0).sum())}")
    open(os.path.join(out, f"markowitz_{ARM}_report.txt"), "w").write("\n".join(lines)); print("\n".join(lines))

if __name__ == "__main__": main()
