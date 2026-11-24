#!/usr/bin/env python3
"""Analysis for row29_explore.py (EXPLORATORY; no confirmatory claim). usage: row29_analyze.py OUTDIR STORED_RAW_CSV > report
Sections: (1) reproduction of stored rows ds 0..199; (2) E0g and sibling coverage on original / new / pooled replicates; (3) bias vs SE;
(4) window selection (data-driven J vs fixed J)."""
import sys, glob, pickle, numpy as np, pandas as pd
from scipy.stats import beta
OUT, STORED = sys.argv[1], sys.argv[2]
rows, ex = [], []
for f in sorted(glob.glob(f"{OUT}/t_*.pkl")):
    r = pickle.load(open(f, "rb")); rows += r[3]; ex += r[4]
R = pd.DataFrame(rows); X = pd.DataFrame(ex); S = pd.read_csv(STORED, low_memory=False); S = S[S.cell.isin([9212, 9213])]
def cp(k, n):
    k, n = int(k), int(n); return (beta.ppf(.025, k, n-k+1) if k > 0 else 0.0, beta.ppf(.975, k+1, n-k) if k < n else 1.0)
def fmt(k, n): lo, hi = cp(k, n); return f"{int(k)}/{int(n)} = {100*k/max(n,1):.1f}% (CP {100*lo:.1f}-{100*hi:.1f})"
print(f"datasets per cell: {R.groupby('cell').ds.nunique().to_dict()}  (ds 0-199 original seeds, 200-{R.ds.max()} new)\n")
print("== 1. REPRODUCTION of stored rows (ds 0-199) ==")
A = R[R.ds < 200].sort_values(["cell", "ds", "est"]).reset_index(drop=True); B = S[S.ds < 200].sort_values(["cell", "ds", "est"]).reset_index(drop=True)
K = A[["cell", "ds"]].drop_duplicates(); B = B.merge(K, on=["cell", "ds"]).sort_values(["cell", "ds", "est"]).reset_index(drop=True)   # compare only the datasets run so far
B = B[A.columns.intersection(B.columns)]; A = A[B.columns]
print(f"rows: new {len(A)} stored(matching datasets) {len(B)}; same keys: {len(A)==len(B) and (A[['cell','ds','est']].values == B[['cell','ds','est']].values).all()}")
num = [c for c in B.columns if c not in ("cell", "ds", "est") and pd.api.types.is_numeric_dtype(B[c])]
d = (A[num].astype(float).values - B[num].astype(float).values); nanm = np.isnan(A[num].astype(float).values) != np.isnan(B[num].astype(float).values)
fin = np.nan_to_num(np.abs(d)); print(f"NaN-pattern mismatches: {int(nanm.sum())}; values compared: {d.size}; exactly equal: {int((fin==0).sum())}; max abs diff: {fin.max():.3e}")
print("\n== 2. COVERAGE (cov_simult, as scored in Amendment 12; bar is >=93% pooled over controls, per-cell CP upper must be >=93%) ==")
for est in ("E0g", "E0", "GKWb", "GKWm", "GKW"):
    for cell in (9212, 9213):
        x = R[(R.cell == cell) & (R.est == est)]; ok = x[x.abstain == 0]
        parts = [("original ds0-199", ok[ok.ds < 200]), ("new ds200-799", ok[ok.ds >= 200]), ("pooled ds0-799", ok)]
        print(f"{est:5s} cell {cell}: abst {len(x)-len(ok)}/{len(x)}; " + "; ".join(f"{nm}: {fmt(p.cov_simult.sum(), len(p))}" for nm, p in parts))
print("\n== 3. BIAS vs SE (E0g-style global fit, data-driven J = 'dataJ'; one number per dataset) ==")
D = X[X.variant == "dataJ"].copy(); D["err"] = D.dhat - D.d_true; D["errpop"] = D.dhat - D.d_pop; D["z"] = D.err / D.se_boot
for cell in (9212, 9213):
    for nm, sel in (("original ds0-199", D.ds < 200), ("new ds200-799", D.ds >= 200), ("pooled", D.ds >= 0)):
        x = D[(D.cell == cell) & sel]; n = len(x)
        if n == 0: continue
        sd = x.dhat.std(ddof=1); se_m = x.se_boot.mean(); bias = x.err.mean(); bse = x.err.std(ddof=1)/np.sqrt(n); bp = x.errpop.mean()
        cov196 = (x.err.abs() <= 1.96*x.se_boot).mean(); covp = (x.errpop.abs() <= 1.96*x.se_boot).mean()
        print(f"cell {cell} {nm} n={n}: mean(dhat-d_true)={bias:+.5f} (se {bse:.5f}, {bias/bse:+.1f} se); mean(dhat-d_pop)={bp:+.5f}; empirical SD(dhat)={sd:.5f} vs mean boot SE={se_m:.5f} (ratio SE/SD={se_m/sd:.3f}); "
              f"z=err/se: mean {x.z.mean():+.3f} SD {x.z.std(ddof=1):.3f}; cover |err|<=1.96se: {100*cov196:.1f}%; vs pseudo-true d_pop: {100*covp:.1f}%")
    x = D[D.cell == cell]
    # what coverage would be with zero bias at the observed SE ratio, and with the observed bias but perfect SE
    from scipy.stats import norm
    sdz = x.err.std(ddof=1) / x.se_boot.mean(); mu = x.err.mean()/x.se_boot.mean()
    print(f"   cell {cell} pooled model check: bias/SE = {mu:+.3f}; SD(dhat)/SE = {sdz:.3f}; normal-theory coverage of 1.96 SE: with bias only {100*(norm.cdf(1.96-mu)-norm.cdf(-1.96-mu)):.1f}%, with SD ratio only {100*(2*norm.cdf(1.96/sdz)-1):.1f}%, both {100*(norm.cdf((1.96-mu)/sdz)-norm.cdf((-1.96-mu)/sdz)):.1f}%")
print("\n== 4. WINDOW SELECTION (same walker data, variants differ only in window length J) ==")
print("data-driven J distribution (counts) and exact_J (the window rule applied to the exact curve) per cell:")
for cell in (9212, 9213):
    x = D[D.cell == cell]; print(f"  cell {cell}: Jdata counts {x.Jdata.value_counts().sort_index().to_dict()}")
    for j, g in x.groupby("Jdata"):
        print(f"     J={j}: n={len(g)} mean err {g.err.mean():+.5f}, cover 1.96se {100*(g.err.abs()<=1.96*g.se_boot).mean():.1f}%, SD(dhat) {g.dhat.std(ddof=1) if len(g)>1 else float('nan'):.5f}, mean se {g.se_boot.mean():.5f}")
print("fixed-J variants on the same data (pooled ds0-799), coverage with +-1.96 boot SE vs d_exact at the fit's log-mean time:")
for cell in (9212, 9213):
    for var in ["dataJ"] + [f"J{j}" for j in (5, 6, 7, 8, 9, 10)]:
        x = X[(X.cell == cell) & (X.variant == var)].copy()
        if x.empty: continue
        x["err"] = x.dhat - x.d_true; sd = x.dhat.std(ddof=1)
        print(f"  cell {cell} {var:6s} n={len(x)} J~{x.J.mean():.2f}: bias {x.err.mean():+.5f}; SD(dhat) {sd:.5f}; mean se {x.se_boot.mean():.5f}; cover {100*(x.err.abs()<=1.96*x.se_boot).mean():.1f}% (CP {cp((x.err.abs()<=1.96*x.se_boot).sum(), len(x))[0]*100:.1f}-{cp((x.err.abs()<=1.96*x.se_boot).sum(), len(x))[1]*100:.1f}); d_true-d_pop (fit bias on the exact curve) {np.mean(x.d_true-x.d_pop):+.5f}")
