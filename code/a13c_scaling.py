# Amendment 13(c), EXPLORATORY (fit only, feeds no gate): two-variable scaling of the minimum detectable size against walker budget W and exact drop Delta-d.
# exec'd by analyze_amend13.py c (names P, RES, HERE, save, report, SFX provided).  Uses stored raw files (Amendment 9b, 11) and the Amendment 13(a) raw files.
import os, numpy as np, pandas as pd
GRID = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768],
        "aniso2_w12": [96, 128, 192, 256, 384], "aniso3_w6": [48, 64, 96, 128, 192]}
XS = {"aniso3_w8": 64, "aniso2_w16": 16, "aniso2_w12": 12, "aniso3_w6": 36}          # fibre cross-section: N = L * XS
DES = {"aniso3_w8": "rep", "aniso2_w16": "rep", "aniso2_w12": "hold", "aniso3_w6": "hold"}
REAL = lambda fam, l: (fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48) or fam in ("aniso2_w12", "aniso3_w6")
rng = np.random.default_rng(1313)


def rdf(name):
    p = os.path.join(RES, name); return pd.read_csv(p, low_memory=False) if os.path.exists(p) else None


def source_frames():
    """list of (source, wpn, DataFrame with est, cell, ds, L, fam, kind, det, md)"""
    out = []
    for w in (8, 32, 128):
        d = rdf(f"raw_amend9_b{w}.csv"); m = rdf(f"cells_meta_amend9_b{w}.csv")
        if d is not None and m is not None: out.append(("9b_stored", w, d.merge(m[["id", "fam", "L", "kind", "max_drop"]], left_on="cell", right_on="id")))
    d = rdf("raw_amend11.csv"); m = rdf("cells_meta_amend11.csv")
    if d is None: d = rdf(os.path.join("..", "results_amend10_11", "raw_amend11.csv")); m = rdf(os.path.join("..", "results_amend10_11", "cells_meta_amend11.csv"))
    if d is not None and m is not None:
        x = d.merge(m[["id", "fam", "L", "kind", "max_drop"]].rename(columns={"id": "cell"}), on="cell")
        for w in (8, 32, 128): out.append(("11_stored", w, x[x.wpn == w]))
    for design in ("rep", "hold"):
        for w in (8, 32, 128):
            d = rdf(f"raw_amend13a_{design}_b{w}.csv"); m = rdf(f"cells_meta_amend13a_{design}_b{w}.csv")
            if d is not None and m is not None: out.append((f"13a_{design}", w, d.merge(m[["id", "fam", "L", "kind", "max_drop"]], left_on="cell", right_on="id")))
    return out


def prep(x):
    x = x[x.kind == "crossover"].copy()
    if "nowin" not in x.columns: x["nowin"] = 0
    x["nowin"] = x.nowin.fillna(0); x["decl_neg"] = x.decl_neg.fillna(0)
    x["det"] = ((x.abstain == 0) & (x.nowin == 0) & (x.decl_neg == 1)).astype(float)
    return x


def minL_from(pw, fam):
    gr = GRID[fam]; idx = len(gr)
    for i in range(len(gr)):
        if all((l in pw) and pw[l] >= .8 for l in gr[i:]): idx = i; break
    return gr[idx] if idx < len(gr) else np.inf


_CACHE = {}


def cached():
    if not _CACHE:
        for src, w, d in source_frames():
            for arm in sorted(d.est.unique()):
                if arm not in ("E4c", "E4o", "E4c_V2", "E4o_V2", "E7"): continue
                x = prep(d[d.est == arm]); x = x[["cell", "fam", "L", "det", "max_drop"]]
                _CACHE[(src, w, arm)] = x
    return _CACHE


def points(boot=False):
    rows = []
    for (src, w, arm), x in cached().items():
        if True:
            for fam in GRID:
                f = x[x.fam == fam]
                if len(f) == 0: continue
                if boot:
                    parts = []
                    for (L_, g) in f.groupby("L"):
                        v = g.det.values; parts.append((L_, v[rng.integers(0, len(v), len(v))].mean()))
                    pw = dict(parts)
                else:
                    pw = f.groupby("L").det.mean().to_dict()
                rows.append(dict(src=src, wpn=w, arm=arm, fam=fam, Lmin=minL_from(pw, fam), dd=float(f[[REAL(fam, l) for l in f.L]].drop_duplicates("cell").max_drop.median())))
    return pd.DataFrame(rows)


def fit(df, ycol, X):
    A = np.column_stack([np.ones(len(df))] + [np.log2(df[c].values) for c in X]); y = np.log2(df[ycol].values)
    co, *_ = np.linalg.lstsq(A, y, rcond=None); res = y - A @ co
    return co, res


P("=" * 20 + " Amendment 13(c) EXPLORATORY: two-variable scaling of the minimum detectable size (fit only; no gate)")
pt = points(); pt["Nmin"] = pt.Lmin * pt.fam.map(XS); pt["design"] = pt.fam.map(DES)
save("amend13c_points", pt)
P(pt.assign(Lmin=pt.Lmin.replace(np.inf, np.nan)).pivot_table(index=["arm", "fam", "src"], columns="wpn", values="Lmin", aggfunc="first").to_string())
fin = pt[np.isfinite(pt.Lmin)].copy(); fin["WNd"] = fin.Nmin * fin.wpn * fin.dd ** (7 / 3)
BOOT = [points(boot=True) for _ in range(200)]
stepmax = max(max(np.array(g[1:]) / np.array(g[:-1])) for g in GRID.values())
res_rows = []
for arm in ("E4c", "E4c_V2", "E4o_V2", "E7"):
    f = fin[fin.arm == arm]
    if len(f) < 4: P(f"{arm}: {len(f)} finite points, no fit"); continue
    for ycol, nm in (("Lmin", "M1 log2 Lmin = a + b log2 W + c log2 dd"), ("Nmin", "M2 log2 Nmin = a + b log2 W + c log2 dd")):
        co, r = fit(f, ycol, ["wpn", "dd"])
        # bootstrap over datasets within cells (recompute the points), B = 200
        bs = []
        for pb in BOOT:
            pb = pb[(pb.arm == arm) & np.isfinite(pb.Lmin)].copy(); pb["Nmin"] = pb.Lmin * pb.fam.map(XS)
            if len(pb) >= 4:
                try: bs.append(fit(pb, ycol, ["wpn", "dd"])[0])
                except Exception: pass
        bs = np.array(bs); ci = (lambda k: (np.percentile(bs[:, k], 2.5), np.percentile(bs[:, k], 97.5)) if len(bs) > 20 else (np.nan, np.nan))
        P(f"{arm} {nm}: n={len(f)} a={co[0]:.2f} b(W)={co[1]:.2f} [{ci(1)[0]:.2f},{ci(1)[1]:.2f}] c(dd)={co[2]:.2f} [{ci(2)[0]:.2f},{ci(2)[1]:.2f}] resid sd {r.std(ddof=3) if len(f) > 3 else np.nan:.2f} log2 units (bootstrap reps {len(bs)})")
        res_rows.append(dict(arm=arm, model=ycol, n=len(f), a=co[0], b_W=co[1], c_dd=co[2], resid_sd=float(r.std(ddof=3)), b_lo=ci(1)[0], b_hi=ci(1)[1], c_lo=ci(2)[0], c_hi=ci(2)[1]))
    r = f.WNd.values; P(f"   {arm} law-implied invariant N*W*dd^(7/3): min {r.min():.0f} max {r.max():.0f} ratio {r.max()/r.min():.2f} (largest adjacent grid ratio {stepmax:.2f})")
    # out-of-sample: fit M2 on replication-design points, predict holdout-design points
    tr = f[f.design == "rep"]; te = f[f.design == "hold"]
    if len(tr) >= 3 and len(te) >= 1:
        co, _ = fit(tr, "Nmin", ["wpn", "dd"]); pred = 2 ** (co[0] + co[1] * np.log2(te.wpn.values) + co[2] * np.log2(te.dd.values)); err = np.log2(te.Nmin.values / pred)
        P(f"   {arm} out-of-sample (fit on rep design, predict hold design, M2): n={len(te)}; |log2 error| median {np.median(np.abs(err)):.2f}, max {np.abs(err).max():.2f}; within one grid step ({np.log2(stepmax):.2f}): {int((np.abs(err) <= np.log2(stepmax)).sum())}/{len(te)}")
    # censoring consistency: every 'not reached' point's predicted Lmin must exceed the largest tested L
    cen = pt[(pt.arm == arm) & ~np.isfinite(pt.Lmin)]
    if len(cen):
        co, _ = fit(f, "Lmin", ["wpn", "dd"]); ok = 0
        for _, q in cen.iterrows():
            pred = 2 ** (co[0] + co[1] * np.log2(q.wpn) + co[2] * np.log2(q.dd)); ok += int(pred > max(GRID[q.fam]))
        P(f"   {arm} censoring consistency (M1 prediction exceeds the largest tested L at 'not reached' points): {ok}/{len(cen)}")
save("amend13c_fits", pd.DataFrame(res_rows))
P("Reading rule (fixed in PLAN): fit only. A coefficient CI that contains the law-implied exponent (b_W = -1 in M2) is reported as 'not rejected', never as 'confirmed'; the data have only four families, two distinct drops and at most three budgets.")
report("amend13c")
