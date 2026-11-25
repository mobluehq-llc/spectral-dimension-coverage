#!/usr/bin/env python3
"""P22 Amendment 13 analysis. Usage: analyze_amend13.py a|b|c [--smoke] [--raw DIR].  Applies the bars of PLAN.md Amendment 13 exactly as written; introduces no new bar.
Reading rules (inherited from Amendments 9-11, verdict_amend10_11.py R1-R7): coverage/FPR denominators are the datasets on which the estimator answered (primary, binding);
the strict reading (abstention = miss) is printed beside; plain numeric bars are tested on the POINT estimate, the plan's own CP rules as written; CP intervals everywhere.
Outputs results/amend13{a,b,c}_*.csv and results/amend13{a,b,c}_report.txt."""
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
args = sys.argv[1:]; PART = args[0]
RES = args[args.index("--raw") + 1] if "--raw" in args else os.path.join(HERE, "results")
SFX = "_smoke" if "--smoke" in args else ""
LINES = []
def P(s=""): print(s); LINES.append(str(s))
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n):
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)
def rd(f):
    p = os.path.join(RES, f); return pd.read_csv(p, low_memory=False) if os.path.exists(p) else None
def save(name, df): df.to_csv(os.path.join(RES, name + SFX + ".csv"), index=False)
def report(name): open(os.path.join(RES, name + SFX + "_report.txt"), "w").write("\n".join(LINES) + "\n")

GRID = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768],
        "aniso2_w12": [96, 128, 192, 256, 384], "aniso3_w6": [48, 64, 96, 128, 192]}
REP_REAL = lambda fam, l: (fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48)
ARMS_A = ("E4o", "E4c", "E4o_V2", "E4c_V2")


# =============================================================================================== (a)
def load_a(design, w):
    d = rd(f"raw_amend13a_{design}_b{w}{SFX}.csv"); m = rd(f"cells_meta_amend13a_{design}_b{w}{SFX}.csv")
    if d is None or m is None: return None
    d = d.merge(m[["id", "fam", "L", "kind"] + (["real"] if "real" in m.columns else [])].rename(columns={"real": "real_meta"}), left_on="cell", right_on="id", how="inner")
    if design == "rep": d["real"] = [bool(REP_REAL(f, l) and k == "crossover") for f, l, k in zip(d.fam, d.L, d.kind)]
    else: d["real"] = d.real_meta.astype(bool) & (d.kind == "crossover")
    d["valid"] = d.abstain == 0
    for c in ("decl_neg", "decl_pos"):
        if c not in d.columns: d[c] = 0
        d[c] = d[c].fillna(0)
    d["decl_any"] = ((d.decl_neg + d.decl_pos) > 0).astype(int)
    return d


def calib(x, ncells_ctrl):
    c = x[x.kind == "control"]; ok = c[c.valid]; k, n = int(ok.cov_simult.sum()), len(ok)
    per = ok.groupby("cell").cov_simult.agg(["sum", "count"]); badc = [int(i) for i, r in per.iterrows() if cp(r["sum"], r["count"])[1] < .93]
    fk = int(ok.decl_any.sum()); pf = ok.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(i) for i, r in pf.iterrows() if cp(r["sum"], r["count"])[0] > .05]
    Cp = n > 0 and k / n >= .93 and not badc; Fp = n > 0 and fk / n <= .05 and not badf
    r = x[x.real]; rv = r[r.valid]; xk, xn = int(rv.cov_simult.sum()), len(rv)
    strict = int(ok.cov_simult.sum()) / max(len(c), 1)
    return dict(k=k, n=n, strict_cov=strict, C=bool(Cp), Cn=fmt(k, n), Cbad=badc, F=bool(Fp), Fn=fmt(fk, n), Fbad=badf, X=bool(xn > 0 and xk / xn >= .90), X93=bool(xn > 0 and xk / xn >= .93),
                Xn=fmt(xk, xn), xk=xk, xn=xn)


def minL(x, fams):
    cx = x[x.kind == "crossover"].copy(); cx["det"] = ((cx.valid) & (cx.decl_neg == 1)).astype(float)
    pw = cx.groupby(["fam", "L"]).det.agg(["sum", "count"]).reset_index(); pw["power"] = pw["sum"] / pw["count"]
    out = {}
    for fam in fams:
        g = pw[pw.fam == fam]; Pw = dict(zip(g.L, g.power)); gr = GRID[fam]; idx = len(gr)
        for i in range(len(gr)):
            if all((l in Pw) and Pw[l] >= .8 for l in gr[i:]): idx = i; break
        out[fam] = gr[idx] if idx < len(gr) else np.inf
        out[fam + "_power"] = " ".join(f"{int(l)}:{p:.2f}" for l, p in zip(g.L, g.power))
    return out


def part_a():
    P("=" * 20 + " Amendment 13(a): independent-split bias-aware variants E4o_V2, E4c_V2 (V0 = frozen E4o, E4c on the same data)")
    rows = []; cells = []; ml = []; V = {}
    for design, fams in (("rep", ["aniso2_w16", "aniso3_w8"]), ("hold", ["aniso2_w12", "aniso3_w6"])):
        for w in (8, 32, 128):
            d = load_a(design, w)
            if d is None: P(f"[{design} wpn {w}] no data"); continue
            for arm in ARMS_A:
                x = d[d.est == arm]
                if len(x) == 0: continue
                r = calib(x, None); r.update(design=design, wpn=w, arm=arm); rows.append(r)
                M = minL(x, fams); ml.append(dict(design=design, wpn=w, arm=arm, **{f: M[f] for f in fams}, **{f + "_power": M[f + "_power"] for f in fams}))
                P(f"[{design} wpn {w:3d}] {arm:7s} C {r['Cn']} -> {'PASS' if r['C'] else 'FAIL'} (cells CP-upper<93: {r['Cbad']}) | F {r['Fn']} -> {'PASS' if r['F'] else 'FAIL'} | X {r['Xn']} -> {'PASS' if r['X'] else 'FAIL'} (>=93: {r['X93']}) | strict ctrl cov {100*r['strict_cov']:.1f}% | minL " + ", ".join(f"{f} {M[f]}" for f in fams))
                for (cid, fam, L_), g in x.groupby(["cell", "fam", "L"]):
                    kind = g.kind.iloc[0]; v = g[g.valid]
                    cells.append(dict(design=design, wpn=w, arm=arm, cell=cid, fam=fam, L=L_, kind=kind, real=bool(g.real.iloc[0]), n=len(v), answered=len(v), cov=int(v.cov_simult.sum()),
                                      decl_any=int(v.decl_any.sum()), decl_neg=int((v.decl_neg == 1).sum()), hw_mean=float(v.hw.mean()) if len(v) else np.nan, M_median=float(v.M.median()) if len(v) else np.nan))
            # paired V2 vs V0 on identical datasets (controls)
            for t in ("E4o", "E4c"):
                a = d[(d.est == t) & (d.kind == "control") & d.valid][["cell", "ds", "cov_simult"]]; b = d[(d.est == t + "_V2") & (d.kind == "control") & d.valid][["cell", "ds", "cov_simult"]]
                m = a.merge(b, on=["cell", "ds"], suffixes=("0", "2"))
                P(f"   paired controls {design} wpn {w} {t}: V0 {int(m.cov_simult0.sum())}/{len(m)} V2 {int(m.cov_simult2.sum())}/{len(m)}; V2-only covers {int(((m.cov_simult2 == 1) & (m.cov_simult0 == 0)).sum())}, V0-only covers {int(((m.cov_simult0 == 1) & (m.cov_simult2 == 0)).sum())}")
    R = pd.DataFrame(rows); save("amend13a_calib", R); save("amend13a_cells", pd.DataFrame(cells)); save("amend13a_minL", pd.DataFrame(ml))
    P("")
    for arm in ("E4o_V2", "E4c_V2"):
        for design in ("rep", "hold"):
            s = R[(R.arm == arm) & (R.design == design)]
            if len(s) < 3: P(f"VERDICT {arm} {design}: incomplete ({len(s)}/3 budgets)"); continue
            Cok = bool(s.C.all()); Xok = bool(s.X.all()); Fok = bool(s[s.wpn == 32].F.all())
            P(f"VERDICT {arm} {design}: C at every wpn {'PASS' if Cok else 'FAIL'}; F at wpn 32 {'PASS' if Fok else 'FAIL'}; X>=90% at every wpn {'PASS' if Xok else 'FAIL'}  =>  {'VALIDATED' if (Cok and Xok and Fok) else 'NOT VALIDATED'}")
    # prediction check recorded in PLAN: E4o_V2 passes C on every ring control cell of the replication design
    c = pd.DataFrame(cells)
    if len(c):
        for w in (8, 32, 128):
            q = c[(c.arm == "E4o_V2") & (c.design == "rep") & (c.wpn == w) & (c.fam == "ring")]
            P(f"   E4o_V2 ring controls rep wpn {w}: " + " ".join(f"L{int(r.L)}:{int(r['cov'])}/{int(r.answered)}" for _, r in q.iterrows()))
    report("amend13a")


# =============================================================================================== (b)
def part_b():
    P("=" * 20 + " Amendment 13(b): random-graph arm (Watts-Strogatz, 2D site-percolation clusters)")
    d = rd(f"raw_amend13b{SFX}.csv"); tab = json.load(open(os.path.join(RES, f"amend13b_cells_exact{SFX}.json")))
    if d is None: P("no data"); report("amend13b"); return
    T = pd.DataFrame([{k: r[k] for k in ("id", "name", "fam", "N_mean", "N_min", "N_max", "exact_J", "max_drop", "max_rise", "max_abs_delta", "d_start", "d_end", "sd_d_between_graphs_mean_over_window", "crossover_like", "flat_like")} for r in tab])
    save("amend13b_cells", T); P(T.to_string())
    d = d.merge(T[["id", "name", "crossover_like", "flat_like"]], left_on="cell", right_on="id"); d["valid"] = (d.abstain == 0) & (d["nowin"].fillna(0) == 0)
    if "infeasible" in d.columns: d.loc[d.infeasible == 1, "valid"] = d.loc[d.infeasible == 1, "valid"]
    for c in ("decl_neg", "decl_pos"):
        d[c] = d[c].fillna(0) if c in d.columns else 0
    d["decl_any"] = ((d.decl_neg + d.decl_pos) > 0).astype(int)
    out = []; ver = []
    for tgt, sfx in (("single-realisation", ""), ("ensemble-mean", "_E")):
        col = "cov_simult" + sfx
        P(f"--- target: {tgt} (column {col})")
        for arm in ("E1q", "E4c", "E4c_V2", "E7"):
            for w in (8, 32, 128):
                x = d[(d.est == arm) & (d.wpn == w)]
                if len(x) == 0 or col not in x.columns: continue
                for grp, flag in (("crossover_like", True), ("flat_like", True)):
                    y = x[x[grp] & x.valid & x[col].notna()]; k, n = int(y[col].sum()), len(y)
                    per = y.groupby("cell")[col].agg(["sum", "count"]); badc = [int(i) for i, r in per.iterrows() if r["count"] > 0 and cp(r["sum"], r["count"])[1] < .93]
                    out.append(dict(target=tgt, arm=arm, wpn=w, group=grp, k=k, n=n, cov=k / n if n else np.nan, cp_lo=cp(k, n)[0], cp_hi=cp(k, n)[1], cells_cpupper_lt93=str(badc),
                                    strict=k / max(int(x[grp].sum()), 1)))
                    P(f"  {arm:7s} wpn {w:3d} {grp:15s} {fmt(k, n)}" + (f" cells with CP-upper<93: {badc}" if grp == "flat_like" else ""))
        P("")
    O = pd.DataFrame(out); save("amend13b_coverage", O)
    # gates
    def cov(tgt, arm, w, grp):
        r = O[(O.target == tgt) & (O.arm == arm) & (O.wpn == w) & (O.group == grp)]; return (float(r["cov"].iloc[0]), int(r.n.iloc[0]), str(r.cells_cpupper_lt93.iloc[0])) if len(r) else (np.nan, 0, "[]")
    for tgt in ("single-realisation", "ensemble-mean"):
        for arm in ("E4c", "E4c_V2", "E7"):
            ok3 = []; ok4 = []
            for w in (8, 32, 128):
                c, n, _ = cov(tgt, arm, w, "crossover_like"); ok3.append(bool(n > 0 and c >= .90))
                f, nf, bad = cov(tgt, arm, w, "flat_like"); ok4.append(bool(nf > 0 and f >= .93 and bad == "[]"))
            P(f"VERDICT B-X ({tgt}) {arm}: crossover-like coverage >= 90% at wpn 8/32/128: {ok3} -> {'PASS' if all(ok3) else 'FAIL'}; B-C flat-like >= 93% and no cell CP-upper<93: {ok4} -> {'PASS' if all(ok4) else 'FAIL'}")
    for tgt in ("single-realisation",):
        c8 = cov(tgt, "E1q", 8, "crossover_like"); c32 = cov(tgt, "E1q", 32, "crossover_like"); c128 = cov(tgt, "E1q", 128, "crossover_like")
        r8 = O[(O.target == tgt) & (O.arm == "E1q") & (O.wpn == 8) & (O.group == "crossover_like")]; r128 = O[(O.target == tgt) & (O.arm == "E1q") & (O.wpn == 128) & (O.group == "crossover_like")]
        g1 = bool(c32[1] > 0 and c32[0] <= .40); g2 = bool(c8[0] > c32[0] > c128[0] and len(r8) and len(r128) and r128.cp_hi.iloc[0] < r8.cp_lo.iloc[0])
        P(f"VERDICT B-G1 E1q crossover-like coverage at wpn 32 <= 40%: {c32[0]:.3f} -> {'PASS' if g1 else 'FAIL'}; B-G2 strictly decreasing 8>32>128 with CP-upper(128) < CP-lower(8): {c8[0]:.3f}, {c32[0]:.3f}, {c128[0]:.3f} -> {'PASS' if g2 else 'FAIL'}")
    # FPR on flat-like cells at wpn 32 for E4c / E4c_V2 (realisation-level exact excursion <= 0.25 as well)
    for arm in ("E4c", "E4c_V2"):
        x = d[(d.est == arm) & (d.wpn == 32) & d.valid & d.flat_like & (d.exc_S <= 0.25)]; k, n = int(x.decl_any.sum()), len(x)
        per = x.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(i) for i, r in per.iterrows() if cp(r["sum"], r["count"])[0] > .05]
        P(f"VERDICT B-F {arm} wpn 32: false-change rate on flat-like cells (realisation exact excursion <= 0.25) {fmt(k, n)}; cells CP-lower > 5%: {badf} -> {'PASS' if (n > 0 and k / n <= .05 and not badf) else 'FAIL'}")
    # between-graph gap diagnostic
    g = d[(d.est == "E4c") & (d.wpn == 32)].groupby("name").agg(gap_SE_mean=("gap_SE", "mean"), exc_S_mean=("exc_S", "mean"), N_r_mean=("N_r", "mean"), eig_s=("eig_s", "mean"))
    P("per-cell mean |d_S - d_E| over the realisation window / exact excursion / N / eig seconds:\n" + g.to_string())
    if "lp_seconds" in d.columns: P("E7 LP seconds: median %.2f, 95th pct %.2f" % (d[d.est == "E7"].lp_seconds.median(), d[d.est == "E7"].lp_seconds.quantile(.95)))
    report("amend13b")


# =============================================================================================== (c)
def part_c():
    exec(open(os.path.join(HERE, "a13c_scaling.py")).read(), {"__name__": "a13c", "RES": RES, "HERE": HERE, "P": P, "LINES": LINES, "report": report, "save": save, "SFX": SFX, "cp": cp, "fmt": fmt})


if __name__ == "__main__":
    {"a": part_a, "b": part_b, "c": part_c}[PART]()
