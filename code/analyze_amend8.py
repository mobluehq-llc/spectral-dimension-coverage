#!/usr/bin/env python3
"""P22 Amendment 8 analysis. Usage: analyze_amend8.py [--smoke]. Reads results/raw_amend8_{core,w8,w128}[_smoke].csv; writes results/amend8_*.csv, prints verdicts."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
SFX = "_smoke" if "--smoke" in sys.argv else ""
G = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768]}
LIM = {"aniso3_w8": 160, "aniso2_w16": 512}          # M bar: min L must be <= this
REAL = lambda fam, l: (fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48)
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n): lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)

def load(group):
    f = os.path.join(RES, f"raw_amend8_{group}{SFX}.csv")
    if not os.path.exists(f): return None, None
    d = pd.read_csv(f); m = pd.read_csv(os.path.join(RES, f"cells_meta_amend8_{group}{SFX}.csv"))
    d = d.merge(m[["id", "fam", "L", "kind", "max_drop", "exact_J"]], left_on="cell", right_on="id")
    d["decl_any"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
    return d, m

def minL(pw, grid, key="power"):
    for i in range(len(grid)):
        if all((l in pw) and pw[l] >= 0.8 for l in grid[i:]): return i
    return len(grid)

def verdicts(d, wpn, gated=True):
    out = {}
    for e in ("E7", "E8"):
        x = d[(d.est == e) & (d.wpn == wpn)]
        ctl = x[x.kind == "control"]
        ok = ctl[(ctl.abstain == 0) & (ctl.nowin.fillna(0) == 0)]
        k, n = int(ok.cov_simult.sum()), len(ok)
        per = ok.groupby("cell").cov_simult.agg(["sum", "count"]); badc = [int(c) for c, r in per.iterrows() if cp(r["sum"], r["count"])[1] < .93]
        Cp = n > 0 and k / n >= .93 and not badc
        fk = int(ok.decl_any.sum()); pf = ok.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(c) for c, r in pf.iterrows() if cp(r["sum"], r["count"])[0] > .05]
        Fp = fk / n <= .05 and not badf
        cx = x[x.kind == "crossover"].copy(); cx["real"] = [REAL(f, l) for f, l in zip(cx.fam, cx.L)]
        rx = cx[cx.real & (cx.abstain == 0) & (cx.nowin.fillna(0) == 0)]
        xk, xn = int(rx.cov_simult.sum()), len(rx); Xp = xn > 0 and xk / xn >= .90
        # M
        cx["det"] = ((cx.decl_neg == 1)).astype(float)
        pw = cx.groupby(["fam", "L"]).det.agg(["sum", "count"]).reset_index(); pw["power"] = pw["sum"] / pw["count"]
        Ms = {}
        for fam, grid in G.items():
            g = pw[pw.fam == fam]; P = dict(zip(g.L, g.power)); lo = {l: cp(r["sum"], r["count"])[0] for l, (_, r) in zip(g.L, g.iterrows())}
            i = minL(P, grid); Lm = grid[i] if i < len(grid) else None
            strict = None
            for j in range(len(grid)):
                if all((l in lo) and lo[l] >= .8 for l in grid[j:]): strict = grid[j]; break
            Ms[fam] = dict(minL=Lm, strict=strict, ok=(Lm is not None and Lm <= LIM[fam]), power=" ".join(f"{int(l)}:{p:.2f}" for l, p in zip(g.L, g.power)))
        out[e] = dict(C=Cp, Cn=fmt(k, n), Cbad=badc, F=Fp, Fn=fmt(fk, n), Fbad=badf, X=Xp, Xn=fmt(xk, xn), M=Ms,
                      passes=bool(Cp and Fp and Xp), beats=bool(Cp and Fp and Xp and all(v["ok"] for v in Ms.values())),
                      infeas=int(((x.infeasible.fillna(0)) == 1).sum()), ntot=len(x))
        if gated:
            xs = rx.groupby("cell").cov_simult.agg(["sum", "count"]); out[e]["xcells"] = " ".join(f"{int(c)}:{int(r['sum'])}/{int(r['count'])}" for c, r in xs.iterrows())
    return out

def main():
    d, m = load("core")
    if d is None: print("no core file"); return
    V = verdicts(d, 32); lines = []
    P = lambda s="": (print(s), lines.append(str(s)))
    for e in ("E7", "E8"):
        v = V[e]
        P(f"=== {e} (wpn 32)")
        P(f"(C) control simultaneous coverage {v['Cn']}; control cells with CP-upper<93%: {v['Cbad']} -> {'PASS' if v['C'] else 'FAIL'}")
        P(f"(F) FPR margin 0.25 any-sign {v['Fn']}; cells with CP-lower>5%: {v['Fbad']} -> {'PASS' if v['F'] else 'FAIL'}")
        P(f"(X) realised-crossover coverage of exact d {v['Xn']} -> {'PASS' if v['X'] else 'FAIL'}  per cell: {v.get('xcells')}")
        for fam, mm in v["M"].items():
            P(f"(M) {fam}: min L (power>=0.8 here and above) = {mm['minL']}; strict(CP-lower>=0.8) = {mm['strict']}; bar <= {LIM[fam]} -> {'PASS' if mm['ok'] else 'FAIL'}; power {mm['power']}")
        P(f"infeasible (E7: no measure inside box / E8 n/a) datasets: {v['infeas']}/{v['ntot']}")
        P(f"VERDICT {e}: {'beats the trade-off' if v['beats'] else ('passes (C,F,X) but not M' if v['passes'] else 'does not pass')}")
    # alarm
    e6 = d[d.est == "E6"].drop(columns=["a1p_alarm", "a1p_maxgap", "a1p_infeasible"], errors="ignore"); e7 = d[d.est == "E7"][["cell", "ds", "a1p_alarm", "a1p_maxgap", "a1p_infeasible"]]
    a = e6.merge(e7, on=["cell", "ds"]); a = a[a.abstain == 0]
    ctl = a[a.kind == "control"]; ck, cn = int(ctl.a1p_alarm.sum()), len(ctl)
    mis = a[a.kind.isin(["twostep", "brush", "gasket"])].copy(); mis["miss"] = 1 - mis.cov_simult
    ms = mis[mis.miss == 1]; mk, mn = int(ms.a1p_alarm.sum()), len(ms)
    P("=== A1' alarm (E6 band vs E7 band disjoint at any window time)")
    P(f"alarm rate on controls: {fmt(ck, cn)} -> {'PASS' if cn and ck / cn <= .05 else 'FAIL'} (bar <= 5%)")
    pc = ctl.groupby("cell").a1p_alarm.agg(["sum", "count"]); P("  per control cell: " + " ".join(f"{int(c)}:{int(r['sum'])}/{int(r['count'])}" for c, r in pc.iterrows()))
    P(f"P(alarm | E6 misses) on misspecified (twostep, brush, gasket): {fmt(mk, mn)} -> {('PASS' if mn and mk / mn >= .5 else 'FAIL') if mn >= 30 else 'too few failures (<30) to say'} (bar >= 50%)")
    for kd in ("twostep", "brush", "gasket"):
        s = mis[mis.kind == kd]; sm = s[s.miss == 1]; P(f"  {kd}: E6 misses {len(sm)}/{len(s)}; alarm among misses {int(sm.a1p_alarm.sum())}/{len(sm)}; alarm overall {int(s.a1p_alarm.sum())}/{len(s)}")
    P(f"E7-infeasible among alarm candidates: controls {int(ctl.a1p_infeasible.sum())}, misspec {int(mis.a1p_infeasible.sum())}")
    # ungated: gate, LP time, widths, misspecified coverage
    x7 = d[d.est == "E7"]
    P("=== ungated")
    P(f"atom gate: datasets using 400->799 only {int((x7.natoms == 799).sum())}, needing 1597 {int((x7.natoms == 1597).sum())}, final gate failure (diff>0.005 at 799->1597) {int(x7.gate_fail.fillna(0).sum())} of {len(x7)}; max gate1 diff (finite) {np.nanmax(np.where(np.isfinite(x7.gate1), x7.gate1, np.nan)):.4f}")
    gf = x7[x7.gate_fail.fillna(0) == 1]
    P(f"gate-failed datasets: {len(gf)}/{len(x7)} ({100*len(gf)/len(x7):.2f}%); their 799->1597 endpoint shift: median {gf.gate2.replace(np.inf, np.nan).median():.4f}, max finite {gf.gate2.replace(np.inf, np.nan).max():.4f}, infinite (feasibility changed) {int(np.isinf(gf.gate2).sum())}; by cell max share {(x7.groupby('cell').gate_fail.mean()).max():.3f}")
    fe = x7[(x7.infeasible.fillna(0) == 0) & (x7.nowin.fillna(0) == 0)]
    for nm, sel in (("controls", fe[fe.kind == "control"]), ("realised crossovers", fe[[REAL(f, l) and k == "crossover" for f, l, k in zip(fe.fam, fe.L, fe.kind)]])):
        P(f"E7 coverage among FEASIBLE datasets only, {nm}: {fmt(int(sel.cov_simult.sum()), len(sel))}")
    nf = x7[(x7.gate_fail.fillna(0) == 0) & (x7.infeasible.fillna(0) == 0)]
    for nm, sel in (("controls", nf[nf.kind == "control"]), ("realised crossovers", nf[[REAL(f, l) and k == "crossover" for f, l, k in zip(nf.fam, nf.L, nf.kind)]])):
        P(f"E7 coverage among gate-passing, feasible datasets only, {nm}: {fmt(int(sel.cov_simult.sum()), len(sel))}")
    P(f"E7 infeasible by kind: " + str(x7[x7.infeasible.fillna(0) == 1].groupby("kind").size().to_dict()))
    P(f"LP solve time per dataset (E7, all window times, 2-3 grids): mean {x7.lp_seconds.mean():.2f}s, median {x7.lp_seconds.median():.2f}s, max {x7.lp_seconds.max():.1f}s")
    cells = d.drop_duplicates("cell").set_index("cell")
    rows = []
    e4 = []
    for f_, flt in (("raw_amend6_orig.csv", lambda r: r.cell < 200), ("raw_amend6_big2d.csv", lambda r: r.cell >= 400), ("raw_amend7.csv", lambda r: (r.cell >= 700))):
        p = os.path.join(RES, f_)
        if os.path.exists(p):
            r = pd.read_csv(p); r = r[(r.est == "E4c") & (r.abstain == 0)]; r = r[flt(r)]; e4.append(r[["cell", "hw", "cov_simult"]])
    e4 = pd.concat(e4).groupby("cell").agg(e4_hw=("hw", "mean"), e4_cov=("cov_simult", "mean")) if e4 else pd.DataFrame(columns=["e4_hw", "e4_cov"])
    for c, g in d[d.wpn == 32].groupby("cell"):
        r = dict(cell=c, fam=g.fam.iloc[0], L=g.L.iloc[0], kind=g.kind.iloc[0])
        for e in ("E7", "E8"):
            s = g[(g.est == e) & (g.nowin.fillna(0) == 0)]
            r[f"{e}_cov"] = round(s.cov_simult.mean(), 3); r[f"{e}_hw"] = round(s.hw.mean(), 3); r[f"{e}_decl_neg"] = round(s.decl_neg.mean(), 3); r[f"{e}_decl_any"] = round(s.decl_any.mean(), 3)
        s6 = g[(g.est == "E6") & (g.abstain == 0)]
        if len(s6): r["E6_cov"] = round(s6.cov_simult.mean(), 3)
        rows.append(r)
    T = pd.DataFrame(rows).merge(e4, left_on="cell", right_index=True, how="left").round(3); T.to_csv(RES + f"/amend8_cell_table{SFX}.csv", index=False)
    P(T.to_string(index=False))
    for gname in ("w8", "w128"):
        dd, _ = load(gname)
        if dd is None: continue
        w = int(gname[1:]); Vw = verdicts(dd, w, gated=False); P(f"=== ungated wpn={w}")
        for e in ("E7", "E8"):
            v = Vw[e]; P(f"{e}: C {v['Cn']} {'PASS' if v['C'] else 'FAIL'}; F {v['Fn']} {'PASS' if v['F'] else 'FAIL'}; X {v['Xn']} {'PASS' if v['X'] else 'FAIL'}; " + "; ".join(f"M {f}: {mm['minL']} (strict {mm['strict']}) power {mm['power']}" for f, mm in v["M"].items()))
    open(RES + f"/amend8_report{SFX}.txt", "w").write("\n".join(lines))

if __name__ == "__main__": main()
