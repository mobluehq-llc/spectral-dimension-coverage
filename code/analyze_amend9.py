#!/usr/bin/env python3
"""P22 Amendment 9 analysis. Usage: analyze_amend9.py [--smoke]. Reads results/raw_amend9_{a,amis,b8,b32,b128,c}.csv (+ earlier raw files for original-design
verdicts) and writes results/amend9_*.csv and amend9_report.txt. Verdict lines are printed literally."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
SFX = "_smoke" if "--smoke" in sys.argv else ""
G = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768]}
LIM = {"aniso3_w8": 160, "aniso2_w16": 512}
REAL = lambda fam, l: (fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48)
LINES = []
def P(s=""): print(s); LINES.append(str(s))
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n):
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)
def rd(f):
    p = os.path.join(RES, f); return pd.read_csv(p) if os.path.exists(p) else None

def norm(d, meta):
    d = d.copy()
    for c in ("nowin", "decl_neg", "decl_pos", "infeasible"):
        if c not in d.columns: d[c] = 0
    d["nowin"] = d.nowin.fillna(0); d["decl_neg"] = d.decl_neg.fillna(0); d["decl_pos"] = d.decl_pos.fillna(0)
    d["decl_any"] = ((d.decl_neg + d.decl_pos) > 0).astype(int)
    d = d.merge(meta[["id", "fam", "L", "kind"]], left_on="cell", right_on="id", how="inner")
    d["real"] = [REAL(f, l) and k == "crossover" for f, l, k in zip(d.fam, d.L, d.kind)]
    d["valid"] = (d.abstain == 0) & (d.nowin == 0)
    return d

def calib(x):
    """C, F on the control rows of x (one arm, one budget). Returns dict."""
    c = x[x.kind == "control"]; ok = c[c.valid]
    k, n = int(ok.cov_simult.sum()), len(ok)
    per = ok.groupby("cell").cov_simult.agg(["sum", "count"]); badc = [int(i) for i, r in per.iterrows() if cp(r["sum"], r["count"])[1] < .93]
    fk = int(ok.decl_any.sum()); pf = ok.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(i) for i, r in pf.iterrows() if cp(r["sum"], r["count"])[0] > .05]
    Cp = n > 0 and k / n >= .93 and not badc; Fp = n > 0 and fk / n <= .05 and not badf
    r = x[x.real]; rv = r[r.valid]; xk, xn = int(rv.cov_simult.sum()), len(rv)
    return dict(k=k, n=n, abst=round(1 - n / max(len(c), 1), 3), C=bool(Cp), Cn=fmt(k, n), Cbad=badc, F=bool(Fp), Fn=fmt(fk, n), Fbad=badf,
                X=bool(xn > 0 and xk / xn >= .90), Xn=fmt(xk, xn), xk=xk, xn=xn)

def power(x):
    cx = x[x.kind == "crossover"].copy(); cx["det"] = ((cx.valid) & (cx.decl_neg == 1)).astype(float)
    return cx.groupby(["fam", "L"]).det.agg(["sum", "count"]).reset_index().assign(power=lambda t: t["sum"] / t["count"])

def minL_of(x, restrict=False):
    pw = power(x); out = {}
    for fam, grid in G.items():
        g = pw[pw.fam == fam]; Pw = dict(zip(g.L, g.power))
        gr = [l for l in grid if l in Pw] if restrict else grid
        idx = len(gr)
        for i in range(len(gr)):
            if all((l in Pw) and Pw[l] >= .8 for l in gr[i:]): idx = i; break
        out[fam] = (gr[idx] if idx < len(gr) else np.inf)
        out[fam + "_power"] = " ".join(f"{int(l)}:{p:.2f}" for l, p in zip(g.L, g.power))
    return out

# ====================================================================== 9a
def part_a():
    d = rd(f"raw_amend9_a{SFX}.csv"); m = rd(f"cells_meta_amend9_a{SFX}.csv")
    if d is None: P("9a: no data"); return
    d = norm(d, m)
    P("=" * 20 + " 9a: E7j joint change test (same cells and data as Amendment 8, wpn 32)")
    c = d[d.kind == "control"]; ok = c[c.valid]
    fk, n = int(ok.decl_any.sum()), len(ok)
    pf = ok.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(i) for i, r in pf.iterrows() if cp(r["sum"], r["count"])[0] > .05]
    Fp = fk / n <= .05 and not badf
    P(f"(F) E7j FPR (margin 0.25, any sign) on controls {fmt(fk, n)}; cells with CP-lower > 5%: {badf} -> {'PASS' if Fp else 'FAIL'}")
    P("  per control cell: " + " ".join(f"{int(i)}:{int(r['sum'])}/{int(r['count'])}" for i, r in pf.iterrows()))
    M = minL_of(d); Mok = {}
    for fam in G:
        v = M[fam]; Mok[fam] = bool(np.isfinite(v) and v <= LIM[fam])
        P(f"(M) {fam}: min L (power>=0.8 here and above) = {'not reached' if not np.isfinite(v) else int(v)}; bar <= {LIM[fam]} -> {'PASS' if Mok[fam] else 'FAIL'}; power {M[fam + '_power']}")
    beats = bool(Fp and all(Mok.values()))
    P(f"VERDICT E7j: {'E7j beats the trade-off' if beats else 'E7j does NOT beat the trade-off (F ' + ('pass' if Fp else 'FAIL') + ', M ' + ('pass' if all(Mok.values()) else 'FAIL') + ')'}")
    # ungated
    P("--- 9a ungated")
    e7 = d.copy(); e7["det7"] = ((e7.valid) & (e7.e7_decl_neg.fillna(0) == 1)).astype(float)
    cx = e7[e7.kind == "crossover"]; pw7 = cx.groupby(["fam", "L"]).det7.mean().reset_index()
    for fam in G:
        g = pw7[pw7.fam == fam]; pj = power(d); pj = pj[pj.fam == fam]
        P(f"{fam} power E7 (interval arithmetic, recomputed) vs E7j: " + " ".join(f"{int(l)}:{a:.2f}/{b:.2f}" for l, a, b in zip(g.L, g.det7, pj.power)))
    P(f"datasets: {len(d)}; no-window {int((d.nowin == 1).sum())}; E7-infeasible (no measure in box; no declaration) {int(d.infeasible.fillna(0).sum())}; E7 gate failures {int(d.e7_gate_fail.fillna(0).sum())}")
    v = d[d.valid & (d.infeasible.fillna(0) == 0)]
    P(f"E7j declared whenever E7 declared (must be 0 violations): {int(v.viol_vs_e7.fillna(0).sum())}")
    P(f"E7j gate: datasets with any decision re-run at twice the atoms {int((v.gate_rechecked > 0).sum())}; with disagreement (decision withheld) {int((v.gate_flag > 0).sum())}; datasets with an undecided (tolerance/unresolved) (j,sign) {int((v.n_amb > 0).sum())}; unresolved LPs {int(v.nbad.sum())}; LPs {int(v.n_lp.sum())}")
    P(f"time per dataset: E7 LP mean {d.lp_seconds.mean():.2f}s; E7j extra mean {d.j_seconds.mean():.2f}s, median {d.j_seconds.median():.2f}s, max {d.j_seconds.max():.1f}s")
    # reproduction of E7 against Amendment 8
    r8 = rd("raw_amend8_core.csv")
    if r8 is not None:
        r8 = r8[r8.est == "E7"][["cell", "ds", "decl_neg", "natoms", "lp_status"]].rename(columns={"decl_neg": "decl_neg8"})
        mm = d.merge(r8, on=["cell", "ds"])
        P(f"reproduction of the Amendment 8 E7 band (same data/seeds): decl_neg equal {int((mm.e7_decl_neg == mm.decl_neg8).sum())}/{len(mm)}; atoms equal {int((mm.e7_natoms == mm.natoms).sum())}/{len(mm)}; status equal {int((mm.e7_status == mm.lp_status).sum())}/{len(mm)}")
    t = d[d.kind == "crossover"].groupby(["fam", "L"]).agg(E7j_drop=("decl_neg", "mean"), E7_drop=("e7_decl_neg", "mean"), E7j_rise=("decl_pos", "mean")).round(3)
    t.to_csv(os.path.join(RES, f"amend9a_power{SFX}.csv"))
    dm = rd(f"raw_amend9_amis{SFX}.csv")
    if dm is not None:
        mm_ = rd(f"cells_meta_amend9_amis{SFX}.csv"); dm = norm(dm, mm_)
        P("misspecified cells (ungated): " + "; ".join(f"{k} {g.fam.iloc[0]} L={g.L.iloc[0]}: drop {g.decl_neg.mean():.3f} rise {g.decl_pos.mean():.3f}" for k, g in dm[dm.valid].groupby("cell")))

# ====================================================================== 9b
ARMS = ["E0", "E1q", "E4c", "E4r", "E6", "E7", "E8"]
def orig_rows(arm, meta):
    """original-design rows (wpn 32) for arm from earlier raw files on cells 0-26, 400-402, 700-706."""
    srcs = {"E0": [("raw.csv", "E0", None), ("raw_amend7.csv", "E0", 32)], "E1q": [("raw_amend.csv", "E1q", None), ("raw_amend7.csv", "E1q", 32)],
            "E4c": [("raw_amend3.csv", "E4c", None), ("raw_amend6_big2d.csv", "E4c", None), ("raw_amend7.csv", "E4c", 32)],
            "E4r": [("raw_amend3.csv", "E4r", None), ("raw_amend6_big2d.csv", "E4r", None), ("raw_amend7.csv", "E4r", 32)],
            "E6": [("raw_amend5.csv", "E6", None), ("raw_amend6_big2d.csv", "E6", None), ("raw_amend7.csv", "E6", 32)],
            "E7": [("raw_amend8_core.csv", "E7", 32)], "E8": [("raw_amend8_core.csv", "E8", 32)]}
    parts = []
    for f, e, w in srcs[arm]:
        r = rd(f)
        if r is None: continue
        r = r[r.est == e]
        if w is not None and "wpn" in r.columns: r = r[r.wpn == w]
        r = r[((r.cell <= 26) | ((r.cell >= 400) & (r.cell <= 402)) | ((r.cell >= 700) & (r.cell <= 706)))]
        parts.append(r)
    if not parts: return None
    r = pd.concat(parts); r["est"] = arm
    return norm(r, meta)

def part_b():
    P("=" * 20 + " 9b: frozen confirmatory replication (code at git 1ba9ca6, fresh seeds [22, id, 2000 + ds])")
    D = {}
    for w in (8, 32, 128):
        d = rd(f"raw_amend9_b{w}{SFX}.csv"); m = rd(f"cells_meta_amend9_b{w}{SFX}.csv")
        if d is not None: D[w] = norm(d, m)
    if not D: P("9b: no data"); return
    rows = []; R = {}
    for w, d in D.items():
        for a in ARMS:
            x = d[d.est == a]
            if len(x) == 0: continue
            c = calib(x); R[(w, a)] = c
            rows.append(dict(wpn=w, arm=a, ctrl_n=c["n"], ctrl_abstain=c["abst"], C=c["C"], C_cov=c["Cn"], C_bad_cells=str(c["Cbad"]), F=c["F"], F_fpr=c["Fn"], F_bad_cells=str(c["Fbad"]), X=c["X"], X_cov=c["Xn"]))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(RES, f"amend9_b_table{SFX}.csv"), index=False)
    P(T.to_string(index=False))
    gates = []
    def G_(name, ok, txt): gates.append(ok); P(f"GATE {name}: {'PASS' if ok else 'FAIL'}  {txt}")
    # (i) E1q realised crossover coverage ordering
    if all((w, "E1q") in R for w in (8, 32, 128)):
        k = {w: (R[(w, "E1q")]["xk"], R[(w, "E1q")]["xn"]) for w in (8, 32, 128)}
        pt = {w: k[w][0] / max(k[w][1], 1) for w in k}
        G_("9b(i) E1q realised-crossover coverage wpn128 < wpn32 < wpn8", pt[128] < pt[32] < pt[8], "; ".join(f"wpn {w}: {fmt(*k[w])}" for w in (128, 32, 8)))
        G_("9b(iii) E1q realised-crossover coverage at wpn 128 <= 20%", pt[128] <= .20, f"wpn 128: {fmt(*k[128])}")
    # (ii)
    for a in ("E4c", "E4r", "E7", "E8"):
        for w in (8, 32, 128):
            if (w, a) in R:
                c = R[(w, a)]; G_(f"9b(ii) {a} controls >= 93% pooled at wpn {w}", c["n"] > 0 and c["k"] / c["n"] >= .93, c["Cn"] + f"; cell rule (no control cell CP-upper<93%): {'ok' if not c['Cbad'] else 'cells ' + str(c['Cbad'])}")
    # (iv) C/F verdicts of Amendments 1-8 reproduce at wpn 32 (recorded: every one of the 7 arms passes C and F)
    meta = rd("cells_meta_amend8_core.csv")
    for a in ARMS:
        if (32, a) not in R: continue
        c = R[(32, a)]
        o = orig_rows(a, meta); oc = calib(o) if o is not None else None
        rec = "recorded C pass / F pass"
        so = f"; same function on the ORIGINAL data: C {'pass' if oc['C'] else 'FAIL'} ({oc['Cn']}), F {'pass' if oc['F'] else 'FAIL'} ({oc['Fn']})" if oc else ""
        G_(f"9b(iv) {a} C at wpn 32 reproduces ({rec})", c["C"], f"fresh C {'pass' if c['C'] else 'FAIL'}: {c['Cn']} cells<93%: {c['Cbad']}{so}")
        G_(f"9b(iv) {a} F at wpn 32 reproduces ({rec})", c["F"], f"fresh F {'pass' if c['F'] else 'FAIL'}: {c['Fn']} cells>5%: {c['Fbad']}")
    P(f"VERDICT 9b: {'ALL gated predictions hold' if all(gates) else 'NOT all gated predictions hold: ' + str(sum(1 for g in gates if not g)) + ' of ' + str(len(gates)) + ' fail'}")
    # ungated min L table
    P("--- 9b ungated: min L (power>=0.8 here and above; abstention = miss) on the realised-crossover cells run here (grid = the 9 sizes per family; sizes below the first are not run)")
    ml = []
    for w, d in D.items():
        for a in ARMS:
            x = d[d.est == a]
            if len(x):
                m_ = minL_of(x, restrict=True); ml.append(dict(wpn=w, arm=a, **{f: ("not reached" if not np.isfinite(m_[f]) else int(m_[f])) for f in G}, xcov=R[(w, a)]["Xn"]))
    ML = pd.DataFrame(ml); ML.to_csv(os.path.join(RES, f"amend9_b_minL{SFX}.csv"), index=False); P(ML.to_string(index=False))

# ====================================================================== 9c
ARMS_C = ["E0", "E1q", "E4c", "E6", "E7", "E8"]
def part_c():
    P("=" * 20 + " 9c: second design (grid 24*2^(j/4), count>=100, N q>=10, J>=8; wpn 32; fresh seeds)")
    d = rd(f"raw_amend9_c{SFX}.csv"); m = rd(f"cells_meta_amend9_c{SFX}.csv")
    if d is None: P("9c: no data"); return
    d = norm(d, m); meta = rd("cells_meta_amend8_core.csv")
    out = []; ok_all = []; ML = {}
    for a in ARMS_C:
        x = d[d.est == a]; c = calib(x); o = orig_rows(a, meta); oc = calib(o)
        # original design restricted to the same cell family (controls + realised crossovers present)
        ML[a] = (minL_of(x, restrict=False), minL_of(o, restrict=True))
        for k_ in ("C", "F", "X"):
            same = c[k_] == oc[k_]; ok_all.append(same)
            P(f"GATE 9c {a} {k_}: new {'PASS' if c[k_] else 'FAIL'} ({ {'C': c['Cn'], 'F': c['Fn'], 'X': c['Xn']}[k_] }) vs original {'PASS' if oc[k_] else 'FAIL'} ({ {'C': oc['Cn'], 'F': oc['Fn'], 'X': oc['Xn']}[k_] }) -> {'EQUAL' if same else 'DIFFERENT'}")
        out.append(dict(arm=a, new_C=c["Cn"], new_F=c["Fn"], new_X=c["Xn"], new_ctrl_abstain=c["abst"], orig_C=oc["Cn"], orig_F=oc["Fn"], orig_X=oc["Xn"],
                        minL3_new=ML[a][0]["aniso3_w8"], minL3_orig=ML[a][1]["aniso3_w8"], minL2_new=ML[a][0]["aniso2_w16"], minL2_orig=ML[a][1]["aniso2_w16"]))
    pd.DataFrame(out).to_csv(os.path.join(RES, f"amend9_c_table{SFX}.csv"), index=False)
    f_ = lambda v: "not reached" if not np.isfinite(v) else str(int(v))
    for a in ARMS_C:
        P(f"min L {a}: 3D new {f_(ML[a][0]['aniso3_w8'])} orig {f_(ML[a][1]['aniso3_w8'])}; 2D new {f_(ML[a][0]['aniso2_w16'])} orig {f_(ML[a][1]['aniso2_w16'])}; power new 3D {ML[a][0]['aniso3_w8_power']}; 2D {ML[a][0]['aniso2_w16_power']}")
    sgn = lambda a, b: 0 if a == b else (-1 if a < b else 1)
    for o in ("E1q", "E6", "E7"):
        for fam in G:
            so = sgn(ML["E4c"][1][fam], ML[o][1][fam]); sn = sgn(ML["E4c"][0][fam], ML[o][0][fam])
            ok = not (so * sn == -1); strict = (so == sn); ok_all.append(ok)
            nm = {-1: "E4c better (smaller L)", 0: "tie", 1: "E4c worse (larger L)"}
            P(f"GATE 9c order E4c vs {o} {fam}: original {nm[so]} (E4c {f_(ML['E4c'][1][fam])}, {o} {f_(ML[o][1][fam])}); new {nm[sn]} (E4c {f_(ML['E4c'][0][fam])}, {o} {f_(ML[o][0][fam])}) -> {'UNCHANGED (ties allowed)' if ok else 'REVERSED'}{'' if strict or not ok else ' [not strictly identical: a tie involved]'}")
    P(f"VERDICT 9c: {'second design reproduces every verdict and ordering (gate PASS)' if all(ok_all) else 'gate FAIL: ' + str(sum(1 for g in ok_all if not g)) + ' of ' + str(len(ok_all)) + ' checks differ'}")
    sub = d[d.valid]; P("--- 9c ungated: window/abstention  " + "; ".join(f"{a} control abstain {calib(d[d.est == a])['abst']}" for a in ARMS_C))

if __name__ == "__main__":
    part_a(); part_b(); part_c()
    open(os.path.join(RES, f"amend9_report{SFX}.txt"), "w").write("\n".join(LINES))
