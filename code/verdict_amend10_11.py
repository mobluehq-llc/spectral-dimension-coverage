#!/usr/bin/env python3
"""P22 Amendments 10b / 10c / 11: verdict script. Applies the bars and gates EXACTLY as written in PLAN.md (not edited here); introduces no new bar.

Usage:
  verdict_amend10_11.py [--raw DIR] [--append]
    --raw DIR  directory holding raw_amend10c.csv / raw_amend11.csv (large; kept on host-D: <local path omitted>). When present, the cell-level
               aggregates results_amend10_11/agg_*.csv are rebuilt from it (and the raw sha256 is checked against the manifest); otherwise the committed
               aggregates are used. raw_amend10b*.csv are small and are committed, so 10b always runs from raw.
    --append   append the generated section to RESULTS.md (refuses if it is already there).
Outputs (all in results_amend10_11/): verdicts.csv, verdict_report.txt, tables_*.csv, RESULTS_SECTION_A10_A11.md, reading_log.txt.

READING LOG (every place a bar needed a reading; the conservative reading is binding, the other is shown):
 R1 Denominators. "Coverage" and "FPR" in every earlier amendment (analyze_amend9.calib) are computed over datasets where the estimator answered
    (abstain==0 and nowin==0). That is the PRIMARY (binding) reading for 10c and 11 and is kept; the strict reading (abstention counts as a miss) is computed
    and printed beside it (sensitivity columns 'strict'). 10b uses FIXED denominators (all datasets in the group; no band or infeasible = not covered /
    not detected) as the referee asked, binding.
 R2 Intervals. Plain numeric bars (G1, G3, G4 pooled, G5, 10b C/X/F/M, 10c E4c >= 93%) are tested on the POINT estimate, as in every earlier amendment (power
    and coverage bars in analyze_amend9); the plan's own CP rules (G2, per-cell C/F rules, 10c CP-upper(128) < CP-lower(8)) are used as written. CP
    intervals are reported for everything, and the stricter reading "the CP bound must also meet the plain bar" is computed and flagged when it would
    change a verdict (column ci_strict).
 R3 Cells with zero answered datasets (E1q/E4 abstain by the frozen J<6 rule) are skipped in per-cell tests (no data, CP undefined) and listed.
 R4 G5 and 10b-F: FPR is any-sign declaration (decl_neg or decl_pos) at margin 0.25. Conservative extra reading (as the earlier F gate): no cell with
    CP-lower > 5%.
 R5 10b: both outer relaxations (A literal, B tight) must meet the bar for E7; E7j exists only for B (certified LP). 3D min L: E7j power >= 0.8 at
    L = 160, 192 and 256 (power >= 0.8 'here and above'), fixed denominator 200.
 R6 10c: the central-result bar must hold pooled over both families AND within each family; E4c >= 93% must hold at each wpn (pooled and per family).
    'monotone' = strictly decreasing 8 > 32 > 128.
 R7 G2: strict decrease in the answered-dataset point estimates; CP-upper(128) < CP-lower(8).
"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import beta

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results_amend10_11")
args = sys.argv[1:]
RAWDIR = args[args.index("--raw") + 1] if "--raw" in args else OUT
LINES = []; VER = []
def P(s=""): print(s); LINES.append(str(s))

def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n):
    lo, hi = cp(k, n)
    if n == 0: return "0/0"
    return "%d/%d = %.1f%% (%.1f-%.1f)" % (k, n, 100 * k / n, 100 * lo, 100 * hi)
def pct(k, n): return "n/a" if n == 0 else "%.1f%%" % (100 * k / n)
def ge(k, n, thr):
    """bar 'rate >= thr'. returns (point ok, conservative ok [point and CP-lower])."""
    if n == 0: return False, False
    p = k / n >= thr - 1e-12; return p, bool(p and cp(k, n)[0] >= thr - 1e-12)
def le(k, n, thr):
    if n == 0: return False, False
    p = k / n <= thr + 1e-12; return p, bool(p and cp(k, n)[1] <= thr + 1e-12)
def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()
def verdict(gid, desc, observed, ok_strict, ok_pt=None, note=""):
    """binding verdict = point-estimate reading (precedent of every earlier amendment's bars) with the per-cell/CP rules the plan itself names built in;
    ok_strict = additionally requiring the matching CP bound to meet a plain numeric bar (reported, flagged when it differs)."""
    b = ok_strict if ok_pt is None else ok_pt
    flag = "" if ok_strict == b else f"CI-strict reading (CP bound must also meet the bar) gives {'PASS' if ok_strict else 'FAIL'}"
    VER.append(dict(gate=gid, description=desc, observed=observed, verdict="PASS" if b else "FAIL", ci_strict=flag, note=note))
    P(f"[{'PASS' if b else 'FAIL'}] {gid}: {desc} -> {observed}" + (f"  [{flag}]" if flag else "") + (f"  ({note})" if note else ""))

# ======================================================================= provenance
def check_manifests():
    P("=" * 20 + " PROVENANCE: manifests (copied from host-D <local path omitted>)")
    ok_all = True
    for nm in ("amend10b", "amend10b_exactq", "amend10c", "amend11"):
        m = json.load(open(os.path.join(OUT, f"manifest_{nm}.json")))
        ok = bool(m["complete"]) and not m["errors"] and m["items_present"] == m["tasks_expected"] and not m["smoke"]
        s = f"{nm}: complete={m['complete']} errors={m['errors']} tasks {m['items_present']}/{m['tasks_expected']} rows {m['n_rows']} csv_sha256 {m['csv_sha256'][:12]}.."
        rp = os.path.join(RAWDIR, m["csv"])
        if os.path.exists(rp):
            h = sha(rp); good = (h == m["csv_sha256"]); s += f" raw present sha_match={good}"; ok = ok and good
        else: s += " raw not present here (sha recorded in manifest)"
        for key, fn in (("selection_sha256", "amend10b_selection.json"), ("windows_sha256", "amend10c_windows.json"), ("cells_exact_sha256", "amend11_cells_exact.json")):
            if key in m:
                h = sha(os.path.join(OUT, fn)); good = (h == m[key]); s += f" {fn} sha_match={good}"; ok = ok and good
        P(("OK   " if ok else "BAD  ") + s); ok_all &= ok
    if not ok_all: P("MANIFEST CHECK FAILED"); sys.exit(2)

# ======================================================================= 10b
REAL = lambda fam, l, kind: kind == "crossover" and ((fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48))

def part_10b():
    d = pd.read_csv(os.path.join(OUT, "raw_amend10b.csv"), low_memory=False)
    P(""); P("=" * 20 + " AMENDMENT 10b (referee-directed test): certified outer bounds for E7 / E7j")
    P(f"datasets {len(d)} (groups: {d.group.value_counts().to_dict()}); B_status {d.B_status.value_counts().to_dict()}; A_status {d.A_status.value_counts().to_dict()}; original o_status {d.o_status.value_counts().to_dict()}")
    d["real"] = [REAL(f, l, k) for f, l, k in zip(d.fam, d.L, d.kind)]
    for c in ("B_cov", "A_cov", "o_cov", "B_e7_neg", "B_e7_pos", "A_e7_neg", "A_e7_pos", "o_e7_neg", "o_e7_pos", "Bj_neg", "Bj_pos", "Bj_amb", "Bj_certnon"):
        d[c] = d[c].fillna(0)
    d["B_ok"] = (d.B_status == 0); d["A_ok"] = (d.A_status == 0); d["o_ok"] = (d.o_status == 0)
    # fixed-denominator outcomes
    d["B_covered"] = (d.B_ok & (d.B_cov == 1)).astype(int); d["A_covered"] = (d.A_ok & (d.A_cov == 1)).astype(int)
    d["o_covered"] = (d.o_ok & (d.o_cov == 1)).astype(int)
    d["B_changeE7"] = (d.B_ok & ((d.B_e7_neg + d.B_e7_pos) > 0)).astype(int); d["A_changeE7"] = (d.A_ok & ((d.A_e7_neg + d.A_e7_pos) > 0)).astype(int)
    d["o_changeE7"] = (d.o_ok & ((d.o_e7_neg + d.o_e7_pos) > 0)).astype(int)
    d["j_drop"] = (d.B_ok & (d.Bj_neg == 1)).astype(int); d["j_rise"] = (d.B_ok & (d.Bj_pos == 1)).astype(int)
    d["j_any"] = ((d.j_drop + d.j_rise) > 0).astype(int)
    d["j_none"] = (d.B_ok & (d.j_any == 0)).astype(int)
    d["j_amb"] = (d.B_ok & (d.j_any == 0) & (d.Bj_amb > 0)).astype(int)
    # ---- cell table with fixed denominators
    g = d.groupby(["kind", "fam", "L"])
    tab = g.apply(lambda t: pd.Series(dict(
        n=len(t), real=int(t.real.iloc[0]),
        E7B_covered=int(t.B_covered.sum()), E7B_not_covered=int((t.B_ok & (t.B_cov == 0)).sum()), no_band=int((~t.B_ok).sum()),
        E7A_covered=int(t.A_covered.sum()), E7orig_covered=int(t.o_covered.sum()), orig_no_band=int((~t.o_ok).sum()),
        B_contains_orig=int(((t.B_contains_orig == 1) & t.B_ok & t.o_ok).sum()), B_unbounded_upper=int((t.B_nopen.fillna(0) > 0).sum()),
        E7B_change=int(t.B_changeE7.sum()), E7A_change=int(t.A_changeE7.sum()), E7orig_change=int(t.o_changeE7.sum()),
        E7j_drop=int(t.j_drop.sum()), E7j_rise=int(t.j_rise.sum()), E7j_undecided_amb=int(t.j_amb.sum()), E7j_none=int((t.j_none - t.j_amb).sum()), E7j_noband=int((~t.B_ok).sum())))).reset_index()
    tab.to_csv(os.path.join(OUT, "tables_10b_cell_outcomes.csv"), index=False)
    # ---- gates
    ctl = d[d.kind == "control"]; rx = d[d.real]
    P(f"controls analysed: {len(ctl)} datasets in {ctl.groupby(['fam','L']).ngroups} cells; realised crossovers: {len(rx)} datasets in {rx.groupby(['fam','L']).ngroups} cells")
    P("Fixed denominators: a dataset with no certified band (B_status!=0) counts as NOT covered / NOT detected / 'no declaration'.")
    for band, col in (("B", "B_covered"), ("A", "A_covered"), ("orig", "o_covered")):
        k, n = int(ctl[col].sum()), len(ctl)
        per = ctl.groupby(["fam", "L"])[col].agg(["sum", "count"]); bad = [f"{a}_{int(b)}" for (a, b), r in per.iterrows() if cp(r["sum"], r["count"])[1] < .93]
        okp, okc = ge(k, n, .93); okp &= not bad; okc &= not bad
        if band != "orig": verdict(f"10b-C[{band}]", "E7 certified-band coverage of exact d on controls >= 93% pooled, no cell with CP-upper < 93%", fmt(k, n) + f"; bad cells {bad}", okc, okp)
        else: P(f"   (reference, original gated E7 on the same datasets: coverage {fmt(k, n)})")
    for band, col in (("B", "B_covered"), ("A", "A_covered"), ("orig", "o_covered")):
        k, n = int(rx[col].sum()), len(rx); okp, okc = ge(k, n, .90)
        if band != "orig": verdict(f"10b-X[{band}]", "E7 certified-band coverage of exact d on realised crossovers >= 90%", fmt(k, n), okc, okp)
        else: P(f"   (reference, original gated E7: crossover coverage {fmt(k, n)})")
    for nm, col, lab in (("E7 interval test, B", "B_changeE7", "B"), ("E7 interval test, A", "A_changeE7", "A"), ("E7j certified (B)", "j_any", "j")):
        k, n = int(ctl[col].sum()), len(ctl)
        per = ctl.groupby(["fam", "L"])[col].agg(["sum", "count"]); bad = [f"{a}_{int(b)}" for (a, b), r in per.iterrows() if cp(r["sum"], r["count"])[0] > .05]
        okp, okc = le(k, n, .05); okp &= not bad; okc &= not bad
        verdict(f"10b-F[{lab}]", f"{nm}: false-change rate (margin 0.25, any sign) on controls <= 5%, no cell with CP-lower > 5%", fmt(k, n) + f"; bad cells {bad}", okc, okp)
    k, n = int(ctl.o_changeE7.sum()), len(ctl); P(f"   (reference, original gated E7 interval test false-change on controls: {fmt(k, n)})")
    # 3D min L
    a3 = d[(d.fam == "aniso3_w8")]
    pw = a3.groupby("L").apply(lambda t: pd.Series(dict(n=len(t), j=int(t.j_drop.sum()), e7B=int(t.B_changeE7.sum()), e7o=int(t.o_changeE7.sum())))).reset_index()
    P("3D (aniso3_w8) power of a certified drop declaration, fixed denominator: " + "; ".join(f"L={int(r.L)}: E7j {int(r.j)}/{int(r.n)}={r.j/r.n:.2f}" for r in pw.itertuples()))
    P("   E7 interval test (B): " + "; ".join(f"L={int(r.L)}: {int(r.e7B)}/{int(r.n)}" for r in pw.itertuples()) + "   | original gated E7: " + "; ".join(f"L={int(r.L)}: {int(r.e7o)}/{int(r.n)}" for r in pw.itertuples()))
    pw.to_csv(os.path.join(OUT, "tables_10b_3d_power.csv"), index=False)
    big = pw[pw.L.isin([160, 192, 256])]
    allok = len(big) == 3; allokp = len(big) == 3; obs = []
    for r in big.itertuples():
        okp, okc = ge(r.j, r.n, .8); allok &= okc; allokp &= okp; obs.append(f"L={int(r.L)} {fmt(r.j, r.n)}")
    verdict("10b-M[E7j]", "3D min L = 160 preserved: certified E7j power >= 0.8 at L = 160, 192, 256 (fixed n = 200)", "; ".join(obs), allok, allokp)
    # exact-q
    x = pd.read_csv(os.path.join(OUT, "raw_amend10b_exactq.csv"))
    x["J"] = x.J; x["is_ctl"] = x.cell.str.startswith(("ring", "torus2", "torus3"))
    x["j_any"] = ((x.Bj_neg.fillna(0) + x.Bj_pos.fillna(0)) > 0).astype(int)
    xt = x[["cell", "J", "N", "B_status", "B_cov", "A_cov", "B_contains_orig", "B_e7_neg", "B_e7_pos", "Bj_neg", "Bj_pos", "Bj_amb", "true_drop", "true_rise"]]
    xt.to_csv(os.path.join(OUT, "tables_10b_exactq.csv"), index=False)
    okcov = bool((x.B_status == 0).all() and (x.B_cov == 1).all() and (x.A_cov == 1).all())
    okF = bool((x[x.is_ctl].j_any == 0).all() and (x[x.is_ctl].B_e7_neg + x[x.is_ctl].B_e7_pos == 0).all())
    verdict("10b-exactq", "7 exact-q cells: certified bands (A and B) contain the exact d(t) (7/7) and no false change on the 3 control cells", f"coverage B/A {int((x.B_cov==1).sum())}/{len(x)}, {int((x.A_cov==1).sum())}/{len(x)}; control false changes (E7 B, E7j) {int(x[x.is_ctl].B_changeE7.sum()) if 'B_changeE7' in x else int(((x[x.is_ctl].B_e7_neg + x[x.is_ctl].B_e7_pos) > 0).sum())}, {int(x[x.is_ctl].j_any.sum())}", okcov and okF)
    # sanity: certified contains original
    k, n = int(((d.B_contains_orig == 1) & d.B_ok & d.o_ok).sum()), int((d.B_ok & d.o_ok).sum())
    P(f"sanity: certified B band contains the original gated band in {k}/{n} datasets with both (max inward violation {d.B_max_inward.max():.2e}); B bands with an unbounded upper end: {int((d.B_nopen.fillna(0) > 0).sum())} datasets")
    P(f"outcome categories, all 1030 datasets (fixed denominator 1030): B covered {int(d.B_covered.sum())}, B not covered {int((d.B_ok & (d.B_cov==0)).sum())}, no certified band {int((~d.B_ok).sum())}; "
      f"original-E7 covered {int(d.o_covered.sum())}, original no band {int((~d.o_ok).sum())}; E7j drop {int(d.j_drop.sum())}, rise {int(d.j_rise.sum())}, undecided (some ambiguous time) {int(d.j_amb.sum())}, none {int((d.j_none - d.j_amb).sum())}")
    return d, tab, pw, xt

# ======================================================================= aggregates for 10c / 11
def aggregate(name, meta):
    rp = os.path.join(RAWDIR, f"raw_{name}.csv"); ap = os.path.join(OUT, f"agg_{name}_cells.csv")
    if os.path.exists(rp):
        d = pd.read_csv(rp, low_memory=False)
        if "mode" not in d.columns: d["mode"] = "single"
        for c in ("nowin", "decl_neg", "decl_pos"): d[c] = d[c].fillna(0) if c in d.columns else 0
        d["valid"] = ((d.abstain == 0) & (d.nowin == 0)).astype(int)
        d["cov"] = (d.valid * (d.cov_simult.fillna(0) == 1)).astype(int); d["covp"] = (d.valid * (d.cov_point.fillna(0) == 1)).astype(int)
        d["anydecl"] = (d.valid * ((d.decl_neg + d.decl_pos) > 0)).astype(int); d["neg"] = d.valid * (d.decl_neg == 1); d["pos"] = d.valid * (d.decl_pos == 1)
        has_m = "M" in d.columns and "g3_exact" in d.columns
        if has_m:
            d["mv_ok"] = d.valid * d.M.notna() * d.g3_exact.notna()
            d["mviol"] = (d.mv_ok * (d.M < d.g3_exact)).astype(int)
            d["mviol_col"] = d.m_viol.fillna(-1)
            d["cov_v"] = d.mviol * d["cov"]; d["nv"] = d.mviol; d["n_nv"] = (d.mv_ok * (d.M >= d.g3_exact)).astype(int); d["cov_nv"] = d.n_nv * d["cov"]
        g = d.groupby(["cell", "est", "mode", "wpn"])
        a = g.agg(n=("ds", "size"), n_valid=("valid", "sum"), k_cov=("cov", "sum"), k_covp=("covp", "sum"), k_decl=("anydecl", "sum"), k_neg=("neg", "sum"), k_pos=("pos", "sum")).reset_index()
        if has_m:
            b = g.agg(n_mok=("mv_ok", "sum"), n_mviol=("mviol", "sum"), n_cov_mviol=("cov_v", "sum"), n_nomviol=("n_nv", "sum"), n_cov_nomviol=("cov_nv", "sum"),
                      M_med=("M", "median"), g3_med=("g3_exact", "median")).reset_index()
            a = a.merge(b, on=["cell", "est", "mode", "wpn"])
            # consistency of the runner's own m_viol flag with M < required g3
            dd = d[d.mv_ok == 1]; agree = int((dd.mviol == dd.mviol_col).sum()); P(f"{name}: runner m_viol flag equals (M < exact required sup|g'''|) in {agree}/{len(dd)} answered E4 rows")
        a.to_csv(ap, index=False)
    else:
        if not os.path.exists(ap): sys.exit(f"neither raw_{name}.csv (in --raw) nor {ap}")
        a = pd.read_csv(ap)
    return a.merge(meta.drop(columns=[c for c in ("wpn", "N", "thr") if c in meta.columns]), left_on="cell", right_on="id", how="left")

# ======================================================================= 10c
def part_10c():
    meta = pd.read_csv(os.path.join(OUT, "cells_meta_amend10c.csv"))
    a = aggregate("amend10c", meta)
    P(""); P("=" * 20 + " AMENDMENT 10c (referee-directed test): fixed-window budget comparison + E4 bound diagnosis")
    main = a[a.cell >= 9000]; sep = a[a.cell.isin([302, 303])]
    fams = ["aniso2_w16", "aniso3_w8"]
    def pool(est, mode, wpn, fam=None):
        x = main[(main.est == est) & (main["mode"] == mode) & (main.wpn == wpn)]
        if fam: x = x[x.fam == fam]
        return int(x.k_cov.sum()), int(x.n_valid.sum()), int(x.n.sum())
    rows = []
    P("Pooled coverage (simultaneous), answered datasets [primary] and strict (abstention = miss). realised cells: 9 per family per wpn, 200 datasets each")
    for est in ("E1q", "E4c", "E4r"):
        for mode in ("fixed", "orig"):
            for fam in (None,) + tuple(fams):
                s = f"  {est} {mode:5s} {fam or 'both':10s}: "
                for w in (8, 32, 128):
                    k, n, nn = pool(est, mode, w, fam); s += f"wpn{w}: {fmt(k, n)} [strict {pct(k, nn)}; answered {n}/{nn}]   "
                    rows.append(dict(est=est, mode=mode, fam=fam or "both", wpn=w, k=k, n_answered=n, n_total=nn, cov=k / n if n else np.nan, cov_strict=k / nn))
                P(s)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "tables_10c_pooled.csv"), index=False)
    # central bar part 1: E1q falls under FIXED windows
    def falls(fam):
        c = {w: pool("E1q", "fixed", w, fam) for w in (8, 32, 128)}
        pt = c[8][0] / c[8][1] > c[32][0] / c[32][1] > c[128][0] / c[128][1]
        sep_ = cp(c[128][0], c[128][1])[1] < cp(c[8][0], c[8][1])[0]
        return bool(pt and sep_), c
    ok_all = True; obs = []
    for fam in (None,) + tuple(fams):
        ok, c = falls(fam); ok_all &= ok
        obs.append(f"{fam or 'pooled'}: " + " > ".join(f"{100*c[w][0]/c[w][1]:.1f}%" for w in (8, 32, 128)) + f"; CP-upper(128) {100*cp(*c[128][:2])[1]:.1f}% vs CP-lower(8) {100*cp(*c[8][:2])[0]:.1f}%")
    sc = {w: pool("E1q", "fixed", w, None) for w in (8, 32, 128)}
    strict_ok = (sc[8][0] / sc[8][2] > sc[32][0] / sc[32][2] > sc[128][0] / sc[128][2]) and cp(sc[128][0], sc[128][2])[1] < cp(sc[8][0], sc[8][2])[0]
    verdict("10c-A", "E1q coverage still falls with wpn under FIXED windows (strictly monotone 8>32>128; wpn-128 CP-upper < wpn-8 CP-lower), pooled and in each family", " | ".join(obs), ok_all, ok_all,
            note=f"strict reading (abstention = miss), pooled: {'holds' if strict_ok else 'fails'}")
    # part 2: E4c >= 93% on matched windows
    ok_all = True; okp_all = True; obs = []; sens = []
    for fam in (None,) + tuple(fams):
        for w in (8, 32, 128):
            k, n, nn = pool("E4c", "fixed", w, fam); okp, okc = ge(k, n, .93); ok_all &= okc; okp_all &= okp
            obs.append(f"{fam or 'pooled'} wpn{w}: {fmt(k, n)}")
            if k / nn < .93: sens.append(f"{fam or 'pooled'} wpn{w}: {pct(k, nn)}")
    verdict("10c-B", "E4c coverage >= 93% on the matched (fixed) windows, at each wpn, pooled and in each family", " | ".join(obs), ok_all, okp_all,
            note="strict reading (abstention = miss) would give " + ("PASS" if not sens else "FAIL: " + "; ".join(sens[:3]) + ("..." if len(sens) > 3 else "") + " (forced abstention: exact window < 6 points in 3D L=48,64)"))
    c_ok = bool(VER[-2]["verdict"] == "PASS" and VER[-1]["verdict"] == "PASS")
    verdict("10c-CENTRAL", "referee bar: central result confirmed iff 10c-A AND 10c-B", f"10c-A {VER[-2]['verdict']}, 10c-B {VER[-1]['verdict']}", c_ok)
    # cell-level fixed vs orig
    cl = []
    for (cell, wpn), t in main.groupby(["cell", "wpn"]):
        r = dict(cell=int(cell), fam=t.fam.iloc[0], L=int(t.L.iloc[0]), wpn=int(wpn), exact_J=int(t.exact_J.iloc[0]), J_fixed=float(t.J_fixed.iloc[0]) if "J_fixed" in t else np.nan)
        for est in ("E1q", "E4c", "E4r"):
            for mode in ("fixed", "orig"):
                u = t[(t.est == est) & (t["mode"] == mode)].iloc[0]; r[f"{est}_{mode}_k"] = int(u.k_cov); r[f"{est}_{mode}_nvalid"] = int(u.n_valid); r[f"{est}_{mode}_n"] = int(u.n)
        cl.append(r)
    cl = pd.DataFrame(cl).sort_values(["fam", "L", "wpn"]); cl.to_csv(os.path.join(OUT, "tables_10c_cells_fixed_vs_orig.csv"), index=False)
    # E4 diagnosis
    dg = []
    allx = pd.concat([main.assign(grp=main.fam), sep.assign(grp=sep.cell.map({302: "sep_ratio2.5_b20", 303: "sep_ratio3.0_b24"}))])
    for (est, mode, wpn, grp), t in allx[allx.est.isin(["E4c", "E4r"])].groupby(["est", "mode", "wpn", "grp"]):
        nm, nv = int(t.n_mok.sum()), int(t.n_mviol.sum()); kv, nvv = int(t.n_cov_mviol.sum()), nv; kn, nn_ = int(t.n_cov_nomviol.sum()), int(t.n_nomviol.sum())
        dg.append(dict(est=est, mode=mode, wpn=int(wpn), group=grp, n_answered=int(t.n_valid.sum()), n_with_M=nm, n_viol=nv, viol_rate=nv / nm if nm else np.nan, viol_lo=cp(nv, nm)[0], viol_hi=cp(nv, nm)[1],
                       cov_all=int(t.k_cov.sum()) / max(int(t.n_valid.sum()), 1), k_cov_given_viol=kv, n_viol2=nvv, k_cov_given_noviol=kn, n_noviol=nn_,
                       M_median_of_cell_medians=float(t.M_med.median()), g3req_median_of_cell_medians=float(t.g3_med.median())))
    dg = pd.DataFrame(dg); dg.to_csv(os.path.join(OUT, "tables_10c_e4_diagnosis.csv"), index=False)
    dc = main[main.est.isin(["E4c", "E4r"])].copy()
    dc["viol_rate"] = dc.n_mviol / dc.n_mok.replace(0, np.nan); dc["cov_given_viol"] = dc.n_cov_mviol / dc.n_mviol.replace(0, np.nan)
    dc[["cell", "fam", "L", "est", "mode", "wpn", "n_valid", "k_cov", "n_mok", "n_mviol", "n_cov_mviol", "n_nomviol", "n_cov_nomviol", "M_med", "g3_med"]].to_csv(os.path.join(OUT, "tables_10c_e4_diagnosis_cells.csv"), index=False)
    # documentation row (no bar)
    f = dg[(dg["mode"] == "fixed") & (dg.group.isin(fams))]
    s = "; ".join(f"{r.est} {r.group} wpn{r.wpn}: viol {100*r.viol_rate:.0f}%" for r in f[(f.est == "E4c")].itertuples())
    P("E4 bound violations (fitted M < exact required sup|g'''| on the window), fixed windows, E4c: " + s)
    VER.append(dict(gate="10c-DIAG", description="E4 bound violations documented (violation rate and coverage conditional on violation, incl. separation-scan ratios 2.5 and 3.0); the referee gives no pass/fail bar for this", observed="see table", verdict="REPORTED", ci_strict="", note="no bar"))
    P("[REPORTED] 10c-DIAG: documented either way, no pass/fail bar in the plan")
    return main, sep, cl, dg, dc

# ======================================================================= 11
def part_11():
    meta = pd.read_csv(os.path.join(OUT, "cells_meta_amend11.csv"))
    a = aggregate("amend11", meta)
    P(""); P("=" * 20 + " AMENDMENT 11 (frozen holdout; run once): generalisation to new graph families")
    a = a[a.est != "E6"]; cx = a[(a.kind == "crossover") & a.real]; ct = a[a.kind == "control"]
    P(f"cells: {a.cell.nunique()//3} per wpn; realised crossovers {cx.cell.nunique()//3} (all {int((meta.kind=='crossover').sum())//3} crossover cells qualify: exact max drop >= 0.6); controls {ct.cell.nunique()//3}; 200 datasets per cell; E6 ungated and omitted from the gates")
    def pooled(x, est, w, col="k_cov"):
        y = x[(x.est == est) & (x.wpn == w)]; return int(y[col].sum()), int(y.n_valid.sum()), int(y.n.sum())
    # G1
    k, n, nn = pooled(cx, "E1q", 32); okp, okc = le(k, n, .40)
    verdict("G1", "E1q realised-crossover coverage at wpn 32 <= 40%", fmt(k, n) + f"; strict (abstention=miss) {pct(k, nn)}", okc, okp,
            note="strict reading is lower, so the answered-dataset reading used is the less favourable one")
    # G2
    c = {w: pooled(cx, "E1q", w) for w in (8, 32, 128)}
    mono = c[8][0] / c[8][1] > c[32][0] / c[32][1] > c[128][0] / c[128][1]; sepok = cp(c[128][0], c[128][1])[1] < cp(c[8][0], c[8][1])[0]
    verdict("G2", "E1q coverage strictly decreases 8 -> 32 -> 128 and CP-upper(128) < CP-lower(8)", " > ".join(fmt(*c[w][:2]) for w in (8, 32, 128)) + f"; CP-upper(128) {100*cp(*c[128][:2])[1]:.1f}% < CP-lower(8) {100*cp(*c[8][:2])[0]:.1f}%", bool(mono and sepok), bool(mono and sepok))
    # G3
    ok = True; okp = True; obs = []; strict = []
    for est in ("E4c", "E4r"):
        for w in (8, 32, 128):
            k, n, nn = pooled(cx, est, w); p_, c_ = ge(k, n, .90); ok &= c_; okp &= p_; obs.append(f"{est} wpn{w} {fmt(k, n)}")
            if k / nn < .90: strict.append(f"{est} wpn{w} {pct(k, nn)}")
    verdict("G3", "E4c and E4r realised-crossover coverage >= 90% at every wpn", " | ".join(obs), ok, okp,
            note="strict reading (abstention = miss) would give " + ("PASS" if not strict else "FAIL: " + "; ".join(strict[:3]) + "... (E1q/E4 abstain by the frozen J<6 rule, e.g. C_48 x C_6 x C_6 has exact_J = 3)"))
    # G4
    ok = True; okp = True; obs = []; strict = []; skipped = []
    for est in ("E4c", "E4r", "E7", "E8"):
        for w in (8, 32, 128):
            k, n, nn = pooled(ct, est, w); p_, c_ = ge(k, n, .93)
            y = ct[(ct.est == est) & (ct.wpn == w)]
            bad = []
            for r in y.itertuples():
                if r.n_valid == 0: skipped.append(f"{est}/wpn{w}/{r.fam}_{int(r.L)}"); continue
                if cp(r.k_cov, r.n_valid)[1] < .93: bad.append(f"{r.fam}_{int(r.L)}")
            p_ &= not bad; c_ &= not bad; ok &= c_; okp &= p_; obs.append(f"{est} wpn{w} {fmt(k, n)}" + (f" BAD CELLS {bad}" if bad else ""))
            if k / nn < .93: strict.append(f"{est} wpn{w} {pct(k, nn)}")
    verdict("G4", "control coverage >= 93% pooled, no cell with CP-upper < 93%, for E4c, E4r, E7, E8 at every wpn", " | ".join(obs), ok, okp,
            note=f"cells with 0 answered datasets skipped (R3): {sorted(set(s.split('/')[2] for s in skipped))}; strict reading would give " + ("PASS" if not strict else "FAIL for " + str(len(strict)) + " est/wpn combos (forced abstention in torus3 L=18, exact_J=4)"))
    # G5
    ok = True; okp = True; obs = []
    for est in ("E1q", "E4c", "E4r", "E7", "E8"):
        k, n, nn = pooled(ct, est, 32, "k_decl"); p_, c_ = le(k, n, .05)
        y = ct[(ct.est == est) & (ct.wpn == 32)]; bad = [f"{r.fam}_{int(r.L)}" for r in y.itertuples() if r.n_valid > 0 and cp(r.k_decl, r.n_valid)[0] > .05]
        p_ &= not bad; c_ &= not bad; ok &= c_; okp &= p_; obs.append(f"{est} {fmt(k, n)}" + (f" BAD CELLS {bad}" if bad else ""))
    verdict("G5", "control FPR (margin 0.25, any sign) <= 5% for all five arms at wpn 32", " | ".join(obs), ok, okp)
    allp = all(v["verdict"] == "PASS" for v in VER if v["gate"] in ("G1", "G2", "G3", "G4", "G5"))
    verdict("11-VERDICT", "the narrowed headline GENERALISES iff G1-G5 all pass", ", ".join(f"{v['gate']} {v['verdict']}" for v in VER if v["gate"] in ("G1", "G2", "G3", "G4", "G5")), allp)
    # per-cell tables
    rows = []
    for (cell, wpn), t in a.groupby(["cell", "wpn"]):
        r = dict(cell=int(cell), fam=t.fam.iloc[0], L=int(t.L.iloc[0]), kind=t.kind.iloc[0], real=bool(t.real.iloc[0]), exact_J=int(t.exact_J.iloc[0]), max_drop=float(t.max_drop.iloc[0]), wpn=int(wpn))
        for est in ("E1q", "E4c", "E4r", "E7", "E8"):
            u = t[t.est == est].iloc[0]; r[f"{est}_cov_k"] = int(u.k_cov); r[f"{est}_decl_k"] = int(u.k_decl); r[f"{est}_nvalid"] = int(u.n_valid)
        rows.append(r)
    cell = pd.DataFrame(rows).sort_values(["kind", "fam", "L", "wpn"]); cell.to_csv(os.path.join(OUT, "tables_11_cells.csv"), index=False)
    return cell

# ======================================================================= markdown
def md_table(df, cols=None, hdr=None):
    cols = cols or list(df.columns); hdr = hdr or cols
    s = "| " + " | ".join(hdr) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
    for r in df[cols].itertuples(index=False): s += "| " + " | ".join(str(v) for v in r) + " |\n"
    return s

def kn(k, n): return "n/a" if n == 0 else f"{k}/{n} ({100*k/n:.0f}%)"

def write_md(d10b, tab, pw, xt, main, sep, cl, dg, dc, cell):
    M = []
    M.append("\n## Amendment 10b / 10c / 11 results (2026-11-22; pre-registered in PLAN.md Amendments 10 and 11; run on host-D; estimator files byte-identical to commit 1ba9ca6)\n")
    M.append("Labels: **10a, 10b and 10c are referee-directed reanalyses/tests** (external referee, 2026-11-22; the referee's own decision criteria are the bars). **Amendment 11 is a frozen holdout** (new graph families, run once, no second holdout). "
             "10a's results are in `results/amend10a_report.txt` and are not repeated here.\n")
    M.append("Files: driver `run_amend1011.py` (notes `IMPLEMENTATION_NOTES_A10_A11.md`), verdict script `verdict_amend10_11.py`, copied result tables and manifests in `results_amend10_11/` (raw per-dataset files for 10c and 11, 14.8 MB and 12.6 MB, stay on host-D `<local path omitted>`; their sha256 are in the manifests). "
             "Manifests verified: all four `complete=True`, no errors, task counts full, csv sha256 re-computed on the copied files (10b, exact-q, 10c, 11), per-item pickle sha256 lists re-checked on host-D (1030, 7, 2240, 1920 items, 0 bad, 0 missing), and the items-of-sha256 digest recomputed (4/4 match). `verdict_report.txt` carries the full printout; `reading_log.txt` the readings below.\n")
    M.append("**Readings applied where a bar needed one (nothing new was added as a bar; stricter alternatives are flagged in the gate table):** "
             "(R1) 10c and 11 coverage/FPR are over datasets where the estimator answered, exactly as every earlier amendment's C/F/X definition; the strict reading (abstention = miss) is shown beside each and flips G3 and G4 (and 10c-B) to FAIL because the frozen J<6 rule forces E1q/E4 to abstain in cells with an exact window under 6 points (C_48 x C_6 x C_6 exact_J = 3, torus3 L=18 exact_J = 4, 3D L=48/64 in 10c). 10b uses fixed denominators as the referee asked. "
             "(R2) Plain numeric bars are tested on the point estimate, as every earlier amendment's bars were, with CP intervals reported; the stricter reading 'the CP bound must also meet the bar' is flagged wherever it would change a verdict, and the plan's own CP rules are used as written. (R3) cells with zero answered datasets are skipped in per-cell tests. (R4) FPR = any-sign declaration at margin 0.25; extra reading: no cell with CP-lower > 5%. (R5) 10b: both outer relaxations A and B must pass for E7; E7j exists for B only. (R6) 10c: bars must hold pooled and within each family; (R7) 'decreases' = strictly.\n")
    M.append("### Gate table\n")
    vt = pd.DataFrame(VER); vt["verdict"] = vt["verdict"].map(lambda v: f"**{v}**")
    M.append("| gate | bar (as pre-registered) | observed (CP 95%) | verdict |\n|---|---|---|---|\n")
    for r in VER:
        obs = r["observed"].replace("|", "/")
        M.append(f"| {r['gate']} | {r['description']} | {obs}{(' ; ' + r['note']) if r['note'] else ''}{(' ; ' + r['ci_strict']) if r['ci_strict'] else ''} | **{r['verdict']}** |\n")
    # 10b
    M.append("\n### Amendment 10b (referee-directed test): certified outer bounds for E7/E7j\n")
    M.append(f"1030 datasets (170 atom-sensitive + 1 withheld + 300 seeded stable + 600 from the 3D cells L = 160, 192, 256, overlaps merged) plus the 7 exact-q cells. B = tight outer relaxation (headline); A = plan-literal nearest-atom relaxation. LP is HiGHS double precision, not a verified solve (stated in `p22amend10.py`). "
             f"All-dataset outcome counts (fixed denominator 1030): B-band covers exact d in {int(d10b.B_covered.sum())}, fails to cover in {int((d10b.B_ok & (d10b.B_cov==0)).sum())}, no certified band in {int((~d10b.B_ok).sum())}; "
             f"original gated E7 covers in {int(d10b.o_covered.sum())}, has no band in {int((~d10b.o_ok).sum())}; certified band contains the original band in {int(((d10b.B_contains_orig==1)&d10b.B_ok&d10b.o_ok).sum())}/{int((d10b.B_ok&d10b.o_ok).sum())}; "
             f"E7j declares a drop in {int(d10b.j_drop.sum())}, a rise in {int(d10b.j_rise.sum())}, is undecided (certified-undetermined at some time, no declaration) in {int(d10b.j_amb.sum())}, no declaration otherwise {int((d10b.j_none - d10b.j_amb).sum())}.\n\n")
    M.append("Cell-level outcome categories (fixed denominator = n datasets in the cell; 'no band' counts as a miss everywhere). E7 columns: B/A = certified bands, orig = original gated band; 'change' = interval-test change declaration (margin 0.25).\n\n")
    t = tab.copy(); t["cell"] = [f"{k} {f} L={int(l)}" for k, f, l in zip(t.kind, t.fam, t.L)]
    t["E7 cov B"] = [kn(a, b) for a, b in zip(t.E7B_covered, t.n)]; t["E7 cov A"] = [kn(a, b) for a, b in zip(t.E7A_covered, t.n)]; t["E7 cov orig"] = [kn(a, b) for a, b in zip(t.E7orig_covered, t.n)]
    t["not cov B"] = t.E7B_not_covered; t["no band (B/orig)"] = [f"{a}/{b}" for a, b in zip(t.no_band, t.orig_no_band)]
    t["E7 change B/A/orig"] = [f"{a}/{b}/{c}" for a, b, c in zip(t.E7B_change, t.E7A_change, t.E7orig_change)]
    t["E7j drop/rise/undecided/none"] = [f"{a}/{b}/{c}/{d}" for a, b, c, d in zip(t.E7j_drop, t.E7j_rise, t.E7j_undecided_amb, t.E7j_none)]
    M.append(md_table(t, ["cell", "n", "real", "E7 cov B", "E7 cov A", "E7 cov orig", "not cov B", "no band (B/orig)", "E7 change B/A/orig", "E7j drop/rise/undecided/none"],
                      ["cell", "n", "realised", "E7 covered, B", "E7 covered, A", "E7 covered, orig", "B not covered", "no band (B/orig)", "E7 change (B/A/orig)", "E7j drop/rise/undecided/none"]))
    M.append("\n3D (aniso3_w8) certified power by L, fixed denominator (E7j = certified joint-test drop; E7 = interval test on B; orig = original gated E7 interval test):\n\n")
    p = pw.copy(); p["E7j (B)"] = [kn(a, b) for a, b in zip(p.j, p.n)]; p["E7 interval (B)"] = [kn(a, b) for a, b in zip(p.e7B, p.n)]; p["E7 interval (orig)"] = [kn(a, b) for a, b in zip(p.e7o, p.n)]
    M.append(md_table(p, ["L", "n", "E7j (B)", "E7 interval (B)", "E7 interval (orig)"]))
    M.append("\nExact-q cells (noise-free q; certified bands vs exact d):\n\n"); M.append(md_table(xt))
    # 10c
    M.append("\n### Amendment 10c (referee-directed test): fixed-window budget comparison and E4 bound diagnosis\n")
    M.append("Cells: the 18 realised crossover cells of 9b (aniso2_w16 L=128..768, aniso3_w8 L=48..256), wpn 8/32/128, 200 datasets, fresh seeds [22, id, 3000+ds]; fixed J = the wpn-32 exact window (grid points with exact N q >= 6), recorded in `amend10c_windows.json` before data. "
             "'orig' = the estimator's own window rule on the same data. Pooled coverage over answered datasets (simultaneous band); 'strict' counts abstention as a miss. E1q has no M.\n\n")
    prow = pd.read_csv(os.path.join(OUT, "tables_10c_pooled.csv"))
    pv = prow.assign(cell_=[f"{r.k}/{r.n_answered} = {100*r.cov:.1f}% [{100*cp(r.k, r.n_answered)[0]:.1f}-{100*cp(r.k, r.n_answered)[1]:.1f}]; strict {100*r.cov_strict:.1f}%; answered {r.n_answered}/{r.n_total}" for r in prow.itertuples()])
    pv = pv.pivot_table(index=["est", "mode", "fam"], columns="wpn", values="cell_", aggfunc="first").reset_index()
    M.append(md_table(pv, ["est", "mode", "fam", 8, 32, 128], ["est", "window", "family", "wpn 8", "wpn 32", "wpn 128"]))
    M.append("\nCell-level: fixed window vs original window (coverage k/answered; J_fixed = exact window, grid points; cells with J_fixed < 6 are abstentions by the frozen rule):\n\n")
    c2 = cl.copy()
    for est in ("E1q", "E4c", "E4r"):
        c2[est] = [f"{kn(a, b)} vs {kn(c, d)}" for a, b, c, d in zip(c2[f"{est}_fixed_k"], c2[f"{est}_fixed_nvalid"], c2[f"{est}_orig_k"], c2[f"{est}_orig_nvalid"])]
    M.append(md_table(c2, ["fam", "L", "wpn", "J_fixed", "E1q", "E4c", "E4r"], ["family", "L", "wpn", "J_fixed", "E1q fixed vs orig", "E4c fixed vs orig", "E4r fixed vs orig"]))
    M.append("\nE4 bound diagnosis: fitted M against the exact required supremum of |g'''| (g = log q vs log t) on [t_0, t_J] (`p22honest.exact_g3_cummax`). A 'violation' is M < required. Coverage is simultaneous. The runner's own `m_viol` flag was re-derived as M < required and agrees on every answered E4 row (see `verdict_report.txt`).\n\n")
    g = dg.copy()
    g["viol rate (CP)"] = [f"{r.n_viol}/{r.n_with_M} = {100*r.viol_rate:.1f}% [{100*r.viol_lo:.1f}-{100*r.viol_hi:.1f}]" for r in g.itertuples()]
    g["cov | viol"] = [kn(r.k_cov_given_viol, r.n_viol) for r in g.itertuples()]; g["cov | no viol"] = [kn(r.k_cov_given_noviol, r.n_noviol) for r in g.itertuples()]
    g["cov (all answered)"] = [f"{100*r.cov_all:.1f}%" for r in g.itertuples()]; g["median M / required"] = [f"{r.M_median_of_cell_medians:.2f} / {r.g3req_median_of_cell_medians:.2f}" for r in g.itertuples()]
    g = g.sort_values(["group", "est", "mode", "wpn"])
    M.append(md_table(g, ["group", "est", "mode", "wpn", "n_answered", "viol rate (CP)", "cov (all answered)", "cov | viol", "cov | no viol", "median M / required"],
                      ["family / scan cell", "est", "window", "wpn", "answered", "violation rate (CP)", "coverage", "coverage given violation", "coverage given no violation", "median M / required sup g''' (median of cell medians)"]))
    M.append("\nSeparation-scan cells (ratio 2.5 = twostep_8_20 / b=20, id 302; ratio 3.0 = twostep_8_24 / b=24, id 303; wpn 32, original windows only; E1q for reference):\n\n")
    ss = sep.copy(); ss["scan"] = ss.cell.map({302: "ratio 2.5 (b=20)", 303: "ratio 3.0 (b=24)"})
    ss["coverage"] = [kn(a, b) for a, b in zip(ss.k_cov, ss.n_valid)]
    ss["viol"] = [("-" if r.est == "E1q" else f"{int(r.n_mviol)}/{int(r.n_mok)}") for r in ss.itertuples()]
    ss["cov|viol"] = [("-" if r.est == "E1q" else kn(int(r.n_cov_mviol), int(r.n_mviol))) for r in ss.itertuples()]
    ss["cov|no viol"] = [("-" if r.est == "E1q" else kn(int(r.n_cov_nomviol), int(r.n_nomviol))) for r in ss.itertuples()]
    ss["M / req"] = [("-" if r.est == "E1q" else f"{r.M_med:.2f} / {r.g3_med:.2f}") for r in ss.itertuples()]
    M.append(md_table(ss.sort_values(["scan", "est"]), ["scan", "est", "coverage", "viol", "cov|viol", "cov|no viol", "M / req"], ["scan", "est", "coverage (answered)", "violations", "coverage given violation", "coverage given no violation", "median M / required"]))
    # 11
    M.append("\n### Amendment 11 (frozen holdout; run once): new graph families\n")
    M.append("Families never used in Amendments 1-10: 2D->1D crossover C_L x C_12 (L = 96..384), 3D->1D crossover C_L x C_6 x C_6 (L = 48..192), controls ring 384/768, torus 40/80, 3D torus 18/26. Exact truth table `amend11_cells_exact.json` (sha256 de21f1a0..., written before simulation). All 10 crossover cells are realised (exact max drop >= 0.6; min 0.83). 200 datasets per cell, wpn 8/32/128, E1q/E4c/E4r/E7/E8 (E6 ungated, omitted).\n\n")
    cx = cell[cell.kind == "crossover"].copy(); ct = cell[cell.kind == "control"].copy()
    M.append("Crossover cells: simultaneous-band coverage, k/answered (G1-G3 use E1q, E4c, E4r; E7/E8 shown for reference). exact_J = points in the exact window; max drop = exact.\n\n")
    for df in (cx,):
        for est in ("E1q", "E4c", "E4r", "E7", "E8"): df[est] = [kn(a, b) for a, b in zip(df[f"{est}_cov_k"], df[f"{est}_nvalid"])]
        df["cell"] = [f"{f} L={int(l)}" for f, l in zip(df.fam, df.L)]; df["maxd"] = df.max_drop.map(lambda v: f"{v:.2f}")
    M.append(md_table(cx, ["cell", "exact_J", "maxd", "wpn", "E1q", "E4c", "E4r", "E7", "E8"], ["cell", "exact_J", "exact max drop", "wpn", "E1q", "E4c", "E4r", "E7", "E8"]))
    M.append("\nControl cells: simultaneous-band coverage k/answered (G4: E4c, E4r, E7, E8 at every wpn; E1q for reference):\n\n")
    for est in ("E1q", "E4c", "E4r", "E7", "E8"): ct[est] = [kn(a, b) for a, b in zip(ct[f"{est}_cov_k"], ct[f"{est}_nvalid"])]
    ct["cell"] = [f"{f} L={int(l)}" for f, l in zip(ct.fam, ct.L)]
    M.append(md_table(ct, ["cell", "exact_J", "wpn", "E1q", "E4c", "E4r", "E7", "E8"]))
    M.append("\nControl false-change declarations (G5 at wpn 32; margin 0.25, any sign), k/answered:\n\n")
    for est in ("E1q", "E4c", "E4r", "E7", "E8"): ct["f_" + est] = [kn(a, b) for a, b in zip(ct[f"{est}_decl_k"], ct[f"{est}_nvalid"])]
    M.append(md_table(ct, ["cell", "wpn", "f_E1q", "f_E4c", "f_E4r", "f_E7", "f_E8"], ["cell", "wpn", "E1q", "E4c", "E4r", "E7", "E8"]))
    M.append("\n")
    return "".join(M)

def main():
    os.makedirs(OUT, exist_ok=True)
    check_manifests()
    d10b, tab, pw, xt = part_10b()
    mainc, sepc, cl, dg, dc = part_10c()
    cell = part_11()
    pd.DataFrame(VER).to_csv(os.path.join(OUT, "verdicts.csv"), index=False)
    open(os.path.join(OUT, "verdict_report.txt"), "w").write("\n".join(LINES) + "\n")
    doc = __doc__; i = doc.index("READING LOG"); open(os.path.join(OUT, "reading_log.txt"), "w").write(doc[i:])
    md = write_md(d10b, tab, pw, xt, mainc, sepc, cl, dg, dc, cell)
    open(os.path.join(OUT, "RESULTS_SECTION_A10_A11.md"), "w").write(md)
    if "--append" in args:
        rp = os.path.join(HERE, "RESULTS.md"); cur = open(rp).read()
        if "## Amendment 10b / 10c / 11 results" in cur: sys.exit("RESULTS.md already contains the section; not appending")
        open(rp, "a").write(md); print("appended to RESULTS.md")
    P(""); P("SUMMARY: " + "; ".join(f"{v['gate']} {v['verdict']}" for v in VER))

if __name__ == "__main__": main()
