#!/usr/bin/env python3
"""P22 Amendment 6 analysis. Writes results/amend6_*.csv and amend6_tables.md. Usage: analyze_amend6.py"""
import os, sys
import numpy as np, pandas as pd
from p22lib import cp
import p22recipe as RC
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
out = []
def P(s=""): out.append(str(s))
def md(df):
    df = df.fillna("")
    L_ = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    for _, r in df.iterrows(): L_.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(L_)
def ci(k, n): lo, hi = cp(int(k), int(n)); return "%.1f (%.1f-%.1f)" % (100 * k / n, 100 * lo, 100 * hi) if n else "n/a"
def rd(name): return pd.read_csv(os.path.join(RES, name))

# ------------------------------------------------------------------ meta / groups (as analyze_amend5)
m_old = pd.read_csv(RES + "/cells_meta.csv", keep_default_na=False); m_old["note"] = m_old["note"].astype(str)
m_old = m_old[(m_old.note == "") & m_old.kind.isin(["control", "crossover"])][["id", "fam", "L", "N", "kind", "max_drop", "max_rise", "exact_J"]]
m_new = rd("cells_meta_amend5.csv")[["id", "fam", "L", "N", "kind", "max_drop", "max_rise", "exact_J"]]
meta = pd.concat([m_old, m_new[m_new.id >= 100]], ignore_index=True)
for c_ in ("max_drop", "max_rise", "exact_J"): meta[c_] = pd.to_numeric(meta[c_], errors="coerce")
def grp(r):
    if r.kind == "control": return "const"
    if r.kind == "twostep": return "twostep"
    if r.kind == "brush": return "brush"
    if r.kind == "gasket": return "gasket"
    return "xover_real" if r.max_drop >= 0.6 else "xover_notreal"
meta["group"] = meta.apply(grp, axis=1); meta["name"] = meta.fam + "_" + meta.L.astype(str)
def prep(df, mt=meta):
    d = df.merge(mt, left_on="cell", right_on="id"); d["decl_any"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
    d["correct"] = np.where(d.max_rise >= 0.6, d.decl_pos, d.decl_neg).astype(float); return d

done = {g: os.path.exists(os.path.join(RES, f"raw_amend6_{g}.csv")) for g in ("big2d", "sep", "fresh", "orig")}
P(f"groups present: {done}")

# ------------------------------------------------------------------ ITEM 1
if done["big2d"]:
    d = RC.load(); big = rd("raw_amend6_big2d.csv"); mb = rd("cells_meta_amend6_big2d.csv")
    big = big.merge(mb[["id", "fam", "L", "N", "max_drop", "exact_J"]], left_on="cell", right_on="id")
    d2 = pd.concat([d[~d.est.isin(["E4c", "E4r"]) | (d.fam != "aniso2_w16") | True], big[big.est.isin(["E4c", "E4r"])]], ignore_index=True)
    g = RC.cell_table(d2)
    cells = RC.exact_cells()
    pred = rd("amend6_recipe_predictions.csv").set_index("est")
    rows = []
    for est in RC.FIRST5 + RC.E4:
        grid = RC.G2 if est in RC.FIRST5 else RC.G2X
        t2 = g[(g.est == est) & (g.fam == "aniso2_w16")]; meas = {int(r.L): r.power for _, r in t2.iterrows()}
        mi = RC.minL(meas, grid, 0.8); pi = int(pred.loc[est, "pred_idx"])
        lab = lambda i: grid[i] if i < len(grid) else f">{grid[-1]}"
        rows.append(dict(est=est, predicted_minL=lab(pi), measured_minL=lab(mi), index_distance=abs(pi - mi), within_one_step=abs(pi - mi) <= 1,
                         power_by_L=" ".join(f"{l}:{meas[l]:.2f}" for l in grid if l in meas)))
    R1 = pd.DataFrame(rows); R1.to_csv(RES + "/amend6_recipe_scores.csv", index=False)
    P("## Item 1: minimum-size recipe (fit on 3D only; power = declaration rate, abstention counted as a miss)"); P(md(R1)); P()
    P(f"First five (consistency): {int(R1[R1.est.isin(RC.FIRST5)].within_one_step.sum())} of 5 within one grid step. E4c/E4r (out of sample): {int(R1[R1.est.isin(RC.E4)].within_one_step.sum())} of 2.")
    pr = pd.read_csv(RES + "/amend6_recipe_params.csv"); P(); P(md(pr.round(3)))
    # predicted vs measured hw for the new cells
    hwrows = []
    for est in RC.E4:
        p = pr[pr.est == est].iloc[0]
        for l in (192, 256, 384, 512, 768):
            r = g[(g.est == est) & (g.fam == "aniso2_w16") & (g.L == l)]
            if len(r): hwrows.append(dict(est=est, L=l, J=int(r.J.iloc[0]), hw_pred=round(p.a * (r.J.iloc[0] / 7) ** (-p.delta), 3), hw_meas=round(r.hw.iloc[0], 3), nonabs=round(r.nonabs.iloc[0], 3), power=round(r.power.iloc[0], 3)))
    P(); P("E4 half-width and power in the large 2D cells (hw_pred is the recipe's J-law, ungated):"); P(md(pd.DataFrame(hwrows)))
    cov = big[big.abstain == 0].groupby(["est", "L"]).agg(cov=("cov_simult", "mean"), n=("cov_simult", "size")).round(3).reset_index()
    P(); P("Coverage of exact d(t) in the large 2D cells (non-abstained):"); P(md(cov)); P()

# ------------------------------------------------------------------ ITEM 2
if done["orig"]:
    o = prep(rd("raw_amend6_orig.csv")); e6 = o[o.est == "E6"].copy(); e6na = e6[e6.abstain == 0].copy(); e6na["miss"] = 1 - e6na.cov_simult
    e6na["a1"] = e6na.a1_alarm.fillna(0); e6na["a2"] = e6na.a2_alarm; e6na["lof"] = e6na.lof_alarm
    e4 = o[o.est == "E4c"]
    P("## Item 2: E6 misfit alarms (original seeds, 37 cells x 200)")
    P(f"E4c abstained where E6 did not: {int(e6na.get('a1_e4_abstain', pd.Series(dtype=float)).fillna(0).sum())} datasets (no alarm).")
    # reproduction vs raw_amend5
    r5 = rd("raw_amend5.csv"); r5 = r5[r5.est == "E6"]; cmp = e6.merge(r5, on=["cell", "ds"], suffixes=("", "_o"))
    P(f"Reproduction of E6 vs raw_amend5: {len(cmp)} datasets, abstain identical {(cmp.abstain == cmp.abstain_o).mean():.4f}, cov_simult identical {(cmp.cov_simult.fillna(-1) == cmp.cov_simult_o.fillna(-1)).mean():.4f}, max |hw diff| {np.nanmax(np.abs(cmp.hw - cmp.hw_o)):.2e}")
    res = []
    for alarm, lab in (("a1", "A1 E6-vs-E4c band disagreement"), ("a2", "A2 one-vs-two-step (fallback)"), ("lof", "LoF (Amendment 5, for reference)")):
        mis = e6na[e6na.group.isin(["twostep", "brush", "gasket"])]; k = int(mis.miss.sum())
        ts320 = e6na[(e6na.fam == "twostep_8_32") & (e6na.L == 320)]; k320 = int(ts320.miss.sum())
        ct = e6na[e6na.group == "const"]
        xr = e6[e6.group == "xover_real"]
        pw = xr.decl_neg.fillna(0).mean(); pwg = ((xr.decl_neg.fillna(0) == 1) & (xr[alarm if alarm != "lof" else "lof_alarm"].fillna(0) == 0) if False else None)
        col = {"a1": "a1_alarm", "a2": "a2_alarm", "lof": "lof_alarm"}[alarm]
        pwg = ((xr.decl_neg.fillna(0) == 1) & (xr[col].fillna(0) == 0)).mean()
        res.append(dict(alarm=lab,
                        a_pooled_P_alarm_given_miss=ci(mis[mis.miss == 1][alarm].sum(), k) + f"  [misses {k}]", a_pass=(mis[mis.miss == 1][alarm].mean() >= 0.5) if k else False,
                        b_twostep320_P_alarm_given_miss=ci(ts320[ts320.miss == 1][alarm].sum(), k320) + f"  [misses {k320}]" if k320 else "no misses", b_pass=(ts320[ts320.miss == 1][alarm].mean() >= 0.8) if k320 else False,
                        c_control_alarm=ci(ct[alarm].sum(), len(ct)), c_pass=ct[alarm].mean() <= 0.05,
                        d_xover_power_E6=round(pw, 3), d_power_E6g=round(pwg, 3), d_loss_points=round(100 * (pw - pwg), 1), d_pass=(pw - pwg) <= 0.05))
        res[-1]["ALL_PASS"] = bool(res[-1]["a_pass"] and res[-1]["b_pass"] and res[-1]["c_pass"] and res[-1]["d_pass"])
    R2 = pd.DataFrame(res); R2.to_csv(RES + "/amend6_alarm_bars.csv", index=False); P(md(R2.T.reset_index().rename(columns={"index": "bar"}).rename(columns={0: "A1", 1: "A2", 2: "LoF"}))); P()
    rows = []
    for g_, x in e6na.groupby("name"):
        mm = meta[meta.name == g_].iloc[0]
        rows.append(dict(cell=g_, group=mm.group, n=len(x), cov=round(x.cov_simult.mean(), 3), misses=int(x.miss.sum()), a1=round(x.a1.mean(), 3), a2=round(x.a2.mean(), 3), lof=round(x.lof.mean(), 3),
                         a1_given_miss=round(x[x.miss == 1].a1.mean(), 3) if x.miss.sum() else np.nan, a2_given_miss=round(x[x.miss == 1].a2.mean(), 3) if x.miss.sum() else np.nan))
    C2 = pd.DataFrame(rows); C2.to_csv(RES + "/amend6_alarm_cells.csv", index=False); P("Per cell (non-abstained E6; a1/a2/lof = alarm rates):"); P(md(C2)); P()
    # E4c coverage reproduced
    P("E4c (re-run) simultaneous coverage by group: " + "; ".join(f"{g_}: {ci(x.cov_simult.sum(), x.cov_simult.notna().sum())}" for g_, x in e4[e4.abstain == 0].groupby("group")))
    P()

# ------------------------------------------------------------------ ITEM 3
def origset(est):
    """original numbers for CP intervals: E6 from raw_amend5; E4c from raw_amend3 (old cells) + raw_amend5 (new cells)."""
    r5 = rd("raw_amend5.csv"); r3 = rd("raw_amend3.csv")
    if est == "E6": return prep(r5[r5.est == "E6"])
    return prep(pd.concat([r5[(r5.est == "E4c") & (r5.cell >= 100)], r3[(r3.est == "E4c") & (r3.cell < 100)]], ignore_index=True))
if done["fresh"]:
    f = prep(rd("raw_amend6_fresh.csv"))
    P("## Item 3: fresh-seed replication (seeds [22, cell, 1000+ds])")
    rows = []; ninside = {}
    for est in ("E6", "E4c"):
        orig = origset(est); fr = f[f.est == est]
        def chk(name, sel_o, sel_f, kind):
            xo = orig[sel_o(orig)]; xf = fr[sel_f(fr)]
            if kind == "cov":
                xo = xo[(xo.abstain == 0) & xo.cov_simult.notna()]; xf = xf[(xf.abstain == 0) & xf.cov_simult.notna()]; ko, no, kf, nf = xo.cov_simult.sum(), len(xo), xf.cov_simult.sum(), len(xf)
            else:
                ko, no, kf, nf = xo.decl_neg.fillna(0).sum(), len(xo), xf.decl_neg.fillna(0).sum(), len(xf)
            lo, hi = cp(int(ko), int(no)); v = kf / nf
            rows.append(dict(est=est, check=name, original=ci(ko, no), replicate=ci(kf, nf), inside=bool(lo <= v <= hi)))
        chk("coverage, realised crossover pooled", lambda x: x.group == "xover_real", lambda x: x.group == "xover_real", "cov")
        for L_ in (128, 192, 256, 320):
            chk(f"coverage, twostep_8_32 L={L_}", lambda x, L_=L_: (x.fam == "twostep_8_32") & (x.L == L_), lambda x, L_=L_: (x.fam == "twostep_8_32") & (x.L == L_), "cov")
        chk("coverage, twostep pooled", lambda x: x.group == "twostep", lambda x: x.group == "twostep", "cov")
        for fam in ("brush_f24", "brush_f16"):
            for L_ in sorted(meta[meta.fam == fam].L.unique()):
                chk(f"coverage, {fam} Lb={L_}", lambda x, fam=fam, L_=L_: (x.fam == fam) & (x.L == L_), lambda x, fam=fam, L_=L_: (x.fam == fam) & (x.L == L_), "cov")
        chk("coverage, brush pooled", lambda x: x.group == "brush", lambda x: x.group == "brush", "cov")
        chk("correct-sign power, realised crossover (abstain=miss)", lambda x: x.group == "xover_real", lambda x: x.group == "xover_real", "pow")
    R3 = pd.DataFrame(rows); R3.to_csv(RES + "/amend6_replication.csv", index=False); P(md(R3)); P()
    e6r = R3[R3.est == "E6"]; P(f"E6 gated checks inside the original CP interval: {int(e6r.inside.sum())} of {len(e6r)}.  E4c (ungated): {int(R3[R3.est=='E4c'].inside.sum())} of {len(R3[R3.est=='E4c'])}."); P()
    # fresh alarm replication (descriptive)
    fe = f[(f.est == "E6") & (f.abstain == 0)].copy(); fe["miss"] = 1 - fe.cov_simult
    mis = fe[fe.group.isin(["twostep", "brush"])]; k = int(mis.miss.sum())
    P(f"Fresh seeds, descriptive alarm check (twostep+brush only, no gasket/controls): misses {k}; P(A1|miss) = " + (ci(mis[mis.miss == 1].a1_alarm.fillna(0).sum(), k) if k else "n/a") + "; P(A2|miss) = " + (ci(mis[mis.miss == 1].a2_alarm.sum(), k) if k else "n/a"))
    ts = fe[(fe.fam == "twostep_8_32") & (fe.L == 320)]; k2 = int(ts.miss.sum())
    if k2: P(f"  twostep_8_32 L=320: misses {k2}, P(A1|miss) = {ci(ts[ts.miss==1].a1_alarm.fillna(0).sum(), k2)}, P(A2|miss) = {ci(ts[ts.miss==1].a2_alarm.sum(), k2)}")
    P()
if done["sep"]:
    sm = rd("cells_meta_amend6_sep.csv"); s = rd("raw_amend6_sep.csv").merge(sm[["id", "fam", "L", "max_drop", "exact_J"]], left_on="cell", right_on="id")
    s["b"] = s.fam.str.split("_").str[-1].astype(int)
    rows = []
    for b, x in s.groupby("b"):
        for est in ("E6", "E4c"):
            y = x[(x.est == est)]; yn = y[y.abstain == 0]
            rows.append(dict(b=b, ratio=b / 8, est=est, n=len(yn), abstain=round(1 - len(yn) / len(y), 3), cov=ci(yn.cov_simult.sum(), len(yn)), cov_lo=cp(int(yn.cov_simult.sum()), len(yn))[0], cov_hi=cp(int(yn.cov_simult.sum()), len(yn))[1],
                             alarm_A1=round(yn.a1_alarm.fillna(0).mean(), 3) if est == "E6" else np.nan, alarm_A2=round(yn.a2_alarm.mean(), 3) if est == "E6" else np.nan, hw=round(yn.hw.mean(), 3)))
    # b = 32 at L = 256: original + fresh pooled
    if done["fresh"] and done["orig"]:
        for est in ("E6", "E4c"):
            parts = [prep(rd("raw_amend6_orig.csv")), prep(rd("raw_amend6_fresh.csv"))]
            y = pd.concat([p_[(p_.est == est) & (p_.fam == "twostep_8_32") & (p_.L == 256)] for p_ in parts]); yn = y[y.abstain == 0]
            rows.append(dict(b=32, ratio=4.0, est=est, n=len(yn), abstain=round(1 - len(yn) / len(y), 3), cov=ci(yn.cov_simult.sum(), len(yn)), cov_lo=cp(int(yn.cov_simult.sum()), len(yn))[0], cov_hi=cp(int(yn.cov_simult.sum()), len(yn))[1],
                             alarm_A1=round(yn.a1_alarm.fillna(0).mean(), 3) if est == "E6" else np.nan, alarm_A2=round(yn.a2_alarm.mean(), 3) if est == "E6" else np.nan, hw=round(yn.hw.mean(), 3)))
    S3 = pd.DataFrame(rows).sort_values(["est", "b"]); S3.to_csv(RES + "/amend6_separation_scan.csv", index=False)
    P("## Step-separation scan (twostep C_256 x C_8 x C_b; E6 simultaneous coverage of exact d(t), non-abstained; b=32 is original+fresh pooled)"); P(md(S3.drop(columns=["cov_lo", "cov_hi"]))); P()
    e = S3[S3.est == "E6"].sort_values("b"); okb = e[e.cov_lo >= 0.93].b.tolist(); badb = e[e.cov_hi < 0.93].b.tolist()
    first_pt = e[e.cov_lo.notna()].assign(pt=lambda x: x["cov"].str.split(" ").str[0].astype(float)); belowpt = first_pt[first_pt.pt < 93].b.tolist()
    P(f"E6: b with CP-lower >= 93%: {okb}; b with point estimate < 93%: {belowpt}; b with CP-upper < 93%: {badb}")
open(os.path.join(RES, "amend6_tables.md"), "w").write("\n".join(out)); print("\n".join(out))
