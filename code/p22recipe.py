#!/usr/bin/env python3
"""P22 Amendment 6, Item 1: minimum-size recipe. Fit on 3D cells only; predict 2D min L. Usage: p22recipe.py predict | p22recipe.py score"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
G2 = [32, 48, 64, 96, 128, 192, 256]; G2X = G2 + [384, 512, 768]
FIRST5 = ["E0", "E1q", "E6", "Mb", "CL"]; E4 = ["E4c", "E4r"]
MEAS_PREV = dict(E0=192, E1q=192, E6=256, Mb=256, CL=128)   # as recorded in RESULTS.md Amendment 5 (check only)

def load():
    m = pd.concat([pd.read_csv(RES + "/cells_meta.csv", keep_default_na=False).query("note=='' and kind in ['control','crossover']")[["id", "fam", "L", "N", "max_drop", "exact_J"]],
                   pd.read_csv(RES + "/cells_meta_amend5.csv").query("id>=100")[["id", "fam", "L", "N", "max_drop", "exact_J"]]])
    for c in ("max_drop", "exact_J"): m[c] = pd.to_numeric(m[c])
    r5 = pd.read_csv(RES + "/raw_amend5.csv"); r0 = pd.read_csv(RES + "/raw.csv"); r0 = r0[r0.est == "E0"]
    r1 = pd.read_csv(RES + "/raw_amend.csv"); r1 = r1[r1.est == "E1q"]; r3 = pd.read_csv(RES + "/raw_amend3.csv"); r3 = r3[r3.est.isin(E4)]
    d = pd.concat([r5, r0[r0.cell < 100], r1[r1.cell < 100], r3[r3.cell < 100]], ignore_index=True)
    return d.merge(m, left_on="cell", right_on="id")

def cell_table(d):
    d = d.copy(); d["na"] = 1 - d.abstain; d["det"] = ((d.abstain == 0) & (d.decl_neg == 1)).astype(float)
    g = d.groupby(["est", "fam", "L"]).agg(nonabs=("na", "mean"), hw=("hw", lambda s: s.dropna().mean()), power=("det", "mean"), J=("exact_J", "first"), N=("N", "first"), S=("max_drop", "first")).reset_index()
    return g

def minL(power_by_L, grid, thr=0.8):
    """smallest L such that it and all larger grid values have power >= thr; index len(grid) = not reached."""
    pw = [power_by_L.get(l, np.nan) for l in grid]
    for i in range(len(grid)):
        if all((not np.isnan(p)) and p >= thr for p in pw[i:]): return i
    return len(grid)

def fit(g, est):
    t3 = g[(g.est == est) & g.fam.isin(["aniso3_w8", "torus3"])]
    low = t3[t3.nonabs < 0.8]
    Jmin = int(low.J.max()) + 1 if len(low) else int(t3.J.min())
    a3 = t3[(t3.fam == "aniso3_w8") & (t3.J >= Jmin) & (t3.nonabs >= 0.8)].sort_values("J")
    out = dict(est=est, Jmin=Jmin, n_fit=len(a3))
    if len(a3) < 2: out["status"] = "recipe not fitted"; return out
    X = np.log(a3.J.values / 7.0); y = np.log(a3.hw.values)
    A = np.vstack([np.ones(len(X)), -X]).T; (la, dl), *_ = np.linalg.lstsq(A, y, rcond=None)
    r = (a3.S.values - 0.25) / a3.hw.values; det = a3.power.values >= 0.8
    Tlo = max(1.0, r[~det].max()) if (~det).any() else 1.0
    Thi = r[det].min() if det.any() else np.inf
    if not np.isfinite(Thi): T, flag = Tlo, "no detected cell"
    elif Tlo < Thi: T, flag = float(np.sqrt(Tlo * Thi)), ""
    else: T, flag = float((Tlo + Thi) / 2), "inconsistent bracket"
    out.update(status="ok", a=float(np.exp(la)), delta=float(dl), Tlo=float(Tlo), Thi=float(Thi), T=T, flag=flag, Jfit=",".join(map(str, a3.J.astype(int))))
    return out

def predict_cells(par, cells):
    res = {}
    for L_, (J, S) in cells.items():
        if par["status"] != "ok": res[L_] = (np.nan, np.nan, False); continue
        hw = par["a"] * (J / 7.0) ** (-par["delta"])
        res[L_] = (hw, (S - 0.25) / hw, bool(J >= par["Jmin"] and (S - 0.25) / hw >= par["T"]))
    return res

def exact_cells():
    """exact J and S for the 2D family (existing from cells_meta; new cells from exact curves)."""
    m = pd.read_csv(RES + "/cells_meta.csv", keep_default_na=False)
    a = m[(m.fam == "aniso2_w16") & (m.note == "")]
    cells = {int(r.L): (int(float(r.exact_J)), float(r.max_drop)) for _, r in a.iterrows()}
    cells.update({384: (14, 1.004), 512: (16, 1.004), 768: (18, 1.004)})   # filled/verified by exact_new() below
    return cells

def exact_new():
    sys.path.insert(0, HERE)
    import p22lib as L, run_p22 as R
    rr = np.random.default_rng(2200); out = {}
    for l in (384, 512, 768):
        s = dict(fam="aniso2_w16", L=l, build=lambda l=l: L.Graph("torus", 0.5, dims=(l, 16)), sizes=None, plateau=None, kind="crossover", grid=(24, 34), lazy=.5, note="", id=900)
        R.prep(s, rr); out[l] = (int(s["exact_J"]), float(s["max_drop"]))
    return out

def main():
    mode = sys.argv[1]
    d = load(); g = cell_table(d)
    cells = exact_cells(); cells.update(exact_new())
    rows = []; par_rows = []
    for est in FIRST5 + E4:
        p = fit(g, est); par_rows.append(p)
        grid = G2 if est in FIRST5 else G2X
        pr = predict_cells(p, {l: cells[l] for l in grid})
        pred_idx = minL({l: float(v[2]) for l, v in pr.items()}, grid, 0.5) if p["status"] == "ok" else np.nan
        # measured on 2D family
        t2 = g[(g.est == est) & (g.fam == "aniso2_w16")]
        meas = {int(r.L): r.power for _, r in t2.iterrows()}
        meas_idx = minL(meas, [l for l in grid if l in meas] if mode == "predict" else grid, 0.8) if len(meas) else np.nan
        # 3D in-sample
        t3 = g[(g.est == est) & (g.fam == "aniso3_w8")]; m3 = {int(r.L): r.power for _, r in t3.iterrows()}
        g3 = [32, 48, 64, 96, 128, 192]
        c3 = {int(r.L): (int(r.J), float(r.S)) for _, r in t3.iterrows()}
        p3 = predict_cells(p, c3) if p["status"] == "ok" else {}
        i3p = minL({l: float(v[2]) for l, v in p3.items()}, g3, 0.5) if p3 else np.nan
        i3m = minL(m3, g3, 0.8)
        lab = lambda i, gr: (gr[int(i)] if int(i) < len(gr) else f">{gr[-1]}") if not np.isnan(i) else "n/a"
        rows.append(dict(est=est, pred_minL_2D=lab(pred_idx, grid), meas_minL_2D=lab(meas_idx, grid) if mode == "score" else "(not looked at yet)" if est in E4 else lab(meas_idx, grid),
                         pred_idx=pred_idx, meas_idx=meas_idx, in_sample_3D_pred=lab(i3p, g3), in_sample_3D_meas=lab(i3m, g3),
                         pred_detail=";".join(f"{l}:hw{v[0]:.2f},r{v[1]:.2f},{'Y' if v[2] else 'n'}" for l, v in pr.items())))
    P = pd.DataFrame(par_rows); R_ = pd.DataFrame(rows)
    if mode == "predict":
        # E4 predictions only: measured side withheld for E4 grid beyond existing data
        P.to_csv(RES + "/amend6_recipe_params.csv", index=False)
        R_["meas_idx"] = np.where(R_.est.isin(E4), np.nan, R_.meas_idx)
        R_.to_csv(RES + "/amend6_recipe_predictions.csv", index=False)
    print(P.round(3).to_string()); print(R_.drop(columns=["pred_detail"]).to_string())
    for _, r in R_.iterrows(): print(r.est, r.pred_detail)

if __name__ == "__main__": main()
