#!/usr/bin/env python3
"""P22 Amendment 7 analysis. Usage: analyze_amend7.py [--smoke]. Writes results/amend7_{calib,minL,scores}.csv and prints verdicts."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
SFX = "_smoke" if "--smoke" in sys.argv else ""
G = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768]}
ESTS = ["E0", "E1q", "E4c", "E4r", "E6"]; GATED = ["E0", "E1q", "E4c", "E6"]

def cp(k, n):
    lo = beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(.975, k + 1, n - k) if k < n else 1.0
    return lo, hi

def load_new():
    d = pd.read_csv(os.path.join(RES, f"raw_amend7{SFX}.csv")); m = pd.read_csv(os.path.join(RES, f"cells_meta_amend7{SFX}.csv"))
    d = d.merge(m[["id", "fam", "L", "kind", "max_drop"]], left_on="cell", right_on="id")
    return d[d.est.isin(ESTS)]   # E0 rows = change test (E0g excluded); abstained rows kept as misses

def load_old32():
    """wpn = 32 crossover rows on the old grid from earlier raw files (same seeds/code as published)."""
    m = pd.read_csv(RES + "/cells_meta.csv", keep_default_na=False); m = m[(m.note == "") & (m.kind == "crossover")][["id", "fam", "L", "max_drop"]]
    parts = []
    for f, e in (("raw.csv", "E0"), ("raw_amend.csv", "E1q"), ("raw_amend3.csv", "E4c"), ("raw_amend3.csv", "E4r"), ("raw_amend5.csv", "E6")):
        r = pd.read_csv(RES + "/" + f); r = r[(r.est == e) & (r.cell < 100)]
        parts.append(r)
    b = pd.read_csv(RES + "/raw_amend6_big2d.csv"); bm = pd.DataFrame(dict(id=[400, 401, 402], fam="aniso2_w16", L=[384, 512, 768], max_drop=1.004))
    o = pd.concat(parts).merge(m, left_on="cell", right_on="id")
    o = pd.concat([o, b.merge(bm, left_on="cell", right_on="id")]); o["wpn"] = 32; o["kind"] = "crossover"
    return o

def minL(pw, grid):
    for i in range(len(grid)):
        if all((l in pw) and pw[l] >= 0.8 for l in grid[i:]): return i
    return len(grid)

def main():
    d = load_new()
    # (A) calibration on controls
    cal = []
    for (w, e), x in d[d.kind == "control"].groupby(["wpn", "est"]):
        ok = x[x.abstain == 0]; k, n = int(ok.cov_simult.sum()), len(ok)
        per = ok.groupby("cell").cov_simult.agg(["sum", "count"]); bad = [c for c, r in per.iterrows() if cp(int(r["sum"]), int(r["count"]))[1] < .93]
        f = ok[ok.decl_neg.notna()]; dk = int(((f.decl_neg.fillna(0) + f.decl_pos.fillna(0)) > 0).sum())
        pf = f.assign(a=((f.decl_neg.fillna(0) + f.decl_pos.fillna(0)) > 0)).groupby("cell").a.agg(["sum", "count"])
        badf = [c for c, r in pf.iterrows() if cp(int(r["sum"]), int(r["count"]))[0] > .05]
        C = n > 0 and k / n >= .93 and not bad; F = len(f) > 0 and dk / len(f) <= .05 and not badf
        cal.append(dict(wpn=w, est=e, n=n, abstain=round(1 - n / len(x), 3), cov=round(k / max(n, 1), 4), cov_cp=tuple(round(v, 3) for v in cp(k, n)),
                        bad_cells=len(bad), fpr=round(dk / max(len(f), 1), 4), bad_fpr_cells=len(badf), PASS=bool(C and F)))
    cal = pd.DataFrame(cal); cal.to_csv(RES + f"/amend7_calib{SFX}.csv", index=False)
    # min L per (est, fam, wpn)
    x = pd.concat([d[d.kind == "crossover"], load_old32() if not SFX else d.iloc[:0]])
    x = x.assign(det=((x.abstain == 0) & (x.decl_neg == 1)).astype(float), cv=np.where(x.abstain == 0, x.cov_simult, np.nan))
    pw = x.groupby(["est", "fam", "wpn", "L"]).agg(power=("det", "mean"), cov=("cv", "mean"), n=("det", "size")).reset_index()
    rows = []
    for (e, f, w), g in pw.groupby(["est", "fam", "wpn"]):
        gr = G[f] if w != 32 else [l for l in G[f] if l in set(g.L)]   # wpn 32: only sizes that exist
        i = minL(dict(zip(g.L, g.power)), gr)
        rows.append(dict(est=e, fam=f, wpn=w, meas_idx=i, meas_minL=gr[i] if i < len(gr) else f">{gr[-1]}",
                         power=";".join(f"{int(l)}:{p:.2f}" for l, p in zip(g.L, g.power)), xcov=round(g["cov"].mean(), 3)))
    M = pd.DataFrame(rows); M.to_csv(RES + f"/amend7_minL{SFX}.csv", index=False); pw.to_csv(RES + f"/amend7_power{SFX}.csv", index=False)
    # (B) scores
    P = pd.read_csv(RES + "/amend7_predictions.csv")
    S = P.merge(M, on=["est", "fam", "wpn"], how="left"); S["dist"] = (S.pred_idx - S.meas_idx).abs(); S["ok"] = S.dist <= 1
    S.to_csv(RES + f"/amend7_scores{SFX}.csv", index=False)
    print("=== (A) calibration on controls\n", cal.to_string(index=False))
    print("=== min L\n", M.pivot_table(index=["fam", "est"], columns="wpn", values="meas_minL", aggfunc="first").to_string())
    for mod in ("P1", "P2"):
        s = S[(S.model == mod) & S.est.isin(GATED) & S.wpn.isin([8, 128])]
        print(f"=== (B) {mod}: {int(s.ok.sum())}/{len(s)} within one fine-grid step -> {'PASS' if s.ok.sum() >= 14 else 'FAIL'}")
        print(s[["fam", "est", "wpn", "pred_minL", "meas_minL", "dist"]].to_string(index=False))

if __name__ == "__main__": main()
