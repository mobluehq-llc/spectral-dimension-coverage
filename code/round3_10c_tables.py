#!/usr/bin/env python3
"""Round-3 response, item 3: complete Amendment 10c reporting + explicit denominators for Amendment 11 gates (G3, G4, G5).
No estimator is re-run; this only counts stored per-dataset rows.  Usage: round3_10c_tables.py --raw DIR --out DIR
DIR must hold raw_amend10c.csv and raw_amend11.csv (sha256 checked against the committed manifests when present).

Denominator vocabulary (used verbatim in every column):
  n          all datasets generated for the (cell, arm, budget, window rule)
  abstain    datasets on which the arm returned no answer (abstain==1: frozen J<6 rule / solver) or no window (nowin==1)
  answered   n - abstain.  For E7/E8 'abstain' is always 0 by construction; an E7 dataset with no band (infeasible) is ANSWERED, scored
             not-covered and not-declared (fixed denominator).
  cov_cond   covered / answered      (the primary 'coverage' of every amendment)
  cov_strict covered / n             (abstentions count as misses)
  power_ans  correct-sign (drop) declarations / answered ;  power_all = / n   (all 10c cells are realised drops, exact drop >= 0.6)
  wrong_sign declarations of a rise on those cells
  fp (Amendment 11 controls) any-sign declaration at margin 0.25 / answered ; fp_all = / n
"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
RAW = args[args.index("--raw") + 1]; OUT = args[args.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
CMT = os.path.join(HERE, "results_amend10_11")


def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load(name):
    p = os.path.join(RAW, f"raw_{name}.csv")
    m = json.load(open(os.path.join(CMT, f"manifest_{name}.json")))
    h = sha(p); assert h == m["csv_sha256"], f"{name}: sha256 MISMATCH want {m['csv_sha256']} got {h}"
    d = pd.read_csv(p, low_memory=False)
    if "mode" not in d.columns: d["mode"] = "single"
    for c in ("nowin", "decl_neg", "decl_pos", "infeasible"): d[c] = d[c].fillna(0) if c in d.columns else 0
    d["answered"] = ((d.abstain == 0) & (d.nowin == 0)).astype(int)
    d["abst"] = 1 - d.answered
    d["cov"] = d.answered * (d.cov_simult.fillna(0) == 1)
    d["neg"] = d.answered * (d.decl_neg == 1); d["pos"] = d.answered * (d.decl_pos == 1)
    d["anydecl"] = d.answered * ((d.decl_neg + d.decl_pos) > 0)
    d["nobandE7"] = ((d.est == "E7") & (d.infeasible == 1)).astype(int)
    return d


def agg(d, keys):
    g = d.groupby(keys).agg(n=("ds", "size"), abstain=("abst", "sum"), answered=("answered", "sum"), k_cov=("cov", "sum"), k_drop=("neg", "sum"), k_rise=("pos", "sum"),
                            k_anydecl=("anydecl", "sum"), n_no_band_E7=("nobandE7", "sum")).reset_index()
    g["cov_cond"] = g.k_cov / g.answered.replace(0, np.nan); g["cov_strict"] = g.k_cov / g.n
    lo_hi = [cp(k, n) for k, n in zip(g.k_cov, g.answered)]; g["cov_cp_lo"] = [x[0] for x in lo_hi]; g["cov_cp_hi"] = [x[1] for x in lo_hi]
    g["power_ans"] = g.k_drop / g.answered.replace(0, np.nan); g["power_all"] = g.k_drop / g.n
    lo_hi = [cp(k, n) for k, n in zip(g.k_drop, g.answered)]; g["power_cp_lo"] = [x[0] for x in lo_hi]; g["power_cp_hi"] = [x[1] for x in lo_hi]
    g["fp_ans"] = g.k_anydecl / g.answered.replace(0, np.nan); g["fp_all"] = g.k_anydecl / g.n
    return g


def main():
    d = load("amend10c"); meta = pd.read_csv(os.path.join(CMT, "cells_meta_amend10c.csv")).set_index("id")
    d["fam"] = d.cell.map(meta.fam); d["L"] = d.cell.map(meta.L); d["real"] = d.cell.map(meta.real); d["exact_J"] = d.cell.map(meta.exact_J)
    sep = d.cell.isin([302, 303]); main_ = d[~sep]
    # (a) per cell x arm x budget x window rule
    cellt = agg(main_, ["fam", "L", "cell", "est", "wpn", "mode"]).sort_values(["fam", "L", "est", "wpn", "mode"])
    cellt = cellt.merge(meta[["exact_J", "J_fixed", "max_drop"]], left_on="cell", right_index=True, how="left")
    cellt.to_csv(os.path.join(OUT, "r3_10c_cell_arm_budget_window.csv"), index=False)
    # (b) pooled: both families, and per family
    pf = agg(main_, ["fam", "est", "wpn", "mode"]); pb = agg(main_.assign(fam="both"), ["fam", "est", "wpn", "mode"])
    pooled = pd.concat([pb, pf]).sort_values(["est", "mode", "fam", "wpn"])
    pooled.to_csv(os.path.join(OUT, "r3_10c_pooled.csv"), index=False)
    # (c) separation-scan cells (ratios 2.5 and 3.0), reported for completeness
    sept = agg(d[sep].assign(fam=d[sep].cell.map({302: "sep_ratio2.5", 303: "sep_ratio3.0"})), ["fam", "est", "wpn", "mode"])
    sept.to_csv(os.path.join(OUT, "r3_10c_separation_scan.csv"), index=False)
    # sanity: counts equal the committed aggregates used by the verdict script
    ag = pd.read_csv(os.path.join(CMT, "agg_amend10c_cells.csv")); chk = cellt.merge(ag, on=["cell", "est", "mode", "wpn"], suffixes=("", "_c"))
    assert len(chk) == len(cellt) and (chk.k_cov == chk.k_cov_c).all() and (chk.answered == chk.n_valid).all() and (chk.n == chk.n_c).all(), "disagrees with committed agg_amend10c_cells.csv"
    print("10c: per-cell table", len(cellt), "rows; agrees with committed aggregates (k_cov, n_valid, n)")
    # (d) Amendment 11
    e = load("amend11"); m11 = pd.read_csv(os.path.join(CMT, "cells_meta_amend11.csv")).set_index("id")
    for c in ("fam", "L", "kind", "real", "exact_J"): e[c] = e.cell.map(m11[c])
    e = e[e.est != "E6"]
    t11 = agg(e, ["kind", "fam", "L", "cell", "est", "wpn"]).sort_values(["kind", "fam", "L", "est", "wpn"])
    t11.to_csv(os.path.join(OUT, "r3_11_cell_arm_budget.csv"), index=False)
    p11 = pd.concat([agg(e.assign(fam="pooled"), ["kind", "fam", "est", "wpn"])])
    p11["real"] = ""
    p11.to_csv(os.path.join(OUT, "r3_11_pooled.csv"), index=False)
    px = agg(e[(e.kind == "crossover") & e.real].assign(fam="realised crossovers"), ["kind", "fam", "est", "wpn"])
    px.to_csv(os.path.join(OUT, "r3_11_pooled_realised.csv"), index=False)
    # G5 denominators spelled out
    g5 = agg(e[(e.kind == "control") & (e.wpn == 32)].assign(fam="controls"), ["fam", "est", "wpn"])
    g5["fp_lo"] = [cp(k, n)[0] for k, n in zip(g5.k_anydecl, g5.answered)]; g5["fp_hi"] = [cp(k, n)[1] for k, n in zip(g5.k_anydecl, g5.answered)]
    g5.to_csv(os.path.join(OUT, "r3_11_G5_denominators.csv"), index=False)
    # E7 solver-failure accounting on controls and realised cells, all budgets
    sf = e[e.est == "E7"].groupby(["kind", "wpn"]).agg(n=("ds", "size"), no_band=("nobandE7", "sum"), lp_status_ne0=("lp_status", lambda s: int((s.fillna(0) != 0).sum())), gate_fail=("gate_fail", lambda s: int(s.fillna(0).sum()))).reset_index()
    sf.to_csv(os.path.join(OUT, "r3_11_E7_solver_accounting.csv"), index=False)
    print(g5.round(4).to_string()); print(sf.to_string())


if __name__ == "__main__": main()
