#!/usr/bin/env python3
"""Round-3 response, item 1: per-dataset crosswalk of the E7/E7j counts between
  VERSION T  = Table IV scoring: results/raw_amend8_core.csv (E7 rows; Amendment 8; run on host-T, numpy 2.4.6); withheld E7j decision from results/raw_amend9_a.csv (host-T)
  VERSION D  = Amendment 10b 'original' columns o_* and certified columns B_*/A_*/Bj_*: results_amend10_11/raw_amend10b.csv (run on host-D, numpy 2.5.5)
Datasets are identified by (cell id, dataset index ds); the seed of every one is SeedSequence([22, cell, ds]) in both versions.
Usage: round3_e7_crosswalk.py --out DIR"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]; OUT = args[args.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
R = os.path.join(HERE, "results"); C = os.path.join(HERE, "results_amend10_11")

a8 = pd.read_csv(os.path.join(R, "raw_amend8_core.csv"), low_memory=False)
for c in ("cell", "ds"): a8[c] = a8[c].astype(int)
e7 = a8[a8.est == "E7"].copy(); e6 = a8[a8.est == "E6"][["cell", "ds", "abstain", "J"]].rename(columns={"abstain": "E6_abstain", "J": "E6_J"})
e7 = e7.merge(e6, on=["cell", "ds"], how="left")
meta = pd.read_csv(os.path.join(R, "cells_meta_amend8_core.csv")).set_index("id")
e7["kind"] = e7.cell.map(meta.kind); e7["fam"] = e7.cell.map(meta.fam); e7["L"] = e7.cell.map(meta.L); e7["max_drop"] = e7.cell.map(meta.max_drop)
e7["group"] = np.where(e7.kind == "crossover", np.where(e7.max_drop >= 0.6, "crossover_realised", "crossover_not_realised"), e7.kind)
# ---- Table IV regeneration (counts the paper prints)
t4 = dict(datasets=len(e7), no_band=int((e7.infeasible == 1).sum()), atom_sensitive=int((e7.gate_fail == 1).sum()),
          atom_sensitive_and_no_band=int(((e7.gate_fail == 1) & (e7.infeasible == 1)).sum()))
print("Table IV (regenerated from raw_amend8_core):", t4)
assert (t4["datasets"], t4["no_band"], t4["atom_sensitive"], t4["atom_sensitive_and_no_band"]) == (9400, 66, 170, 9)

b = pd.read_csv(os.path.join(C, "raw_amend10b.csv"), low_memory=False)
sel = json.load(open(os.path.join(C, "amend10b_selection.json")))
S170 = {tuple(x) for x in sel["sens"]}; WH = {tuple(x) for x in sel["withheld"]}
x = e7[["cell", "ds", "group", "fam", "L", "J", "E6_J", "E6_abstain", "infeasible", "gate_fail", "natoms", "lp_status", "cov_simult", "n_boot", "n_cp", "c_box"]].rename(
    columns={"J": "T_J", "E6_J": "T_E6_J", "E6_abstain": "T_E6_abstain", "infeasible": "T_no_band", "gate_fail": "T_atom_sensitive", "natoms": "T_final_natoms", "lp_status": "T_lp_status",
             "cov_simult": "T_covered", "n_boot": "T_n_boot", "n_cp": "T_n_cp", "c_box": "T_c_box"})
x["in_170_sens"] = [(c, d) in S170 for c, d in zip(x.cell, x.ds)]
bb = b[["cell", "ds", "group", "J", "o_status", "o_natoms", "o_gate_fail", "o_cov", "B_status", "B_cov", "A_status", "A_cov", "B_contains_orig", "Bj_neg", "Bj_pos", "Bj_amb", "Bj_certnon", "Bj_minU_drop"]].rename(
    columns={"group": "D_group_10b", "J": "D_J", "o_status": "D_orig_status", "o_natoms": "D_orig_natoms", "o_gate_fail": "D_orig_atom_sensitive", "o_cov": "D_orig_covered", "B_status": "D_B_status",
             "B_cov": "D_B_covered", "A_status": "D_A_status", "A_cov": "D_A_covered", "B_contains_orig": "D_B_contains_orig", "Bj_neg": "D_E7j_drop", "Bj_pos": "D_E7j_rise", "Bj_amb": "D_E7j_ambiguous_times",
             "Bj_certnon": "D_E7j_certified_nondeclarable_times", "Bj_minU_drop": "D_E7j_minU_drop"})
cw = x.merge(bb, on=["cell", "ds"], how="left")
cw["in_10b_1030"] = cw.D_group_10b.notna()
cw["D_orig_no_band"] = (cw.D_orig_status != 0).where(cw.in_10b_1030)
cw["D_B_no_band"] = (cw.D_B_status != 0).where(cw.in_10b_1030)
cw["J_equal"] = (cw.T_J == cw.D_J).where(cw.in_10b_1030)
# withheld E7j decision (Amendment 9a)
a9 = pd.read_csv(os.path.join(R, "raw_amend9_a.csv"), low_memory=False)
a9 = a9[a9.gate_flag.fillna(0) > 0][["cell", "ds", "decl_neg", "minU_drop", "gate_flag", "gate_rechecked", "n_amb", "nbad"]].rename(columns=lambda c: c if c in ("cell", "ds") else "T9a_" + c)
a9["cell"] = a9.cell.astype(int); a9["ds"] = a9.ds.astype(int)
cw = cw.merge(a9, on=["cell", "ds"], how="left"); cw["withheld_E7j"] = [(c, d) in WH for c, d in zip(cw.cell, cw.ds)]
# host-T re-evaluation of the original E7 LP on the 10b-code path (11 ids)
xp = pd.DataFrame(json.load(open(os.path.join(HERE, "results_round3", "xplat_hostT_E7_status.json")))).rename(
    columns={"J": "X_hostT_J", "status": "X_hostT_orig_status", "natoms": "X_hostT_natoms", "gate_fail": "X_hostT_atom_sensitive"})
cw = cw.merge(xp, on=["cell", "ds"], how="left")
keep = cw[cw.in_170_sens | cw.in_10b_1030 | cw.withheld_E7j | (cw.T_no_band == 1)].sort_values(["cell", "ds"])
keep.to_csv(os.path.join(OUT, "r3_e7_crosswalk_per_dataset.csv"), index=False)

# ---- summary tables
rows = []
def add(label, df, ver_cols):
    rows.append(dict(set=label, n=len(df), **ver_cols(df)))
T = lambda d: dict(version="T (Table IV, host-T)", no_band=int((d.T_no_band == 1).sum()), atom_sensitive=int((d.T_atom_sensitive == 1).sum()), covered_with_band=int(((d.T_no_band == 0) & (d.T_covered == 1)).sum()),
                   band_not_covered=int(((d.T_no_band == 0) & (d.T_covered == 0)).sum()))
D = lambda d: dict(version="D (10b original, host-D)", no_band=int((d.D_orig_no_band == True).sum()), atom_sensitive=int((d.D_orig_atom_sensitive == 1).sum()), covered_with_band=int(((d.D_orig_status == 0) & (d.D_orig_covered == 1)).sum()),
                   band_not_covered=int(((d.D_orig_status == 0) & (d.D_orig_covered == 0)).sum()))
DB = lambda d: dict(version="D (10b certified B, host-D)", no_band=int((d.D_B_status != 0).sum()), atom_sensitive=np.nan, covered_with_band=int(((d.D_B_status == 0) & (d.D_B_covered == 1)).sum()),
                    band_not_covered=int(((d.D_B_status == 0) & (d.D_B_covered == 0)).sum()))
s170 = cw[cw.in_170_sens]; s1030 = cw[cw.in_10b_1030]
for lab, df in (("170 atom-sensitive ids (Table IV selection rule)", s170), ("1030 ids of Amendment 10b", s1030), ("all 9400 datasets", cw)):
    for f in (T, D, DB):
        if lab.startswith("all") and f is not T: continue
        add(lab, df, f)
summ = pd.DataFrame(rows); summ.to_csv(os.path.join(OUT, "r3_e7_crosswalk_summary.csv"), index=False); print(summ.to_string())
# concordance between versions on the 1030 ids
cc = pd.crosstab([s1030.T_no_band.map({0.0: "T: band", 1.0: "T: no band"})], s1030.D_orig_no_band.map({False: "D: band", True: "D: no band"}))
cc.to_csv(os.path.join(OUT, "r3_e7_crosswalk_nobands_matrix_1030.csv")); print(cc)
cc2 = pd.crosstab(s170.T_no_band.map({0.0: "T: band", 1.0: "T: no band"}), s170.D_orig_no_band.map({False: "D: band", True: "D: no band"})); cc2.to_csv(os.path.join(OUT, "r3_e7_crosswalk_nobands_matrix_170.csv")); print(cc2)
jc = pd.crosstab(s1030.cell.map(lambda c: "3D min-L cells" if c in (702, 26, 703) else "other"), s1030.J_equal); print("window length J equal between T and D:\n", jc)
jc.to_csv(os.path.join(OUT, "r3_e7_crosswalk_J_concordance.csv"))
print("atom-sensitivity flag concordance (170 ids): T flag = 1 for all; D flag = 1 for", int((s170.D_orig_atom_sensitive == 1).sum()), "of", len(s170))
print("withheld:", keep[keep.withheld_E7j][["cell", "ds", "T9a_decl_neg", "T9a_minU_drop", "T9a_gate_flag", "D_E7j_drop", "D_E7j_minU_drop", "D_E7j_ambiguous_times", "D_B_status"]].to_string())
