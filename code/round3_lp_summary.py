#!/usr/bin/env python3
"""Summarises results_round3/lp_verify_results_final.json (exact-rational LP verification, round 3 item 2). Usage: round3_lp_summary.py --out DIR"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); args = sys.argv[1:]; OUT = args[args.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
src = os.path.join(HERE, "results_round3", "lp_verify_results_final.json")
if not os.path.exists(src): src = os.path.join(HERE, "results_round3", "lp_verify_results.json")
d = json.load(open(src)); sel = {(s["cell"], s["ds"]): s["why"] for s in json.load(open(os.path.join(HERE, "results_round3", "lp_verify_selection.json")))}
rows = []
for r in d:
    rows.append(dict(cell=r["cell"], ds=r["ds"], why="; ".join(sel[(r["cell"], r["ds"])]), J=r["J"], B_status=r["B_status"], n_lp=r.get("n_lp", 0), n_exact=r.get("n_exact", 0), n_failed=r.get("n_failed", 0),
                     n_needed_exact_pivots=r.get("n_pivoted", 0), max_solver_basis_dual_violation=r.get("max_solver_basis_dual_viol", np.nan), max_primal_violation_exact=r.get("max_primal_viol", np.nan),
                     max_dual_violation_exact=r.get("max_dual_viol", np.nan), max_abs_gap_dualsimplex_vs_ipm=r.get("max_ds_ipm_gap", np.nan), max_rel_gap_reported_vs_exact=r.get("max_rel_gap_exact_vs_solver", np.nan),
                     empty_set_certified_exactly=r.get("empty_set_certified_exactly", ""), retried=r.get("retried", False)))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "lp_verify_summary.csv"), index=False)
tot = dict(datasets=len(T), empty_set_datasets=int((T.B_status != 0).sum()), empty_certified=int((T.empty_set_certified_exactly == True).sum()), lps=int(T.n_lp.sum()), lps_exact=int(T.n_exact.sum()), lps_not_exact=int(T.n_failed.sum()),
           lps_needing_exact_pivots=int(T.n_needed_exact_pivots.sum()), max_solver_basis_dual_violation=float(T.max_solver_basis_dual_violation.max()), max_abs_gap_ds_vs_ipm=float(T.max_abs_gap_dualsimplex_vs_ipm.max()),
           max_rel_gap_reported_vs_exact=float(T.max_rel_gap_reported_vs_exact.max()), datasets_with_any_unverified_lp=int((T.n_failed > 0).sum()))
json.dump(tot, open(os.path.join(OUT, "lp_verify_totals.json"), "w"), indent=1); print(tot)
