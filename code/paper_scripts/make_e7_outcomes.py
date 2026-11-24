#!/usr/bin/env python3
"""Exhaustive E7 outcome categories at 32 walkers per node (Amendment 8 raw data), fixed denominators.
Stdlib only. Reads ../results/raw_amend8_core.csv and ../results/cells_meta_amend8_core.csv,
writes ../results/amend8_e7_outcomes.csv. No estimator is re-run; this only counts stored rows.
Categories (mutually exclusive, per dataset): infeasible (no measure in the box, so no band; a miss),
atom-sensitive (feasibility changed when atoms were doubled), covered (band contains the exact d(t))."""
import csv, collections, pathlib
R = pathlib.Path(__file__).resolve().parent.parent / "results"
meta = {r["id"]: r for r in csv.DictReader(open(R / "cells_meta_amend8_core.csv"))}
def f(v):
    try: return float(v)
    except ValueError: return 0.0
cnt = collections.Counter()
for r in csv.DictReader(open(R / "raw_amend8_core.csv")):
    if r["est"] != "E7": continue
    m = meta[str(int(float(r["cell"])))]
    grp = m["kind"]
    if grp == "crossover": grp = "crossover_realised" if float(m["max_drop"]) >= 0.6 else "crossover_not_realised"
    key = (grp, int(f(r["infeasible"]) == 1), int(f(r["gate_fail"]) == 1), int(f(r["cov_simult"]) == 1))
    cnt[key] += 1; cnt[("ALL",) + key[1:]] += 1
groups = sorted({k[0] for k in cnt})
with open(R / "amend8_e7_outcomes.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["group", "n", "covered", "not_covered_with_band", "infeasible_no_band", "atom_sensitive",
                "infeasible_and_atom_sensitive", "atom_sensitive_with_band_covered", "atom_sensitive_with_band_not_covered",
                "stable_feasible_covered", "stable_feasible_not_covered"])
    for g in groups:
        c = lambda i, a, v: sum(n for k, n in cnt.items() if k[0] == g and (i is None or k[1] == i) and (a is None or k[2] == a) and (v is None or k[3] == v))
        row = [g, c(None, None, None), c(None, None, 1), c(0, None, 0), c(1, None, None), c(None, 1, None),
               c(1, 1, None), c(0, 1, 1), c(0, 1, 0), c(0, 0, 1), c(0, 0, 0)]
        assert c(1, None, 1) == 0, "an infeasible dataset is scored as covered"
        w.writerow(row); print(row)
