# P22 Amendment 13: results (2026-11-24). Bars exactly as registered in PLAN.md Amendment 13; nothing added.

Runs: 13a on host-D (WSL, `<local path omitted> six budget/design runs, all manifests complete=True, 12,800 rows in the last); 13b on host-T (`~/p22_a13`, ens then data, manifest complete=True, 36,000 rows). Both finished without restarts. Outputs copied to `results_amend13/`; sha256 of every raw CSV, manifest and the 13b cell table re-checked locally against the host: all match. Analysis: `analyze_amend13.py a|b|c` run on host-T in the pinned venv (reports in `results_amend13/amend13{a,b,c}_report.txt`).

## Harness fixes (analysis script only; no bar, estimator or data changed)
Two bugs in `analyze_amend13.py` surfaced on the first real 13b data: (1) line 140 used `r.cov` (a DataFrame method) instead of `r["cov"]`, crashing the report; (2) line 117 required `nowin == 0`, but `nowin` is NaN for E1q/E4c/E4c_V2 rows (only E7 writes it), so those arms had zero valid datasets (a false "FAIL" on every arm but E7, printed before the fix). Both are fixed; numbers below are from the fixed script.

## 13a: V2 variants (verdict = VALIDATED in both designs)
- E4o_V2 and E4c_V2: VALIDATED in the replication design and in the holdout design (C pooled control coverage >= 93% at every wpn with no cell CP-upper < 93%; F at wpn 32 <= 5%; X >= 90% at every wpn).
- E4o_V2 control coverage 95.8-97.0% (V0 E4o: 92.7-95.4%, FAILS C at wpn 8 and 32 in rep and at wpn 8 and 32 in hold; the ring L=512 wpn 32 defect is repaired). E4c_V2 99.2-99.8% on controls, crossover X 97.1-99.5%, F 0.0% (corrected 2026-11-24 to the report file amend13a_report.txt; the earlier prose said 98.8-99.8% and 96.9-99.5%). Paired controls: V2-only covers outnumber V0-only covers by 3-6x for E4o and 2-10x for E4c.
- Prediction 1 holds (E4o_V2 passes C at every ring control cell: L=128 at wpn 8 is 7/9, pooled ring cells at L >= 256 are 190-198/200). Prediction 2 holds (E4c_V2 within 2 points of E4c). Predictions 3-4 (half-width ratio, minimum detectable L) are in the report; minL equals V0's within one grid step for E4c (rep wpn 32 aniso3_w8 160 -> 192 is one step) and equals V0 for E4o_V2 at wpn 128.
- Paper reading: V2 is validated on both designs and may be adopted; V0 E4c results stay as published.

## 13b: random graphs (WS and 2D percolation; 300 datasets per cell, 10 cells, 3 wpn)
Cell classes (ensemble truth): crossover-like = the five WS cells (max excursion 0.93-2.28); flat-like = the five percolation cells (0.06-0.20).

Single-realisation target (PRIMARY):
- B-X crossover-like coverage >= 90% at wpn 8/32/128: E4c 98.9/98.8/98.3, E4c_V2 99.4/99.4/99.3, E7 99.9/100/100: all PASS.
- B-C flat-like >= 93% and no cell CP-upper < 93%: E4c 99.3/98.5/98.7, E4c_V2 99.3/99.1/99.1, E7 about 100: all PASS.
- B-G1 E1q crossover-like coverage at wpn 32 <= 40%: 56.1% -> FAIL (prediction 1 half wrong; E1q still degrades).
- B-G2 E1q strictly decreasing 8 > 32 > 128 with CP-upper(128) < CP-lower(8): 76.3 > 56.1 > 38.0, 40.5 < 74.1 -> PASS.
- B-F false-change rate on flat-like cells, wpn 32: E4c 0/1341 (CP 0-0.3%), E4c_V2 0/1341 -> PASS.

Ensemble-mean target (reported as a property of the target, not an estimator failure, per registered prediction 3):
- E7 PASS everywhere (crossover 99.2-100%, flat 99.9-100%).
- E4c and E4c_V2 PASS at wpn 8 and 32 but FAIL at wpn 128: crossover-like 80.4% (E4c) and 82.7% (E4c_V2); flat-like 85.2% / 87.9%, percolation cells 4, 5, 9 have CP-upper < 93. Coverage falls as budget grows (98 -> 94-95 -> 80-83), the signature of a band that tightens around one graph while the target is the ensemble average. E1q on this target: crossover 68.3/39.6/22.2%, flat 87.6/65.2/36.1%.
- Prediction 3 holds (ensemble coverage <= single-realisation for every arm). Prediction 4 (V2 minus V0 within 2 points) FAILS on the ensemble-mean target at wpn 128: +2.3 points on crossover-like cells (82.7 vs 80.4) and +2.7 on flat-like cells (87.9 vs 85.2); the flat-like cells at wpn 32 also exceed 2 points (+2.1: 1462/1500 vs 1430/1500); it holds on the single-realisation target at every budget (largest +1.0) and on the ensemble target at wpn 8 (corrected 2026-11-24 to the report file amend13b_report.txt; the earlier prose said it held). No verdict changes.

## 13c (EXPLORATORY, fit only)
`amend13c_report.txt`: E7 has 0 finite points so no fit; E4o_V2 censoring consistency 3/7. No gate.

## Unblocks
A v8 may state: V2 validated on both designs; the random-graph arm passes B-X, B-C, B-F for E4c, E4c_V2, E7 on the single-realisation target; E1q fails B-G1 (56% not <= 40%) though B-G2 holds; ensemble-mean coverage of E4c/E4c_V2 fails at wpn 128 (target mismatch, said as such).
