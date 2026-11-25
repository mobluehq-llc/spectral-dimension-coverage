# Round-3 referee response: data and analyses (paper A, p22)

2026-11-23. Answers the analysis items of `ideas/A_spectral_dimension_coverage/REVIEWS/external_referee_chatgpt_round3_2026-11-23.md` (sections 3, 6, 7).
Every number below is produced by a script in this directory from per-dataset files; the one command `archive/regenerate_all.py` regenerates all of the
tables (item 6). Only item 4 and the exact LP checks of item 2 are new compute; both are labelled. No confirmatory run was made, no bar, gate, seed or
estimator was changed.

**Compute log.** E4 diagnosis: host-D (WSL, 4 workers, pinned venv: numpy 2.5.5 / scipy 1.18.1 / pandas 3.1.1 / python-flint 0.9.0 / highspy 1.16.0), after checking
load (17 on 16 threads, another job running; my 4 workers were the only addition, 8 at the peak with the LP job). LP verification: host-D, 4 workers. Cross-platform
check of item 1: 3 and 11 datasets on host-T (1 core, system python, numpy 2.4.6). All other items are pure counting on the Mac in the pinned environment of `archive/`.

## 1. E7/E7j counts: per-dataset crosswalk (`results_round3/r3_e7_crosswalk_per_dataset.csv`)

**The two sets of counts are two computations of the same seeds, not two scorings of the same data.** Datasets are identified by (cell id, ds); the seed is
`SeedSequence([22, cell, ds])` in both. Table IV (Amendment 8) was run on host-T (numpy 2.4.6); Amendment 10b was run on host-D (numpy 2.5.5). The walk simulation (`rng.binomial`) does not give the same walks in different environments
(a Mac arm64 run with numpy 2.5.5 gives a third answer, J = 13, 13, 14 on the three test ids below), so the 10b run regenerated different realisations of the same 1030 seeds. Library version versus CPU was not isolated.

| evidence | result |
|---|---|
| window length J (host-T realisation vs 10b/host-D) equal | 752 of 1030 ids (3D minimum-size cells 384/600; other cells 368/430) |
| cells 102/103, ds 16, 46, 27: J on host-T / host-D / the Table IV file | host-T 13, 14, 14 and host-D 14, 13, 13; Table IV file 13, 14, 14 (the host-T run reproduces the stored file; host-D does not) |
| the 9 no-band atom-sensitive ids of Table IV, rerun on host-T with the 10b code path | all 9 return LP status 4 (no band), atom-sensitive, 1597 atoms (9/9 reproduce) |
| the 2 ids that are no-band only in 10b (402/85, 400/119), rerun on host-T | both feasible (status 0), atom-sensitive |

Counts, by version (T = Table IV, host-T; D = 10b "original" and certified columns, host-D):

| set | version | no band | atom-sensitive | covered with band | band, not covered |
|---|---|---|---|---|---|
| the 170 atom-sensitive ids | T Table IV | 9 | 170 | 161 | 0 |
| the 170 ids | D original E7 | 3 | 14 | 167 | 0 |
| the 170 ids | D certified (B) | 0 | | 170 | 0 |
| all 1030 ids of 10b | T | 18 | 170 | 1012 | 0 |
| all 1030 ids | D original E7 | **6** | 44 | 1024 | 0 |
| all 1030 ids | D certified (B) | **1** (id 26/60) | | 1029 | 0 |
| all 9400 | T | 66 | 170 | 9332 | 2 |

Reconciliation of the referee's "nine vs six": Table IV's nine (the 66 no-band datasets include nine of the 170) are host-T realisations. The six in 10b are 3 of
the 170 ids plus 3 ids from the 3D minimum-size cells, all on host-D. Overlap on the 1030 ids: both versions no band 1 (8/96); T only 17; D only 5 (matrix in
`r3_e7_crosswalk_nobands_matrix_1030.csv`; 170 ids: 1 / 8 / 2). Only 14 of the 170 ids are atom-sensitive on host-D (and 30 of the 860 other ids are), so atom
sensitivity and infeasibility are properties of a realisation near a feasibility boundary, not of a seed. The 10b statement that it "includes all 170" is true of
the ids and false of the realisations; the sentence should say "the 170 seeds flagged in the host-T run, regenerated on host-D".

**Withheld E7j decision (Amendment 9a): id 26/192.** host-T: E7j declares a drop (minU = -0.365), gate flag 1 (rechecked 1), so it was withheld. host-D: the original E7
band exists (status 0, not atom-sensitive); the certified E7j declares a drop with minU = -0.274 (margin 0.024 below -0.25), one ambiguous time, no solver trouble,
certified band contains the original. So the withheld decision is a drop on both realisations.

Suggested text: Table IV counts are the host-T run; 10b counts are a host-D regeneration of the same seeds; the six original no-band cases are from that regeneration;
the certified result (1 no band in 1030, id 26/60, emptiness certified exactly in item 2) is the one on which the structure-based claims rest.

## 2. Verification of the 10b outer-relaxation LPs (`round3_lp_verify.py`; results `results_round3/lp_verify_*`)

**Tail treatment (as implemented, `p22amend10.cell_constants`).** Relaxation B partitions the support of the spectral measure in s = -ln(lambda), s in [0, inf]
(lazy walk, lambda in [0,1]), into cells: the atom Z = {0}; 1597 geometric cells C_k tiling (0, 1], whose first cell (0, e_1] reaches s = 0 so nothing
between 0 and the smallest grid atom 0.02/t_max is lost; and a tail cell INF = [1, inf]. For each time t_j, a measure with mass m_k on cell k has
q_j in [sum m_k e^{-b_k t_j}, sum m_k e^{-a_k t_j}] and numerator sum m_k s e^{-s t_j} in [sum m_k Gm_kj, sum m_k Gp_kj] (cell minima/maxima of s e^{-st}, using the
peak at s = 1/t when inside the cell). The tail INF therefore contributes e^{-t_j} (at most 4e-11 at t = 24) to the upper q-enclosure and 0 to the lower one, and at
most the cell maximum of s e^{-st} (at s = 1, e^{-t}) to the numerator: the tail is enclosed, not truncated. All e^{-st} and s e^{-st} constants come from python-flint arb balls (80 bits) and are
rounded outward to doubles; the box constraints are relaxed independently per time (an outer set). The upper endpoint is the minimum of two LPs (denominator normalised; denominator replaced by the box bound), the lower endpoint the maximum.

**What was verified and how.** For each selected dataset, the box is regenerated from its seed and every endpoint LP of `bound_B` (both variants, both sides, all
window times) is rebuilt and solved with HiGHS dual simplex (presolve off, to obtain a basis). The basis is then checked **in exact rational arithmetic**
(every double read as the exact rational it is; python-flint `fmpq_mat` solves): primal feasibility of the basic solution and sign of every reduced cost. If the
solver's basis is exactly primal feasible but dual infeasible by a solver-tolerance amount, exact rational simplex pivots (Bland) continue to optimality. If a basis
is not exactly primal feasible, or a pivot limit is hit, the LP is retried once with tolerance 1e-10 and no scaling and a 1000-pivot limit. A passing LP has a
rational optimal value, so no solver tolerance enters its verdict. A second algorithm (HiGHS interior point) is compared on every LP. Emptiness of the relaxed set (the one 'no certified band' dataset)
is certified by the always-feasible phase-1 LP min sum of artificials, whose exact optimum must be strictly positive.
Sample (fixed by rule and seed 3101 and committed before any LP was solved: `lp_verify_selection.json`): all datasets with any non-zero B/A/original status (7), 8 random of each of the three 3D
minimum-size cells (L = 160, 192, 256), 8 random atom-sensitive ids, 4 random stable ids, and the withheld dataset: 44 datasets.

**Result (44 datasets, 2016 endpoint LPs plus one emptiness certificate).**

| quantity | value |
|---|---|
| endpoint LPs verified exactly (rational optimum, primal and dual feasibility exact) | **2016 of 2016** |
| of these, solver basis already exact | 1646 |
| of these, needed exact rational pivots after HiGHS stopped (its basis was dual infeasible by up to 9.9e-8, its own 1e-7 tolerance) | 370 |
| 14 datasets with a first-pass failure (pivot limit 60, one solver error, one primal infeasibility of 2.5e-8) | all 14 pass on the retry (1000 pivots, tolerance 1e-10, no scaling) |
| largest relative difference between the reported endpoint (HiGHS, double) and the exact optimum, all endpoints | 8.2e-8 |
| largest absolute difference between HiGHS dual simplex and HiGHS interior point on any LP | 2.9e-9 |
| dataset with an empty relaxed set (26/60): phase-1 LP exact optimum strictly positive | **certified** |

So every sampled band endpoint is the exact optimum of the stated relaxed problem to within 8.2e-8 relative, and the one 'no certified band' dataset (26/60) is empty exactly. The first-pass failures were solver stalls, not wrong values (dual violations of 1e-16 to 1e-14 at the pivot limit); all are resolved.
Per-dataset detail: `lp_verify_summary.csv`, raw: `lp_verify_results.json` (first pass), `lp_verify_results_final.json` (after retry). Compute: host-D, about 20 minutes on 4 workers plus about 15 minutes of retry on 5.

**Joint-test decisions against solver tolerance (all 1030 datasets, no recomputation; `lp_tolerance_margins.csv`).** The E7j decision "drop" requires the certified upper bound U
on the largest decrease to be below -0.25. Of the 570 datasets with a declared drop, the smallest margin (-0.25 - U) is 2.4e-5 (id 26/143, L = 192); 2 are within 1e-4 (2.4e-5 and 7.1e-5, ids 26/143 and 702/186) and 8 within 1e-3 of the threshold;
none within 1e-5. Every margin exceeds the largest solver-versus-exact endpoint discrepancy found above (8.2e-8 relative, so at most about 4e-7 absolute for d up to 5) by a factor of at least 60, and by more than 1000 for all but 5 of the 570 declarations (those with margin below 4e-4). The nearest NOT-declared finite U is 2.7e-3 above the threshold. Stress test at 3D L = 160: if every drop declared within 1e-3 of the threshold (4 datasets) were lost, power would be
166/200 = 83.0% against 170/200 = 85.0%, still above 0.8 on the point estimate; the CP interval problem noted in the text (lower end 79.3%) is unchanged. No dataset has solver-trouble flags (Bj_nbad = 0 for all 1029 with a band); 428 have at least one
ambiguous time (counted as no declaration, as in the text). Limit: the exact check covers the band LPs; the branch-and-bound U values of E7j (370,259 LPs) are checked by the margin analysis, not LP by LP.

## 3. Amendment 10c complete reporting (`results_round3/r3_10c_*.csv`; script `round3_10c_tables.py`)

* `r3_10c_cell_arm_budget_window.csv`: 324 rows = 18 realised cells x 3 arms (E1q, E4c, E4r) x 3 budgets (8, 32, 128) x 2 window rules (fixed, orig). Columns: n, abstain, answered, covered,
  conditional coverage with CP, strict coverage (abstentions = misses), drops declared (power, correct sign), power over answered and over all n with CP, rises declared, plus the cell's exact J, fixed J and exact max drop.
  Counts agree with the committed aggregates the verdict script uses (asserted).
* `r3_10c_pooled.csv`: the same columns pooled over both families and per family; `r3_10c_separation_scan.csv` for the ratio 2.5 and 3.0 cells.
* **Denominators.** n = all datasets generated for the cell/arm/budget/rule. abstain = no answer (frozen J < 6 rule; forced when the exact window has fewer than six points). answered = n - abstain.
  Conditional coverage = covered/answered (the primary figure of every amendment). Strict coverage = covered/n. Power = declared drops/answered, and /n; all 10c cells are realised drops (exact drop >= 0.6) so a declared drop is the correct sign. No false-positive column applies to 10c (no control cells; the 10c cells are all realised crossovers).

Pooled over both families (full table in the csv):

| arm | rule | wpn | n | abstain | answered | covered | cov_cond | cov_strict | drops | power_answered | power_all |
|---|---|---|---|---|---|---|---|---|---|---|---|
| E1q | fixed | 8 | 3600 | 400 | 3200 | 1292 | 1292/3200 = 40.4% (38.7-42.1) | 35.9% | 2306 | 72.1% | 64.1% |
| E1q | fixed | 32 | 3600 | 400 | 3200 | 627 | 627/3200 = 19.6% (18.2-21.0) | 17.4% | 2871 | 89.7% | 79.8% |
| E1q | fixed | 128 | 3600 | 400 | 3200 | 114 | 114/3200 = 3.6% (2.9-4.3) | 3.2% | 3106 | 97.1% | 86.3% |
| E1q | orig | 8 | 3600 | 685 | 2915 | 1395 | 1395/2915 = 47.9% (46.0-49.7) | 38.8% | 1841 | 63.2% | 51.1% |
| E1q | orig | 32 | 3600 | 403 | 3197 | 613 | 613/3197 = 19.2% (17.8-20.6) | 17.0% | 2852 | 89.2% | 79.2% |
| E1q | orig | 128 | 3600 | 395 | 3205 | 120 | 120/3205 = 3.7% (3.1-4.5) | 3.3% | 3088 | 96.3% | 85.8% |
| E4c | fixed | 8 | 3600 | 402 | 3198 | 3123 | 3123/3198 = 97.7% (97.1-98.2) | 86.8% | 413 | 12.9% | 11.5% |
| E4c | fixed | 32 | 3600 | 402 | 3198 | 3101 | 3101/3198 = 97.0% (96.3-97.5) | 86.1% | 1240 | 38.8% | 34.4% |
| E4c | fixed | 128 | 3600 | 409 | 3191 | 3070 | 3070/3191 = 96.2% (95.5-96.8) | 85.3% | 2212 | 69.3% | 61.4% |
| E4c | orig | 8 | 3600 | 686 | 2914 | 2863 | 2863/2914 = 98.2% (97.7-98.7) | 79.5% | 326 | 11.2% | 9.1% |
| E4c | orig | 32 | 3600 | 404 | 3196 | 3091 | 3091/3196 = 96.7% (96.0-97.3) | 85.9% | 1277 | 40.0% | 35.5% |
| E4c | orig | 128 | 3600 | 397 | 3203 | 3085 | 3085/3203 = 96.3% (95.6-96.9) | 85.7% | 2240 | 69.9% | 62.2% |
| E4r | fixed | 8 | 3600 | 400 | 3200 | 3082 | 3082/3200 = 96.3% (95.6-96.9) | 85.6% | 563 | 17.6% | 15.6% |
| E4r | fixed | 32 | 3600 | 400 | 3200 | 3061 | 3061/3200 = 95.7% (94.9-96.3) | 85.0% | 1432 | 44.8% | 39.8% |
| E4r | fixed | 128 | 3600 | 404 | 3196 | 3043 | 3043/3196 = 95.2% (94.4-95.9) | 84.5% | 2347 | 73.4% | 65.2% |
| E4r | orig | 8 | 3600 | 685 | 2915 | 2826 | 2826/2915 = 96.9% (96.3-97.5) | 78.5% | 450 | 15.4% | 12.5% |
| E4r | orig | 32 | 3600 | 403 | 3197 | 3035 | 3035/3197 = 94.9% (94.1-95.7) | 84.3% | 1496 | 46.8% | 41.6% |
| E4r | orig | 128 | 3600 | 398 | 3202 | 3035 | 3035/3202 = 94.8% (94.0-95.5) | 84.3% | 2365 | 73.9% | 65.7% |

**Holdout G5 denominator (`r3_11_G5_denominators.csv`, `r3_11_E7_solver_accounting.csv`).** G5 is a false-positive rate on the 6 control cells at 32 walkers per node, 200 datasets each (n = 1200 per arm):
FP = datasets with a declared change of either sign (margin 0.25) divided by **answered** datasets. E1q, E4c, E4r abstain on 202 of the 1200 (forced J < 6 in torus3 L = 18 and similar), so their denominator is 998; E7 and E8 never abstain, so theirs is 1200.
Results: E1q 3/998 = 0.30% (CP 0.06-0.88), E4c 0/998, E4r 0/998 (CP upper 0.37%), E7 0/1200, E8 0/1200 (CP upper 0.31%). Convention for solver failures: an E7 dataset with no band (4 of 1200 on controls at 32 walkers, 13 at 128, 1 at 8) is ANSWERED and counted not-declared in
the G5 denominator and not-covered in the coverage denominators (fixed denominators, as in 10b). That is the conservative direction for coverage (G4: E7 1196/1200 = 99.67% instead of 1196/1196) and the favourable direction for FP;
dropping them from the FP denominator changes E7 from 0/1200 to 0/1196 (CP upper 0.31%, unchanged to two figures). The all-n reading for the three abstaining arms gives FP 3/1200 = 0.25% (E1q); the answered-only reading used is the larger.
Coverage in the same table: E1q 926/998 (92.8%), E4c 993/998 (99.5%), E4r 990/998 (99.2%) conditional; strict 77.2%, 82.8%, 82.5%; E7, E8 1196/1200 under both readings.

## 4. E4 oracle-arm diagnosis (EXPLORATORY; spec `E4_ORACLE_DIAGNOSIS_SPEC.md`, committed before any run: 74d7535)

New compute: `e4_oracle_diag.py` re-executes the frozen E4o branch with M fixed at the exact window supremum (the arm's own M) and keeps b_j (exact smoothing bias), B_j (the allowance), n_j (stochastic part), sd_j and c for every dataset.
Replay check, stated in advance: on the original seeds the rerun reproduces the stored raw_amend3 E4o rows (abstain, J, hw within 1e-6, cov_simult) in all four ring cells: 200/200 abstain flags, 127 + 200 + 200 + 200 answered rows, every J, hw and coverage flag equal.
Every miss is attributed (bound / calibration); the replicate arm used 4000 fresh seeds per cell (seed block 13000 + id) for coverage and noise-scale diagnostics.

| cell | arm | answered | coverage (CP) | misses | bound | calibration | max (b_j / B_j) | P(max_j abs n_j/sd_j <= c) | median sd_boot / sd_emp |
|---|---|---|---|---|---|---|---|---|---|
| ring 128 | replay | 127 (73 abstain) | 89.8% (83.1-94.4) | 13 | 0 | 13 | 0.50 | 89.0% | 0.91 |
| ring 256 | replay | 200 | 91.5% (86.7-95.0) | 17 | 0 | 17 | 0.48 | 91.0% | 0.89 |
| ring 512 | replay | 200 | 88.5% (83.2-92.6) | 23 | 0 | 23 | 0.35 | 83.5% | 0.85 |
| ring 1024 | replay | 200 | 95.5% (91.6-97.9) | 9 | 0 | 9 | 0.35 | 85.0% | 0.83 |
| ring 128 | fresh x4000 | 2398 (1602 abstain) | 92.3% (91.2-93.4) | 184 | 0 | 184 | 0.50 | 92.2% | 0.95 |
| ring 256 | fresh | 4000 | 92.2% (91.4-93.0) | 311 | 0 | 311 | 0.49 | 91.1% | 0.92 |
| ring 512 | fresh | 4000 | 92.0% (91.1-92.8) | 320 | 0 | 320 | 0.36 | 87.8% | 0.87 |
| ring 1024 | fresh | 4000 | 92.6% (91.8-93.4) | 295 | 0 | 295 | 0.36 | 84.1% | 0.85 |

Reading (exploratory, no gate): (1) **no miss is a bound miss.** In the 14,398 fresh and 727 replay answered datasets the exact smoothing bias never exceeded the allowance; the largest ratio is 0.50, so the worst-case allowance is
never the cause and is on average generous. All 62 original misses (13 + 17 + 23 + 9) are calibration misses: the stochastic part exceeded its bootstrap critical band. (2) The calibration shortfall is general, not specific to ring 512:
with 4000 datasets per cell the oracle arm covers 92.0-92.6% in every ring cell, below the nominal 95%; the 88.5% of ring 512 in the 200 original datasets (CP 83.2-92.6) is compatible with the 92.0% of the 4000-dataset run (it lies inside that interval), so the one failing cell is consistent with a true level near 92% plus sampling noise. (3) The measured cause is the scale of the bootstrap standard error: the median bootstrap sd is
83-95% of the empirical sd of the same stochastic part across datasets (minimum 0.77), and the stochastic part alone is covered only 84-92% of the time at the bootstrap critical value. Centring is not the problem: the mean stochastic part is 0.03-0.13 empirical sd
from zero. A plausible mechanism, **not tested here**: weights are optimised against a covariance estimated from 200 bootstrap replicates of half B, so the optimiser exploits noise in that estimate and the realised variance exceeds the
estimated one. What a pre-registered follow-up would test: coverage of E4o with the weights computed from the true covariance (estimated from the 4000 replicates) against the bootstrap covariance, and with 1000 rather than 200 bootstrap replicates.
Caveat: the fitted-M arms E4c/E4r are unaffected by this diagnosis (their allowance is violated in 17-53% of datasets, 10c); the E4o result explains only why the benchmark itself misses the 95% target.

## 5. Native-target numbers for E0 and E2 (`results_round3/r3_native_E0_E2_*.csv`; script `round3_native_targets.py`)

**Targets, stated exactly.** E2 estimates, on each grid pair (t_k, t_{k+2}), d_hat_k = -2 ln[q_hat(t_{k+2})/q_hat(t_k)] / ln(t_{k+2}/t_k); its exact counterpart D2_k, with q for q_hat, equals the
average of d over the log-time interval [ln t_k, ln t_{k+2}] with respect to d(ln t), because d = -2 d ln q / d ln t. The manuscript scored E2 against d at the geometric midpoint. E0 estimates two weighted least-squares slopes of ln q_hat on ln t, over the lower and upper half of the window
(weights = return counts); its exact counterpart D0 is the same weighted slope of ln q (the stored `cov_pop` uses each dataset's observed counts; the descriptive curves below use expected counts 32 N q, which differ negligibly). The manuscript scored E0 against d at exp(mean ln t) of each half.

Counts at the tested sizes bracketing the reported thresholds (E0: 2D L = 192, 3D L = 96; E2: 2D L = 192, 3D L = 128). Power = declared drops / 200 (abstention = miss), CP 95%. "native min change" is the most negative
change of the exact native target across the window relative to its first value, "pointwise" the same for the pointwise d the original coverage used; both exceed the 0.25 margin once the window is long enough.

| est | family | L | thr | answered | cov_pointwise | cov_native | power | native_first_last | native_min_change | pointwise_min_change |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 | aniso2_w16 | 96 | 192 | 200/200 | 183/200 = 91.5% | 183/200 = 91.5% (86.7-95.0) | 10/200 = 5.0% (2.4-9.0) | 2.00 -> 1.88 | -0.12 | -0.11 |
| E0 | aniso2_w16 | 128 | 192 | 200/200 | 194/200 = 97.0% | 194/200 = 97.0% (93.6-98.9) | 50/200 = 25.0% (19.2-31.6) | 2.00 -> 1.54 | -0.46 | -0.49 |
| E0 | aniso2_w16 | 192 | 192 | 200/200 | 178/200 = 89.0% | 184/200 = 92.0% (87.3-95.4) | 198/200 = 99.0% (96.4-99.9) | 1.99 -> 1.21 | -0.78 | -0.89 |
| E0 | aniso2_w16 | 256 | 192 | 200/200 | 178/200 = 89.0% | 188/200 = 94.0% (89.8-96.9) | 200/200 = 100.0% (98.2-100.0) | 1.97 -> 1.06 | -0.92 | -0.99 |
| E0 | aniso3_w8 | 48 | 96 | 0/200 | 0/0 | - | 0/200 = 0.0% (0.0-1.8) | - | - | - |
| E0 | aniso3_w8 | 64 | 96 | 26/200 | 21/26 = 80.8% | 21/26 = 80.8% (60.6-93.4) | 25/200 = 12.5% (8.3-17.9) | - | - | - |
| E0 | aniso3_w8 | 96 | 96 | 200/200 | 186/200 = 93.0% | 190/200 = 95.0% (91.0-97.6) | 200/200 = 100.0% (98.2-100.0) | 2.67 -> 1.10 | -1.57 | -1.62 |
| E0 | aniso3_w8 | 128 | 96 | 200/200 | 184/200 = 92.0% | 187/200 = 93.5% (89.1-96.5) | 200/200 = 100.0% (98.2-100.0) | 2.47 -> 1.01 | -1.46 | -1.44 |
| E0 | aniso3_w8 | 192 | 96 | 200/200 | 41/200 = 20.5% | 187/200 = 93.5% (89.1-96.5) | 200/200 = 100.0% (98.2-100.0) | 2.28 -> 1.00 | -1.28 | -1.15 |
| E2 | aniso2_w16 | 96 | 192 | 1/200 | 0/1 = 0.0% | 0/1 = 0.0% (0.0-97.5) | 1/200 = 0.5% (0.0-2.8) | 2.00 -> 1.36 | -0.64 | -0.65 |
| E2 | aniso2_w16 | 128 | 192 | 148/200 | 141/148 = 95.3% | 144/148 = 97.3% (93.2-99.3) | 148/200 = 74.0% (67.3-79.9) | 2.00 -> 1.36 | -0.64 | -0.65 |
| E2 | aniso2_w16 | 192 | 192 | 200/200 | 179/200 = 89.5% | 187/200 = 93.5% (89.1-96.5) | 200/200 = 100.0% (98.2-100.0) | 2.00 -> 1.03 | -0.97 | -0.98 |
| E2 | aniso2_w16 | 256 | 192 | 200/200 | 181/200 = 90.5% | 184/200 = 92.0% (87.3-95.4) | 200/200 = 100.0% (98.2-100.0) | 2.00 -> 1.00 | -1.00 | -1.00 |
| E2 | aniso3_w8 | 64 | 128 | 0/200 | 0/0 | - | 0/200 = 0.0% (0.0-1.8) | - | - | - |
| E2 | aniso3_w8 | 96 | 128 | 48/200 | 34/48 = 70.8% | 39/48 = 81.2% (67.4-91.1) | 48/200 = 24.0% (18.3-30.5) | 2.81 -> 1.02 | -1.80 | -1.84 |
| E2 | aniso3_w8 | 128 | 128 | 200/200 | 170/200 = 85.0% | 186/200 = 93.0% (88.5-96.1) | 200/200 = 100.0% (98.2-100.0) | 2.81 -> 1.00 | -1.81 | -1.85 |
| E2 | aniso3_w8 | 192 | 128 | 200/200 | 164/200 = 82.0% | 184/200 = 92.0% (87.3-95.4) | 200/200 = 100.0% (98.2-100.0) | 2.81 -> 1.00 | -1.81 | -1.85 |

Reading: coverage on the native target is at 91.5-97.3% where the estimators answer in 2D and 92-95% in 3D at the threshold sizes; on the pointwise target E0 is 89% in 2D at L = 192, 256 and **20.5% in 3D at L = 192**, and E2's pointwise
coverage falls to 82-85% in 3D at L = 128-192 whereas native coverage stays at 92-93%: the secant bias is what the pointwise score penalises. The power thresholds are sharp: E2 in 2D has power 74.0% (CP 67.3-79.9) at L = 128 and 100% (98.2-100) at 192; E0 in 2D 25.0% (19.2-31.6) at 128 and 99.0% (96.4-99.9) at 192; in 3D E0 12.5% (8.3-17.9) at 64
with 26 answered and 100% (98.2-100) at 96; E2 24.0% (18.3-30.5) at 96 with 48 answered and 100% at 128. At L = 192 (2D) E0's CP lower bound (96.4%) and E2's (98.2%) clear 0.8; so do all 3D threshold sizes.
At the sizes just below the thresholds the arms answer on few datasets (E0 3D L = 64: 26; E2 3D L = 96: 48; E2 2D L = 96: 1), so the threshold also marks where the estimator starts to answer at all (frozen J < 6 rule).

**Fig. 3 / Table I annotation.** CL's plotted L is the largest of three graphs (L, 0.75 L, 0.5 L). In both crossover families the node count is proportional to the long-axis length, so the three graphs have 1 + 0.75 + 0.5 = 2.25 times the nodes of the single graph of size L, hence 2.25 times the walkers at fixed walkers per node. Caption text: "CL: L is the largest member of a three-graph run (L, 0.75 L, 0.5 L; 2.25 times the walkers of one graph of size L)."

## 6. The archive (local, NOT published, no DOI)

`archive/`: `setup_env.sh` (Python 3.12, exact versions in `requirements.lock`), `regenerate_all.py` (one command), `data_manifest.json` (sha256 of the 72 per-dataset inputs: 54 in git, 18 in `archive/data/`, not in git, 54 MB),
`README.md`. Run: `archive/.venv/bin/python -s archive/regenerate_all.py`. It checks the environment and every input hash, copies code and inputs to a work directory (the committed tree is not touched), runs 19 analysis steps (replication, Amendments 1, 3-9, 10a, 10b/10c/11, 12, and the round-3 tables),
deletes stale outputs before each step, and compares every regenerated table with the committed one. Result of the final run: all 19 steps ran with exit 0, all 72 input hashes matched, and all 61 regenerated tables are byte-identical to the committed ones (61 IDENTICAL, 0 DIFFERS, 0 NOT-WRITTEN), run on the Mac (arm64, Python 3.12.13) in about 60 s. The gate verdicts it prints (G1-G5, 10b-C/X/F/M, 10c-A/B/CENTRAL, the Amendment 1-9 pass bars) are collected in `archive/out/gate_verdicts.txt` and `gate_verdicts_10b_10c_11.csv`.
Limits: it regenerates tables and verdicts from the stored per-dataset results; it does not re-simulate (item 1 shows why a re-simulation on another numpy gives different realisations). Not covered: causal-set and track-audit analyses. Publication and a DOI are the author's call.

## 7. What is not done or not claimed

* Exact LP verification covers 44 datasets, not 1030, and relaxation B only (relaxation A, the plan's literal one, and the E7j branch-and-bound bounds are not LP-by-LP verified; the margin analysis covers all 1030 E7j decisions).
* The crosswalk explains the count difference by cross-version realisations (demonstrated on 14 ids by rerun on host-T and by J concordance on all 1030); the exact numpy change responsible was not isolated.
* The E4 calibration mechanism (optimiser exploiting noise in the estimated covariance) is a hypothesis, not tested; the measured facts are the bound never failing and the sd ratio of 0.77-0.95.
* Item 5 gives native-target numbers on the original replication cells (wpn 32); the secant arms in later amendments (GK, M, Mb, E3, CL) have no native-target table here.
* Nothing here changes a verdict; items 1, 2, 3 and 5 are corrective reporting, item 4 is exploratory.

brick-criteria completion and product-acceptance footer: not applicable (research seat; no brick/gate state of the blueMonster canon was touched).
