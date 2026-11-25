# ROUND 4 cheap checks: EXPLORATORY specification (written and committed BEFORE any of these checks is run)

Date: 2026-11-23. Status: **EXPLORATORY diagnostics.** They answer the cheap checks named by external referees Gemini and ChatGPT in
round 4 (`ideas/A_spectral_dimension_coverage/REVIEWS/*round4*`). They feed no gate, bar or verdict. No estimator, bar, seed or cell
definition is changed (the six frozen estimator files stay byte-identical to 1ba9ca6; the new scripts only call their functions).
Any fix suggested by a result would be a new pre-registered amendment. Results go to `ROUND4_CHECKS.md`; raw output to `round4_checks/`.
Environment: host-D, the pinned venv of REPRO_CHECK.md (`python -I`, numpy 2.5.5), at most 6 cores (OMP/OPENBLAS/MKL = 1 per process).

## 1. Gemini 7.1: E7 bootstrap RNG and the 3.5176 vs 3.5541 box-constant drift (`round4_e7_rng.py`)
Cell: the REPRO_CHECK 9b cell (ring L=128, wpn 32, control, id as in `repro_check.r9b`), dataset 0 (stored c_box 3.5541, n_boot 16; re-run 3.5176, n_boot 15).
Prediction written now (so it can fail): E7's bootstrap draws come only from `rngs[7]`, a child of `SeedSequence([22,id,2000+ds]).spawn(9)`, so the
bootstrap is deterministic and independent of other estimators and of global state; the drift comes from the point set `big = cnt >= 30`
(`cnt = U.mean(0)*Wt`, float) having one more point in the stored run, i.e. a floating-point near-tie at the integer count 30.
Procedure, all on dataset 0 unless stated:
 a. Static check: grep the E7 code path (`p22amend8.make_box`, `p22lib.boot_means`, `p22lib.sample_blocks`) for `np.random.` global use, `default_rng()` without seed, set/dict-order dependence.
 b. Run `RI.run_dataset` twice in one process and in two separate processes; log for the E7 call: the `rngs[7]` bit-generator state before, the first 5 bootstrap multinomial
    count vectors (the "index arrays": `rng.multinomial(n, 1/n)` rows 0-4), the first 5 studentised maxima T[0:5], c_box, n_boot, n_cp. Compare all.
 c. Call `make_box` in isolation with a freshly built `rngs[7]` (no other estimator run first) and compare with (b) (tests inheritance of mutated state).
    Repeat with `np.random.seed(0)` and `np.random.seed(1)` set beforehand (tests global-state dependence).
 d. Near-tie test: for every time point print `cnt - 30` (float rule used by the code) next to the exact integer total count `sum(C)` (from the raw counts C before division by m);
    list points with |cnt-30| < 1e-6 and whether the float rule and the integer rule (`C.sum(0) >= 30`) disagree. Recompute c_box and n_boot under both rules with the same Ps
    draws, and compare with the stored 3.5541/16 and the re-run 3.5176/15.
 e. If (d) confirms a float-rule tie: scan, with sampling only (no estimators), every dataset of the 9b b32 cells (and the 10b, 10c, 11 cells if time allows, at most 6 cores) and
    count datasets where the float rule and the integer rule disagree, per amendment; those are the only datasets whose E7/E7j/E0-window numbers can depend on the platform. The scan
    uses each runner's own seed formula and `sample_blocks` call.
Reported: yes/no seeded deterministically; the 5 draws and maxima; the explanation of 3.5176 vs 3.5541; exactly which reported numbers can be affected (by amendment, with counts).

## 2. Gemini 7.2 / 3: E4o in-sample covariance (Markowitz bias) (`round4_e4_markowitz.py`)
Cell: ring control L=512 (cell id 2), wpn 32, E4o (oracle M), datasets 0..199 with the replay seeds `SeedSequence([22,2,ds])` exactly as `e4_oracle_diag.py`.
The three variants differ ONLY in how S (which sets the weights and the sd) and the critical value are obtained:
 - V0 baseline: S from 200 bootstrap draws of half B (the frozen design; must reproduce the stored/diagnosis E4o rows, including coverage 177/200).
 - V1 B=1000: S and the max-t critical value from 1000 bootstrap draws of half B (own child rng, `SeedSequence([22,2,ds,1000])`; V0's stream is untouched).
 - V2 independent split: S (hence weights and sd) from 200 bootstrap draws of half A (independent of the half-B data that gives the estimate; own child rng
   `SeedSequence([22,2,ds,2000])`); the critical value from the same 200 half-B draws as V0, applied to those weights. Window J, M, estimate yB unchanged.
Extra (labelled as extra, because 200 datasets give coverage only to about +-3 points): the same three variants on fresh seeds `SeedSequence([22,13002,ds])`, ds 0..999.
Reported per variant: answered/abstained, coverage with CP, `sd ratio` as defined in the earlier diagnosis (median over j of median bootstrap sd_j / empirical sd of n_j across
datasets within a J stratum of at least 30 datasets for the 200 arm, at least 100 for the 1000 arm) AND the pooled ratio 1/std(n_j/sd_j) over all datasets and j, the stochastic-part
coverage P(max_j |n_j|/sd_j <= c), and the bound-violation count (must be 0 under M exact). Referee's prediction to test: sd ratio rises from ~0.85 to ~1.0 and coverage from ~92% to ~95%.

## 3. ChatGPT 7: fixed-window E4 abstentions (`round4_fixedwindow_abstain.py`)
Source: stored `raw_amend10c.csv` on host-D (sha256 asserted against `manifest_amend10c.json`), mode == "fixed", est in E4c, E4r, realised cells (cell >= 9000).
Tally every abstain==1 row by (family, cell, wpn, J_fixed < 6 vs >= 6, recorded `why`, any other recorded column that explains it). List every row with J_fixed >= 6 and abstain == 1
(e.g. the C768xC16 cell at 128 wpn with 195/200). Then regenerate bar 10c-B (E4c coverage >= 93% at each wpn, pooled and per family) three ways: (i) as published (answered only);
(ii) eligible failures (abstentions with J_fixed >= 6) counted as misses, structural abstentions (J_fixed < 6, forced by the frozen J<6 rule) excluded; (iii) all abstentions as misses.
Same point-estimate rule and CP intervals as `verdict_amend10_11.py`. Also the same for E4r and for E1q's fixed-window abstentions for completeness.

## 4. ChatGPT 3: E7j minimum size (`round4_e7j_minsize.py`, reads only stored files)
From stored `raw_amend10b.csv` (committed) and `results/amend9a_power.csv`/`raw_amend9_a.csv`: for each tested 3D (aniso3_w8) size state, under the 10b outer-B certified joint test and
under the 9a discretised E7j: n datasets, number declaring a drop, fraction, CP interval, and how those datasets were chosen (full 200, or the 10b selected subsets). State what the stored
results support for "minimum size", and whether smaller sizes failed under the same scoring (same-scoring means the 10b certified outer-B column on the same datasets; sizes with only a
selected subset are flagged as not a power estimate).

## Addendum A (2026-11-23, committed BEFORE the scan below was run; after the single-dataset step 1b-1d was run)
Step 1 (b-d) was run first. Result that changes step (e): the near-tie prediction FAILED (no point within 1e-4 of count 30; float and integer rules agree, n_boot 15 under both). The drift is not an RNG or
rounding effect: the walk counts at the last two time points differ (stored/numpy-2.4.6 vs numpy-2.5.5). Step (e) is therefore replaced by: **numpy cross-version walk scan.** For the 9b groups b8, b32, b128
(all cells, datasets 0..199, the runner's own seeds `SeedSequence([22,id,2000+ds])`, `sample_blocks` only, no estimator) save the raw block-count matrix C in the pinned numpy 2.5.5 venv and in an
otherwise identical venv with numpy 2.4.6 (all other pins equal; the 9b stored results came from host-T with unrecorded versions, and dataset 0 matches 2.4.6 exactly, which is the evidence for using it as proxy).
Report, per group: datasets with any C difference, first differing time index distribution, and how many differing datasets have the first difference inside the E0-style window rule's reach
(count >= 30 and N*p >= thr, initial run), since only those can change an estimator that uses the window (E1q, E4*, E6, E0); E7/E7j/E8 use every time point and are exposed whenever any point differs and
that point is in the grid. Not claimed: that host-T ran 2.4.6.

## Addendum B (2026-11-23, committed BEFORE the reduced scan was run)
The addendum-A scan over 200 datasets of groups b8, b32, b128 had not finished its first group after about 25 minutes (host-D was at load average about 25, shared with other jobs; 6 cores used by this
task). It was stopped with no results read. Reduced scope, fixed now: groups **b8 and b32, datasets 0..29 of every cell** (`--nd 30`), both numpy versions, same comparison. b128 is not scanned.
So the scan estimates the fraction of datasets whose walks differ across numpy versions from 30 datasets per cell (about 960 datasets per group in b32), not an exhaustive count.
