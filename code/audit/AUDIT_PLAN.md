# P22 real-data audit: PRE-REGISTERED PLAN (written 2026-11-21, committed BEFORE any estimator was run on AnDi or CET data)
What was looked at before this plan: (i) file listings and label structure of the AnDi-2 benchmark (track/experiment names, model names, track lengths,
counts of tracks with a given number of label changes, the label-implied exponent excursion S_lab); NO tracks were passed through any estimator;
(ii) HadCET: first/last rows, row count, sha256; no statistics. (iii) one smoke run of the DFA code on a synthetic fGn record (H = 0.65) to check the code runs.
Estimator code is the p22 code UNCHANGED (E0 `p22lib`, E1q `p22amend`, E4c `p22honest`, E6 `p22shape`/`p22amend5`/`p22amend6`); `audit_lib.py` only adapts inputs.
Adaptation: pass y = -0.5 log F (F = MSD or DFA fluctuation), x = log scale; then d = -2 dy/dx = dlogF/dlogs = alpha (MSD exponent) or H (DFA slope).
Window weights (the role of return counts) = 1/bootstrap variance of the split-A curve. The equivalence margin DELTA is patched per call (default 0.25).

## Data (official sources only, downloaded into fresh directories, read with `python -I`, nothing executed)
* **A. AnDi-2 benchmark dataset**, Zenodo record 14281479, doi 10.5281/zenodo.14281479 (Munoz-Gil et al., CC-BY-4.0; the 2nd AnDi Challenge benchmark, cites
  Nat Commun 2025 doi 10.1038/s41467-025-61949-x), file `andi2_benchmark_dataset.zip`, 749,404,326 bytes, md5 877ae83269e54ed083481845be3c07b5 (verified equal to the Zenodo checksum).
  Only .csv/.txt members extracted (videos skipped) to `~/dev/.p22audit_data/andi2_txt` (outside the repo). `track_2` is byte-identical to `track_1` (checked on 3 files) so only `track_1` is used.
  **THIS IS SIMULATED BENCHMARK DATA WITH LABELS (fBm-type diffusion with state switches, traps, confinement, dimerization), NOT real measurements.**
  Per-frame labels (alpha, D, state) are in `trajs_fov_*.csv`. Tracks are at most 200 frames. (The 1st-AnDi challenge public data carry no labels, see Zenodo 23126182; the AnDi
  generator is code and is not executed here, so the 1st-challenge single-track task is not used.)
* **B. HadCET** daily mean central England temperature, Met Office Hadley Centre, https://www.metoffice.gov.uk/hadobs/hadcet/data/meantemp_daily_totals.txt,
  1,768,845 bytes, 93,096 daily values 1772-01-01 .. 2026-11-20, sha256 da61cbb432af14a1b9c9856d0c43231178f248dff34eab2876d2fb69f5885f45 (file as served 2026-11-21).

## Part A: validation on labelled tracks (AnDi-2)
**Unit of analysis = one track.** First T frames, T in {100, 150, 200}. Lags Delta_k = unique round(2^(k/2)) <= floor(T/12): T=100 {1,2,3,4,6,8} (J=6), T=150 {..,11} (J=7), T=200 {..,16} (J=8).
**Statistic.** Time-averaged MSD per coordinate, TAMSD_c(Delta) = mean_t (u_c(t+Delta)-u_c(t))^2, c = x, y. The estimators' two independent halves are the two coordinates
(A = x, B = y; independent for the single-state models by construction; this replaces the time/walker split of p22 and works also when states switch). All = pooled x+y (used by E0, which has no split).
Bootstrap (B = 200): origin blocks of length b = max lag (T=200: 16 -> 12 blocks per coordinate) resampled with replacement within each coordinate, sums and counts of squared displacements pooled by block.
Dependence between neighbouring blocks (windows overlap) is not removed, so the bootstrap is expected to understate variance at the longest lags; this is a stated weakness whose effect is measured by the coverage numbers.
alpha_hat(Delta) = dlog TAMSD/dlog Delta (E0: two halves of the lag range; E1q, E4c, E6 pointwise/bands as in p22). Truth: constant alpha for single-state tracks (E[TAMSD]=MSD for stationary-increment fBm, so no secant bias).
For switching tracks the label-implied expected local exponent d_lab(Delta) = dlogE[TAMSD]/dlogDelta with E[TAMSD](Delta)=sum_i n_i(Delta) D_i Delta^alpha_i / sum_i n_i(Delta), n_i = max(len_i - Delta, 0)
(within-segment pairs only; pairs that straddle a change are ignored: an approximation, flagged).
**Sets (fixed from labels only).**
* SS (single state, false-change and coverage): length >= T frames, (alpha, D) constant over the first T frames, from experiments with model `single_state` (exp_8) or `multi_state` (exp_1, exp_7). Capped at 400 tracks per T (seeded sample if more). At T=200 the labels give 252 + 26 = 278 tracks.
* SW (switching, power): exp_7 `multi_state` tracks with length >= T, 1 to 5 label changes in the first T frames, every segment >= 0.1 T... for T = 200: >= 20 frames (T=150: 15, T=100: 10), and label-implied S_lab = max_j |d_lab(Delta_j) - d_lab(Delta_1)| >= 0.30.
  At T=200 the labels give 29 tracks (exp_1 has none: it switches too fast), so power estimates will carry wide Clopper-Pearson intervals; this is reported, not hidden.
**Scores (pre-registered).** Per estimator and T: non-abstention rate (J<6 -> abstain; E0 needs J>=6 too); on SS: simultaneous coverage of true alpha at all lags (E0: coverage of the two segment values) among non-abstained, with CP 95%; false-change rate = P(any-sign declaration at DELTA=0.25) with CP;
median band half-width hw (in alpha units). **C** passes iff pooled SS coverage >= 93% and CP upper >= 93%; **F** passes iff false-change rate <= 5% and CP lower <= 5% (the p22 rules).
On SW: power = P(declaration of the correct sign at DELTA=0.25), abstention counted as a miss, where the correct sign is that of the largest |d_lab(Delta_j)-d_lab(Delta_1)|; also coverage of d_lab (E1q, E4c, E6) as a diagnostic.
**Minimum-size recipe.** `recipe_view` uses the frozen constants in `results/amend6_recipe_params.csv` (a, delta, T, Jmin per estimator). Per stratum (T, estimator) we compute the recipe's predicted hw_law(J) = a (J/7)^-delta and the verdict
"detect iff J >= Jmin and (S-0.25)/hw_law >= T" with S = median S_lab of the SW stratum; compare with measured power >= 0.8 (match / mismatch), and compare measured median hw with hw_law (ratio).
The recipe was fitted on J = 7..11 grid points of walker data; J = 6..8 here is at its edge and the data are a different kind (single TAMSD tracks), so any agreement is a test of transfer, not a fit. We also report the J (and hence track length T) the recipe says is needed for S = 0.5: T such that J(T) >= J_needed (J(T) from the lag rule above; T-needed reported as a lower extrapolation of the lag rule).
**Replication.** None (one pass, one seed stream SeedSequence([2201, track_uid, T])); the unit sizes give the sampling error through CP intervals.

## Part B: genuine real-data audit, HadCET
**Claim under test (stated before looking):** the daily central-England temperature series has a weather-to-climate persistence crossover near 10 days, i.e. the local DFA slope H(s) changes
between short scales (s of about 7-14 d) and long scales (s of about 40-900 d) by at least a margin, with the transition located between 5 and 20 days. Expected direction (literature, from memory, not re-verified): the slope FALLS with s
(steep for the smooth short-memory weather scale, near 0.6-0.7 at long scales; claimed size S_claim = 0.30).
**Preprocessing (fixed).** Deseasonalise on the TRAINING period 1772-1850: mean seasonal cycle = harmonic regression (4 harmonics of the day of year / 365.25) of the daily values; seasonal standard deviation = sqrt of a 2-harmonic regression of the squared training residuals; anomaly z = (T - mean_doy)/sd_doy applied to all years. Missing days: none expected (assert).
**Analysis periods.** P1 = 1851-1950 (primary, 100 y), P2 = 1951-2025 (out-of-period confirmation, 75 y; 2026 partial year excluded). Each period analysed separately and identically; nothing from P2 enters P1 choices. 
**DFA.** DFA-2 (quadratic detrending of the profile), profile built inside 5-year blocks (so boxes never cross blocks), boxes of s days with 50% overlap, scales s_k = round(7 * 2^(k/2)), k = 0..14 (7 ... 896 d, J = 15); F(s)^2 = mean box mean-square residual.
y = -0.5 log F (so d = H(s) = dlogF/dlogs). Halves: A = even 5-year blocks, B = odd 5-year blocks (inference on B, tuning on A; E0 uses all). Bootstrap: calendar-year clusters (a box belongs to the year of its first day), resampled with replacement within each half, B = 300. (Year blocks understate dependence at s > 1 y; stated.)
**Estimators.** E4c primary (data-driven M); E6 secondary (carries the logistic location m and step size; its misfit flag is reported; E6 is known to fail silently on two separated steps, p22 Amendment 5/6); E0 and E1q as comparators. 
**Margins.** Primary DELTA_H = 0.10 (H units); also 0.25 (the recipe/p22 calibration margin). 
**Verdict rule (fixed).** SUPPORTED iff, in BOTH P1 and P2, E4c declares a change of the expected (negative) sign at DELTA_H = 0.10 AND the E6 logistic-location confidence set (exp(m) in days, 99% chi2 RSS set on the half-A data) in P1 intersects 5-20 d.
NOT SUPPORTED iff (a) in both periods every E4c contrast interval lies inside (-0.10, +0.10) (equivalence: no change of that size), or (b) a change is declared in P1 and P2 but the E6 location set lies wholly outside 5-20 d in both. UNDECIDABLE otherwise.
**Recipe view.** For E4c and E6: J=15 (extrapolation beyond the fitted J 7-11, flagged), hw_law, measured hw, and J_needed for S_claim = 0.30 at margin 0.25 (and the smallest S the record can certify, S_min = 0.25 + T*hw_law); and whether the 100 y / 75 y records are "long enough" (predicted detect, and measured agreement).
**Pipeline calibration on synthetic records with the SAME layout** (100 y x 365 d, same code, 100 replicates each, SeedSequence([2202, rep, kind])): N0 = fGn H = 0.65 (no crossover; false-change rate and coverage of the population DFA-2 local slope, computed from 200 long fGn records, finite-size curvature included); N1 = AR(1) phi = 0.85 plus fGn H = 0.65 (known crossover; power of E4c/E6). Reported with CP.
**Figures (PNG).** A1 coverage/false-change by estimator and T; A2 power vs track length with recipe prediction; A3 measured half-width vs J with the recipe law; B1 DFA F(s) and local slope with E4c and E6 bands, P1 vs P2; B2 synthetic calibration. 
**Heavy runs** on host-T (`~/p22audit`, system python3, NCPU <= 16, nohup); plots made locally; nothing is killed that this session did not start.
## What this plan cannot do
AnDi tracks are short (<= 200 frames), so most estimators will abstain or be uninformative at the smallest T; the audit then measures how much, via the recipe. The AnDi labels describe fBm-type models with state switches; confinement, traps and dimerization effects on TAMSD are not label-derived, so those experiments are not used for truth. 
CET is one record; the literature claim is paraphrased from memory and its numbers are not verified here. A negative verdict on this record is not a verdict on the physics.

---
## AMENDMENT 1 (2026-11-21, written after the synthetic population curves were computed and BEFORE `cet_real.csv` or any AnDi output was opened)
Observed on the SYNTHETIC null only (200 fGn H = 0.65 records of the CET layout, population local DFA-2 slope on the 15 scales): 0.82, 0.74, 0.69, 0.67, 0.66, 0.65, ... i.e. the DFA-2 finite-size bias alone produces a
DECLINE of about 0.17 between s = 7 d and s ~ 50 d for a process with NO crossover. So a declared decline at the pre-registered margin 0.10 cannot, by itself, be read as a weather-to-climate crossover.
Change to the verdict rule (fixed now): the PRIMARY verdict uses margin DELTA_H = 0.25 (above the 0.17 null drift; also the p22/recipe calibration margin) in the rule exactly as worded above (both periods E4c-neg declaration AND E6 location set intersecting 5-20 d => SUPPORTED;
equivalence inside +-0.25 in both periods or declared-but-wrong-location => NOT SUPPORTED; else UNDECIDABLE). The margin-0.10 verdict is still computed and reported, labelled "confounded by DFA small-scale bias, cannot discriminate". Also reported: the observed contrast vs the N0 population contrast at each scale (a null-referenced reading).
For the synthetic kinds: the N0 false-change rate is a valid false-positive rate at margin 0.25 only; at margin 0.10 N0 contains a real drift of 0.17 so its declaration rate is a power figure for a small true drift. N1 weights: AR(1) phi = 0.85 (unit variance) + 0.5 x fGn(H = 0.65), fixed in code before any run.
