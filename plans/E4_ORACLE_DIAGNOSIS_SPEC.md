# E4 oracle-arm diagnosis: EXPLORATORY specification (written and committed BEFORE any diagnostic run)

Date: 2026-11-23. Status: **EXPLORATORY.** Answers external referee round 3, major item 3 ("finish the E4 diagnosis"). It
feeds no gate, bar or verdict, and no earlier verdict, bar, seed or code is changed. If the result suggests a fix to E4,
that fix would be a new pre-registered amendment, not an edit to this one.

## The question
E4o (oracle M) failed the cell-level coverage rule C on one control cell: ring L=512, 177/200 = 88.5% simultaneous coverage
(raw_amend3.csv; also ring L=128 89.8% and L=256 91.5% among answered datasets). With M equal to the exact
supremum of |g'''| over the window the worst-case smoothing-bias allowance is, by the Taylor remainder, a guaranteed bound
on the deterministic part of the error. So a miss should come from the stochastic part, yet the manuscript calls the
failure undiagnosed. We separate the two causes dataset by dataset.

## Definitions (per dataset, per window point j; g = ln q exact, u = ln t, beta = g')
E4o's weights a_j (row j of W) are those the frozen code returns (`p22honest.honest_weights`) with M = max(Mo_cum[J-1], 1e-3),
the exact window supremum already used by E4o. The 'rerun with M fixed at its exact value' is therefore the E4o arm itself,
re-executed with every intermediate quantity retained:
- b_j = sum_i a_ji g(u_i) - beta(u_j)            exact smoothing bias (no noise; exact g from the spectrum)
- B_j = (M/6) sum_i |a_ji| |u_i-u_j|^3             the allowance the estimator adds (as `Bj`)
- n_j = sum_i a_ji (ln Phat_B,i - g(u_i))         stochastic part (Phat_B the half-B sample mean, floored as in the code)
- sd_j = the estimator's bootstrap sd; c = its bootstrap max-t 0.95 critical value; hw_j = c sd_j + B_j
- err_j = b_j + n_j (identity; asserted to 1e-9 against the code's own `err`). Miss = any j with |err_j| > hw_j.

## Attribution of each miss (mutually exclusive, exhaustive; stated before the data)
1. **BOUND miss**: some j has |b_j| > B_j (+1e-12). The deterministic allowance was violated although M is the exact sup
   (would indicate a flaw in the bound as implemented: range of u, quadratic annihilation, or the weights).
2. **CALIBRATION miss**: |b_j| <= B_j for every j. Then |err_j| <= B_j + |n_j|, so a miss forces |n_j| > c sd_j for some j:
   the stochastic critical value was too small for the realised noise. Sub-split by cause, using the replicate
   population below: (a) *centring*: the cell mean of n_j/sd_emp,j is non-zero (Jensen effect of ln of a sample mean, or
   weight dependence on the data); (b) *scale*: the bootstrap sd_j or c is below the empirical spread of n_j across replicates.
A miss is never reported as 'both'; BOUND takes precedence (and the count of misses that also have a calibration excess is given).

## Compute design (new compute; EXPLORATORY labels on every output file)
- Cells: the four ring controls (ids 0-3: L=128, 256, 512, 1024), `wpn`=32 as in Amendment 3 (the primary object is L=512).
- **Replay arm:** seeds SeedSequence([22, id, ds]), ds = 0..199, exactly as run_p22c. Acceptance check, stated now: the
  rerun must reproduce the stored raw_amend3 E4o rows (abstain, J, cov_simult, hw within 1e-6) on this platform; any
  difference is reported and the attribution of the original 200 is then marked 'platform-differs'.
- **Replicate arm:** fresh seeds SeedSequence([22, 13000 + id, ds]), ds = 0..3999 per cell (so coverage is estimated to
  about +-0.5 point and the noise-scale diagnostics have usable power). Never mixed with the replay arm when computing coverage.
- Quantities per dataset retained: J, M, c, per-j (b_j, B_j, n_j, sd_j) and the miss flags above.
- Strata by J (the window length is data-dependent); diagnostics are computed within J strata with at least 100 datasets.

## Reported (no bars)
1. Coverage (answered datasets) and CP interval, per cell and arm; abstention count. Denominator = non-abstained datasets, plus the
   abstentions-as-misses figure beside it.
2. Misses: BOUND vs CALIBRATION counts per cell and arm; for the 200 replay datasets, the table per dataset id.
3. Bound check: fraction of datasets with any |b_j| > B_j, and max b_j/B_j.
4. Noise calibration: P(max_j |n_j|/sd_j <= c) (nominal 0.95; this is the coverage of the stochastic part alone);
   median sd_j/sd_emp,j; mean n_j/sd_emp,j per J stratum.
5. One sentence stating which cause the data support, and what a pre-registered follow-up would test.

## Not claimed
No claim about E4c/E4r (fitted M), whose bound violations are reported in Amendment 10c. No claim beyond ring controls at
wpn 32. Nothing here can change the verdict that E4o 'fails C' in Amendment 3; it can only explain it.
