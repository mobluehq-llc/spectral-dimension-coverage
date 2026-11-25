# Checks of the Gemini and ChatGPT final reviews of A v9 (2026-11-24). Registered in PREREG_REDTEAM_FINAL.md. Stored data only; nothing new simulated. Python: pinned venv numpy 2.5.5.

## 1. Spectral-variable domain (Gemini 2.2)
- **The paper's Eq. (1) was ambiguous, not wrong in the integral.** v9 main Sec. III.F said the spectrum of S lies in [0,1], so q(t) = integral_0^inf e^{-st} dnu(s): a non sequitur as written (it names the eigenvalue range and then integrates over s). SM S7 said the same. SM S3 already had s in [0, inf], a tail cell [1, inf], and the bound e^{-t_j}. Fixed in main.tex and supplement.tex: q(t) = integral_[0,1] lambda^t dmu, s = -ln lambda in [0, inf] (lambda = 0 is s = inf). The representation in lambda is now stated; the code stays in s.
- **Code:** `p22amend8.atoms_grid` = atom at s=0, n log-spaced atoms on [0.02/tmax, 1], atom at s=inf. So lambda in (0, e^-1) is NOT represented in the discretised programmes (only lambda = 0). The outer relaxation (SM S3) covers all s with the tail cell.
- **Size of the effect.** Mass at s > 1 changes q(t) by at most e^-t <= 3.8e-11 at t >= 24 and the numerator s e^{-st} by less. Measured on stored paired runs (`results/rt_monolp_p2.csv`, variants B0 vs B0x with 60 extra atoms on (1, 30], 2D crossover cells, 1300 datasets): 1280 datasets answered by both; largest change in mean band width 7.4e-5 (relative 5.1e-5); changes in simultaneous coverage 0, in the endpoint-test decision 0, in the joint-test decision 0. **No reported bound or decision changes.** 7 datasets that failed under B0 are answered under B0x.
- **A new, real wording error found:** the SM said a dataset gets no band "when no measure on the grid fits the box". Status codes in `raw_amend9_a.csv` (7400) and `raw_amend9_amis.csv` (2000): 66 non-zero, ALL status 4 (HiGHS numerical difficulty); none status 2 (infeasible). The 20 to 13 improvement from the extended grid is in solver failures, not infeasibility. Fixed in the SM MonoLP section and the structural-constraints paragraph. (ChatGPT's "infeasibility 20 to 13" reading inherits the error.)
- Gemini's claim that the method "breaks at t = 2, 4, 6" is correct as a scope limit of the frozen grid and is not tested here; the outer relaxation would cover it.

## 2. Denominators (Gemini 2)
No mismatch; five different runs, all "answered" counts, all seeds different. Table II sources (`paper/make_figs_v7.py`): Frozen replication = Amendment 9b tables; Fixed windows = Amendment 10c fixed mode; Holdout = Amendment 11 stored cell file. Table S2 (tab_r3_10c_pooled) = Amendment 10c, own-window mode (E4c own 32: 404 abstentions, 3196 answered of 3600). Tables S3 = Amendment 11 gates (1327 = answered of 2000 on the ten crossover cells at 8 walkers). Table S7 = Amendment 13(a), fresh seeds, own rules per arm (holdout 8 walkers: 1291 of 2000; replication 32: 3200/3202; replication 128: 3193/3196; the oracle and CV arms abstain differently). Gemini's "3193 at 32 walkers in Table II" is wrong: 3193 is the CV count at 32 in the stored replication (Table II) and the oracle count at 128 in S7; his three-way "3193/3196/3200" is three different runs. 1327 vs 1291: abstention rate 33.7% vs 35.5%; the SD of the difference of two such counts is 29.9, so 36 is 1.2 SD. Nothing was recomputed or changed except one explanatory sentence added to the S7 caption. Also: Gemini's "Table II 1310/1327" is the E4c holdout count in tab_r5_gates.

## 3. h = 0.7 floor (Gemini 1)
Recomputed from `archive/data/raw_amend10c.csv` (E1q, fixed mode, answered): 23.7%, 40.0%, 62.9% at 8, 32, 128 walkers per node, 3200 answered each. Already in the paper. Gemini's "ensures RBC cannot shrink" claim is the untested hypothesis of P1.
Fixed-window answered counts (E4c, E4r, fixed): E4c 3198, 3198, 3191; E4r 3200, 3200, 3196: matches Table II (3198/3198/3191).

## 4. Fair RBC on integer steps (P1)
NOT RUN on host-W. Stored walks do not keep every step: `sample_blocks` returns block counts at the 11 grid times only. Fresh walks needed: about 5e8 steps per dataset at 128 walkers per node, 100 datasets x 3 budgets, plus a new nprobust-style estimator. The host-W machine was at load 8-35. Plan and bar are in PREREG_REDTEAM_FINAL.md P1 for host-D after L. The v9 claims are scoped accordingly (RBC is "our local-cubic correction at the original CV bandwidth on 6 cells").

## 5. ChatGPT items
- (a) certification of the 37: registered (P2), not run. One of the 37 DROP datasets sits among the 3 datasets where the endpoint/joint comparator hit a solver failure (stored `e7_status`, `j_infeasible`); the other 36 had a completed joint test that declared no drop. The stored fact is added to the P3 table.
- (b) same-box comparison, stored version (same regenerated 98.125% box; 100 datasets each; secant = log-convex secant bound through the box with sampled times, in `rt_monolp`; widths = pointwise band mean width over the window):

| cell | secant >0.25 | pointwise endpoint | joint | global DROP | median pointwise width | median seconds: pointwise / joint / global |
|---|---|---|---|---|---|---|
| C768xC16 (2D) | 0 | 0 | 0 | 37 | 1.16 | 3.3 / 4.3 / 2.7 |
| C160xC8^2 (3D) | 0 | 27 | 87 | 100 | 1.76 | 1.3 / 3.0 / 1.3 |

Secant certificates fail even in 3D because the box is wide at t >= 4352 (median certified drop -0.81). The extension to a 95% CP box is P3, for host-D.
- (c) near-null controls: registered (P4), for host-D. The three present controls have exact drops 0.000, 0.004, 0.053.
- (d) wording, all applied in v9 draft text: "spectral dimension of the mean return curve" (with E[q_g d_g]/E[q_g]) replaces "ensemble mean"; RBC is "our local-cubic correction at the original CV bandwidth", six cells; "smoothing bias explains most of the loss"; the 37/100 is "rejected on both atom grids, to numerical tolerance, without a continuum certificate"; the secant sentence clarified; eligibility rate and worst cell were already beside the pooled BiasAware figures (Sec. III.C: 86.8/86.1/85.3% of all, holdout 65.5/77.4/78.5%, worst cell 170/200).
- (e) physical caveat added to the no-drop Limits: finite-graph, walk-specific result and the factor [1 - 1/(Nq)]. ChatGPT's 5/6 at Nq = 6 and 0.33 are algebra not re-derived here beyond the identity d = [1 - 1/(Nq)] d_excess, which follows from q - 1/N differentiated.
- Not done: ChatGPT's "paired grid-by-bandwidth" experiment is P1's design on a dense grid; SLQ baseline; hero figure and master table (presentation, not done).

## 6. Build
main 8 pages, supplement 10 pages (v9: 8 and 10). To stay in PROCESS 20 the following v9 text was cut: the abstract prior-art sentence (kept in Sec. I), the A13(a) predictions detail (kept as "three of four held"), the holdout/regression scope sentence, the literature-arms sentence, the minimum-size remark of S3, the A13(a) reproducibility sentence, the A13(c) bar and per-3D widths, and the "search not exhaustive" sentence. The scope sentence on the holdout and the single check in a second design was cut from the Conclusion limits (the random-graph result stays in its own subsection).
