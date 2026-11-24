# P22 - scale-dependent spectral dimension d_s(t) on finite graphs: estimator simulation (PRE-REGISTERED PLAN)
Written 2026-11-20 BEFORE any simulation run. Only the exact (spectral, noise-free) d(t) curves of a few tori were inspected, to choose
the time grid and the anisotropic widths (see "Design choices made from exact curves"). Nothing below is tuned on simulated data.
Any post-hoc change is listed under "Amendments" in RESULTS.md. Context: ~/dev/brainstorm_111426/RESEARCH_TRACKS_AND_BM_LESSONS_111926.md,
DIMFLOW_LITCHECK_111926.md (literature gap: no d_s estimator with a coverage table and a minimum-size criterion).

## Question
Estimate d_s(t) = -2 dlogP(t)/dlogt from random-walk return data on finite graphs with honest error bars; find the smallest graph on which a
real dimensional crossover is detectable at power >= 0.8 with false-positive rate <= 5%.

## Walk, estimand, truth
- Lazy walk S = l*I + (1-l)*P, P = D^-1 A, l = 1/2 (ablation: l = 0). Walkers start uniformly over nodes. P(t) := q(t) = tr(S^t)/N (return probability
  averaged over start nodes). Exact truth from the spectrum: q(t) = mean(lambda^t), d(t) = -2 t q'(t)/q(t), q' = mean(lambda^t ln lambda), continuous t.
  Tori / products: closed-form eigenvalues l + (1-l)(1/m) sum_a cos(2 pi k_a/n_a). Gasket and random regular: dense eigh of the graph.
- Time grid: union {24*2^j} U {34*2^j} (ratios 1.417/1.412, so 2*t_k = t_{k+2} exactly), up to a generation cap T_cap(graph) = last grid time with exact N q >= 4.5
  (also hard cap 4096 for the step-simulated graphs). t_lo = 24 (random regular: grid {2*2^j} U {3*2^j}, see below).
- True d(t) used for coverage is the EXACT POINTWISE DERIVATIVE at the estimator's own reference scale (E1: t_k; E2: sqrt(t_k * 2 t_k); E0/E3: segment centres).
  Smoothing / secant bias is therefore PART of the error being scored ("bandwidth leakage" is measured, not excluded). A second, diagnostic coverage vs the
  estimator's own population functional (same estimator on exact q) is reported to separate noise from bias.

## Graphs (cells: graph x size; >= 200 datasets per cell)
Controls (no crossover; scored for coverage and FPR): ring C_L, L in {128,256,512,1024}; torus C_L x C_L, L in {24,32,48,64,96}; 3D torus C_L^3, L in {16,20,24,28,32}.
Known crossover (scored for power): C_L x C_16 (2D->1D, exact curve crosses between t~136 and ~768), L in {32,48,64,96,128,192,256};
C_L x C_8 x C_8 (3D->1D, crossover t~34-100), L in {32,48,64,96,128,192}.
Fractal: Sierpinski gasket levels n in {4,5,6,7} (N=123,366,1095,3282); asymptotic d_s = 2 ln3/ln5 = 1.365. The finite-graph exact d(t) (log-periodic, boundary-affected) is
the truth; coverage only, not scored in the pass bar (its d(t) legitimately varies), plus a report of |d_hat - 1.365| in mid-window.
Random regular (degree 3, N in {512,2048}, 4 independent graphs per dataset => genuine graph level of the two-level bootstrap; truth = mean exact q over 24 reference graphs):
degree/expander stress control, coverage only (d(t) is not constant; grid {2,3}*2^j; window threshold N q >= 1.3), NOT scored.
Parity/laziness ablation (E1, E2 only, reported separately): torus L=48 and ring L=256 with l = 0 (non-lazy; bipartite parity zeros).
Deterministic graphs have a trivial graph level: the two-level bootstrap there is walkers only (stated, not hidden).

## Data per dataset
Walker budget W = 32 N walkers (E1/E0/E3) and P = 32 Rademacher probes (E2) per dataset. Matched in information: one probe ~ N walkers. Walker sampling is exact: tori
via independent multinomial/binomial increments between grid times (verified against step simulation in checks), gasket / random regular by step simulation.
Walkers are grouped into 200 equal blocks per graph; block counts of returns at grid times are the bootstrap unit (cluster bootstrap over iid walker blocks).
Seeds: numpy SeedSequence([22, cell_id, dataset_id]); every estimator sees the same walker blocks / probes (common random numbers).

## Common window rule (all estimators, data-driven, uses only known N)
Window = longest initial run of grid times with (return count >= 30 [E1: in the selection half]) and N*P_hat >= 6 (N*q_hat >= 6 for E2), except random regular (1.3).
Threshold 6 chosen because exact d(t) of tori deviates < 0.1 from its plateau while N q >= 6 (read from the exact curves). Fewer than 6 grid points => the estimator ABSTAINS
(counted separately, never counted as covered or as a detection; reported per cell). E2 additionally loses the last two grid points (needs t and 2t).

## Estimators
- E0 baseline: global straight-line log-log slope over the window (weighted by counts) => single d; plus its natural change test, two straight-line slopes on the lower and
  upper half of the window, contrast = d_up - d_low (single contrast, bootstrap-sd band, +-1.96 sd). CI by the same two-level block bootstrap (B=200).
- E1 (blueMonster hero): local-linear slope b(t) of log P_hat on x = log t, Gaussian kernel, weights = counts; d = -2 b. SPLIT SAMPLE: blocks alternate into halves A/B.
  On A only: window (rule above, counts >= 30 on A) and bandwidth h in {0.5,0.7,1.0,1.4,2.0} (units of log t) by leave-one-out weighted CV of the local-linear fit of log P_A; kernel weights
  from A counts. Inference only on B (log P_B through the A-chosen smoother). Two-level bootstrap (graphs x walker blocks), B=200, on B. 95% simultaneous max-t band: c = 0.95 quantile of
  max_j |D*_j - d_hat_j|/se_j over the window; pointwise band +-1.96 se. CHANGE = band for Delta_j = d(t_j) - d(t_0), j >= 1 (own max-t band) excludes the null region (below).
  Ablation E1x: identical but no split (window, bandwidth and inference on the same data) to quantify bandwidth/selection leakage.
- E2 runner-up: Hutchinson/Rademacher trace q_hat(t) = mean_p z_p' S^t z_p / N with exact sparse powers (symmetrised S), d_k = -2 log[q(2t_k)/q(t_k)]/log 2. Bands: probe bootstrap (resample
  the 32 probes, B=400, same max-t construction) and a Hoeffding variant (per-probe terms lie in [0,1]; union bound over grid; interval arithmetic for the ratio; guaranteed but expected to be vacuous).
  No split (no tuning parameter beyond the window rule).
- E3 published baseline (Reitz & Bianconi, arXiv:2006.10421, spirit): FAITHFUL SIMPLE VERSION - finite-size scaling collapse of N P(t) against t/N^{2/d} across three sizes of the same family
  (L, 0.75L, 0.5L of the long axis; gasket/random regular not applicable). d chosen by minimising the squared difference of log(N P_hat) between interpolated curves on overlapping u = log t - (2/d) log N
  (grid 0.5..3.5 step 0.1, refine +-0.1 step 0.02), curves taken to the equilibration regime (count >= 30, N P_hat >= 1.05). To obtain a scale-resolved output, the collapse is done separately for
  the early and late halves (in t, split at the geometric mean of 24 and the largest window end). Output d_early, d_late; CI from a block bootstrap (B=80, +-1.96 sd; simultaneous over the two scales
  via max-t); contrast = d_late - d_early. It estimates the SCALING dimension, so its truth is the thermodynamic value (plateau) and it is scored on controls only for coverage.
  This is not their full pipeline (no network-ensemble averaging, no fit of corrections to scaling).

## Change declaration and truth for change tests
Convention Delta = d(late) - d(early); a real flow has Delta < 0. Primary (pre-registered) CHANGE rule: band excludes the equivalence interval [-0.25, +0.25].
Reason: with these data volumes, tiny true nonzero Delta from lattice discreteness (|Delta| up to ~0.1) would be "detected" by any consistent estimator and make FPR meaningless.
Secondary, reported but not scored: band excludes 0 (the literal rule).
Power counts a declaration only if the sign is correct (band entirely below -0.25). A crossover is REAL at (graph, L) iff exact max over window grid of [d(t_0) - d(t)] >= 0.6 where the
window is the exact-truth analogue of the rule (grid times >= 24 up to exact N q >= 6). Cells where it is not real are reported as "not realised" (the detection rate there is a finite-size / short-window artefact rate, not power).

## Metrics
Per estimator x cell: abstention rate; pointwise coverage (fraction of (dataset, scale) pairs); simultaneous coverage (fraction of datasets with all scales covered); FPR on controls; power on real
crossover cells; each with Clopper-Pearson 95% interval; mean absolute error and mean signed bias (d_hat - d_true); mean band half-width. Bootstrap (cluster) CIs where ratios are used.

## PRE-REGISTERED PASS BAR
An estimator PASSES iff all hold:
 (C) simultaneous coverage pooled over non-abstained control datasets (ring, torus2, torus3) >= 93%, AND no control cell has a Clopper-Pearson 95% UPPER bound below 93%;
 (F) FPR (primary margin rule, above) pooled over control datasets <= 5%, AND no control cell has a Clopper-Pearson 95% LOWER bound above 5%;
 (M) a minimum detectable L is reported for each crossover family: smallest tested L (real-crossover cell) such that power >= 0.8 at that L and every larger tested L; reported with its CP 95% interval.
     Strict variant: CP lower bound >= 0.8. If no L qualifies, "not reached within the sweep" (<= max L tested).
Estimators failing C or F get no valid minimum size (their M is shown for information only, marked INVALID). E0's pass is assessed on its two-segment version.
Gasket and random-regular results are descriptive (not scored).

## Convergence / sanity checks (run before the grid; failures are reported)
(a) direct-increment walker sampler vs exact q and vs step-simulated walkers (3 graphs, |z|<4 per time); (b) E2 trace estimator unbiased vs exact q; (c) bootstrap B convergence: sd at B=200 vs 1000 within 10% on 20 datasets;
(d) first-100 vs last-100 dataset coverage agree within 2 SE per estimator (reported in RESULTS); (e) exact spectral truth vs sparse-matrix-power q for gasket n=5 and ring (1e-10); (f) gasket exact mid-window d(t) vs 1.365 versus level.
Runtime: target < ~2 h on 10 cores; if exceeded, CELLS are cut (not replicates), and the cut is stated.

## Design choices made from exact curves (before any simulation)
Base grid 24/34 (not 12/17): lazy 3D torus exact d overshoots to 3.19 at t=8-12 (lattice artefact), 3.02 by t=24. Width w=16 for 2D->1D (exact crossover t~136-768, task's example); w=8 for 3D->1D
(w=16 puts the 3D->1D crossover near t~100-400 but N=256 L is too large for 200 datasets x many L); threshold N q >= 6 from the 3D torus exact curve (d=3.01 at N q=6.9, 2.76 at 2.5).

---
## AMENDMENT 1 (2026-11-20, written BEFORE the amendment runs; the original pass bar above is UNCHANGED)
Motivation: in the first run E1 (hero) failed C at 91.9% and covered the exact d(t) on realised crossovers only 0-38% (smoothing bias, not noise). Amendment tests whether standard bias control repairs it,
and quantifies a gap the first report left open (E3 premature change declarations). Same pass bar: simultaneous coverage pooled over non-abstained controls >= 93% (and no control cell with CP upper bound < 93%), FPR (margin rule) <= 5%,
minimum-L reported. Nothing is tuned on amendment data; every variant below is fixed here.

**A1. E1b bias-corrected family** (all keep: split sample A/B, window + bandwidth chosen on A, inference on B, two-level block bootstrap B=200 that RE-RUNS the whole bias-corrected estimator on each resample, so the
correction's estimation variance is inside the band; max-t simultaneous band; identical change rule and margin 0.25). Slope operator = linear in log P_B, exact same Gaussian kernel and count weights.
 - **E1u** (undersmoothing): local-linear, h = 0.5 * h_cv where h_cv is the original leave-one-out CV choice on A (same grid H_GRID). Rate below MSE-optimal.
 - **E1q** (robust bias-corrected, CCT spirit): local-QUADRATIC fit, slope read from the linear coefficient; h = leave-one-out CV of the local-quadratic fit on A over grid (0.7,1.0,1.4,2.0,2.8,4.0). The quadratic term removes the leading O(h^2 d'') bias of the local-linear slope; the bootstrap supplies the RBC variance.
 - **E1qu** (RBC + undersmoothing): local-quadratic with h = 0.7 * h_cv(quadratic).
 - PRE-REGISTERED HERO = **E1q**; E1u and E1qu are reported as arms of the same family, not selected post hoc. A family-level statement ("some bias-correction arm passes") is reported as such, flagged as 3 looks.
 - Reported additionally: simultaneous coverage vs the EXACT d(t) at each estimator's own reference scale on the REALISED crossover cells (aniso2_w16 L>=128, aniso3_w8 L>=48), pooled and per cell; first run values E1 0-38% are the comparator.
 - Local-quadratic needs >= 6 window points (same abstention rule); no other change to window rule.
**A2. Premature-change rate.** Defined as the primary-rule declaration rate (band for Delta entirely beyond -0.25, i.e. correct-sign "flow" declaration; also any-sign rate) among NON-abstained datasets on (i) crossover cells with real == False
(exact crossover not realised in the observed window: aniso2_w16 L in {32,48,64,96}, aniso3_w8 L in {32}) and (ii) controls (= FPR). Computed for ALL estimators: E0, E1, E1x, E2 (and E2h) from the original raw.csv (no rerun needed; data are unchanged), plus E1u/E1q/E1qu/E3s from the new run.
 Note the exact d(t) in non-realised cells still drifts by < 0.6 inside the window; a declaration there is "premature" with respect to the pre-registered 'real' definition (max drop >= 0.6), and the rate is reported next to that cell's max_drop so a small genuine drift is visible.
 **E3s (single-size E3)**: one graph only; collapse across sizes is degenerate, so E3s = the early/late scaling-exponent estimate of E3 (split at the same geometric mean tsplit, E3's own equilibration-regime window rule count >= 30 and N P >= 1.05, B=80 block bootstrap,
 same contrast rule +-1.96 sd), implemented as a count-weighted log-log straight line per half. Coverage is scored against exact d(t) at the two segment centres (like E0) so it is comparable with E0-E2. This is E0 with E3's window and bootstrap size, stated honestly: with one size the
 finite-size-scaling collapse reduces to a power-law slope. The 3-size E3 results from the first run are re-reported unchanged for comparison.
**A3. E0b not run.** Condition was "if E0's band is pointwise". Checked in code: E0's band is the max-t quantile over its two segment slopes (`score()` with J=2) and its change band is the single contrast, so it is already simultaneous over the two segments. Only >2-segment structure would need more; not added.
**A4. Design of the amendment run.** Same cells needed: 14 control cells + 13 crossover cells (27 of 39), 200 datasets, same seeds (SeedSequence([22, cell, ds]); walker blocks bit-identical to the first run because rngs[1] is untouched; verified by reproducing the original E1 row). Gasket / RRG / parity ablation / E2 probes / E3 3-size are not rerun.
 Runtime target < ~2 h on 9 workers; if exceeded, cells (not replicates) are cut and stated.
**A5. Reading rule.** The amendment can only ADD estimators with their own verdicts; it cannot change the verdicts of E0-E3 from the first run. The final table is the union.

---
## AMENDMENT 2 (2026-11-20, written BEFORE the causal-set run; original pass bar and Amendment 1 UNCHANGED)
Purpose: apply the p22 estimators to a discrete-spacetime case where published spectral-dimension results disagree: causal sets. Only code smoke-tested (3 sprinklings at N=500/1000, plus N=500/2000 timings; no estimates inspected beyond "it runs").

**What is the published disagreement (read from the arXiv texts, 2026-11-20).** (i) Eichhorn & Mizera, arXiv:1311.2530 (CQG 31, 125007): random walk on the *Hasse diagram* (transitive reduction; links are edges; the walk can step to future and past neighbours, i.e. the
link graph made UNDIRECTED), p(x,x') = 1/n(x) with n(x) = number of links at x, return probability P(x,x,sigma) exact (not simulated), averaged over starting points and ~100 sprinklings; Poisson sprinklings in finite regions of Minkowski 1+1, 2+1, 3+1 (N up to 1e8 in 1+1). They report d_s that INCREASES at small scales
(d_s > d, "superdiffusion"), approaches d at intermediate scales with an IR cutoff on link length, and falls to 0 at large sigma from boundary effects; they attribute it to nonlocality (infinitely many nearest neighbours).
(ii) Belenchia, Benincasa, Marciano, Modesto, arXiv:1507.00330 (PRD 93, 044017): NOT a graph random walk. d_s from the heat kernel P(s) of the Euclideanised, regularised *continuum* non-local d'Alembertian (kernel of the causal-set d'Alembertian, all layers), giving reduction to d_s = 2 at short scales in every dimension. **The task text suggested a "nearest-neighbour proxy used by Belenchia et al."; I could not find a graph/nearest-neighbour proxy in that paper and do not implement one under their name.**
The disagreement is therefore between a graph-walk on the link structure and a continuum-operator heat kernel. What our graph experiment CAN test is whether enlarging the edge set from links toward the support of the non-local operator changes what the walk reports (below). It cannot reproduce (ii).

**Constructions implemented (all: lazy walk S = I/2 + P/2, P = D^-1 A, A = undirected adjacency, exactly p22's walk; EM use the non-lazy walk, so t_here ~ 2 sigma_EM and the parity wiggle EM mention is removed; deviation stated).**
 - **k0 = Hasse/link graph made undirected (Eichhorn-Mizera construction).** i~j iff i precedes j with NO element strictly between.
 - **k3 = "BDG-support" graph:** i~j iff causally related with <= 3 elements strictly between (the layers 0..3 on which the Benincasa-Dowker-Glaser d'Alembertian has non-zero coefficients in 2D and 4D, ignoring the signed coefficients - an UNWEIGHTED positive proxy of the operator's support, my construction, not a published one).
 - **k15:** same with <= 15 intermediates (a "non-locality dial", my construction).
 - **Full comparability graph** (all causal relations) is NOT analysed: mean degree ~N/4 (1+1) so N*q(2) is below any window threshold; only its mean degree / ordering fraction are recorded as a diagnostic.
 Largest connected component is kept (isolated pieces are recorded in the graph meta). Sprinkling regions: Alexandrov diamond of Minkowski 1+1 (unit square in light-cone coordinates u,v, relation = both coordinates increase, as in p14_lorentz_sprinkling.py) and 3+1 (double cone |x| < 1 - |t|, rejection sampling). Poisson sprinkling replaced by fixed-N uniform sampling (the Poisson N-fluctuation is irrelevant at this precision; stated). No IR cutoff on link length (EM's L) is implemented; boundary effects stay in.
 **Not implemented:** EM's cutoff, other manifolds, non-manifoldlike causets, any path-integral weighting.

**Controls (null for "flow").** For every causet graph (dim, construction, N, sprinkling) a random geometric graph (RGG) in the non-periodic unit cube of dimension D = 2 (for 1+1) / 4 (for 3+1) with the radius bisected so that its mean degree equals that sprinkling's graph (same N, boundary effects included, ~same degree, but local metric structure).
Also a regular-lattice reference (torus L^2 / m^4, N nearest) is NOT run (its exact curve has no sprinkling variability and adds nothing the RGG does not); stated. **No exact d_s(t) exists for any of these graphs** (random, finite, boundary-affected): coverage of a "true d(t)" is NOT scored; p22's cov_*/bias columns are meaningless here and not reported. The RGG controls are the
null: an estimator "reports flow" only to the extent it differs from what it reports on the control at the same N and degree. The honest reference for the control is its own sprinkled-ensemble estimate, not a known dimension (a finite RGG with high degree has d_s(t) of its own).

**Data and estimators.** Per (dim, N): R independent sprinklings; per sprinkling and construction the EXACT q_g(t) = tr(S^t)/N_eff from the dense symmetrised spectrum on the time set (union of the two grids below), then the ensemble of R curves is the dataset (units = sprinklings; the bootstrap resamples sprinklings; no walker noise). To feed p22's estimators unchanged, the pseudo walker budget is
set to 32 N per dataset (as in p22), so the p22 window rule (return count >= 30 and N*q_hat >= thr) is identical in form. Window threshold: **thr = 6 (p22 primary) and thr = 1.3 (p22's expander setting)**, both reported. Grids: **G2 = {2,3}*2^j** (p22's expander grid) and **G4 = {4,5,6,7}*2^j** (denser, E2 step 4), t <= 1024. A dataset with < 6 window points abstains (as in p22).
 Estimators (p22 code, unchanged except that E2 takes the grid step and E3's split uses t0): **E0** two-segment line, **E1q** (Amendment 1 hero: split sample, local-quadratic, CV bandwidth, B=200 max-t band), **E2** secant of q(2t)/q(t) from per-sprinkling values (Hutchinson role played by per-sprinkling exact q; Hoeffding variant not run), **E3** Millan/Reitz-Bianconi-style collapse over triplets of consecutive sizes of the ladder (N, N/2, N/4 for N>=2000, early and late halves, split at sqrt(2 * window end)), thr 1.05, B=80.
 Change convention as p22: Delta = d(late) - d(early); declaration if the max-t band excludes [-0.25, +0.25]; BOTH signs are reported (decl_neg = flow to lower d, decl_pos = rise; EM predict a small-scale rise, i.e. a negative Delta with this convention if d_s falls to d, positive if d_s rises with t - this is recorded, not assumed).
 Sizes **N in {500, 1000, 2000, 4000, 8000}** for dim 1 and 3; replicates **R = 100 (N <= 2000), 40 (N = 4000), 12 (N = 8000)** - the 8000 and 4000 cuts are because the dense eigendecomposition of 6 graphs per sprinkling is ~7 min per sprinkling at N=8000 (timed in the smoke); R is cut, N-range kept. 12 and 40 are below the ">= 100 where feasible" target and are flagged.
 Seeds: SeedSequence([22, 900+dim, N, rep]) for sprinkling and RGG, [22,777,...] / [22,778,...] for bootstraps. Runtime target < ~1.5 h on 9 workers.

**PRE-REGISTERED QUESTIONS**
 (a) For each estimator, construction, dim: does it report a monotone d(t) change (declared Delta outside +-0.25 with correct sign) and is the *same sign and approximately the same size* reported as N grows 500 -> 8000? Stability criterion: Delta band for N=8000 overlaps that for N=2000 (sign and band); report d_hat(t) curves per N. Also report the fraction of (dim, construction, N, grid, thr) cells where the estimator abstains - ABSTENTION IS A RESULT: these graphs have N*q(t) -> 1 within a few tens of lazy steps (measured on the smoke ensemble: N*q(24) ~ 1-2), so windows are short.
 (b) Does the conclusion (sign/size of Delta; d_hat(t) shape) depend on the construction (k0 vs k3 vs k15) beyond what the matched-degree RGG control does? A construction dependence "explains" the published disagreement only if it is (i) larger than the control's construction dependence and (ii) shows a monotone trend in k. The reading limit is stated in advance: since (ii) of the published pair is a different object (continuum operator heat kernel), at most this shows that the walk's answer is sensitive to what counts as a neighbour; it cannot show the continuum result is wrong or right.
 (c) Minimum N to detect a real change of the claimed size (2 -> 4 or 4 -> 2, |Delta| = 2), using p22 calibration. Operational rule fixed here: from p22 raw.csv / amend results take, for each estimator, the realised-crossover cells at the minimum detectable L (power >= 0.8): record the median window J (grid points) and the span log10(t_last / t_first) of the window there. For each causet construction the estimator is said to be "able to detect a p22-sized change" at the smallest N whose ensemble window has J >= J* and span >= span*; if never in the ladder, extrapolate J and span linearly in log2(N) and report the extrapolated N, FLAGGED as extrapolation. This is a necessary-condition (window sufficiency) calibration, NOT a power curve: p22 tested only drops (3->1, 2->1) at 32N walkers, with a smooth deterministic crossover, and no equivalent power experiment was run on causets (no injection of a known flow). Stated limit.
 Nothing is tuned on causet data; all choices above are fixed. Reading rule: this amendment adds a section; it changes no earlier verdict.

---
## AMENDMENT 3 (2026-11-20, written BEFORE the run; original pass bar, Amendment 1 and 2 UNCHANGED)
Purpose: test whether bias-aware "honest" fixed-length intervals (Armstrong & Kolesar) break the coverage-vs-power trade-off seen in Amendment 1 (E1u: crossover coverage 86% but min L not reached at 2D; E1q: passes the bar, crossover coverage 27%).
References (verified 2026-11-20 by search of publisher/arXiv/RePEc records): T. B. Armstrong & M. Kolesar, "Optimal inference in a class of regression models", Econometrica 86(2), 2018, 655-683 (doi 10.3982/ECTA14434; arXiv 1511.06028);
T. B. Armstrong & M. Kolesar, "Simple and honest confidence intervals in nonparametric regression", Quantitative Economics 11(1), 2020, 1-39 (doi 10.3982/QE1199; arXiv 1606.01200), whose Sec. 3.3 proposes a global-polynomial rule of thumb for M (after Fan & Gijbels 1996). I did not read the formula in their Sec. 3.3; E4r below is MY adaptation (global quartic, see below), not a claim about their exact rule.

**E4 estimator (all variants share everything except M).** Same data as E1: split sample (blocks alternate into halves A/B), window J chosen on A by the common rule, inference on B, same bootstrap draws as E1 (rngs[2] stream, two-level, B=200, on B).
 y_i = log p_hat_B(t_i), u_i = log t_i, g = log q (truth), beta_j = g'(u_j), d = -2 beta. Class: |g'''| <= M on the window range. For each grid point j the weights a^(j) (length J, linear in y) solve the convex programme
   min_a 1.96*sqrt(a' S a) + (M/6) sum_i |a_i| |u_i-u_j|^3   s.t. sum a_i = 0, sum a_i (u_i-u_j) = 1, sum a_i (u_i-u_j)^2 = 0,
 S = bootstrap covariance of log p* (J x J, from the B-half bootstrap). **Deviation from the task text:** the quadratic constraint is NOT optional here: with a third-derivative bound the cubic remainder bound (M/6)|u-u_j|^3 is only valid if the linear and quadratic Taylor terms are annihilated; without it the bias is unbounded in g''. A second-derivative-class variant is not run.
 Estimate beta_hat_j = a^(j)'y; sd_j = sqrt(a'Sa); worst-case bias B_j = (M/6) sum|a_i||u_i-u_j|^3. Pointwise interval beta_hat_j +- (1.96 sd_j + B_j). Simultaneous interval over the window grid: beta_hat_j +- (c sd_j + B_j), c = 0.95 quantile of max_j |a^(j)'(y*-y_hat)|/sd_j over the bootstrap (same max-type construction as score()); bias enters only through B_j (worst case), noise through the bootstrap. Intervals converted to d by x(-2). Solver: SLSQP on a = p - n, p,n >= 0 (cost c(p+n) so p*n=0 at optimum), started at the sd-minimal feasible a; non-convergence counted and reported; if constraint residual > 1e-6 after a retry the dataset abstains (counted).
 Change test: Delta_j = d_j - d_0 (j >= 1, d_0 = first window point, as in E1), contrast weights c_j = a^(j) - a^(0), sd from the bootstrap, simultaneous max-t critical value from the bootstrap, half-width cv*sd_cj + B_j + B_0. **Deviation from the task text:** B(c) = (M/6) sum |c_i||u_i-u0|^3 is not a valid bias bound for a contrast of weights centred at two different points; I use the triangle bound B_j + B_0 (valid, slightly conservative). Declaration rule: interval for Delta excludes [-0.25,+0.25] (primary, same as all estimators); correct-sign = entirely below -0.25.
**Choosing M (pre-registered, three variants, no other variant looked at):**
 - **E4o oracle**: M = max |g'''(u)| over the window range [log t_0, log t_{J-1}] computed from the exact spectrum (analytic: g_uuu = t g_t + 3 t^2 g_tt + t^3 g_ttt with g_t, g_tt, g_ttt the first three cumulants of ln(lambda) under weights lambda^t/sum). Upper benchmark only, not a feasible estimator.
 - **E4c held-out**: on half A only, local-QUARTIC fit of log p_hat_A (Gaussian kernel, count weights, bandwidth by leave-one-out CV over {1.0,1.4,2.0,2.8,4.0} in log t units), M = max over window grid points of |6 x cubic coefficient|; floor 1e-3.
 - **E4r rule of thumb**: on half A, global weighted quartic polynomial fit of log p_hat_A on u over the window (count weights), M = max over the window range of |third derivative of the fit|; floor 1e-3. (Fan-Gijbels/AK style; quartic because the third derivative of a quartic is the lowest-order non-degenerate case.)
 Caveat stated in advance: weights depend on S estimated from the same B half they are applied to (second-order; no further split).
**Pass bar for E4o/E4c/E4r (unchanged):** simultaneous coverage pooled over non-abstained controls >= 93% (no control cell with CP upper bound < 93%), FPR (margin 0.25 rule) <= 5% (no control cell with CP lower bound > 5%), minimum detectable L per crossover family. Additionally reported as in Amendment 1: simultaneous coverage of the EXACT d(t) at the grid points on REALISED crossover cells (aniso2_w16 L>=128, aniso3_w8 L>=48), premature-change rate on non-realised cells, mean half-width.
**Coverage-vs-power question (decided now):** an estimator "beats the trade-off" iff it passes C and F, has realised-crossover coverage >= 93% pooled, AND min detectable L <= that of E1q (256 / 96) in both families, i.e. it is not dominated by E1u/E1q. Otherwise the statement is "no".
**Cells/run:** the same 27 cells (14 controls + 13 crossovers, note == ""), 200 datasets, same seeds and bit-identical walker blocks; E1 is re-run in the check path to verify reproduction on a smoke subset. Runtime target < 2 h; if exceeded, cells (not replicates) are cut and stated. **Causal-set cells are NOT run**: they have no exact truth so coverage cannot be scored, and the ensemble-level datasets would not yield a calibrated comparison; stated as not done.
Reading rule: adds estimators with their own verdicts; changes no earlier verdict.

---
## AMENDMENT 4 (2026-11-21, written BEFORE the run; pass bar UNCHANGED: simultaneous coverage pooled over non-abstained controls >= 93% with no control cell CP-upper < 93%; FPR (margin 0.25) <= 5% with no cell CP-lower > 5%)
Added reporting: coverage of the EXACT d(t) at grid points on REALISED crossover cells (aniso2_w16 L>=128, aniso3_w8 L>=48) with Clopper-Pearson bounds; min detectable L (power >= 0.8, abstention = miss, CP-lower variant also reported) per family. Reference verified: L. Duembgen (2003), "Optimal confidence bands for shape-restricted curves", Bernoulli 9(3), 423-449, doi 10.3150/bj/1065444812. I use only its basic device (interval-average envelopes under monotonicity); I do NOT use its scale-dependent penalty, so no optimality claim is made.
Same 27 cells, same 200 datasets/seeds/walker blocks, same split (blocks alternate A/B), same window rule and abstention, B=200 two-level bootstrap on half B (rngs[2], same draws as E1). Runtime target < 2 h; cells (not replicates) cut if exceeded.

**E5 monotone multiscale band.** Slope estimates dhat_j (j<J window grid, on half B) from a linear local-polynomial operator chosen on A. Two arms fixed now: **E5qu (pre-registered hero)** = E1qu slope operator (local-quadratic, h = 0.7 h_cv(quad)); E5q = E1q operator (h_cv). Direction declared a priori: d NON-INCREASING in t for every graph family in the run (the hypothesis is 2D->1D / 3D->1D flow down; controls are tested under the same declaration, so a control whose exact d(t) is non-monotone, e.g. the lazy torus overshoot, can legitimately break coverage and that counts against E5).
 For every contiguous index interval I: m_I = mean of dhat over I, sd_I = bootstrap sd of the interval mean (uses the full bootstrap covariance, i.e. the correct dependence). c = 0.95 quantile over bootstrap draws of max_I |m*_I - m_I|/sd_I. For non-increasing f: U(x_j) = min over I with max(I) <= j of (m_I + c sd_I); L(x_j) = max over J with min(J) >= j of (m_J - c sd_J). (The task text's L/U assignment is for an increasing curve; swapped here for the declared direction.) Band covers d(t_j) simultaneously if f is non-increasing and the slope estimator is unbiased.
 Change rule: drop declared iff L(x_0) - U(x_{J-1}) > 0.25 (L, U are non-increasing so this is the maximum over pairs); rise never declared (excluded by the shape). Delta coverage: Delta_j=d_j-d_0 in [L_j-U_0, U_j-L_0].
 Stated limitation in advance: the band covers the smoothed (expected-estimator) curve, so slope-estimator bias is NOT removed; undersmoothed local-quadratic is chosen to limit it.
**Causal-set diagnostic (no truth, no coverage scoring):** E5qu on the same Hasse k0/k3/k15 and RGG datasets as Amendment 2 (G4 grid, thr 1.3 and 6, units = sprinklings, B=200). Direction is NOT a priori there, so BOTH directions are run (non-increasing and non-decreasing, by sign flip) and reported as two looks: whether the monotone band declares a drop, a rise, or neither.

**E6 sample-split logistic-step inference.** Model d(u) = d1 - H sigma((u-m)/w), u = log t, so g = log q = c0 - (d1/2)u + (H w/2) log(1+exp((u-m)/w)) (linear in c0,d1,H given m,w). Step 1 (half A, window J from A): grid m in 21 points on [u_0-1, u_{J-1}+1], w in {0.1,0.2,0.4,0.8,1.5,2.5}; GLS fit of (c0,d1,H) with covariance S_A = bootstrap covariance of log p_A (B=200, rngs[3]) shrunk 10% to its diagonal; profile RSS(m,w); confidence set C_A = {(m,w): RSS - RSS_min <= chi2_{2,0.99} = 9.21} (level 0.99). Step 2 (half B): for each (m,w) in C_A, the GLS estimate on B; (d1,H) = coefficients after partialling out c0 (orthogonalised matched contrast; sandwich covariance K S_B K' from the half-B bootstrap, no shrinkage), Wald region level 0.96 (chi2_2). Overall level >= 0.95 by union bound (0.01 + 0.04). Curve band at each j = union over C_A and over the Wald ellipse of d1 - H sigma_j (closed form per (m,w)). Change: Delta_j intervals the same way; drop declared iff the union upper bound of Delta_j is < -0.25 for some j (so every (m,w) in C_A must be compatible with a drop of more than 0.25; if C_A contains a step outside the window the band contains 0). E6 abstains with E1's window rule. Misspecification (logistic vs true crossover shape, lattice overshoot) is not corrected and is part of what is tested.
**Verdict rule (fixed now):** E5 or E6 "beats the trade-off" iff it passes C and F, has realised-crossover coverage >= 93% pooled, AND min L <= E1q's (256 for 2D->1D, 96 for 3D->1D) in both families. Else "no". The Amendment 3 table is extended with E5qu, E5q, E6 beside E0, E1q, E1u, E2, E3, E4c/E4r.

---
## AMENDMENT 5 (2026-11-21, written BEFORE the run; pass bar UNCHANGED: simultaneous coverage pooled over non-abstained controls >= 93% with no control cell CP-upper < 93%; FPR (margin 0.25, any-sign on cells whose exact curve has no real change) <= 5% with no cell CP-lower > 5%; minimum detectable L reported)
Only exact (noise-free) spectral d(t) curves of candidate graphs were inspected to choose cells (listed below); no simulated data. Candidate (b) "sharper single step" was examined on exact curves (ring x ring with w = 6..32, hypercube x ring cutoff products): every single-step product curve inside the p22 window (J = 6..14 grid points) is fitted by a logistic in log t to max abs error <= 0.04 (ring x ring <= 0.01), so none is a materially non-logistic single step; I did NOT build the hypercube product (cutoff steps are also logistic-fitted, max error 0.036, and need a custom sampler). The non-logistic cases are therefore: the two-step graph (exact logistic misfit 0.05-0.24), the brush (rise then end-effect, misfit 0.04-0.14, 0.34 at Lb16/Lf8), and the gasket (log-periodic wiggle). This is stated rather than hidden: the "sharper/asymmetric single step" request is covered only by the brush and by the existing ring x ring cells (which are asymmetric in log t: slow onset, super-exponential end).

**Part 1: misspecification of E6.**
Cells (all lazy walk l = 1/2, 200 datasets, W = 32N, same window rule/abstention, same split A/B, grid {24,34}*2^j):
 - (a) TWO-STEP tori C_L x C_8 x C_32, L in {128,192,256,320} (exact 3D -> 2D plateau at d = 2.00 for t ~ 150-400 -> 1D, window ends before the 1D plateau for L <= 256). Chosen from exact curves.
 - (b/iii) BRUSH of Grimaldi et al. (arXiv 2510.20520, read 2026-11-21 from the PDF text): Dirac brush = 2D periodic base (Lb x Lb nodes) with a 1D periodic ring fibre of Lf nodes through every base node (the base node is fibre site 0; base edges join fibre-site-0 nodes of adjacent base nodes; "periodic boundary conditions"); graph walk P = D^-1 A, lazy l = 1/2. **Assumption:** their paper's Fig. 3 / SM does not specify fibre attachment beyond the "fibre with a localised impurity of strength l^b" Laplacian equation (3); I take it as the ring-through-base-node graph above. Cells (Lb, Lf) in {(8,24),(12,16),(16,16)}. Exact curve = -2 dlogq/dlogt of tr(S^t)/N from the dense spectrum. NOTE this is p22's d(t), not their C(tau) = tau^2 Var(lambda) (different object, flagged in P22_COMPETITOR_READ). Exact d(t) RISES from 1.2-1.3 to 3.3-5.0 inside the window (no plateau at 2), so "real change" here is a RISE; drop declarations are the wrong sign.
 - (c) SIERPINSKI GASKET levels 5, 6, 7 (level 4 abstains: exact J = 1; recorded not run). Exact d(t) from the dense spectrum (window values 1.24-1.37, log-periodic wiggle of amplitude ~0.03-0.05). No real change expected.
 - (d) CONSTANT-DIMENSION controls: the 14 existing control cells (ring 128..1024, torus2 24..96, torus3 16..32), same seeds, rerun for the new arms; E6 reproduced against raw_amend4.csv as a check.
 - The 13 existing crossover cells (aniso2_w16, aniso3_w8) are rerun for the new arms (they are the minimum-L reference).
Truth labels from exact curves in the exact window (as in `prep`): `real_drop` iff max_drop >= 0.6, `real_rise` iff max_rise >= 0.6, else `null` (no real change). A "false change" is any band-excludes-+-0.25 declaration (either sign) on a null cell; wrong-sign declaration on a real cell is reported separately.
Estimators on every cell: E0 (two-segment), E1q, E4c, E6 (+ lack-of-fit below). Old cells: E0, E1q, E4c reused UNCHANGED from raw.csv / raw_amend.csv / raw_amend3.csv; new cells: computed with the same code, seeds SeedSequence([22, 100+k, ds]).
**E6 lack-of-fit (LoF) diagnostic, fixed now.** For each (m,w) in the half-A confidence set C_A: GLS fit on half B with the half-B bootstrap covariance S_B (same as E6), residual r_B = y_B - X theta_hat_B(m,w), Q(m,w) = r_B' S_B^-1 r_B. Statistic Q_min = min over C_A of Q(m,w) (best-fitting model compatible with the A half; if the truth is in C_A (prob >= 0.99) this is <= the Q at the truth, so conservative). Alarm iff Q_min > chi2_{J-3, 0.95} * (B-1)/(B-J) with B = 200 (the second factor is a first-order Hotelling-type correction for the bootstrap covariance, fixed beforehand). Nothing is tuned. Calibration (alarm rate on correctly specified constant controls) is reported, not gated.
Reported per cell and pooled by cell class: E6 simultaneous coverage of exact d(t) (CP), pointwise coverage, false-change rate (null cells), power (real cells, correct sign), LoF alarm rate, and the 2x2 (miss vs cover) x (alarm vs no alarm). **Loud/silent rule, fixed now:** among E6 datasets that MISS the exact curve (simultaneous coverage fails) pooled over the non-constant misspecified cells, E6 fails LOUDLY iff P(alarm | miss) >= 0.8, SILENTLY iff < 0.5, partially otherwise; if fewer than 30 misses are available the loud/silent question is reported as "too few failures to say" and the rates of alarm in the misspecified cells (two-step, brush, gasket) vs constant controls is the descriptive answer. Guarded E6 (E6g = E6 declaration suppressed when the alarm fires) is reported with its power and false-change rate. E0, E1q, E4c are compared on the same cells, same metrics (no LoF defined for them).
**Part 2: published comparison arms** (all use p22's common window, split at tsplit = sqrt(24 t_last) into an early and a late half exactly as E3/E3s so every arm yields d_early, d_late, the contrast, and the same declaration rule).
 - **M (Millan et al., arXiv 2006.10421, Sec. IV A and App. B, read 2026-11-21): ML power-law fit of the return probability, window optimised.** Their text: power-law fit of P0(t) by a maximum-likelihood algorithm adapted from Clauset-Shalizi-Newman (their ref. 73) that finds optimal t_min and t_max for an initial t_max; error as "d_s = 1.569 +- 0.004" (inverse-information type). The paper does not state the likelihood, the distance, or the error formula. **My assumptions, all fixed here:** (1) data = return counts c_k at the grid times, Poisson with mean mu_k = W A t_k^(-d/2), ML in (A,d) (profile in d solved by Newton); (2) per half, scan all contiguous sub-windows of >= 4 grid points (initial window = the half) and choose the one minimising the Kolmogorov-Smirnov distance max_k |F_emp(k) - F_mod(k)| between the normalised cumulative counts and the fitted cumulative means (the CSN criterion); the minimum 4 points is mine; a half with < 4 points => abstain; (3) error: sd(d) = 2 / sqrt( sum mu * Var_mu(log t) ) (Fisher information, treats counts at different times as independent and ignores the window choice, as the paper's +- does); band = +-1.96 sd per half, contrast = difference with sd = sqrt(sd1^2 + sd2^2); truth for coverage = exact d at the log-mean time of the SELECTED window. **M** = this analytic-error arm (the procedure as described). **Mb** = same estimate with block bootstrap (B = 80, window re-optimised on every resample, max-t band over the two halves, contrast +-1.96 sd of the bootstrap difference): not in the paper, added to separate the point estimator from its error bar.
 - **CL (true optimised finite-size collapse, Millan Eq. 11: P0^N(t) = (1/N) f(N t^(-d_s/2))).** Same three sizes (L, 0.75L, 0.5L of the long axis; brush: Lb scaled, Lf fixed) and curves (count >= 30, N P_hat >= 1.05) as E3. Differences from the simplified E3: (i) cost is a variance-weighted Houdayer-Hartmann-type statistic S(d) = sum over curve pairs and points of (y_i - interp y_j)^2 / (1/c_i + 1/c_j,interp) / n over the overlap (y = log N P, u = log t - (2/d) log N), (ii) minimised by continuous optimisation (coarse grid of 40 values on [0.4, 8], then bounded Brent in the best bracket), no 3.5 ceiling; (iii) B = 80 block bootstrap with the full re-optimisation. Early/late halves as E3. **Assumption:** the amplitude exponent is fixed at 1 (the ansatz's 1/N), so one exponent is fitted; I do not fit a second free exponent (would make d unidentified on short windows). CL coverage target = exact d at the log-mean time of the half; on constant controls additionally reported vs the plateau (E3's convention).
 - **Gaspar and Kukal (Physica A 602, 127624)** is paywalled (ScienceDirect 403, no arXiv) and is NOT implemented; abstract/slides only were seen elsewhere.
 - Arms M, Mb, CL are scored on the same pass bar (C, F, M) over the 14 controls and 13 crossover cells, plus premature-change rate (non-realised crossover cells and the new null cells), plus abstention rate, and the minimum L (power >= 0.8, abstention = miss) per family. Where an arm's window rule needs more points than the window provides, abstention is a result and is reported.
**Run design:** 14 + 13 old + 4 + 3 + 3 new cells = 37 cells (CL not run on gasket: no size family), 200 datasets per cell, same walker data per dataset for all estimators. Runtime target < 2 h on 9 workers; if exceeded, CELLS (not replicates) are cut and stated. **Reading rule:** adds estimators/cells with their own verdicts; changes no earlier verdict.

## AMENDMENT 6 (2026-11-21, written BEFORE the run and BEFORE any prediction is computed; pass bar UNCHANGED: simultaneous coverage pooled over non-abstained controls >= 93% with no control cell CP-upper < 93%; FPR (margin 0.25) <= 5%; minimum detectable L reported)
Only exact (noise-free) spectral curves were inspected to choose cells (2D aniso L = 384, 512, 768: exact window J = 14, 16, 18, max drop 1.004; twostep tori C_256 x C_8 x C_b: b = 12, 16, 20, 24 have window J = 13, max drop 2.04, a visible 2D plateau only from b >= 20). Heavy jobs run on host-T (code rsynced to ~/p22bench excluding .venv and results/causet*, system python3, NCPU = 16, nohup). Three items, each with a fixed pass bar. Nothing may be re-tuned after the data are seen; a failure is reported as a failure.

**Item 1: minimum-size recipe.** Rule: an estimator e detects a drop in a cell iff (i) the cell is FEASIBLE, J >= Jmin_e, and (ii) (S - 0.25) / hw_e(J) >= T_e, with S = exact max drop inside the exact window (`max_drop` in cells_meta), J = exact window length (number of grid points; it is set by walkers x volume: the window ends where N q falls below the threshold, walkers per node fixed at 32, so J is the carrier of "walkers and L"), hw_e(J) = the estimator's mean band half-width (for E4c/E4r this already contains noise + worst-case bias bound B_j; for the others bias sits in T_e), modelled as hw_e(J) = a_e (J/7)^(-delta_e). Predicted min L = smallest L of the family such that it and every larger L of the grid detect (grid for the first five: 32,48,64,96,128,192,256; for E4c/E4r: the same plus 384,512,768; "not reached" = one index past the last grid value).
 Constants fitted on 3D cells ONLY (aniso3_w8 L = 32..192 and torus3 L = 16..32, existing raw files; no 2D datum enters the fit):
 - Jmin_e = 1 + (largest exact J among 3D cells whose non-abstention rate for e is < 0.8).
 - (a_e, delta_e): OLS of ln(mean hw over non-abstained datasets) on ln(J/7) over the feasible aniso3_w8 cells (J >= Jmin_e). Fewer than 2 such cells => the estimator is reported "recipe not fitted".
 - T_e: over the feasible aniso3_w8 cells, ratio r = (S - 0.25)/hw_observed; cells with power >= 0.8 (abstention counted as a miss) are "detected", the rest "not detected". T_lo = max(1, max r over not-detected feasible cells), T_hi = min r over detected cells; T_e = sqrt(T_lo * T_hi) if T_lo < T_hi, else the arithmetic mean (flagged).
 Predictions: 2D (aniso2_w16) min L for E0, E1q, E6, Mb, CL; and for E4c and E4r over the extended grid. The first five are compared with the already-measured 2D values (E0 192, E1q 192, E6 256, Mb 256, CL 128, recomputed from the raw files by the same rule); this part is therefore a post-hoc consistency check of a pre-fitted rule, not out-of-sample. The E4c/E4r prediction is out-of-sample: predictions are committed (results/amend6_recipe_predictions.csv) BEFORE the new 2D cells L = 384, 512, 768 (aniso2_w16, 200 datasets, ids 400-402, E4c and E4r) are run.
 **Pass bar:** every prediction within ONE grid step (index distance <= 1 on the grid above) of the measured value. Reported separately: the five consistency predictions (5 of 5 needed), the E4c/E4r out-of-sample predictions (2 of 2 needed), plus the in-sample 3D fit. Also reported, ungated: predicted vs measured hw.

**Item 2: one attempt at a better E6 misfit detector.** E6 (with its LoF diagnostic) and E4c are re-run on the SAME datasets (original seeds SeedSequence([22, cell, ds]), ds 0..199, all 37 cells of Amendment 5; E6 streams rngs[2], rngs[3] and E4c stream rngs[4] as in run_p22e) because the raw CSVs kept only band summaries (dmin, dmax, hw), not per-point bands. Reproduction check before any claim: re-run E6 coverage/hw must equal raw_amend5 exactly on a smoke subset. E4c may differ from raw_amend3 on old cells (the old run used stream rngs[2]); stated, not hidden.
 - **Primary alarm A1:** E6 simultaneous band [cLo_j, cHi_j] on d at the J window points and the E4c simultaneous point band [d_hat_j -+ 2 hwj_j] (same window J, same split A/B) are DISJOINT at at least one grid point. Alarm requires both non-abstained; if E4c abstains and E6 does not, no alarm (count reported).
 - **Fallback A2 (used only if A1 fails the bar below; computed in the same run, reported either way, verdict follows the rule):** one-step vs two-step test on split A: Q1 = min over the E6 grid (21 m x 6 w, linear + one softplus step) of the GLS residual sum of squares with the shrunk half-A bootstrap covariance, Q2 = min over (Q1, all pairs m1 < m2, w1, w2 on the same grid, linear + two softplus steps), Lambda = Q1 - Q2 >= 0; alarm iff Lambda > chi2_{3,0.95} (B-1)/(B-J), B = 200. Nothing tuned.
 - **Pass bar (all four, fixed now):** (a) P(alarm | E6 misses exact curve) >= 0.5 pooled over the misspecified cells (twostep 4, brush 3, gasket 3); (b) P(alarm | miss) >= 0.8 on twostep_8_32 L = 320; (c) alarm rate on constant controls (14 cells, non-abstained) <= 5%; (d) realised-crossover power loss <= 5 percentage points: correct-sign declaration rate (abstention = miss) of E6 minus that of E6g (declaration suppressed on alarm), pooled over aniso2_w16 L >= 128 and aniso3_w8 L >= 48. If A1 passes: verdict "E6 guarded by E4c". Else if A2 passes: "E6 guarded by the one-vs-two-step test". Else E6 becomes a caveat. No third detector.

**Item 3: fresh-seed replication and the step-separation scan.** Fresh seeds SeedSequence([22, cell, 1000 + ds]), ds 0..199, for: the realised crossover cells (aniso2_w16 L = 128, 192, 256; aniso3_w8 L = 48, 64, 96, 128, 192), the four twostep_8_32 cells, the three brush cells; estimators E6 and E4c (same code path as Item 2; A1/A2 alarms are also reported on these seeds as a descriptive replication). New separation cells twostep C_256 x C_8 x C_b, b = 12, 16, 20, 24 (ids 300-303, seeds [22, id, ds]), E6 and E4c, with b = 32 (L = 256, original + fresh) completing the scan; separation ratio b/8 = 1.5, 2, 2.5, 3, 4. **Pass bar:** every replicated number lies inside the Clopper-Pearson 95% interval of the corresponding original (Amendment 5 / 3) number: E6 simultaneous coverage pooled over the realised crossover cells, each of the four twostep cells, pooled twostep, each brush cell, pooled brush, and E6 correct-sign power on the realised crossover cells (11 checks; all 11 must be inside; with 11 independent 95% intervals about 0.5 outside is expected by chance, which is stated rather than corrected). E4c checks are reported the same way, ungated. The 93% crossing of E6 coverage is reported as the bracket of b/8 between the last b with pooled coverage CP-lower >= 93% and the first with CP-upper < 93%.
**Reading rule:** adds results and verdicts; changes no earlier verdict. Cell/dataset counts are not cut; if runtime forces a cut it is stated.

## AMENDMENT 7 (2026-11-21, written BEFORE the run and BEFORE any prediction is computed; pass bar UNCHANGED: simultaneous coverage pooled over non-abstained controls >= 93% with no control cell CP-upper < 93%; FPR (margin 0.25) <= 5% with no cell CP-lower > 5%; minimum detectable L reported)
Purpose: test whether the minimum sizes and the Amendment 6 recipe are specific to the one sampling design (32 walkers per node), as the paper's Limitations section concedes. Only one thing changes: the walker budget. Window rule, thresholds (count >= 30, N q_hat >= 6), lazy walk, time grid, 200 blocks, B, margins and estimator code are unchanged.

**Observation from arithmetic only, recorded before any run:** the expected count at a grid point is wpn * N * q(t). With wpn = 32 the count rule (>= 30) binds only where N q < 0.94, so the N q >= 6 rule always ends the window first; the same holds at wpn = 8 (N q < 3.75) and wpn = 128. So the exact window length J does **not** depend on the walker budget for E0, E4c, E4r, E6 — RESULTS.md Amendment 6 calling J "the carrier of walkers x volume" is wrong for these estimators: J carries volume only. Exception: E1q's selection half (wpn/2 walkers per node) at wpn = 8 has count 4 N q, so its count rule binds at N q < 7.5 and its window is shorter. This means the frozen recipe (hw a function of J only) predicts the SAME minimum L at every walker budget (except E1q at wpn = 8). That is the falsifiable prediction tested here.

**Cells** (lazy l = 1/2, 200 datasets each, seeds SeedSequence([22, id, ds]), ids 500+ for wpn = 8, 600+ for wpn = 128, 700+ for wpn = 32 fill-in):
 - wpn in {8, 128}: controls ring L = 256, 1024; torus2 L = 24, 48, 96; torus3 L = 16, 24, 32 (8 cells); crossover aniso3_w8 L in {32, 48, 64, 80, 96, 112, 128, 160, 192, 256} and aniso2_w16 L in {128, 160, 192, 224, 256, 320, 384, 512, 768} (19 cells). 27 cells per budget.
 - wpn = 32 fill-in on the finer grid: aniso3_w8 L = 80, 112, 160, 256 and aniso2_w16 L = 160, 224, 320 (7 cells), so the baseline minimum L is re-stated on the same fine grid. Existing wpn = 32 cells are reused unchanged from earlier raw files.
 - Estimators on every cell: E0, E1q, E4c, E4r, E6 (with its A2 alarm, reported, not gated). Same walker data for all estimators within a dataset.
 - Fine grids: 3D {32,48,64,80,96,112,128,160,192,256}; 2D {128,160,192,224,256,320,384,512,768}. Steps are 1.14-1.5x, finer than Amendment 6, so "within one grid step" is a STRICTER bar than before (stated).

**Predictions, committed to results/amend7_predictions.csv before the run** by a script using only exact curves and the frozen Amendment 6 constants (amend6_recipe_params.csv; nothing refitted): (P1) the frozen J-only recipe's predicted min L per (estimator, family, wpn) on the fine grids; (P2) an alternative noise-scaled model, fixed now: hw_e(J, wpn) = a_e (J/7)^(-delta_e) (wpn/32)^(-1/2), same Jmin_e and T_e (pure walker noise scaling; for E4c/E4r the bias bound is not separately modelled — a stated simplification).

**Pass bars (fixed now; failures reported as failures):**
 - (A) Calibration robustness, per (estimator, wpn) for E4c, E4r, E6, E1q, E0: the original pass bar on the 8 controls of that budget. Verdict per pair; the headline is whether E4c passes at both budgets.
 - (B) Recipe transfer: P1 within one fine-grid step of the measured min L (power >= 0.8, abstention = miss; "not reached" = one index past the grid) for E4c, E6, E1q, E0 on both families at both budgets (16 predictions). PASS iff >= 14 of 16. P2 is scored the same way and reported; if P1 fails and P2 passes (>= 14 of 16), the paper says the recipe must carry the walker budget explicitly. If both fail, the recipe is reported as design-specific and not transferable.
 - (C) Ungated, reported: min L vs wpn per estimator; whether the honest-coverage-vs-power ordering of Amendment 3-5 (E4c costs more L than E1q/E6) holds at each budget; crossover coverage of exact d(t) per estimator and budget.
**Runtime:** host-T, NCPU = 16, nohup; target < 6 h. If exceeded, cells are cut largest-L first within each family and the cut is stated; replicates are never cut.
**Reading rule:** adds results and verdicts; changes no earlier verdict. Amendment 6's "J carries walkers x volume" sentence is corrected in RESULTS.md only after this run, with a pointer here.
**Amendment 7 erratum (2026-11-21, before any prediction or run):** the arithmetic note above is wrong for the split-sample estimators. E4c, E4r and E6 (like E1q) take the window from the selection half's counts (wpn/2 walkers per node; `est_all`: cntA = pA*WA). At wpn = 8 that count is 4 N q, so the count rule binds at N q < 7.5 for E1q, E4c, E4r and E6; only E0 (full counts) keeps the volume-only window. Consequence for P1: the frozen recipe predicts larger min L at wpn = 8 for the split-sample estimators (through a shorter J) and the same min L as wpn = 32 at wpn = 128. J for each (estimator, cell, wpn) is computed from the exact curve with the estimator's own count factor (E0: wpn; others: wpn/2). Nothing else changes.

## AMENDMENT 8 (2026-11-21, written BEFORE any run; earlier pass bars, verdicts and code are UNCHANGED)
Source: Fable's paper-A leads (`~/dev/brainstorm_111426/FABLE_PAPER_A_LEADS_112126.md`).

**Structure used (a theorem, not an assumption).** For the lazy walk (l = 1/2), the spectrum of S lies in [0, 1], so q(t) = ∫ e^{−st} dν(s) with ν ≥ 0, q(0) = 1, and an atom of mass 1/N at s = 0. Hence d(t) = −2t q′/q = 2t E_t[s] (the mean under the tilted measure). ln q is convex in t, so for t₋ < t < t₊:

  2t[ln q(t) − ln q(t₊)]/(t₊ − t) ≤ d(t) ≤ 2t[ln q(t₋) − ln q(t)]/(t − t₋).

(Checked by the coordinator from convexity: (ln q)′ lies between the adjacent secant slopes, and multiplying by −2t flips the inequalities.)

**Arms** (same walker data and seeds as the existing cells; wpn = 32; all 14 control cells and the 13 + 3 crossover cells of Amendments 1–6, plus the Amendment-5 misspecified cells for the alarm):
- **E7, moment-cone band.**
  - Inputs: return counts at every recorded grid time. A simultaneous 97.5% box for q(t_j): block-bootstrap max-t over the 200 walker blocks where counts ≥ 30; exact Clopper–Pearson with a Bonferroni split elsewhere. The two parts share the 2.5% equally.
  - For each window time t, the band is [min, max] of d(t) = −2t Σ w_k λ_k^t ln λ_k / Σ w_k λ_k^t over w ≥ 0, Σ w = 1, w(λ = 1) ≥ 1/N, with every q(t_j) inside its box. This is linear-fractional; solve it via Charnes–Cooper as two LPs per t (scipy HiGHS).
  - Atoms: 400 values log-spaced in s = −ln λ, plus λ = 1 and λ = 0.
  - Validity gate: doubling the atoms must move no band endpoint by more than 0.005. If it does, refine (800, 1600) and report; a band that never stabilises is reported as a failure of the arm.
  - E7 never abstains.
- **E8, secant bracket.**
  - At each window point t_j: lower = max over the next two grid times t₊ of the left secant expression; upper = min over the previous two grid times t₋ of the right one, using the estimated ln q.
  - Noise: the band is widened by c·se, with c the block-bootstrap max-t quantile (95%) over all secant statistics used.
  - Existing grid ratio is √2, so the noise-free bracket is [0.84d, 1.18d] (stated, not tuned).
- **A1′ alarm:** E6 band and E7 band are disjoint at any window time.

**Pass bars** (E7 and E8 separately, each the paper's bar):
- (C) pooled simultaneous coverage on non-abstained controls ≥ 93%, with no control cell CP-upper < 93%;
- (F) FPR (margin 0.25) ≤ 5%;
- (X) coverage of the exact d(t) on realised crossovers ≥ 90%;
- (M) min L at wpn = 32 at least one grid step below E4c in both families: ≤ 512 for 2D→1D and ≤ 160 for 3D→1D, on the Amendment-7 fine grids, with fill-in cells run as needed.
- **Verdict:** "passes" iff C, F and X all pass; "beats the trade-off" iff C, F, X and M all pass.
- **Alarm bar:** rate on controls ≤ 5%, and P(alarm | E6 misses) ≥ 50% on the misspecified cells (two-step, brush, gasket).

**Reported, ungated:** E7 and E8 at wpn = 8 and 128 on the Amendment-7 cells; band width against E4c; the LP solve time.

**Reading rule:** either outcome is reported. If E7 passes C/F/X but fails M, the paper says even the exact constraint set does not buy size. If it beats the trade-off, the title and claims are revised and the result is stated plainly. E7 does not transfer to the audit data (MSD/DFA are not completely monotone); say so.

## AMENDMENT 9 (2026-11-21, written BEFORE any run; all earlier verdicts UNCHANGED)
Three items. Seeds for fresh data: SeedSequence([22, id, 2000 + ds]).

**9a. Joint change test for E7 (closes the Amendment 8 caveat).**
- For each window time t_j: Δ_j^min / Δ_j^max = min / max of d(t_j) − d(t_0) over the same measure set as E7 (same box, atoms and gate). Solve with Charnes–Cooper as a fractional program: for fixed j it is the difference of two linear-fractional terms, so use a bisection on the value with LP feasibility, to tolerance 0.002.
- Declare a drop iff Δ_j^max < −0.25 for some j (a rise symmetrically).
- Same cells and data as Amendment 8 at wpn = 32.
- **Bars:** F ≤ 5% (no cell CP-lower > 5%) on controls; M as in Amendment 8 (≤ 512 in 2D, ≤ 160 in 3D).
- **Verdict:** "E7j beats the trade-off" iff F and M both pass (E7's C and X are unchanged by the test).

**9b. Frozen confirmatory replication.**
- Every published arm (E0, E1q, E4c, E4r, E6, E7, E8) with code frozen at commit 1ba9ca6, on fresh seeds.
- Cells: the 14 controls plus the realised crossover cells of Amendments 1–7, at wpn ∈ {8, 32, 128}, 200 datasets each.
- **Gated predictions** (all must hold):
  - (i) E1q realised-crossover coverage at wpn 128 < at wpn 32 < at wpn 8;
  - (ii) E4c, E4r, E7, E8 controls ≥ 93% pooled at each budget;
  - (iii) E1q realised-crossover coverage at wpn 128 ≤ 20%;
  - (iv) every C/F verdict of Amendments 1–8 reproduces at wpn 32.
- Report each check with CP intervals.

**9c. Second design.**
- Time grid ratio 2^(1/4) (grid {24·2^(j/4)}) and a stricter window: count ≥ 100 and N q ≥ 10, with J ≥ 8 (else abstain).
- Arms E0, E1q, E4c, E6, E7, E8 at wpn = 32 on the 14 controls and the aniso2/aniso3 crossover families (fine grids, Amendment 7).
- **Gate:** for each arm, the C/F/X verdicts equal those in the original design, AND the ordering of min L between E4c and each of E1q/E6/E7 is unchanged (ties allowed).

## AMENDMENT 10 (2026-11-22, written BEFORE any run). REVIEW RESPONSE to the external referee (external referee report, round 1; not distributed)
All three checks are the referee's own, with the referee's decision criteria adopted as the bars. They are tagged CONFIRMATORY with respect to the referee's stated outcomes. The paper's headline may only be narrowed by them.

**10a. Margin-consistent false-change scoring (referee check 1; analysis of existing data only)**
- For every non-realised crossover cell in Table I, compute the exact excursion Δ_j = d(t_j) − d(t_0) at each estimator's tested times from the exact curve.
- Re-score each declaration as follows:
  - **true change:** the sign matches and |exact excursion| > 0.25 at a time the declaration covers;
  - **false change:** declared, but the exact excursion is ≤ 0.25 in that direction over the window;
  - **no declaration**;
  - **abstain**.
- Cell by cell, with CP intervals, for E3, CL, E0 and E1q; document the mapping from each collapse output to a declaration.
- **Bar (referee):** the "premature change" headline survives only if false-change rates stay high (say ≥ 20%, CP-lower > 10%) on the near-flat cells (exact drop ≤ 0.02) with ≥ 30 non-abstained datasets. Otherwise the headline is narrowed or withdrawn.

**10b. Certified outer bounds for E7/E7j (referee check 2)**
- Use a rigorous discretisation-error bound. For any measure ν on s ∈ [0, s_max], replacing each mass by its nearest grid atom changes e^{−s t_j} by at most δ_j = max over cells of |e^{−s t_j} − e^{−s_k t_j}|. So widening every box constraint by δ_j·(total mass) gives an LP whose feasible set contains every continuous measure: an OUTER relaxation. Also widen the objective's d(t) to a certified enclosure, by bounding numerator and denominator separately with interval arithmetic (python-flint arb).
- Re-run E7 and E7j on all 170 atom-sensitive datasets, the withheld decision, a fixed seeded sample of 300 stable datasets, and exact-q versions of the 7 checked cells.
- **Bar (referee):** coverage, controls, F, and the 3D min L = 160 result are preserved under certified bounds. Otherwise report the widened results as the corrected ones.

**10c. Fixed-window budget comparison + E4 bound diagnosis (referee check 3)**
- On the two crossover families from 9b (realised cells), run E1q, E4c and E4r at wpn 8/32/128 with IDENTICAL deterministic windows: J fixed per cell as the wpn-32 exact window. The window is chosen from exact curves before any simulation and recorded in the run log before data generation. Fresh seeds [22, id, 3000 + ds]; 200 datasets; original-window results reported alongside.
- **E4 diagnosis:** for each fitted M, compare it with the exact required third-derivative supremum on the window (from the exact curve). Report the violation rate and coverage conditional on violation. Include the separation-scan cells at ratios 2.5 and 3.0.
- **Bar (referee):** the central result is confirmed if E1q coverage still falls with wpn under fixed windows (monotone, wpn-128 CP-upper < wpn-8 CP-lower) AND E4c stays ≥ 93% on matched windows. Bound violations are documented either way.

**Not run here (blocked or separate):**
- referee item 5 (Gašpar–Kukal full text: PDF not available);
- referee item 6 (temperature dependence: a separate audit amendment).

## AMENDMENT 11 (2026-11-22, written BEFORE any run). GENERALISATION HOLDOUT for the narrowed paper (Adam approved narrowing)
**Purpose:** answer the referees' objection that the 9b holdout only re-seeded the same graphs. This holdout uses graph families NEVER used in Amendments 1–10. It is run ONCE; a failure is reported as a failure, and no second holdout is run.

**New families** (lazy walk l = 1/2; the same window rule and time grid as the main benchmark; exact truth from the spectrum):
- Crossover 2D→1D with a new fibre: C_L × C_12, L ∈ {96, 128, 192, 256, 384}.
- Crossover 3D→1D with a new cross-section: C_L × C_6 × C_6, L ∈ {48, 64, 96, 128, 192}.
- Controls of a new size: ring C_L, L ∈ {384, 768}; torus C_L × C_L, L ∈ {40, 80}; 3D torus C_L³, L ∈ {18, 26}.
- "Realised crossover" uses the same definition as before (exact max drop ≥ 0.6 inside the exact window), computed from exact curves BEFORE simulation. Cells that do not qualify are reported, not dropped.

**Budgets and seeds:** wpn ∈ {8, 32, 128}; 200 datasets per cell; seeds SeedSequence([22, 11000 + id, ds]).

**Arms:** E1q, E4c, E4r, E7 and E8, with code frozen at the commit of this amendment. The estimator files must be byte-identical to 1ba9ca6 (E1q/E4) and to Amendment 8 (E7/E8).

**Gates** (numeric, fixed now; each reported with CP intervals):
- G1: E1q realised-crossover coverage at wpn 32 ≤ 40%.
- G2: E1q coverage strictly decreases from wpn 8 → 32 → 128, with CP-upper(128) < CP-lower(8).
- G3: E4c and E4r realised-crossover coverage ≥ 90% at every wpn.
- G4: control coverage ≥ 93% pooled (no cell with CP-upper < 93%) for E4c, E4r, E7 and E8 at every wpn.
- G5: control FPR (margin 0.25) ≤ 5% for all five arms at wpn 32.

**Verdict:**
- the narrowed headline GENERALISES iff G1–G5 all pass;
- otherwise the paper states exactly which gate failed, and the claim is restricted to the families where it held.

## AMENDMENT 12 (2026-11-23, written BEFORE any run; all earlier verdicts, bars and code UNCHANGED). GASPAR-KUKAL ESTIMATOR (closes referee item 5 of Amendment 10)
Source: the author-shared full text of Gaspar and Kukal, Physica A 602 (2022) 127624 (private copy, read in full 2026-11-23; analysis in ideas/A_spectral_dimension_coverage/GASPAR_KUKAL_COMPARISON.md). The earlier line "is NOT implemented" in Amendment 5 is superseded by this amendment.
**What is implemented (their estimator as published).** Spectral dimension from the return probability: regress ln p_hat(t_k) on ln t_k with weights w_k = p_hat/(1 - p_hat) (their Eqs 17, 23, 24, WSSQ), their classical reference without weights (Eq 17, LSSQ), d_s = -2 x slope. Their uncertainty is the regression standard error ("point estimates ... together with standard errors", Sec 3.1) used in a normal Z test at alpha = 0.05 (Eq 37); they never report interval coverage.
**Fixed adaptations, because our data differ (all fixed now):** (i) time points are our geometric grid, not consecutive even steps; (ii) the window is the common data-driven rule of this plan (J >= 6, else abstain), NOT their window (t0, m), which they choose to minimise MSE among windows unbiased against the KNOWN truth (Sec 4, Tables 4-7) and so is not available to a user; (iii) standard error = residual-scaled weighted least squares (df = J - 2), band +-1.96 se; (iv) walk dimension and d_f are not estimated.
**Arms (same walker data per dataset as the reference arms):** GKW (WSSQ, analytic se; THE primary arm), GKL (LSSQ, analytic se), GKWm (WSSQ, model-based se with sigma^2 = 1/n, n = walkers, no residual scaling), GKWb (WSSQ point estimate, our two-level block bootstrap, B = 200, weights recomputed on each resample, same max-t construction as E0g). Two-segment extension (OURS, not in their paper, to obtain a change test exactly as E0): GK2W, GK2L (analytic se per segment, contrast se = sqrt(se1^2 + se2^2), band +-1.96), GK2Wb (bootstrap). Reference arms E0g and E0 (frozen p22lib) are recomputed on the same data. Truth for coverage: exact d at the log-mean time of each segment (as E0g/E0). Additionally reported (ungated): coverage of the NOMINAL constant dimension (1, 2, 3) on controls, which is the target their Z test addresses.
**Cells, budgets, seeds:** exactly the Amendment 9b cell set (14 controls + the realised crossover cells; ids 9000/9100/9200 + index), wpn in {32 (gated), 8, 128 (reported)}; 200 datasets per cell; fresh seeds SeedSequence([22, 12000 + id, ds]). Code: p22amend12.py, run_amend12.py, analyze_amend12.py at the commit of this amendment; p22lib.py and run_p22.py must be byte-identical to Amendment 11 (sha256 asserted at start).
**Bars (Amendment 9 scoring, unchanged):** (C) simultaneous coverage pooled over non-abstained controls >= 93% with no control cell CP-upper < 93%; (F) FPR (margin 0.25) <= 5% with no cell CP-lower > 5% (two-segment arms only; one-scale arms have no change test); (M) minimum detectable L at power >= 0.8 for that L and all larger tested L (abstention = miss), per crossover family, two-segment arms only. Coverage on realised crossovers (X) is reported descriptively for one-scale arms. Gated at wpn 32; wpn 8 and 128 reported.
**Predictions recorded before data (not gated):** (1) GKW point estimates are nearly identical to E0g (their weight p/(1-p) differs from count weighting only by the factor 1/(1-p), at most about 20% at the first grid point); (2) analytic se (GKW) is narrower than the bootstrap se (GKWb), because time points share walkers and the regression se treats them as independent, so GKW coverage <= GKWb coverage; (3) GKW/GKL fail C at one or more budgets; (4) the weighting gives at most a small change in coverage or minimum L relative to E0g/E0.
**Reading rule:** either outcome is reported as found. If GKW (their estimator with their se) passes C, the paper says that on these graphs the published estimator with its regression se is calibrated for a constant exponent, and that what remains unmeasured is the crossover/minimum-size question; if it fails C, the paper reports that and the budget at which it fails. No earlier verdict is changed. A failure is a failure; no second run.
**Compute:** host-D (WSL), at most 6 worker processes, order wpn 32, 8, 128, resumable with a manifest. If wpn 128 does not finish in a reasonable time its unfinished cells are cut and stated.
