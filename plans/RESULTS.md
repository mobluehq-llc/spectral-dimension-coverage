# P22 RESULTS - scale-dependent spectral dimension d_s(t) on finite graphs
Reproduce: `uv venv .venv && uv pip install -r requirements.txt && .venv/bin/python checks.py && .venv/bin/python run_p22.py && .venv/bin/python analyze.py`
(run took 88.5 min on 10 cores; 39 cells x 200 datasets, all replicates kept, no cells cut; the 20-dataset-per-cell smoke run is `raw_smoke.csv`).
Pre-registered plan: `PLAN.md` (committed before any simulation). **No amendments to the pass bar.** Deviations/clarifications are listed at the end.
Files: `results/raw.csv` (52,200 dataset-estimator rows), `cell_summary.csv`, `verdicts.csv`, `min_detectable_L.csv`, `power_table.csv`, `ablation_parity.csv`,
`descriptive_gasket_rrg.csv`, `checks.csv`, `convergence_halves.csv`, `fig1..fig4*.png`.

## Plain verdict (pre-registered bar: C = simultaneous coverage >= 93% pooled on controls and no cell significantly below; F = FPR <= 5%; M = minimum detectable L reported)
| estimator | simult. coverage, controls (95% CI) | worst cells | FPR (margin rule) | min L 2D->1D (C_L x C_16) | min L 3D->1D (C_L x C_8^2) | verdict |
|---|---|---|---|---|---|---|
| E0 two-segment straight line | 94.3% (93.2-95.3) | none | 0.26% | 192 (N=3072) | 96 (N=6144) | **PASS** |
| E1 local-linear, split-sample (blueMonster hero) | **91.9%** (90.6-93.1) | torus3 L=28, 32 | 0% | 256 (N=4096) | 96 (N=6144) | **FAIL (C)** |
| E1x same, no split (leakage ablation) | **87.9%** (86.4-89.3) | ring 256, torus3 L=24, 28, 32 | 0% | 192 | 96 | FAIL (C) |
| E2 Hutchinson trace, probe bootstrap | **91.7%** (90.0-93.2) | ring 1024 | 0% | 192 (N=3072) | 128 (N=8192) | **FAIL (C)** |
| E2h Hutchinson, Hoeffding band | 100% | none | 0% | not reached | not reached | passes C,F vacuously; **fails M** (half-width ~76 in d) |
| E3 finite-size collapse (Reitz-Bianconi style, simple version) | 94.1% (93.1-95.0) *vs the thermodynamic value, an easier target* | none | 0% | 128 | 48 | **PASS** (with caveats below) |

Min-L = smallest tested L with power >= 0.8 at it and every larger tested L (power = correct-sign change declared; abstention counts as a miss). The strict Clopper-Pearson-lower-bound
version gives identical L in every row. Power at the reported minima is 0.95-1.00 except E3 3D->1D (0.89). Real crossovers exist in the exact window only for L>=128 (2D->1D) and L>=48 (3D->1D);
the pre-registered "realised" flag excludes smaller cells from power.

**So: of the four requested estimators, E0 and E3 pass the pre-registered bar; the blueMonster hero E1 and the runner-up E2 narrowly fail on simultaneous coverage (91.9% and 91.7% vs 93%).**
Both are good "noise" estimators and bad "bias" estimators; see below. E2h is valid but useless.

## What the numbers actually say
1. **E1's failure is smoothing bias, not noise.** Against its own population functional E1 covers 92-96% in every cell; against the exact d(t) it falls to 84-88% on 3D tori (L=28,32) and to **0-38% on realised crossovers**
   (3D->1D: 0-6%; 2D->1D: 26-38%) and 21-55% on the gasket (n=6,7), 0-6% on random regular graphs. The data are so precise (32 walkers/node) that the local-linear bias (CV picks h ~1.5 in log t on controls, ~0.75 on crossovers)
   exceeds the standard error where d(t) bends. A change test is still sensitive (power 0.97-1.0 at L>=96 for 3D->1D) because bias does not hide a 1-2 unit drop; the band is just not a calibrated statement about d(t).
2. **Bandwidth leakage is real and measurable.** Choosing window and bandwidth on the same data as inference (E1x) lowers control coverage from 91.9% to 87.9% (torus3 L=32: 84% -> 74%) with narrower bands (mean half-width 0.096 vs 0.149).
   The split-sample version spends half the data to buy ~4 coverage points; it costs ~1 grid step of minimum L on 2D->1D (E1 256 vs E1x 192).
3. **E0's pass is cheap.** The two-segment line only reports two numbers (low-t and high-t half of the window), cannot localise the crossover, and its exact-d coverage on crossovers is 20-97% (the segments straddle the bend);
   it passes because on constant-d controls any unbiased linear summary is calibrated and the equivalence margin (0.25) absorbs the rest. It is nonetheless a serious baseline: it beats E1 on minimum L for 2D->1D (192 vs 256).
4. **E2 (Hutchinson) has the best bias but the shortest window.** Its secant needs t and 2t, and its N*q>=6 window rule leaves < 6 usable points on small graphs: **abstention on controls is 57% (E0/E1 31%, E3 10%)**.
   Where it runs it is accurate (mean absolute error 0.044) and powerful (power 1.00 at L>=192 / 128). Its coverage shortfall comes from ring 1024 (86%) and a few cells with 12-18% literal-rule FPR; the first-100 vs last-100 coverage differ by 2.3 SE (0.90 vs 0.94), so its 91.7% is itself noisy.
   Literal (equivalence-free) rule "band excludes 0": FPR 11.3% for E2, ~4.7-5.1% for E0/E1/E1x/E3 (column `fpr_zero_rule`).
5. **E3's pass needs three caveats.** (i) It is handed three sizes of each graph (~1.75x the walkers) and a different, easier target (the thermodynamic exponent, constant on controls). (ii) On non-realised cells it "detects" a change
   (51% at C_64xC_16, 93% at C_96xC_16, 78% at C_32xC_8^2): collapse failure flags anisotropy / non-scaling early, which is not the same as detecting the 2D->1D crossover and cannot say where. Its apparent min L (128 / 48) is therefore an anisotropy-detection size, not a localisation size.
   (iii) It is our simple version (collapse of N P(t) vs t/N^(2/d) over L, 0.75L, 0.5L, early/late halves), not the published pipeline.
6. **Minimum size.** For a genuine 2D->1D crossover with w=16 the crossover sits at t~140-770; a detector needs the data window (N*P>=6, i.e. t up to ~L^2/(4 pi) scale) to reach well past it: **L~192 (N~3000)** for E0/E1x/E2, 256 for E1, i.e. L ~ 12 w. For 3D->1D with w=8 (crossover t~35-100) **L~96 (N~6000)**, L ~ 12 w again (E2: 128). Below these, the window ends inside the crossover. This is a statement about this design (32 walkers or probes per node, lazy walk, threshold 6), not a universal law.
7. **Parity / laziness (ablation, E1/E2).** With the non-lazy walk and a time grid containing odd times, the return probability is exactly 0 at odd t on even-length bipartite graphs: **every estimator abstains in 100% of datasets** (ring 256, torus 48).
   With an all-even grid the non-lazy walk works with the same coverage as the lazy walk (E1 92-96%). So parity is a grid-design trap, and laziness removes it; it did not change calibration otherwise.
   (Truth in the non-lazy ablation is the plateau, since d(t) has no analytic continuation for negative eigenvalues.)
8. **Fractal and expander stress (descriptive, not scored).** Exact finite-graph d(t) on the gasket in the window is 1.30-1.34 (n=5-7), i.e. 0.03-0.06 BELOW the asymptotic 2ln3/ln5=1.365 even with no noise, and it oscillates +-0.05 (log-periodic).
   E2 covers the exact curve at ~91% (n=6,7); E1 covers only 55% and 21% (bias vs log-periodic wiggles) while covering its own functional 93-96%. On random 3-regular graphs d(t) is not constant (1.33 at t=2 then rising/falling with equilibration): E1 coverage 5% / 0%, E2 92% / 90%. Degree effects are essentially a bias-curvature effect; the secant (E2) tolerates it better than a log-t local-linear fit.
9. **Convergence checks all passed** (`checks.csv`): sampler vs exact (|z|<2.2), probes unbiased (|z|<1.5), spectral q vs dense matrix powers (1e-16), d(t) vs finite differences (1e-9), bootstrap sd B=200 vs 1000 within 6.6%.
   Replicate halves agree within 2 SE for all estimators except E2 (z=-2.3; see 4).

## Honest limits
- Abstention is large for small graphs (window rule needs >= 6 grid points), so "pooled coverage" is conditional on non-abstention; 3D tori with L<=20 and 2D tori with L<=32 are essentially unanalysable at N*P>=6.
- Coverage targets the exact pointwise derivative at each estimator's reference scale, as pre-registered; a gentler target (own population functional) would put E1 at 93-96% everywhere (shown, not scored).
- Margin 0.25 was pre-registered to keep FPR meaningful; the literal "excludes 0" rates are reported. Deterministic graphs have a trivial graph level in the "graphs x walkers" bootstrap (only the random-regular cells have a real graph level; there E1 two-level coverage is poor for bias reasons).
- Walker blocks (200) are the bootstrap unit, not single walkers; E1 bootstraps 100 blocks (half sample).
- E3 is a simple re-implementation; no claim about the published method's actual performance. Novelty claims need a hand literature check (DIMFLOW_LITCHECK was snippet-level).

## Deviations / clarifications (none changes the pass bar)
- `E0g` (single global line) is reported but not scored; the scored E0 is the two-segment version (a single slope cannot test for change).
- E2h min-L listed "not reached" and flagged as a vacuous pass of C/F; verdict column `PASS` in `verdicts.csv` therefore reads True for E2h but the table above overrides it (fails M).
- Non-lazy ablation cells use plateau truth, no E3, E1/E2/E0 only.

---
# Amendment results (PLAN.md AMENDMENT 1, pre-registered and committed `cd56021` before the run)
Run: `run_p22b.py` (27 cells = 14 controls + 13 crossovers x 200 datasets, same seeds; walker blocks verified bit-identical to the first run, E1 reproduced exactly on 216 smoke rows), 15.2 min on 9 workers; `analyze_amend.py`.
Files: `results/raw_amend.csv`, `amend_final_table.csv`, `amend_cell_summary.csv`. The first-run verdicts above are unchanged (A5); this section adds estimators. Three bias-correction arms were looked at (3 looks); the pre-registered hero is E1q.

## Final table (all estimators; controls = ring/torus2/torus3, 7,000 datasets pooled, abstention excluded)
| est | coverage on controls, simult. % (CP 95%) | coverage vs EXACT d(t), realised crossovers (cell range) | FPR % (margin rule) | premature-change rate, non-realised crossover cells (n non-abstained) | min L 2D->1D | min L 3D->1D | C / F |
|---|---|---|---|---|---|---|---|
| E0 two-segment line | 94.3 (93.2-95.3) | 80.1 (20-97) | 0.3 | 5.5% (201) | 192 | 96 | pass / pass |
| E1 local-linear split (orig. hero) | 91.9 (90.6-93.1) | 16.2 (0-38) | 0.0 | 0.0% (203) | 256 | 96 | FAIL / pass |
| E1x no split | 87.9 (86.4-89.3) | 12.9 (0-32) | 0.0 | 1.0% (201) | 192 | 96 | FAIL / pass |
| E2 Hutchinson | 91.7 (90.0-93.2) | 87.2 (71-95) | 0.0 | n=1, uninformative (abstains) | 192 | 128 | FAIL / pass |
| E2h Hoeffding | 100 | 100 | 0.0 | n=1 | not reached | not reached | vacuous (fails M) |
| E3 3-size collapse | 94.1 (93.1-95.0) (vs thermodynamic d) | n/a (targets thermodynamic d) | 0.0 | **47.4%** (999) | 128 | 48 | pass / pass |
| **E3s single-size (new)** | 89.4 (88.2-90.5) | 23.1 (0-64) | 3.2 (CP 2.6-4.0; torus3 L=16 cell > 5%) | **88.8%** (1000) | 128 | 48 | FAIL / FAIL (cell) |
| E1u undersmoothed (h = 0.5 h_cv) | 93.7 (92.5-94.7) | **85.7 (64-94)** | 0.1 | 0.0% (203) | not reached (>256; power 0.72 at 256) | 96 | pass / pass |
| **E1q bias-corrected, pre-reg. hero** | **94.4 (93.2-95.3)** | 27.2 (0-88) | 0.6 | 1.5% (203) | 256 | 96 | **pass / pass** |
| E1qu quadratic + undersmooth | 94.1 (93.0-95.2) | 43.4 (20-88) | 0.8 | 1.5% (203) | 256 | 96 | pass / pass |
Min L is the first L with power >= 0.8 at it and all larger tested L (CP-lower-bound strict variant identical except E2 3D, where both agree, and E1u 2D "not reached"). Premature rate = primary-rule correct-sign declaration, among non-abstained datasets; any-sign rate within 1 point everywhere.

## What the amendment shows
1. **Standard bias control fixes the pre-registered failure, partly.** E1q (quadratic + same CV scheme, bootstrap re-running the correction) lifts control coverage 91.9 -> 94.4% with FPR 0.6%, so the hero now PASSES C and F (min L 256 / 96, same as E1). Every arm passes C and F. E1x-style leakage is not involved (all arms keep the split).
2. **But the pre-registered bar does not measure what E1 was failing at on crossovers, and E1q does not repair it there.** Coverage of the exact d(t) on realised crossovers is 27% for E1q (0-88%; worst L=192/256 2D and 3D L>=128 where it is 0-3%), versus 16% for E1. The local-quadratic fit removes the O(h^2) term but CV picks h = 1.7-2.6 (log-t units, vs 0.5-1.2 for linear) so on the sharp 3D->1D step the quadratic is still too smooth; undersmoothing is what works: **E1u (h = 0.5 h_cv) covers 85.7% (64-94%)** at cost of wider bands (mean half-width 0.3-0.9 vs 0.2-0.4) and lower power (2D->1D power 0.72 at L=256 so "not reached" in this sweep). E1qu sits between (43%).
   So: noise is calibrated for all arms (own-functional coverage 90-97%); **bias on crossovers is controlled only by undersmoothing here, and it costs power.** Coverage-of-exact-d and power trade off; no tested arm has both >= 93% crossover coverage and min-L 192.
3. **Premature change (A2).** The two-segment/local-polynomial estimators (E0, E1 family) are essentially never premature (0-5.5%), but that is mostly **protection by abstention**: below L=96 (2D) and L=48 (3D) they abstain in ~100% of datasets (window rule), and the few that run (L=96, ~200 datasets) declare 0-5%. They do not "see" a crossover that is not in the window. E3 (3-size collapse) declares a flow in **47%** of non-abstained non-realised datasets (13.5% at C_48xC_16, 51% at C_64, 93% at C_96, 78.5% at C_32xC_8^2), and the single-size E3s in **89%** (50-100% per cell). Reading: collapse-based change tests respond to finite-size anisotropy/equilibration and not to a realised crossover; E3's apparent min L (128/48) is an anisotropy detection size. Caveats: the exact drop at C_96xC_16 is 0.34 and at C_32xC_8^2 is 0.57 (below the 0.6 'real' threshold), so a sliver of E3's declarations there track genuine drift; at C_32/48xC_16 (drop <= 0.02) they do not. E3s is E0 with E3's equilibration-regime window (N P >= 1.05): the window then runs into saturation, where the log-log slope falls toward 0 for any estimator, which is what drives its 89% and its 3.2% FPR (torus3 L=16 cell significantly > 5%).
4. **E3s verdict:** the single-size version of the collapse estimator FAILS both C (89.4%) and F; E3's pass depends on the cross-size collapse and the extra ~1.75x walkers, as flagged in caveat (i) of the first run.
5. **E0b not run** (A3): E0's band was already max-t over both segments.
## Limits of the amendment
Three bias-correction arms were examined; E1q was pre-designated hero but the best crossover-coverage arm (E1u) was not, and calling the family "passes" is a 3-look statement. The premature-rate comparison across estimators confounds detection behaviour with abstention (denominators shown). E3s is not a published method. Truth for E3/E3s controls differs (thermodynamic vs exact d at segment centres, as pre-registered).

---
# Causal-set results (PLAN.md AMENDMENT 2, committed `4f7fe74` before the run; 2026-11-20)
Code: `causet_lib.py` (sprinkling, graphs, exact q(t)), `run_causets.py`, `analyze_causets.py`, `report_causets.py`, `calib_from_p22.py`, `ref_peak.py` (post hoc), `e3_wide.py` (exploratory). Data: `results/causet/*.pkl` (704 task files, 21 of them N=8000 reps; 3 empty marker files).
Outputs: `causet_estimates(_annotated).csv`, `causet_curves.json`, `causet_delta_table.csv`, `causet_abstention.csv`, `causet_stability_t8.csv`, `causet_reference_peak.csv`, `causet_window_span_extrapolation.csv`, `causet_p22_calibration.csv`, `causet_graph_diag.csv`, `causet_E3_wide_exploratory.csv`, `fig5..fig7*.png`.
**No exact d_s(t) exists for any of these graphs** (random, finite, boundary-affected); nothing below is "coverage of the truth". Null = matched-degree random geometric graphs (RGG) in a box of dimension 2 (for 1+1) / 4 (for 3+1), same N and mean degree, same boundary effects.

## What was implemented (exactly)
Poisson-like (fixed-N uniform) sprinkling into the Alexandrov diamond of flat Minkowski 1+1 (light-cone square, as in `p14_lorentz_sprinkling.py`) and 3+1; N = 500, 1000, 2000, 4000, 8000; **R = 100 / 100 / 100 / 40 / 12 sprinklings in 1+1 and 100 / 100 / 100 / 40 / 9 in 3+1** (3+1 N=8000 was cut from 12 to 9 for time: three of its tasks were replaced by skip markers after the first wave measured ~35 min per sprinkling under load; ">= 100 where feasible" is met only for N <= 2000). Total run 130 min on 9 workers (over the ~2 h target).
Graphs (all with the p22 lazy walk S = I/2 + P/2, exact q(t) = tr S^t / N from the dense spectrum, largest connected component):
 * **k0: Hasse diagram (links) made undirected, p = 1/degree - the Eichhorn-Mizera (arXiv:1311.2530) construction** (their walk is non-lazy; here lazy, so t_here ~ 2 sigma_EM; no IR cutoff on link length; boundary effects kept).
 * k3 (relations with <= 3 elements between: support of the BDG operator layers) and k15: my unweighted "non-locality dial" constructions, NOT published graph walks.
 * **The Belenchia-Benincasa-Marciano-Modesto study (arXiv:1507.00330) is a continuum heat kernel of the regularised non-local d'Alembertian, not a graph walk; no nearest-neighbour proxy was found in it and none is implemented under that name.** So the "two published studies" differ in the object studied, not only in a graph choice (and its d_s -> 2 small-scale result is outside what a walk at t >= 2 lazy steps can see).
Estimators = p22 code: E0 (two-segment line), E1q (Amendment-1 hero), E2 (secant of per-sprinkling q; units = sprinklings, Hoeffding variant not run), E3 (collapse across consecutive size triplets). Dataset = the ensemble of R sprinklings at one N (bootstrap over sprinklings, B = 200 / 400 / 80); pseudo walker budget 32 N to keep the p22 window rule. Grids G2 = {2,3}*2^j and G4 = {4,5,6,7}*2^j; window thresholds N q >= 6 and >= 1.3.

## Headline numbers (Hasse/link graph k0 vs its matched RGG; G4 grid, window N q >= 1.3)
Reference column = numerical d_ref(t) = -2 dln qbar/dln t of the ensemble-mean q (post hoc diagnostic, see below; not a truth). Delta = extremal d_j - d_0 with the declared sign; r = rise declared, d = drop declared, n = neither (margin 0.25).
| | N=500 | 1000 | 2000 | 4000 | 8000 |
|---|---|---|---|---|---|
| 1+1 causet: d_ref peak (at t = 10 -> 14) | 4.57 | 5.19 | 5.84 | 6.49 | 7.20 |
| 1+1 causet: **E1q** Delta | +0.70 r | +1.49 r | +1.97 r | +3.04 r | +3.70 r |
| 1+1 causet: **E0** Delta | -0.77 d | -0.32 d | -0.49 d | +0.12 n | +0.37 r |
| 1+1 causet: **E2** Delta | -1.16 d | -0.78 d | -0.93 d | -0.41 d | -1.25 d |
| 1+1 RGG control: d_ref peak (t = 4) | 2.41 | 2.62 | 2.81 | 2.97 | 3.11 |
| 1+1 RGG control: E1q / E0 / E2 Delta (all drop) | -1.0 / -0.23 n / -0.9 | -1.2 / -0.25 n / -1.0 | -1.4 / -0.27 / -1.05 | -1.5 / -0.28 / -1.1 | -1.65 / -0.31 / -1.2 |
| 3+1 causet: d_ref peak (t = 7 -> 10) | 6.15 | 7.36 | 8.52 | 9.71 | 10.89 |
| 3+1 causet: **E1q** Delta | +1.96 r | +2.98 r | +4.63 r | +6.11 r | +7.79 r |
| 3+1 causet: **E0** Delta | -0.13 n | +0.69 r | +1.77 r | +1.24 r | +2.92 r |
| 3+1 causet: **E2** | abstains at every N | | | | |
| 3+1 RGG control: d_ref peak | 4.02 | 4.84 | 5.63 | 6.48 | 7.22 |
| 3+1 RGG control: E1q Delta (drop) | -2.4 | -3.4 | -4.0 | -4.8 | -5.3 |
Abstention (pooled over all cells; share of datasets with < 6 window points): thr 1.3 - causets E0 0%, E1q 13%, E2 72%; RGG 0% / 0% / 0%. thr 6 - causets 65% / 65% / 100%; RGG 25% / 25% / 62%. E3 never abstains.
Figures: `fig5_causet_dcurves.png` (E1q curves, causet vs RGG, all k), `fig6_causet_delta_vs_N.png` (Delta vs N), `fig7_causet_reference_dt.png` (reference d(t): rises to a peak, then collapses as N q -> 1).

## What the estimators say (answers to the pre-registered questions)
**(a) Dimension flow and stability in N.** The graphs do NOT have a plateau dimension: the ensemble-mean d(t) of every causet graph rises with t to a peak of 4.6-11 (more than the topological 2 or 4) at t ~ 7-14 lazy steps, then collapses as the walker reaches the boundary / equilibrates (N q -> 1; EM report the same "increasing spectral dimension" and boundary collapse). The peak grows with N (exponent d_peak ~ N^0.13-0.16 in 1+1, 0.19-0.21 in 3+1, `causet_reference_peak.csv`) and so is NOT stable as N grows; but the matched RGG controls also grow (exponent 0.07-0.09 in 1+1, 0.20-0.21 in 3+1). The level difference to the control is large in 1+1 (7.2 vs 3.1 at N = 8000) and small in 3+1 (10.9 vs 7.2).
 * **E1q reports a monotone RISE on every causet graph, in all 73 non-abstained causet datasets (both grids and thresholds) (Delta +0.5 ... +9.4), growing with N** - but most of that is an edge artefact: the reference curve peaks at t ~ 12 and falls to ~3.3 by t = 28 (1+1, N = 4000) whereas E1q's d keeps rising to the last window point (7.1), because the CV-chosen bandwidth (0.7, the smallest allowed for the quadratic) cannot follow the sharp downturn at the window end. Its band is a statement about the wrong curve shape there.
 * **E2 (the secant) runs on only 17 causet datasets (all 1+1, thr 1.3)** and its sign depends on the grid: on G4 (window reaches past the turn-over, 10 datasets) it reports a DROP (Delta -0.4 ... -2.4), on G2 (shorter window, ends near the peak, 7 datasets) a RISE (+1.7 ... +4.3); 3+1 never reaches it. So E1q (rise) and E2-on-G4 (drop) give OPPOSITE signs on the same graph, and E2 itself disagrees with E2 across grids. E0's sign also depends on the window start (G4: drop at N <= 2000, rise at N = 8000 in 1+1; G2/G4 sign agreement 61% on causets, 88% on RGG; E2 0% on causets (7 pairs), 100% on RGG; E1q 100% on causets (consistent artefact) and 62% on RGG). On the RGG controls E1q and E2 agree on a drop.
 * **E3 (collapse, p22 settings) returns the grid ceiling 3.6 for both halves on every causet** (search range 0.5-3.5): no collapse exponent exists in its range; the "no change" it reports is a failure, not a measurement. An EXPLORATORY widening of the range to 14 (`e3_wide.py`, not pre-registered) gives early/late exponents of 6-14 that depend on the triplet and pin at 14.1 in 3+1: collapse over N is not meaningful for these graphs (the walk is not diffusive in N^(2/d)).
**(b) Does the construction matter?** Sign and shape do not: k0, k3, k15 all give rise-then-collapse in every estimator that runs. Magnitude does: the peak d at N = 8000 is 7.2 / 8.4 / 9.8 (1+1) and 10.9 / 11.2 / 11.4 (3+1) for k0 / k3 / k15, but the matched RGG rises with k too (3.1 / 5.0 / 7.0 in 1+1), so the causet-over-control excess SHRINKS with k (x2.3, x1.7, x1.4); the window shortens with k (span 0.9 -> 0.6 decades). So in this experiment the link/undirected-Hasse choice versus a bigger edge set changes the size of the transient, not the qualitative answer, and there is no construction-driven reversal that would reconcile "d_s rises" (walk on links) with "d_s -> 2" (continuum non-local operator). That second result is a different object and is not reproduced or refuted here. Answer: **no, not from this family of graph constructions; the published disagreement is more plausibly about the object (walk on a graph vs heat kernel of a continuum operator) and about N and boundaries (EM use N up to 1e8; we use <= 8000).**
**(c) Minimum N for a change of the claimed size, p22 calibration.** At p22's minimum detectable L the estimator data windows spanned (decades of t): E0/E1x 1.36 (2D->1D) and 0.90 (3D->1D), E1q 1.66 / 0.90, E2 1.36 / 1.20 (`causet_p22_calibration.csv`; E0 via E1x since E0's own J = 2 segments; E3 has no comparable window). The causet windows (E1q rule, thr 1.3) span only 0.78-1.08 decades (G2) / 0.40-0.85 (G4) in 1+1 and 0.40-0.78 in 3+1 at N = 500-8000, growing ~0.04-0.06 decades per doubling of N in 1+1 and ~0.00-0.04 in 3+1 (`causet_window_span_extrapolation.csv`). Extrapolating linearly in log2 N (FLAGGED as extrapolation, grid- and threshold-dependent): the 0.9-decade p22 3D->1D requirement is met at N ~ 5e2 (1+1, k0, G2) to ~1e4-4e5 (k3/k15, G4); 3+1 needs N ~ 1e7 or more (G4) and never reaches it on G2 (flat); the 1.66-decade 2D->1D requirement needs N ~ 1e7 (best case) up to > 1e12. E2 (needs t and 2t) abstains at thr 6 everywhere and at thr 1.3 in all of 3+1. This is a window-sufficiency (necessary) condition, not a power curve: since the dimension itself is not constant here, "a change of size 2 -> 4 or 4 -> 2" is not a well-posed target on these graphs at N <= 8000 (d_ref moves by 2-7 units rise and fall inside the window by itself, and RGG nulls move by 1-5 units).

## Verdict per estimator on this case
| estimator | what it reports on causet graphs | usable? |
|---|---|---|
| E0 two-segment | sign depends on grid/N (drop at small N, rise at large N, "none" at the crossing) | no (unstable); abstains 65% at thr 6 |
| E1q (hero) | monotone rise, Delta +0.5..+9, grows with N | **misleading at the window edge** (reference curve turns over); calibrated noise band, wrong shape |
| E2 | drop when its window passes the turn-over (G4), rise when it does not (G2) | follows the curve shape best where it runs, but abstains 72-100% (needs t and 2t) and its sign depends on the window |
| E3 | pinned at the search ceiling 3.6, "no change" | no collapse exists; not informative |
## Caveats and deviations (none changes the p22 pass bar or earlier verdicts)
* Sizes are 4-5 orders of magnitude below EM's (1e6-1e8); R = 40 (N = 4000) and 12 / 9 (N = 8000) are below the target of 100; the 130 min run exceeded the 2 h target.
* Fixed-N uniform sprinkling, not Poisson-N; lazy walk (EM: non-lazy); no IR cutoff; sprinkling region is a diamond, boundary effects dominate the large-t end (that is where the collapse comes from). The k3/k15 graphs and the RGG nulls are my constructions; the RGG null has a boundary but a different geometry (Euclidean box), so "excess over the control" is suggestive, not a clean null.
* Unit of analysis is the sprinkling ensemble (one dataset per cell): no coverage/power/FPR was computed; abstention is pooled over cells, not a rate with intervals. Changes of `delta` convention: the table uses dmin when a drop is declared and dmax otherwise (the first annotated table used dmax only; fixed before any result was written up).
* Post hoc, labelled: `ref_peak.py` (ensemble-mean central-difference d_ref and its peak, added after E1q and E2 disagreed in sign; has secant/ensemble bias, is NOT an exact truth) and `e3_wide.py` (E3 search range widened after the pinned result). The span extrapolations assume linearity in log N and hold the window rule fixed. p22 calibration was built on drops of 1-2 units in deterministic lattices at 32 N walkers; the causet changes here are rises and falls inside one window.
* p22lib changes: `est_E2(step=)`, `est_E3(t0=, dmax=)` with unchanged defaults.

---
# Amendment 3: honest intervals (PLAN.md AMENDMENT 3, committed before the run; 2026-11-20)
Estimator E4 = Armstrong & Kolesar bias-aware fixed-length intervals for beta(u)=g'(u), u=log t, g=log P, d=-2beta (refs verified: Econometrica 86(2) 2018, 655-683; Quantitative Economics 11(1) 2020, 1-39). Weights per grid point from the convex programme in PLAN (third-derivative class, linear+quadratic Taylor terms annihilated; the quadratic constraint is NOT optional), S from the two-level bootstrap on half B, simultaneous critical value by bootstrap max-t, worst-case bias added linearly. M variants: **E4o** oracle sup|g'''| from the exact spectrum over the window range (benchmark), **E4c** local-quartic CV fit on half A, **E4r** global weighted quartic on half A (my adaptation of the AK/Fan-Gijbels rule; I did not read their Sec. 3.3 formula). Change test: contrasts a_j - a_0 with bias bound B_j + B_0 (the task's B(c) centred at one point is not valid for a difference of weights centred at two points).
Code `p22honest.py`, `run_p22c.py` (27 cells x 200 datasets, same seeds/walker blocks and the same bootstrap draws as E1; 15.8 min on 9 workers), `analyze_amend3.py`. Files: `results/raw_amend3.csv`, `amend3_final_table.csv`, `amend3_cell_summary.csv`, `amend3_run.log`. Abstention is identical to E1's (window rule; 31% on controls); one E4c dataset abstained for solver failure. Causal-set cells were not run (no exact truth).

## Final comparison (controls = ring/torus2/torus3, 7,000 datasets, abstention excluded; premature = non-realised crossover cells)
| est | control coverage, simult. % (CP 95%) | coverage vs EXACT d(t), realised crossovers (cell range) | FPR % | premature | min L 2D->1D | min L 3D->1D | C / F |
|---|---|---|---|---|---|---|---|
| E0 two-segment | 94.3 (93.2-95.3) | 80.1 (20-97) | 0.3 | 5.5% | 192 | 96 | pass / pass |
| E1 local-linear split | 91.9 (90.6-93.1) | 16.2 (0-38) | 0.0 | 0.0% | 256 | 96 | FAIL / pass |
| E1u undersmoothed | 93.7 (92.5-94.7) | 85.7 (64-94) | 0.1 | 0.0% | not reached (>256) | 96 | pass / pass |
| E1q bias-corrected | 94.4 (93.2-95.3) | 27.2 (0-88) | 0.6 | 1.5% | 256 | 96 | pass / pass |
| E2 Hutchinson | 91.7 (90.0-93.2) | 87.2 (71-95) | 0.0 | n=1 | 192 | 128 | FAIL / pass |
| E3 3-size collapse | 94.1 (93.1-95.0), thermodynamic d | n/a | 0.0 | 47.4% | 128 | 48 | pass / pass (caveats above) |
| **E4o oracle M** | 94.0 (92.8-95.0) | **98.5** (98-100) | 0.4 | 0.5% | not reached (>256) | 128 | **FAIL C** (ring 512 cell CP upper < 93%) / pass |
| **E4c held-out M** | 99.6 (99.2-99.8) | **96.4** (92-100) | 0.0 | 0.0% | not reached (>256) | 192 | pass / pass |
| **E4r rule-of-thumb M** | 99.2 (98.7-99.6) | **95.0** (90-100) | 0.0 | 0.0% | not reached (>256) | 192 | pass / pass |
(E1x, E2h, E3s as in Amendment 1/original runs; see `amend3_final_table.csv`.) Median M used: E4o 0.12, E4c 1.9, E4r 1.4 (data-driven M is ~12-16x the oracle). Median half-width in d: E4o 0.47, E4c 0.94, E4r 0.84 (E1q about 0.3). Power of E4 at 2D L=256: 0.25 (o), 0.24 (c), 0.32 (r) vs E1q 0.95; at 3D L=128: 0.89 / 0.58 / 0.67 vs 1.00.

## Reading
1. **Honest intervals deliver what they promise on coverage**, including through realised crossovers: 95-98% simultaneous coverage of the exact d(t) (E1q 27%, E1u 86%), 0% / <=0.5% premature and false positives. The data-driven M is conservative (~12-16x oracle) so control coverage is 99%+, i.e. the bands are wider than needed; with the oracle M, control coverage is on target (94.0%) but one cell (ring 512) fails the cell-level C rule, so even the benchmark does not formally pass.
2. **The price is power, and it is large.** No E4 variant reaches power 0.8 at any tested 2D L (<=256); on 3D->1D min L is 128 (oracle, infeasible) / 192 (E4c, E4r) versus 96 for E1u and E1q. Worst-case bias is about as large as the sampling noise in the half-width (bias bound 0.23-0.27 vs noise 0.23-0.25 for E4c/E4r), so the change band must clear the margin 0.25 by a wide gap.
3. **Does any estimator now beat the coverage-vs-power trade-off? No.** The pre-registered criterion (pass C and F, crossover coverage >= 93%, min L no worse than E1q's 256/96) is met by none: E4c/E4r have the coverage but min L 192 in 3D and not reached in 2D; E1q/E1u keep power but not both coverages. E4 moves along the same frontier (high coverage, low power) rather than shifting it; the closest thing to a free lunch is E2 (87% crossover coverage, min L 192/128) which abstains 57% of the time, and E0, which is coarse. The oracle benchmark shows part of E4's loss is M uncertainty (3D min L 128 vs 192) but even with the true M the 2D result is not reached, so the loss is mostly intrinsic to certifying a slope bound with this much noise.
## Limits
Weights are computed from S estimated on the same half B they are applied to (second-order). Simultaneous interval adds worst-case bias to a bootstrap max-t quantile (not the AK noncentral critical value; 1.96*sd+B in the objective as specified). Quartic rules for M are my adaptation. The "beats the trade-off" criterion was fixed in PLAN before the run. 3 M variants looked at, all reported.

---
# Amendment 4: shape-restricted band (E5) and split-sample logistic step (E6) (PLAN.md AMENDMENT 4, committed before the run; 2026-11-21)
Code `p22shape.py`, `run_p22d.py` (same 27 cells x 200 datasets, same seeds/blocks; 12.5 min on 9 workers), `analyze_amend4.py`, `e5_causets.py`. Files: `results/raw_amend4.csv`, `amend4_final_table.csv`, `amend4_cell_summary.csv`, `causet_E5.csv`, `amend4_run.log`. Dümbgen reference verified: Bernoulli 9(3), 423-449 (2003), doi 10.3150/bj/1065444812. Only its interval-average envelope is used (no scale penalty, so no optimality claim). E5 direction declared non-increasing for ALL families including controls. Abstention identical to E1 (31% of control datasets).

## Full comparison (controls pooled 7,000 datasets; crossover coverage = simultaneous coverage of EXACT d(t) on realised cells; min L = power >= 0.8, abstention = miss)
| est | control coverage % (CP 95%) | FPR % | coverage through realised crossover % (CP) | min L 2D->1D | min L 3D->1D | C / F |
|---|---|---|---|---|---|---|
| E0 | 94.3 (93.2-95.3) | 0.3 | 80.1 (77.8-82.3) | 192 | 96 | pass / pass |
| E1q | 94.4 (93.2-95.3) | 0.6 | 27.2 (24.7-29.8) | 192 (strict 256) | 96 | pass / pass |
| E1u | 93.7 (92.5-94.7) | 0.1 | 85.7 (83.6-87.6) | not reached | 96 | pass / pass |
| E2 | 91.7 (90.0-93.2) | 0.0 | 87.2 (85.0-89.3) | 192 | 128 | FAIL / pass |
| E3 | 94.1 (93.1-95.0) (thermodynamic d) | 0.0 | n/a | 128 | 48 | pass / pass (premature 47%) |
| E4c | 99.6 (99.2-99.8) | 0.0 | 96.4 (95.2-97.4) | not reached | 192 | pass / pass |
| E4r | 99.2 (98.7-99.6) | 0.0 | 95.0 (93.7-96.2) | not reached | 192 | pass / pass |
| **E5qu** (hero) | 92.6 (91.4-93.8) | 0.1 | 46.8 (44.0-49.6) | 256 | 96 | **FAIL C** (ring 512, 1024) / pass |
| E5q | 92.4 (91.1-93.5) | 0.1 | 30.1 (27.5-32.7) | 256 | 96 | FAIL C (ring 512, 1024) / pass |
| **E6** | 99.6 (99.3-99.9) | 0.0 | **99.1 (98.4-99.6)** (cells 98-100) | **256** | **96** | **pass / pass** |
Premature declarations on non-realised crossover cells: E5qu/E5q 1.5%, E6 0.5%. Median half-width of the d band: E5qu ~0.3-0.4, E6 0.5-2.2 (wide when the split-A set is large, narrow once the crossover is identified). Power E6 vs E1q: 2D L=192 0.56 vs 0.83; L=256 0.96 vs 0.95; 3D L=96 0.89 vs 0.995, L=128 0.995 vs 1.0. E6 uses on average 107 of 126 (m,w) grid points in the split-A set (median 126 on controls, where the step is not identified).

## Reading
1. **E5 does not beat the trade-off.** Its power equals E1q's (it has the same slope estimator, so the same 256/96) but it fails C on the ring (flat exactly-constant d: coverage 88-92% at ring 128-1024) and covers the exact crossover d(t) only 47% (E5qu) / 30% (E5q). The monotone envelope removes shape/noise error but NOT the slope estimator's smoothing bias at the crossover kink, and the lattice's exact d(t) is not monotone near the start of the torus controls (torus 2D L=32: 67%). The shape restriction therefore did not fix coverage; undersmoothing the slopes (qu) helped crossover coverage only partly (47% vs 30%).
2. **E6 passes the pre-registered criterion** (C and F pass, realised-crossover coverage 99.1% >= 93%, min L 256 / 96 = E1q's): it is the first estimator here with honest coverage through the crossover AND E1q-level minimum L. Qualifications, none of which changes the verdict but all of which limit it: (i) the tie is at the grid resolution — at 2D L=192 E6 power is 0.56 vs 0.83 for E1q, so E6 is weaker than E1q in the middle of the range and only equal at the pass point; E1u (86% crossover coverage) is the only comparator with the full 3D power; (ii) E6 is a parametric (logistic-in-log t) model: the crossover truth in these cells is a smooth single step, which is the favourable case; no non-logistic truth (two steps, ramp, non-monotone) was tested, and control coverage of 99.6% shows the union is conservative; (iii) the conclusion rests on one run of 200 datasets/cell with an alpha split of 0.01/0.04 chosen in advance.
3. Pre-registered answer: **E5: no. E6: yes (by the PLAN rule), with the stated qualifications.** E4 stays on the high-coverage/low-power end; E6 sits at E1q's power with E4-level coverage.

## Causal-set diagnostic (E5qu, no truth, G4 grid, both directions as two looks; `causet_E5.csv`)
100 of 120 (dim x N x construction x threshold) datasets are non-abstained. For every one, **both** the non-increasing and the non-decreasing E5 bands declare a change > 0.25. This is not evidence of a drop or a rise: the ensemble-mean curve rises then collapses (Amendment 2), so neither monotone shape is compatible with it and the multiscale envelopes separate (max interval mean minus min interval mean exceeds the margin) whichever direction is imposed. The monotone-band diagnostic therefore only rejects monotonicity of d(t) on these graphs (rise-and-fall), and the matched RGG controls behave identically; it cannot assign a direction or a dimension flow.
## Limits
Slope-estimator bias is not removed in E5; c uses the plain max-|z| (two-sided) rather than a scale-penalised statistic; E6 grid m/w discrete (21 x 6), A-covariance shrunk 10% to its diagonal (fixed beforehand), and its coverage is conditional on the logistic family. 3 new variants looked at (E5qu, E5q, E6), all reported.

---
# Amendment 5: E6 under misspecification, and published comparison arms (PLAN.md Amendment 5, commit 2775f1a, written before the run; 2026-11-21)
Code `p22amend5.py`, `run_p22e.py` (run on host-T with NCPU = 16; the local copy was NCPU = 9, the only difference; 1480 tasks, 35,000 rows, 2469 s), `analyze_amend5.py` (no estimator or pass-bar changes; the script makes no figures, so none are included). Files: `results/raw_amend5.csv`, `cells_meta_amend5.csv`, `amend5_run.log`, `amend5_tables.md`, `amend5_e6_cells.csv`, `amend5_e6_lof_by_group.csv`, `amend5_groups_all_estimators.csv`, `amend5_passbar.csv`. 37 cells x 200 datasets.

## (a) E6 under misspecification (non-abstained datasets; coverage = simultaneous coverage of exact d(t), CP 95%)
| class | n | E6 cov | E0 cov | E1q cov | E4c cov | LoF alarm rate | false change (null cells) |
|---|---|---|---|---|---|---|---|
| constant controls | 1930 | 99.6 (99.2-99.8) | 94.3 | 94.4 | 99.6 | 3.7 (2.9-4.7) | 0.0 |
| realised crossover | 1231 | 99.2 (98.5-99.6) | 80.1 | 27.2 | 96.4 | 3.2 | n/a |
| two-step (4 cells) | 800 | **54.2 (50.7-57.7)** | 22.5 | 3.8 | 95.2 | 11.9 (9.7-14.3) | n/a |
| brush (3 cells) | 600 | 93.3 (91.0-95.2) | 58.0 | 8.3 | 97.7 | 5.0 | n/a |
| gasket wiggle (3 cells) | 558 | 93.7 (91.4-95.6) | 78.1 | 67.7 | 98.7 | 6.5 | 0.0 (0.0-0.7) |

E6 coverage collapses on the two-step as the step sharpens (twostep_8_32: 0.995, 0.83, 0.31, 0.035 cell coverage from L=128 to 320; the pointwise coverage stays 0.71-0.99). No false change on gasket wiggle or constant controls; the sharper steps are declared (power 1.0), always with the right sign. E6 is still above E0 and E1q everywhere, but E4c (the conservative bound) is the only one that stays at about 95-99% on every misspecified class, at the price of power.

**Loud or silent.** Pooled over the misspecified cells (two-step, brush, gasket; n=1958): 441 misses, alarm given miss 18.8% (15.3-22.8), alarm rate 8.2% overall, 3.7% on constant controls (calibration near the nominal 5%). By the pre-registered rule (silent iff P(alarm|miss) < 0.5): **E6 fails SILENTLY.** Two-step alone: 366 misses, 20.2% alarmed. Guarding E6 by suppressing declarations when the alarm fires costs power (two-step declarations 100% to 88.1%, brush 99.8% to 94.8%) and does not fix the coverage.

## (b) Published arms (same pass bar; abstention counted as a miss for power)
| arm | control cov | FPR | cov through realised crossover | premature (non-realised crossover) | abstain (controls) | min L 2D / 3D | C / F |
|---|---|---|---|---|---|---|---|
| M (Millan windowed ML, analytic error) | 86.0 (83.9-87.9) | 0.2 | 87.4 | n/a (100% abstain) | 0.57 | 192 / 128 | FAIL (ring 256-1024, torus296) / pass |
| Mb (same, block bootstrap error) | 97.9 (97.0-98.7) | 0.1 | 98.9 (98.0-99.4) | n/a (100% abstain) | 0.57 | 256 / 128 | pass / pass |
| CL (optimised true collapse, Eq. 11) | 95.9 (95.0-96.6) | 0.2 | 32.2 (30.0-34.6) | **41.9 (38.8-45.0)** | 0.11 | 128 / 64 | pass / pass |
| E6 (reference) | 99.6 | 0.0 | 99.2 | 0.5 | 0.31 | 256 / 96 | pass / pass |
| E1q (reference) | 94.4 | 0.6 | 27.2 | 2.5 | 0.31 | 192 / 96 | pass / pass |

Brush case: M 60.2% coverage, Mb 89.3%, CL 26.2%; CL declares the brush change with power 1.0 but covers poorly. Two-step: M 80.5, Mb 99.5, CL 0.0.

## Reading
1. **E6 fails silently outside its logistic family.** The lack-of-fit diagnostic is calibrated (3.7% on controls) but catches fewer than one in five misses. The E6 claim stays "honest coverage through a smooth logistic crossover", not "honest under arbitrary d(t)". Sharp steps break it; use E4c-level conservatism, or report the claim as conditional.
2. **The published methods do not change the paper's claims.** Optimised true collapse (CL) reproduces the E3 result: good thermodynamic-d coverage but a 42% premature-change rate and 32% coverage of the exact crossover d(t), so it cannot be used to certify a change. Millan's windowed ML with its stated analytic error under-covers controls (86%) and abstains over half the time; with a bootstrap error (Mb, my addition, not in the paper) it is honest (98%) but needs L = 256 in 2D versus E6's 256 and E1q's 192, and abstains heavily. No published arm gives E1q-level power with honest crossover coverage.
## Limits
M, Mb and CL rest on my fixed assumptions about unstated details (PLAN Amendment 5); Gaspar-Kukal not implemented (paywalled). [Superseded 2026-11-23: the authors kindly shared the paper and the estimator was benchmarked in Amendment 12, below.] CL not run on gasket. Windows (>= 4 points) and the logistic LoF threshold were fixed beforehand.

# Amendment 6: minimum-size recipe, E6 misfit alarm, fresh-seed replication (PLAN.md Amendment 6, commit 39b9e27; predictions committed in 6ea633f before the new 2D cells ran; 2026-11-21)
Code `p22recipe.py`, `p22amend6.py`, `run_p22f.py`, `analyze_amend6.py` (heavy runs on host-T, system python3, NCPU = 16, ~2,600 s wall in total; local Mac smoke run used only to debug). Files: `results/raw_amend6_{big2d,sep,fresh,orig}.csv`, `cells_meta_amend6_*.csv`, `amend6_tables.md` (all tables), `amend6_recipe_{params,predictions,scores}.csv`, `amend6_alarm_{bars,cells}.csv`, `amend6_replication.csv`, `amend6_separation_scan.csv`, `amend6_run_*.log`. No pass bar, estimator or earlier verdict changed. Reproduction: the re-run of E6 on the original seeds reproduces raw_amend5 exactly on host-T (7,400 datasets: abstention and simultaneous coverage identical, max |hw diff| 0). On the Mac the same code differed in about 7% of datasets through floating-point ties, which is why the heavy jobs run on host-T.

## Item 1: minimum-size recipe -- PASS (7 of 7 within one grid step; the first five are a consistency check, only E4c/E4r are out of sample)
Rule (fitted on 3D cells only): detect iff J >= Jmin_e and (S - 0.25)/hw_e(J) >= T_e, hw_e(J) = a_e (J/7)^(-delta_e); fitted constants in `amend6_recipe_params.csv` (e.g. E6: Jmin 6, a 1.28, delta 2.3, T 1.19; E4c: Jmin 6, a 1.55, delta 1.46, T 1.90). Power = declaration rate, abstention a miss; min L = smallest L with power >= 0.8 for it and all larger.
| estimator | predicted 2D min L | measured 2D min L | grid-index distance | within one step |
|---|---|---|---|---|
| E0 | 128 | 192 | 1 | yes |
| E1q | 192 | 192 | 0 | yes |
| E6 | 192 | 256 | 1 | yes |
| Mb | >256 | 256 | 1 | yes |
| CL | 192 | 128 | 1 | yes |
| **E4c (out of sample)** | 768 | **>768** (power 0.79 at 768) | 1 | yes |
| **E4r (out of sample)** | 768 | **768** (power 0.92) | 0 | yes |
New 2D cells (aniso2_w16, 200 datasets each): E4c power 0.41 / 0.55 / 0.79 and E4r 0.49 / 0.71 / 0.92 at L = 384 / 512 / 768; predicted hw (J-law) vs measured: E4c 0.56/0.46/0.39 vs 0.56/0.51/0.43, E4r 0.48/0.38/0.32 vs 0.51/0.46/0.38 (the law under-predicts at large J by 10-20%). Coverage of exact d(t): E4c 0.95/0.98/0.95, E4r 0.93/0.96/0.92, **E6 0.96/0.90/0.645** (E6 starts to fail on plain 2D crossovers once L reaches 768).
**Plain verdict:** the recipe predicts the minimum size to within one grid step for every estimator tested, and it correctly said that the honest-bound estimators need L of about 768 (N = 12,288) in 2D, so honest certification of a 2D->1D change is expensive. Caveats: the grid is coarse (steps of 1.33-1.5x), so "within one step" is a lenient bar; for the first five the 2D answers were known before the rule was written; the E4c miss (0.79 vs 0.8) is a razor-edge call that I counted as ">768"; Mb's fitted delta is negative (two fit points), so the recipe never lets Mb detect, which happened to be one step from the truth.

## Item 2: better E6 misfit detector -- FAIL on both attempts; E6 stays a caveat
Original seeds, 37 cells x 200; 441 E6 misses on the misspecified cells (two-step, brush, gasket). E4c abstained where E6 did not on 4 datasets (no alarm).
| bar | A1 (E6 band disjoint from E4c band) | A2 (one- vs two-step on split A) | LoF (Amendment 5) |
|---|---|---|---|
| (a) P(alarm given miss), pooled >= 0.5 | **0.7% (0.1-2.0)** FAIL | **36.5% (32.0-41.2)** FAIL | 18.8% FAIL |
| (b) twostep_8_32 L=320, >= 0.8 | **1.6% (0.3-4.5)** FAIL | **56.0% (48.7-63.1)** FAIL | 30.6% FAIL |
| (c) control alarm rate <= 5% | 0.0% pass | 3.3% (2.6-4.2) pass | 3.7% pass |
| (d) realised-crossover power loss <= 5 points | 0.0 pass | 1.4 pass | 1.6 pass |
**Plain verdict:** E6 is NOT guarded by E4c, and the one-vs-two-step fallback does not rescue it either. A1 is useless because E4c's honest band is so wide that it always overlaps E6's: E6 can be wrong by more than its own half-width and still sit inside E4c. A2 is calibrated (3.3% on controls) and catches 56% of the worst two-step failures (L=320) but only 36% of misses overall. Per the rule there is no third detector; **E6's claim stays "honest through a smooth logistic crossover", with silent failure when two steps are separated by a ratio of about 2.5 or more (next section).** Guarding by A2 would cost 1.4 points of power for a 36% catch rate.

## Item 3: fresh-seed replication -- PASS (11 of 11 gated checks inside the original Clopper-Pearson intervals; E4c 11 of 11 as well)
Fresh seeds [22, cell, 1000+ds], 15 cells x 200. E6 simultaneous coverage: realised crossover pooled 99.0% (98.3-99.5) vs original 99.2; twostep_8_32 L = 128/192/256/320: 99.0 / 85.0 / 26.0 / 1.5 vs 99.5 / 83.0 / 31.0 / 3.5; twostep pooled 52.9 (49.3-56.4) vs 54.2; brush cells 95.0 / 95.5 / 92.5 vs 94.0 / 96.5 / 89.5, pooled 94.3 vs 93.3; E6 correct-sign power on the realised crossovers (abstention = miss) 55.6% vs 56.9%. Full list in `amend6_replication.csv`. The A1/A2 alarm results also replicate on fresh seeds (two-step plus brush: P(A1|miss) 0.0%, P(A2|miss) 39%; L=320: 0.0% and 59%).

**Step-separation scan** (twostep C_256 x C_8 x C_b, L = 256; separation ratio b/8; E6 simultaneous coverage, 200 datasets each, b = 32 pooled original + fresh = 398-400):
| b | ratio | E6 coverage (CP) | E4c coverage (CP) | A2 alarm rate |
|---|---|---|---|---|
| 12 | 1.5 | 95.5 (91.6-97.9) | 94.0 (89.8-96.9) | 1.5% |
| 16 | 2.0 | 92.5 (87.9-95.7) | 97.0 (93.6-98.9) | 3.0% |
| 20 | 2.5 | 42.0 (35.1-49.2) | 89.0 (83.8-93.0) | 11% |
| 24 | 3.0 | 12.0 (7.8-17.3) | 89.5 (84.4-93.4) | 21.5% |
| 32 | 4.0 | 28.5 (24.1-33.2) | 94.0 (91.2-96.1) | 27.5% |
**Plain verdict:** E6 coverage falls through 93% between a separation ratio of 1.5 and 2.0 (point estimates 95.5 then 92.5; the CP lower bound never reaches 93%, so it is "about 2", not a sharp line) and collapses between 2.0 and 2.5 (CP upper bound below 93% from ratio 2.5). Two steps closer than a factor of about 2 in the 2D/1D width (b/a) are safe; farther apart is not. New finding that tempers the earlier E4c summary: **E4c itself dips to 89% at ratios 2.5 and 3.0** (CP upper bounds 93.0 and 93.4, so it is at or just under 93% there); pooled over all two-step cells it stays at 95.2% (93.5-96.6).

## Limits
The recipe uses J as the carrier of walkers x volume, with walkers per node fixed at 32, so it has not been tested at other walker budgets. Item 1 does not test the abstention rule outside the J range seen in 3D (J >= 6). The separation scan holds a = 8 and L = 256 fixed and varies only b; other (a, L) combinations may shift the threshold. A2's threshold chi2_3 with the bootstrap correction is not an exact null for a gridded two-step search; its 3.3% control rate shows it is acceptable in practice, not proven. Gaspar-Kukal remains unimplemented. [Superseded 2026-11-23 by Amendment 12, below.] The earlier mid-run E4c from raw_amend3 used a different bootstrap stream; its re-run here (rngs[4]) agrees to within sampling noise (group coverages 95-99%).

## Amendment 7: walker budget (pre-registered 7339c90, erratum + predictions dc310b5; run on host-T, 61 cells x 200 datasets, 1207 s, NCPU 16)
Data: `results/raw_amend7.csv`, `cells_meta_amend7.csv`; scoring `analyze_amend7.py` -> `amend7_calib.csv`, `amend7_minL.csv`, `amend7_power.csv`, `amend7_scores.csv`.
**Analysis fix made after the data were seen (disclosed):** the first scoring pass dropped E0's abstained datasets (their declaration field is empty) instead of counting them as misses as the plan requires, and it scored the wpn = 32 baseline on sizes that do not exist for E0/E1q (384-768), giving spurious ">768". Both fixed; no verdict changed (P1 7/16, P2 13/16 before and after).

**(A) Calibration on controls -- PASS for every estimator at both budgets.** Pooled simultaneous coverage (CP 95%), FPR at margin 0.25:
| wpn | E0 | E1q | E4c | E4r | E6 |
|---|---|---|---|---|---|
| 8 | 94.9% (93.5-96.1), FPR 0.7% | 94.4% (92.9-95.6), 1.3% | 99.5%, 0% | 99.4%, 0% | 99.5%, 0.1% |
| 128 | 94.0% (92.5-95.3), 0% | 93.4% (91.9-94.8), 0.2% | 99.1%, 0% | 98.9%, 0% | 99.9%, 0% |
No control cell has CP-upper < 93% or CP-lower FPR > 5%. Abstention 25-29% (the small torus2 L = 24 and torus3 L = 16 cells, as before).

**(B) Recipe transfer -- FAIL.** P1 (frozen J-only recipe): 7 of 16 within one fine-grid step (bar 14). P2 (noise-scaled, (wpn/32)^-1/2): 13 of 16 -- also FAIL. Per the pre-registered rule the recipe is reported as **design-specific and not transferable** across walker budgets. P1 misses are all at wpn = 8 (under-predicts min L by 2-4 steps) plus E4c at wpn = 128 (over-predicts); P2's three misses: E1q 2D wpn 8 (predicted >768, measured 384), E4c 2D wpn 128 (256 vs 384), E1q 3D wpn 8 (112 vs 160).
Measured min L (power >= 0.8, abstention = miss; fine grids):
| family | est | wpn 8 | wpn 32 | wpn 128 |
|---|---|---|---|---|
| 2D->1D | E0 | 224 | 192 | 160 |
| | E1q | 384 | 192 | 160 |
| | E4c | >768 | >768 (0.79 at 768) | 384 |
| | E4r | >768 | 768 | 384 |
| | E6 | 512 | 224 | 160 |
| 3D->1D | E0 | 112 | 80 | 80 |
| | E1q | 160 | 96 | 80 |
| | E4c | >256 | 192 | 112 |
| | E4r | >256 | 160 | 112 |
| | E6 | 192 | 96 | 80 |
(wpn 32 = earlier runs on the old grid + 7 fill-in cells; on the fine grid E0 3D is 80, not the coarse-grid 96.) Mean band half-width on realised crossovers falls with budget roughly as wpn^-1/2 for E0 (0.39 -> 0.095 over 16x, ratio 4.1) and slightly slower for E4c (1.82 -> 0.52, ratio 3.5); E6 faster (1.75 -> 0.33). So J (window length) is not a sufficient carrier: half-width depends on walkers at fixed J, which is why P1 fails.

**(C) Ungated.** (1) Ordering holds at every budget: the honest-bound estimators (E4c/E4r) always need the largest graphs. (2) **Coverage of the exact d(t) on realised crossovers falls as walkers increase for every estimator that does not bound bias**: E1q 48.1% (wpn 8) / 20.2% (32, fill-in cells) / 3.5% (128); E0 82.2 / 74.7 / 53.5%; E6 99.1 / 97.5 / 84.4%; E4c 97.6 / 97.2 / 96.5%; E4r 96.3 / 95.6 / 95.0%. More data shrink the noise part of the band but not the smoothing bias, so un-corrected bands tighten around the wrong curve. Only the bias-aware intervals keep nominal coverage at every budget.
**Plain verdict:** calibration on constant-dimension controls does not depend on the walker budget; minimum sizes do, and the Amendment 6 recipe does not transfer (the paper must drop or heavily qualify it). The new and strongest evidence for the paper's headline: buying more data makes biased estimators *less* honest about the true d(t), while honest-bound intervals stay honest -- "no free coverage" holds across designs.
**Correction to Amendment 6 text:** "J is the carrier of walkers x volume" is wrong; J is set by volume except where the split-half count rule binds (wpn 8), and half-width depends on walkers at fixed J. See PLAN.md Amendment 7 and its erratum.

## Amendment 8: structure-using arms E7 (moment-cone LP) and E8 (secant bracket) (pre-registered aa84083; run on host-T, 47 cells × 200, 2312 s)
Files: `p22amend8.py`, `run_p22h.py`, `analyze_amend8.py`, `noisefree_e7.py`; `results/amend8_report.txt`, `amend8_cell_table.csv`, `amend8_noisefree.csv`, `raw_amend8_{core,w8,w128}.csv`. The walker data is bit-identical to the earlier runs (checked on 16 datasets).

| Arm | C (controls) | F (FPR) | X (exact d on crossovers) | M (min L vs E4c) | Verdict |
|---|---|---|---|---|---|
| E7 | 99.8% (CP 99.6–99.9) PASS | 0/2800 PASS | 98.8% (CP 98.4–99.1) PASS | 3D: 256 (bar ≤ 160) FAIL; 2D: power 0.00 up to L = 768 FAIL | **passes C, F, X but not M** |
| E8 | 99.7% PASS | 0/2800 PASS | 99.8% PASS | 3D: not reached (0.61 at 256) FAIL; 2D: ≤ 0.01 FAIL | **passes C, F, X but not M** |
| A1′ alarm (E6 vs E7 disjoint) | controls 0/1930 PASS | | P(alarm \| E6 miss) 0/441 FAIL | | **FAIL: never fires** (E7 band ±0.5–0.7 is wider than E6's errors) |

**Noise-free check:** with the exact q, the E7 band contains the exact d(t) in all 7 cells tested.

**E7 coverage loss** is entirely from 66 of 9400 datasets (0.7%) whose bootstrap box admits no completely monotone measure. Among feasible datasets, coverage is 2795/2796 on controls and 3556/3556 on crossovers.

**Atom gate:** in 170 datasets (1.8%) feasibility still flips between 799 and 1597 atoms. That is reported as an arm failure in those datasets.

At wpn = 128, E7 reaches 3D min L = 80 (vs E4c 112), but still not 2D (0.72 at 768).

**Caveat on M (important, stated for the paper).** The E7/E8 change test is interval arithmetic on the simultaneous band, Δ_j ∈ [lo_j − hi_0, hi_j − lo_0], which is conservative. Part of the zero 2D power may come from this test rather than the band. A joint LP bound on Δ = d(t_j) − d(t_0) over the same measure set is tighter. It was not pre-registered, so it was not run. It is a candidate Amendment 9.

**Plain verdict:** using the walk's guaranteed structure (complete monotonicity / log-convexity) gives honest coverage with no tuning constant at every budget, and E7 is much narrower than E4c on small graphs (ring 128: 0.84 vs 2.44). With this change test it does not buy smaller detectable sizes. As pre-registered: "even the exact constraint set does not buy size", qualified by the change-test caveat above.

## Amendment 9 (pre-registered 5b88a19; run on host-T, and on host-D for 9c; estimator code byte-identical to 1ba9ca6)
Files:
- code: `p22amend9.py`, `run_p22i.py`, `analyze_amend9.py`;
- verdicts: `results/amend9_report.txt`;
- tables: `results/amend9_b_table.csv`, `amend9_b_minL.csv`, `amend9_c_table.csv`, `amend9a_power.csv`.

### 9a. Joint change test for E7 (E7j), same data as Amendment 8
- **F:** 0/2800 = 0.0% → **PASS**.
- **M:** 3D min L = **160** (bar ≤ 160) → **PASS**. 2D: power 0.00 up to L = 768 → **FAIL**.
- **Verdict:** E7j does not beat the trade-off.
- **Ungated:** the joint test turns E7's 3D power from 0.06 to 0.61 at L = 128, and from 0.28 to 0.90 at L = 160. So in 3D a structure-only honest method reaches the size bar. In 2D the certified largest possible drop stays at about −0.2, never below −0.25, so the null there is not an artefact of the change test.
- **Implementation:** the "bisection" in the plan was implemented as an exact branch-and-bound over d(t₀) with one Charnes–Cooper LP per value. Tolerance 0.002, plus an atom-doubling gate; 1 of 740 gated decisions was withheld.

### 9b. Frozen confirmatory replication (fresh seeds, wpn 8/32/128): ALL 28 gated predictions PASS
1. E1q realised-crossover coverage falls with the walker budget: **49.6% → 18.3% → 4.0%** (wpn 8 → 32 → 128).
2. E4c, E4r, E7 and E8 controls are ≥ 93% at every budget (98.8–100.0%).
3. E1q at wpn 128 is ≤ 20% (4.0%).
4. Every C/F verdict from Amendments 1–8 reproduces at wpn 32 (14/14).

**Ungated:** E0's C fails at wpn 8 by the cell rule.

### 9c. Second design (time-grid ratio 2^(1/4), strict window): gate FAIL, 1 of 24 checks
- **The one difference:** E0's C is 92.98% (CP 91.7–94.1), below the 93.0% bar by 0.02 points. The original design gave 94.3%.
- Every other C/F/X verdict is unchanged, and the min-L ordering of E4c against E1q, E6 and E7 is unchanged (6/6, ties allowed).

**Plain reading:** the headline replicates on fresh seeds, all 28 gated predictions hold. The pre-registered second-design gate fails narrowly, on the baseline E0's control coverage, so the paper must report it as a failed robustness gate, not round it. The joint change test makes the structure-only band competitive in 3D, but not in 2D.

## Amendment 10b / 10c / 11 results (2026-11-22; pre-registered in PLAN.md Amendments 10 and 11; run on host-D; estimator files byte-identical to commit 1ba9ca6)
Labels: **10a, 10b and 10c are referee-directed reanalyses/tests** (external referee, 2026-11-22; the referee's own decision criteria are the bars). **Amendment 11 is a frozen holdout** (new graph families, run once, no second holdout). 10a's results are in `results/amend10a_report.txt` and are not repeated here.
Files: driver `run_amend1011.py` (notes `IMPLEMENTATION_NOTES_A10_A11.md`), verdict script `verdict_amend10_11.py`, copied result tables and manifests in `results_amend10_11/` (raw per-dataset files for 10c and 11, 14.8 MB and 12.6 MB, stay on host-D `<local path omitted>`; their sha256 are in the manifests). Manifests verified: all four `complete=True`, no errors, task counts full, csv sha256 re-computed on the copied files (10b, exact-q, 10c, 11), per-item pickle sha256 lists re-checked on host-D (1030, 7, 2240, 1920 items, 0 bad, 0 missing), and the items-of-sha256 digest recomputed (4/4 match). `verdict_report.txt` carries the full printout; `reading_log.txt` the readings below.
**Readings applied where a bar needed one (nothing new was added as a bar; stricter alternatives are flagged in the gate table):** (R1) 10c and 11 coverage/FPR are over datasets where the estimator answered, exactly as every earlier amendment's C/F/X definition; the strict reading (abstention = miss) is shown beside each and flips G3 and G4 (and 10c-B) to FAIL because the frozen J<6 rule forces E1q/E4 to abstain in cells with an exact window under 6 points (C_48 x C_6 x C_6 exact_J = 3, torus3 L=18 exact_J = 4, 3D L=48/64 in 10c). 10b uses fixed denominators as the referee asked. (R2) Plain numeric bars are tested on the point estimate, as every earlier amendment's bars were, with CP intervals reported; the stricter reading 'the CP bound must also meet the bar' is flagged wherever it would change a verdict, and the plan's own CP rules are used as written. (R3) cells with zero answered datasets are skipped in per-cell tests. (R4) FPR = any-sign declaration at margin 0.25; extra reading: no cell with CP-lower > 5%. (R5) 10b: both outer relaxations A and B must pass for E7; E7j exists for B only. (R6) 10c: bars must hold pooled and within each family; (R7) 'decreases' = strictly.
### Gate table
| gate | bar (as pre-registered) | observed (CP 95%) | verdict |
|---|---|---|---|
| 10b-C[B] | E7 certified-band coverage of exact d on controls >= 93% pooled, no cell with CP-upper < 93% | 112/112 = 100.0% (96.8-100.0); bad cells [] | **PASS** |
| 10b-C[A] | E7 certified-band coverage of exact d on controls >= 93% pooled, no cell with CP-upper < 93% | 112/112 = 100.0% (96.8-100.0); bad cells [] | **PASS** |
| 10b-X[B] | E7 certified-band coverage of exact d on realised crossovers >= 90% | 760/761 = 99.9% (99.3-100.0) | **PASS** |
| 10b-X[A] | E7 certified-band coverage of exact d on realised crossovers >= 90% | 760/761 = 99.9% (99.3-100.0) | **PASS** |
| 10b-F[B] | E7 interval test, B: false-change rate (margin 0.25, any sign) on controls <= 5%, no cell with CP-lower > 5% | 0/112 = 0.0% (0.0-3.2); bad cells [] | **PASS** |
| 10b-F[A] | E7 interval test, A: false-change rate (margin 0.25, any sign) on controls <= 5%, no cell with CP-lower > 5% | 0/112 = 0.0% (0.0-3.2); bad cells [] | **PASS** |
| 10b-F[j] | E7j certified (B): false-change rate (margin 0.25, any sign) on controls <= 5%, no cell with CP-lower > 5% | 0/112 = 0.0% (0.0-3.2); bad cells [] | **PASS** |
| 10b-M[E7j] | 3D min L = 160 preserved: certified E7j power >= 0.8 at L = 160, 192, 256 (fixed n = 200) | L=160 170/200 = 85.0% (79.3-89.6); L=192 186/200 = 93.0% (88.5-96.1); L=256 198/200 = 99.0% (96.4-99.9) ; CI-strict reading (CP bound must also meet the bar) gives FAIL | **PASS** |
| 10b-exactq | 7 exact-q cells: certified bands (A and B) contain the exact d(t) (7/7) and no false change on the 3 control cells | coverage B/A 7/7, 7/7; control false changes (E7 B, E7j) 0, 0 | **PASS** |
| 10c-A | E1q coverage still falls with wpn under FIXED windows (strictly monotone 8>32>128; wpn-128 CP-upper < wpn-8 CP-lower), pooled and in each family | pooled: 40.4% > 19.6% > 3.6%; CP-upper(128) 4.3% vs CP-lower(8) 38.7% / aniso2_w16: 38.3% > 16.3% > 2.2%; CP-upper(128) 3.0% vs CP-lower(8) 36.0% / aniso3_w8: 43.1% > 23.8% > 5.4%; CP-upper(128) 6.7% vs CP-lower(8) 40.5% ; strict reading (abstention = miss), pooled: holds | **PASS** |
| 10c-B | E4c coverage >= 93% on the matched (fixed) windows, at each wpn, pooled and in each family | pooled wpn8: 3123/3198 = 97.7% (97.1-98.2) / pooled wpn32: 3101/3198 = 97.0% (96.3-97.5) / pooled wpn128: 3070/3191 = 96.2% (95.5-96.8) / aniso2_w16 wpn8: 1751/1799 = 97.3% (96.5-98.0) / aniso2_w16 wpn32: 1742/1798 = 96.9% (96.0-97.6) / aniso2_w16 wpn128: 1737/1791 = 97.0% (96.1-97.7) / aniso3_w8 wpn8: 1372/1399 = 98.1% (97.2-98.7) / aniso3_w8 wpn32: 1359/1400 = 97.1% (96.0-97.9) / aniso3_w8 wpn128: 1333/1400 = 95.2% (94.0-96.3) ; strict reading (abstention = miss) would give FAIL: pooled wpn8: 86.8%; pooled wpn32: 86.1%; pooled wpn128: 85.3%... (forced abstention: exact window < 6 points in 3D L=48,64) | **PASS** |
| 10c-CENTRAL | referee bar: central result confirmed iff 10c-A AND 10c-B | 10c-A PASS, 10c-B PASS | **PASS** |
| 10c-DIAG | E4 bound violations documented (violation rate and coverage conditional on violation, incl. separation-scan ratios 2.5 and 3.0); the referee gives no pass/fail bar for this | see table ; no bar | **REPORTED** |
| G1 | E1q realised-crossover coverage at wpn 32 <= 40% | 417/1597 = 26.1% (24.0-28.3); strict (abstention=miss) 20.9% ; strict reading is lower, so the answered-dataset reading used is the less favourable one | **PASS** |
| G2 | E1q coverage strictly decreases 8 -> 32 -> 128 and CP-upper(128) < CP-lower(8) | 761/1327 = 57.3% (54.6-60.0) > 417/1597 = 26.1% (24.0-28.3) > 161/1605 = 10.0% (8.6-11.6); CP-upper(128) 11.6% < CP-lower(8) 54.6% | **PASS** |
| G3 | E4c and E4r realised-crossover coverage >= 90% at every wpn | E4c wpn8 1310/1327 = 98.7% (98.0-99.3) / E4c wpn32 1548/1596 = 97.0% (96.0-97.8) / E4c wpn128 1570/1605 = 97.8% (97.0-98.5) / E4r wpn8 1299/1327 = 97.9% (97.0-98.6) / E4r wpn32 1524/1597 = 95.4% (94.3-96.4) / E4r wpn128 1532/1605 = 95.5% (94.3-96.4) ; strict reading (abstention = miss) would give FAIL: E4c wpn8 65.5%; E4c wpn32 77.4%; E4c wpn128 78.5%... (E1q/E4 abstain by the frozen J<6 rule, e.g. C_48 x C_6 x C_6 has exact_J = 3) | **PASS** |
| G4 | control coverage >= 93% pooled, no cell with CP-upper < 93%, for E4c, E4r, E7, E8 at every wpn | E4c wpn8 886/890 = 99.6% (98.9-99.9) / E4c wpn32 993/998 = 99.5% (98.8-99.8) / E4c wpn128 989/1000 = 98.9% (98.0-99.4) / E4r wpn8 887/890 = 99.7% (99.0-99.9) / E4r wpn32 990/998 = 99.2% (98.4-99.7) / E4r wpn128 990/1000 = 99.0% (98.2-99.5) / E7 wpn8 1199/1200 = 99.9% (99.5-100.0) / E7 wpn32 1196/1200 = 99.7% (99.1-99.9) / E7 wpn128 1187/1200 = 98.9% (98.2-99.4) / E8 wpn8 1193/1200 = 99.4% (98.8-99.8) / E8 wpn32 1196/1200 = 99.7% (99.1-99.9) / E8 wpn128 1199/1200 = 99.9% (99.5-100.0) ; cells with 0 answered datasets skipped (R3): ['torus3_18']; strict reading would give FAIL for 6 est/wpn combos (forced abstention in torus3 L=18, exact_J=4) | **PASS** |
| G5 | control FPR (margin 0.25, any sign) <= 5% for all five arms at wpn 32 | E1q 3/998 = 0.3% (0.1-0.9) / E4c 0/998 = 0.0% (0.0-0.4) / E4r 0/998 = 0.0% (0.0-0.4) / E7 0/1200 = 0.0% (0.0-0.3) / E8 0/1200 = 0.0% (0.0-0.3) | **PASS** |
| 11-VERDICT | the narrowed headline GENERALISES iff G1-G5 all pass | G1 PASS, G2 PASS, G3 PASS, G4 PASS, G5 PASS | **PASS** |

### Amendment 10b (referee-directed test): certified outer bounds for E7/E7j
1030 datasets (170 atom-sensitive + 1 withheld + 300 seeded stable + 600 from the 3D cells L = 160, 192, 256, overlaps merged) plus the 7 exact-q cells. B = tight outer relaxation (headline); A = plan-literal nearest-atom relaxation. LP is HiGHS double precision, not a verified solve (stated in `p22amend10.py`). All-dataset outcome counts (fixed denominator 1030): B-band covers exact d in 1029, fails to cover in 0, no certified band in 1; original gated E7 covers in 1024, has no band in 6; certified band contains the original band in 1024/1024; E7j declares a drop in 570, a rise in 12, is undecided (certified-undetermined at some time, no declaration) in 44, no declaration otherwise 403.

Cell-level outcome categories (fixed denominator = n datasets in the cell; 'no band' counts as a miss everywhere). E7 columns: B/A = certified bands, orig = original gated band; 'change' = interval-test change declaration (margin 0.25).

| cell | n | realised | E7 covered, B | E7 covered, A | E7 covered, orig | B not covered | no band (B/orig) | E7 change (B/A/orig) | E7j drop/rise/undecided/none |
|---|---|---|---|---|---|---|---|---|---|
| brush brush_f16 L=12 | 3 | 0 | 3/3 (100%) | 3/3 (100%) | 3/3 (100%) | 0 | 0/0 | 3/0/3 | 0/3/0/0 |
| brush brush_f16 L=16 | 4 | 0 | 4/4 (100%) | 4/4 (100%) | 4/4 (100%) | 0 | 0/0 | 4/0/4 | 0/4/0/0 |
| brush brush_f24 L=8 | 9 | 0 | 9/9 (100%) | 9/9 (100%) | 9/9 (100%) | 0 | 0/0 | 5/0/7 | 0/5/4/0 |
| control ring L=128 | 6 | 0 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/6 |
| control ring L=256 | 9 | 0 | 9/9 (100%) | 9/9 (100%) | 9/9 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/9 |
| control ring L=512 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/7 |
| control ring L=1024 | 12 | 0 | 12/12 (100%) | 12/12 (100%) | 12/12 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/12 |
| control torus2 L=24 | 6 | 0 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/6 |
| control torus2 L=32 | 9 | 0 | 9/9 (100%) | 9/9 (100%) | 9/9 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/9 |
| control torus2 L=48 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/7 |
| control torus2 L=64 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/7 |
| control torus2 L=96 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 6/7 (86%) | 0 | 0/1 | 0/0/0 | 0/0/0/7 |
| control torus3 L=16 | 8 | 0 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| control torus3 L=20 | 8 | 0 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| control torus3 L=24 | 4 | 0 | 4/4 (100%) | 4/4 (100%) | 4/4 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/4 |
| control torus3 L=28 | 8 | 0 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| control torus3 L=32 | 14 | 0 | 14/14 (100%) | 14/14 (100%) | 14/14 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/14 |
| crossover aniso2_w16 L=32 | 4 | 0 | 4/4 (100%) | 4/4 (100%) | 4/4 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/4 |
| crossover aniso2_w16 L=48 | 6 | 0 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/6 |
| crossover aniso2_w16 L=64 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/7 |
| crossover aniso2_w16 L=96 | 7 | 0 | 7/7 (100%) | 7/7 (100%) | 7/7 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/7 |
| crossover aniso2_w16 L=128 | 8 | 1 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| crossover aniso2_w16 L=160 | 6 | 1 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/6 |
| crossover aniso2_w16 L=192 | 6 | 1 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/6 |
| crossover aniso2_w16 L=224 | 9 | 1 | 9/9 (100%) | 9/9 (100%) | 9/9 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/9 |
| crossover aniso2_w16 L=256 | 12 | 1 | 12/12 (100%) | 12/12 (100%) | 12/12 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/12 |
| crossover aniso2_w16 L=320 | 17 | 1 | 17/17 (100%) | 17/17 (100%) | 17/17 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/17 |
| crossover aniso2_w16 L=384 | 16 | 1 | 16/16 (100%) | 16/16 (100%) | 15/16 (94%) | 0 | 0/1 | 0/0/0 | 0/0/0/16 |
| crossover aniso2_w16 L=512 | 15 | 1 | 15/15 (100%) | 15/15 (100%) | 15/15 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/15 |
| crossover aniso2_w16 L=768 | 14 | 1 | 14/14 (100%) | 14/14 (100%) | 13/14 (93%) | 0 | 0/1 | 0/0/0 | 0/0/0/14 |
| crossover aniso3_w8 L=32 | 8 | 0 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| crossover aniso3_w8 L=48 | 9 | 1 | 9/9 (100%) | 9/9 (100%) | 9/9 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/9 |
| crossover aniso3_w8 L=64 | 12 | 1 | 12/12 (100%) | 12/12 (100%) | 12/12 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/12 |
| crossover aniso3_w8 L=80 | 4 | 1 | 4/4 (100%) | 4/4 (100%) | 4/4 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/4 |
| crossover aniso3_w8 L=96 | 13 | 1 | 13/13 (100%) | 13/13 (100%) | 13/13 (100%) | 0 | 0/0 | 0/0/0 | 2/0/2/9 |
| crossover aniso3_w8 L=112 | 6 | 1 | 6/6 (100%) | 6/6 (100%) | 6/6 (100%) | 0 | 0/0 | 0/0/0 | 2/0/2/2 |
| crossover aniso3_w8 L=128 | 14 | 1 | 14/14 (100%) | 14/14 (100%) | 14/14 (100%) | 0 | 0/0 | 0/0/0 | 7/0/3/4 |
| crossover aniso3_w8 L=160 | 200 | 1 | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) | 0 | 0/0 | 35/0/54 | 170/0/14/16 |
| crossover aniso3_w8 L=192 | 200 | 1 | 199/200 (100%) | 199/200 (100%) | 197/200 (98%) | 0 | 1/3 | 46/0/86 | 186/0/7/6 |
| crossover aniso3_w8 L=256 | 200 | 1 | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) | 0 | 0/0 | 147/0/177 | 198/0/2/0 |
| gasket gasket L=5 | 5 | 0 | 5/5 (100%) | 5/5 (100%) | 5/5 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/5 |
| gasket gasket L=6 | 8 | 0 | 8/8 (100%) | 8/8 (100%) | 8/8 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/8 |
| gasket gasket L=7 | 10 | 0 | 10/10 (100%) | 10/10 (100%) | 10/10 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/10 |
| twostep twostep_8_32 L=128 | 28 | 0 | 28/28 (100%) | 28/28 (100%) | 28/28 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/28 |
| twostep twostep_8_32 L=192 | 18 | 0 | 18/18 (100%) | 18/18 (100%) | 18/18 (100%) | 0 | 0/0 | 0/0/0 | 0/0/0/18 |
| twostep twostep_8_32 L=256 | 22 | 0 | 22/22 (100%) | 22/22 (100%) | 22/22 (100%) | 0 | 0/0 | 2/0/5 | 2/0/4/16 |
| twostep twostep_8_32 L=320 | 18 | 0 | 18/18 (100%) | 18/18 (100%) | 18/18 (100%) | 0 | 0/0 | 3/0/8 | 3/0/6/9 |

3D (aniso3_w8) certified power by L, fixed denominator (E7j = certified joint-test drop; E7 = interval test on B; orig = original gated E7 interval test):

| L | n | E7j (B) | E7 interval (B) | E7 interval (orig) |
|---|---|---|---|---|
| 32 | 8 | 0/8 (0%) | 0/8 (0%) | 0/8 (0%) |
| 48 | 9 | 0/9 (0%) | 0/9 (0%) | 0/9 (0%) |
| 64 | 12 | 0/12 (0%) | 0/12 (0%) | 0/12 (0%) |
| 80 | 4 | 0/4 (0%) | 0/4 (0%) | 0/4 (0%) |
| 96 | 13 | 2/13 (15%) | 0/13 (0%) | 0/13 (0%) |
| 112 | 6 | 2/6 (33%) | 0/6 (0%) | 0/6 (0%) |
| 128 | 14 | 7/14 (50%) | 0/14 (0%) | 0/14 (0%) |
| 160 | 200 | 170/200 (85%) | 35/200 (18%) | 54/200 (27%) |
| 192 | 200 | 186/200 (93%) | 46/200 (23%) | 86/200 (43%) |
| 256 | 200 | 198/200 (99%) | 147/200 (74%) | 177/200 (88%) |

Exact-q cells (noise-free q; certified bands vs exact d):

| cell | J | N | B_status | B_cov | A_cov | B_contains_orig | B_e7_neg | B_e7_pos | Bj_neg | Bj_pos | Bj_amb | true_drop | true_rise |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ring_256 | 10 | 256 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0.009954849434608226 |
| torus2_24 | 3 | 576 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | -0.0034127040670899333 | 0.0 |
| torus3_16 | 4 | 4096 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | -0.04777958408304528 | 0.0 |
| aniso2_w16_128 | 8 | 2048 | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 0 | 1 | -0.8933336360709865 | 0.0 |
| aniso3_w8_32 | 3 | 2048 | 0 | 1 | 1 | 1 | 0 | 0 | 1 | 0 | 0 | -0.5725041036076042 | 0.0 |
| aniso2_w16_256 | 12 | 4096 | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 0 | 1 | -1.0041429394301133 | 0.0 |
| aniso3_w8_96 | 7 | 6144 | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 0 | 0 | -2.0126046310571333 | 0.0 |

### Amendment 10c (referee-directed test): fixed-window budget comparison and E4 bound diagnosis
Cells: the 18 realised crossover cells of 9b (aniso2_w16 L=128..768, aniso3_w8 L=48..256), wpn 8/32/128, 200 datasets, fresh seeds [22, id, 3000+ds]; fixed J = the wpn-32 exact window (grid points with exact N q >= 6), recorded in `amend10c_windows.json` before data. 'orig' = the estimator's own window rule on the same data. Pooled coverage over answered datasets (simultaneous band); 'strict' counts abstention as a miss. E1q has no M.

| est | window | family | wpn 8 | wpn 32 | wpn 128 |
|---|---|---|---|---|---|
| E1q | fixed | aniso2_w16 | 689/1800 = 38.3% [36.0-40.6]; strict 38.3%; answered 1800/1800 | 294/1800 = 16.3% [14.7-18.1]; strict 16.3%; answered 1800/1800 | 39/1800 = 2.2% [1.5-3.0]; strict 2.2%; answered 1800/1800 |
| E1q | fixed | aniso3_w8 | 603/1400 = 43.1% [40.5-45.7]; strict 33.5%; answered 1400/1800 | 333/1400 = 23.8% [21.6-26.1]; strict 18.5%; answered 1400/1800 | 75/1400 = 5.4% [4.2-6.7]; strict 4.2%; answered 1400/1800 |
| E1q | fixed | both | 1292/3200 = 40.4% [38.7-42.1]; strict 35.9%; answered 3200/3600 | 627/3200 = 19.6% [18.2-21.0]; strict 17.4%; answered 3200/3600 | 114/3200 = 3.6% [2.9-4.3]; strict 3.2%; answered 3200/3600 |
| E1q | orig | aniso2_w16 | 845/1788 = 47.3% [44.9-49.6]; strict 46.9%; answered 1788/1800 | 307/1800 = 17.1% [15.3-18.9]; strict 17.1%; answered 1800/1800 | 48/1800 = 2.7% [2.0-3.5]; strict 2.7%; answered 1800/1800 |
| E1q | orig | aniso3_w8 | 550/1127 = 48.8% [45.8-51.8]; strict 30.6%; answered 1127/1800 | 306/1397 = 21.9% [19.8-24.2]; strict 17.0%; answered 1397/1800 | 72/1405 = 5.1% [4.0-6.4]; strict 4.0%; answered 1405/1800 |
| E1q | orig | both | 1395/2915 = 47.9% [46.0-49.7]; strict 38.8%; answered 2915/3600 | 613/3197 = 19.2% [17.8-20.6]; strict 17.0%; answered 3197/3600 | 120/3205 = 3.7% [3.1-4.5]; strict 3.3%; answered 3205/3600 |
| E4c | fixed | aniso2_w16 | 1751/1799 = 97.3% [96.5-98.0]; strict 97.3%; answered 1799/1800 | 1742/1798 = 96.9% [96.0-97.6]; strict 96.8%; answered 1798/1800 | 1737/1791 = 97.0% [96.1-97.7]; strict 96.5%; answered 1791/1800 |
| E4c | fixed | aniso3_w8 | 1372/1399 = 98.1% [97.2-98.7]; strict 76.2%; answered 1399/1800 | 1359/1400 = 97.1% [96.0-97.9]; strict 75.5%; answered 1400/1800 | 1333/1400 = 95.2% [94.0-96.3]; strict 74.1%; answered 1400/1800 |
| E4c | fixed | both | 3123/3198 = 97.7% [97.1-98.2]; strict 86.8%; answered 3198/3600 | 3101/3198 = 97.0% [96.3-97.5]; strict 86.1%; answered 3198/3600 | 3070/3191 = 96.2% [95.5-96.8]; strict 85.3%; answered 3191/3600 |
| E4c | orig | aniso2_w16 | 1751/1787 = 98.0% [97.2-98.6]; strict 97.3%; answered 1787/1800 | 1736/1799 = 96.5% [95.5-97.3]; strict 96.4%; answered 1799/1800 | 1747/1798 = 97.2% [96.3-97.9]; strict 97.1%; answered 1798/1800 |
| E4c | orig | aniso3_w8 | 1112/1127 = 98.7% [97.8-99.3]; strict 61.8%; answered 1127/1800 | 1355/1397 = 97.0% [96.0-97.8]; strict 75.3%; answered 1397/1800 | 1338/1405 = 95.2% [94.0-96.3]; strict 74.3%; answered 1405/1800 |
| E4c | orig | both | 2863/2914 = 98.2% [97.7-98.7]; strict 79.5%; answered 2914/3600 | 3091/3196 = 96.7% [96.0-97.3]; strict 85.9%; answered 3196/3600 | 3085/3203 = 96.3% [95.6-96.9]; strict 85.7%; answered 3203/3600 |
| E4r | fixed | aniso2_w16 | 1720/1800 = 95.6% [94.5-96.5]; strict 95.6%; answered 1800/1800 | 1713/1800 = 95.2% [94.1-96.1]; strict 95.2%; answered 1800/1800 | 1733/1796 = 96.5% [95.5-97.3]; strict 96.3%; answered 1796/1800 |
| E4r | fixed | aniso3_w8 | 1362/1400 = 97.3% [96.3-98.1]; strict 75.7%; answered 1400/1800 | 1348/1400 = 96.3% [95.2-97.2]; strict 74.9%; answered 1400/1800 | 1310/1400 = 93.6% [92.2-94.8]; strict 72.8%; answered 1400/1800 |
| E4r | fixed | both | 3082/3200 = 96.3% [95.6-96.9]; strict 85.6%; answered 3200/3600 | 3061/3200 = 95.7% [94.9-96.3]; strict 85.0%; answered 3200/3600 | 3043/3196 = 95.2% [94.4-95.9]; strict 84.5%; answered 3196/3600 |
| E4r | orig | aniso2_w16 | 1721/1788 = 96.3% [95.3-97.1]; strict 95.6%; answered 1788/1800 | 1700/1800 = 94.4% [93.3-95.5]; strict 94.4%; answered 1800/1800 | 1731/1797 = 96.3% [95.4-97.1]; strict 96.2%; answered 1797/1800 |
| E4r | orig | aniso3_w8 | 1105/1127 = 98.0% [97.1-98.8]; strict 61.4%; answered 1127/1800 | 1335/1397 = 95.6% [94.3-96.6]; strict 74.2%; answered 1397/1800 | 1304/1405 = 92.8% [91.3-94.1]; strict 72.4%; answered 1405/1800 |
| E4r | orig | both | 2826/2915 = 96.9% [96.3-97.5]; strict 78.5%; answered 2915/3600 | 3035/3197 = 94.9% [94.1-95.7]; strict 84.3%; answered 3197/3600 | 3035/3202 = 94.8% [94.0-95.5]; strict 84.3%; answered 3202/3600 |

Cell-level: fixed window vs original window (coverage k/answered; J_fixed = exact window, grid points; cells with J_fixed < 6 are abstentions by the frozen rule):

| family | L | wpn | J_fixed | E1q fixed vs orig | E4c fixed vs orig | E4r fixed vs orig |
|---|---|---|---|---|---|---|
| aniso2_w16 | 128 | 8 | 8.0 | 168/200 (84%) vs 171/188 (91%) | 196/200 (98%) vs 187/188 (99%) | 197/200 (98%) vs 187/188 (99%) |
| aniso2_w16 | 128 | 32 | 8.0 | 118/200 (59%) vs 127/200 (64%) | 197/200 (98%) vs 200/200 (100%) | 197/200 (98%) vs 199/200 (100%) |
| aniso2_w16 | 128 | 128 | 8.0 | 34/200 (17%) vs 42/200 (21%) | 198/200 (99%) vs 198/200 (99%) | 196/200 (98%) vs 197/200 (98%) |
| aniso2_w16 | 160 | 8 | 9.0 | 158/200 (79%) vs 177/200 (88%) | 200/200 (100%) vs 198/200 (99%) | 200/200 (100%) vs 198/200 (99%) |
| aniso2_w16 | 160 | 32 | 9.0 | 64/200 (32%) vs 69/200 (34%) | 198/200 (99%) vs 197/200 (98%) | 196/200 (98%) vs 194/200 (97%) |
| aniso2_w16 | 160 | 128 | 9.0 | 5/200 (2%) vs 6/200 (3%) | 197/200 (98%) vs 196/200 (98%) | 197/200 (98%) vs 196/200 (98%) |
| aniso2_w16 | 192 | 8 | 10.0 | 129/200 (64%) vs 154/200 (77%) | 200/200 (100%) vs 200/200 (100%) | 198/200 (99%) vs 198/200 (99%) |
| aniso2_w16 | 192 | 32 | 10.0 | 43/200 (22%) vs 45/200 (22%) | 197/200 (98%) vs 193/200 (96%) | 189/200 (94%) vs 190/200 (95%) |
| aniso2_w16 | 192 | 128 | 10.0 | 0/200 (0%) vs 0/200 (0%) | 192/200 (96%) vs 193/200 (96%) | 191/200 (96%) vs 190/200 (95%) |
| aniso2_w16 | 224 | 8 | 11.0 | 83/200 (42%) vs 129/200 (64%) | 196/200 (98%) vs 198/200 (99%) | 194/200 (97%) vs 197/200 (98%) |
| aniso2_w16 | 224 | 32 | 11.0 | 17/200 (8%) vs 28/200 (14%) | 196/200 (98%) vs 193/200 (96%) | 196/200 (98%) vs 191/200 (96%) |
| aniso2_w16 | 224 | 128 | 11.0 | 0/200 (0%) vs 0/200 (0%) | 182/200 (91%) vs 185/200 (92%) | 184/200 (92%) vs 183/200 (92%) |
| aniso2_w16 | 256 | 8 | 12.0 | 39/200 (20%) vs 84/200 (42%) | 198/200 (99%) vs 198/200 (99%) | 195/200 (98%) vs 196/200 (98%) |
| aniso2_w16 | 256 | 32 | 12.0 | 27/200 (14%) vs 19/200 (10%) | 188/200 (94%) vs 190/200 (95%) | 187/200 (94%) vs 186/200 (93%) |
| aniso2_w16 | 256 | 128 | 12.0 | 0/200 (0%) vs 0/200 (0%) | 197/200 (98%) vs 193/199 (97%) | 193/200 (96%) vs 191/200 (96%) |
| aniso2_w16 | 320 | 8 | 13.0 | 30/200 (15%) vs 49/200 (24%) | 193/200 (96%) vs 195/200 (98%) | 187/200 (94%) vs 191/200 (96%) |
| aniso2_w16 | 320 | 32 | 13.0 | 13/200 (6%) vs 7/200 (4%) | 186/200 (93%) vs 188/200 (94%) | 179/200 (90%) vs 180/200 (90%) |
| aniso2_w16 | 320 | 128 | 13.0 | 0/200 (0%) vs 0/200 (0%) | 196/199 (98%) vs 197/200 (98%) | 197/200 (98%) vs 196/200 (98%) |
| aniso2_w16 | 384 | 8 | 14.0 | 36/200 (18%) vs 42/200 (21%) | 191/200 (96%) vs 193/200 (96%) | 187/200 (94%) vs 191/200 (96%) |
| aniso2_w16 | 384 | 32 | 14.0 | 12/200 (6%) vs 12/200 (6%) | 195/200 (98%) vs 192/200 (96%) | 190/200 (95%) vs 187/200 (94%) |
| aniso2_w16 | 384 | 128 | 14.0 | 0/200 (0%) vs 0/200 (0%) | 196/199 (98%) vs 197/200 (98%) | 197/200 (98%) vs 197/200 (98%) |
| aniso2_w16 | 512 | 8 | 16.0 | 25/200 (12%) vs 22/200 (11%) | 187/200 (94%) vs 191/200 (96%) | 182/200 (91%) vs 181/200 (90%) |
| aniso2_w16 | 512 | 32 | 16.0 | 0/200 (0%) vs 0/200 (0%) | 195/200 (98%) vs 195/200 (98%) | 194/200 (97%) vs 191/200 (96%) |
| aniso2_w16 | 512 | 128 | 16.0 | 0/200 (0%) vs 0/200 (0%) | 197/198 (99%) vs 198/200 (99%) | 198/198 (100%) vs 200/200 (100%) |
| aniso2_w16 | 768 | 8 | 18.0 | 21/200 (10%) vs 17/200 (8%) | 190/199 (95%) vs 191/199 (96%) | 180/200 (90%) vs 182/200 (91%) |
| aniso2_w16 | 768 | 32 | 18.0 | 0/200 (0%) vs 0/200 (0%) | 190/198 (96%) vs 188/199 (94%) | 185/200 (92%) vs 182/200 (91%) |
| aniso2_w16 | 768 | 128 | 18.0 | 0/200 (0%) vs 0/200 (0%) | 182/195 (93%) vs 190/199 (95%) | 180/198 (91%) vs 181/197 (92%) |
| aniso3_w8 | 48 | 8 | 4.0 | n/a vs n/a | n/a vs n/a | n/a vs n/a |
| aniso3_w8 | 48 | 32 | 4.0 | n/a vs n/a | n/a vs n/a | n/a vs n/a |
| aniso3_w8 | 48 | 128 | 4.0 | n/a vs n/a | n/a vs n/a | n/a vs n/a |
| aniso3_w8 | 64 | 8 | 5.0 | n/a vs 2/3 (67%) | n/a vs 3/3 (100%) | n/a vs 3/3 (100%) |
| aniso3_w8 | 64 | 32 | 5.0 | n/a vs 14/20 (70%) | n/a vs 20/20 (100%) | n/a vs 20/20 (100%) |
| aniso3_w8 | 64 | 128 | 5.0 | n/a vs 4/8 (50%) | n/a vs 8/8 (100%) | n/a vs 8/8 (100%) |
| aniso3_w8 | 80 | 8 | 6.0 | 180/200 (90%) vs 38/40 (95%) | 199/200 (100%) vs 40/40 (100%) | 199/200 (100%) vs 40/40 (100%) |
| aniso3_w8 | 80 | 32 | 6.0 | 166/200 (83%) vs 135/178 (76%) | 198/200 (99%) vs 176/178 (99%) | 199/200 (100%) vs 176/178 (99%) |
| aniso3_w8 | 80 | 128 | 6.0 | 68/200 (34%) vs 61/197 (31%) | 198/200 (99%) vs 195/197 (99%) | 199/200 (100%) vs 195/197 (99%) |
| aniso3_w8 | 96 | 8 | 7.0 | 175/200 (88%) vs 93/104 (89%) | 199/200 (100%) vs 103/104 (99%) | 199/200 (100%) vs 103/104 (99%) |
| aniso3_w8 | 96 | 32 | 7.0 | 130/200 (65%) vs 110/199 (55%) | 200/200 (100%) vs 199/199 (100%) | 199/200 (100%) vs 198/199 (99%) |
| aniso3_w8 | 96 | 128 | 7.0 | 7/200 (4%) vs 6/200 (3%) | 195/200 (98%) vs 195/200 (98%) | 193/200 (96%) vs 192/200 (96%) |
| aniso3_w8 | 112 | 8 | 8.0 | 134/200 (67%) vs 152/182 (84%) | 199/200 (100%) vs 181/182 (99%) | 199/200 (100%) vs 181/182 (99%) |
| aniso3_w8 | 112 | 32 | 8.0 | 32/200 (16%) vs 42/200 (21%) | 199/200 (100%) vs 199/200 (100%) | 196/200 (98%) vs 196/200 (98%) |
| aniso3_w8 | 112 | 128 | 8.0 | 0/200 (0%) vs 1/200 (0%) | 196/200 (98%) vs 193/200 (96%) | 192/200 (96%) vs 189/200 (94%) |
| aniso3_w8 | 128 | 8 | 9.0 | 56/200 (28%) vs 141/198 (71%) | 196/200 (98%) vs 196/198 (99%) | 197/200 (98%) vs 196/198 (99%) |
| aniso3_w8 | 128 | 32 | 9.0 | 5/200 (2%) vs 5/200 (2%) | 196/200 (98%) vs 192/200 (96%) | 196/200 (98%) vs 193/200 (96%) |
| aniso3_w8 | 128 | 128 | 9.0 | 0/200 (0%) vs 0/200 (0%) | 199/200 (100%) vs 198/200 (99%) | 195/200 (98%) vs 197/200 (98%) |
| aniso3_w8 | 160 | 8 | 10.0 | 28/200 (14%) vs 80/200 (40%) | 193/200 (96%) vs 196/200 (98%) | 189/200 (94%) vs 193/200 (96%) |
| aniso3_w8 | 160 | 32 | 10.0 | 0/200 (0%) vs 0/200 (0%) | 192/200 (96%) vs 193/200 (96%) | 192/200 (96%) vs 190/200 (95%) |
| aniso3_w8 | 160 | 128 | 10.0 | 0/200 (0%) vs 0/200 (0%) | 193/200 (96%) vs 192/200 (96%) | 195/200 (98%) vs 190/200 (95%) |
| aniso3_w8 | 192 | 8 | 11.0 | 18/200 (9%) vs 35/200 (18%) | 192/200 (96%) vs 197/200 (98%) | 189/200 (94%) vs 196/200 (98%) |
| aniso3_w8 | 192 | 32 | 11.0 | 0/200 (0%) vs 0/200 (0%) | 192/200 (96%) vs 189/200 (94%) | 191/200 (96%) vs 186/200 (93%) |
| aniso3_w8 | 192 | 128 | 11.0 | 0/200 (0%) vs 0/200 (0%) | 182/200 (91%) vs 182/200 (91%) | 179/200 (90%) vs 178/200 (89%) |
| aniso3_w8 | 256 | 8 | 13.0 | 12/200 (6%) vs 9/200 (4%) | 194/199 (97%) vs 196/200 (98%) | 190/200 (95%) vs 193/200 (96%) |
| aniso3_w8 | 256 | 32 | 13.0 | 0/200 (0%) vs 0/200 (0%) | 182/200 (91%) vs 187/200 (94%) | 175/200 (88%) vs 176/200 (88%) |
| aniso3_w8 | 256 | 128 | 13.0 | 0/200 (0%) vs 0/200 (0%) | 170/200 (85%) vs 175/200 (88%) | 157/200 (78%) vs 155/200 (78%) |

E4 bound diagnosis: fitted M against the exact required supremum of |g'''| (g = log q vs log t) on [t_0, t_J] (`p22honest.exact_g3_cummax`). A 'violation' is M < required. Coverage is simultaneous. The runner's own `m_viol` flag was re-derived as M < required and agrees on every answered E4 row (see `verdict_report.txt`).

| family / scan cell | est | window | wpn | answered | violation rate (CP) | coverage | coverage given violation | coverage given no violation | median M / required sup g''' (median of cell medians) |
|---|---|---|---|---|---|---|---|---|---|
| aniso2_w16 | E4c | fixed | 8 | 1799 | 396/1799 = 22.0% [20.1-24.0] | 97.3% | 358/396 (90%) | 1393/1403 (99%) | 1.55 / 0.68 |
| aniso2_w16 | E4c | fixed | 32 | 1798 | 579/1798 = 32.2% [30.0-34.4] | 96.9% | 529/579 (91%) | 1213/1219 (100%) | 1.18 / 0.68 |
| aniso2_w16 | E4c | fixed | 128 | 1791 | 782/1791 = 43.7% [41.3-46.0] | 97.0% | 736/782 (94%) | 1001/1009 (99%) | 0.78 / 0.68 |
| aniso2_w16 | E4c | orig | 8 | 1787 | 323/1787 = 18.1% [16.3-19.9] | 98.0% | 296/323 (92%) | 1455/1464 (99%) | 2.24 / 0.68 |
| aniso2_w16 | E4c | orig | 32 | 1799 | 623/1799 = 34.6% [32.4-36.9] | 96.5% | 568/623 (91%) | 1168/1176 (99%) | 1.01 / 0.68 |
| aniso2_w16 | E4c | orig | 128 | 1798 | 803/1798 = 44.7% [42.3-47.0] | 97.2% | 758/803 (94%) | 989/995 (99%) | 0.75 / 0.68 |
| aniso2_w16 | E4r | fixed | 8 | 1800 | 618/1800 = 34.3% [32.1-36.6] | 95.6% | 548/618 (89%) | 1172/1182 (99%) | 0.90 / 0.68 |
| aniso2_w16 | E4r | fixed | 32 | 1800 | 805/1800 = 44.7% [42.4-47.1] | 95.2% | 724/805 (90%) | 989/995 (99%) | 0.82 / 0.68 |
| aniso2_w16 | E4r | fixed | 128 | 1796 | 951/1796 = 53.0% [50.6-55.3] | 96.5% | 895/951 (94%) | 838/845 (99%) | 0.69 / 0.68 |
| aniso2_w16 | E4r | orig | 8 | 1788 | 513/1788 = 28.7% [26.6-30.8] | 96.3% | 456/513 (89%) | 1265/1275 (99%) | 1.50 / 0.68 |
| aniso2_w16 | E4r | orig | 32 | 1800 | 886/1800 = 49.2% [46.9-51.6] | 94.4% | 793/886 (90%) | 907/914 (99%) | 0.71 / 0.68 |
| aniso2_w16 | E4r | orig | 128 | 1797 | 973/1797 = 54.1% [51.8-56.5] | 96.3% | 913/973 (94%) | 818/824 (99%) | 0.67 / 0.68 |
| aniso3_w8 | E4c | fixed | 8 | 1399 | 242/1399 = 17.3% [15.4-19.4] | 98.1% | 225/242 (93%) | 1147/1157 (99%) | 4.00 / 1.32 |
| aniso3_w8 | E4c | fixed | 32 | 1400 | 421/1400 = 30.1% [27.7-32.5] | 97.1% | 385/421 (91%) | 974/979 (99%) | 2.58 / 1.32 |
| aniso3_w8 | E4c | fixed | 128 | 1400 | 538/1400 = 38.4% [35.9-41.0] | 95.2% | 475/538 (88%) | 858/862 (100%) | 1.87 / 1.32 |
| aniso3_w8 | E4c | orig | 8 | 1127 | 175/1127 = 15.5% [13.5-17.8] | 98.7% | 168/175 (96%) | 944/952 (99%) | 6.99 / 1.32 |
| aniso3_w8 | E4c | orig | 32 | 1397 | 446/1397 = 31.9% [29.5-34.4] | 97.0% | 408/446 (91%) | 947/951 (100%) | 2.99 / 1.32 |
| aniso3_w8 | E4c | orig | 128 | 1405 | 560/1405 = 39.9% [37.3-42.5] | 95.2% | 498/560 (89%) | 840/845 (99%) | 2.18 / 1.32 |
| aniso3_w8 | E4r | fixed | 8 | 1400 | 354/1400 = 25.3% [23.0-27.6] | 97.3% | 326/354 (92%) | 1036/1046 (99%) | 3.17 / 1.32 |
| aniso3_w8 | E4r | fixed | 32 | 1400 | 563/1400 = 40.2% [37.6-42.8] | 96.3% | 515/563 (91%) | 833/837 (100%) | 1.94 / 1.32 |
| aniso3_w8 | E4r | fixed | 128 | 1400 | 654/1400 = 46.7% [44.1-49.4] | 93.6% | 566/654 (87%) | 744/746 (100%) | 1.80 / 1.32 |
| aniso3_w8 | E4r | orig | 8 | 1127 | 247/1127 = 21.9% [19.5-24.4] | 98.0% | 232/247 (94%) | 873/880 (99%) | 5.92 / 1.32 |
| aniso3_w8 | E4r | orig | 32 | 1397 | 609/1397 = 43.6% [41.0-46.2] | 95.6% | 549/609 (90%) | 786/788 (100%) | 2.22 / 1.32 |
| aniso3_w8 | E4r | orig | 128 | 1405 | 689/1405 = 49.0% [46.4-51.7] | 92.8% | 590/689 (86%) | 714/716 (100%) | 2.06 / 1.32 |
| sep_ratio2.5_b20 | E4c | orig | 32 | 200 | 141/200 = 70.5% [63.7-76.7] | 85.5% | 112/141 (79%) | 59/59 (100%) | 0.45 / 0.68 |
| sep_ratio2.5_b20 | E4r | orig | 32 | 200 | 179/200 = 89.5% [84.4-93.4] | 77.0% | 133/179 (74%) | 21/21 (100%) | 0.26 / 0.68 |
| sep_ratio3.0_b24 | E4c | orig | 32 | 200 | 105/200 = 52.5% [45.3-59.6] | 91.5% | 89/105 (85%) | 94/95 (99%) | 0.65 / 0.68 |
| sep_ratio3.0_b24 | E4r | orig | 32 | 200 | 166/200 = 83.0% [77.1-87.9] | 81.0% | 129/166 (78%) | 33/34 (97%) | 0.40 / 0.68 |

Separation-scan cells (ratio 2.5 = twostep_8_20 / b=20, id 302; ratio 3.0 = twostep_8_24 / b=24, id 303; wpn 32, original windows only; E1q for reference):

| scan | est | coverage (answered) | violations | coverage given violation | coverage given no violation | median M / required |
|---|---|---|---|---|---|---|
| ratio 2.5 (b=20) | E1q | 10/200 (5%) | - | - | - | - |
| ratio 2.5 (b=20) | E4c | 171/200 (86%) | 141/200 | 112/141 (79%) | 59/59 (100%) | 0.45 / 0.68 |
| ratio 2.5 (b=20) | E4r | 154/200 (77%) | 179/200 | 133/179 (74%) | 21/21 (100%) | 0.26 / 0.68 |
| ratio 3.0 (b=24) | E1q | 4/200 (2%) | - | - | - | - |
| ratio 3.0 (b=24) | E4c | 183/200 (92%) | 105/200 | 89/105 (85%) | 94/95 (99%) | 0.65 / 0.68 |
| ratio 3.0 (b=24) | E4r | 162/200 (81%) | 166/200 | 129/166 (78%) | 33/34 (97%) | 0.40 / 0.68 |

### Amendment 11 (frozen holdout; run once): new graph families
Families never used in Amendments 1-10: 2D->1D crossover C_L x C_12 (L = 96..384), 3D->1D crossover C_L x C_6 x C_6 (L = 48..192), controls ring 384/768, torus 40/80, 3D torus 18/26. Exact truth table `amend11_cells_exact.json` (sha256 de21f1a0..., written before simulation). All 10 crossover cells are realised (exact max drop >= 0.6; min 0.83). 200 datasets per cell, wpn 8/32/128, E1q/E4c/E4r/E7/E8 (E6 ungated, omitted).

Crossover cells: simultaneous-band coverage, k/answered (G1-G3 use E1q, E4c, E4r; E7/E8 shown for reference). exact_J = points in the exact window; max drop = exact.

| cell | exact_J | exact max drop | wpn | E1q | E4c | E4r | E7 | E8 |
|---|---|---|---|---|---|---|---|---|
| aniso2_w12 L=96 | 6 | 0.83 | 8 | 37/38 (97%) | 38/38 (100%) | 38/38 (100%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=96 | 6 | 0.83 | 32 | 148/175 (85%) | 175/175 (100%) | 175/175 (100%) | 200/200 (100%) | 198/200 (99%) |
| aniso2_w12 L=96 | 6 | 0.83 | 128 | 124/198 (63%) | 198/198 (100%) | 198/198 (100%) | 200/200 (100%) | 200/200 (100%) |
| aniso2_w12 L=128 | 8 | 1.00 | 8 | 154/169 (91%) | 169/169 (100%) | 168/169 (99%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=128 | 8 | 1.00 | 32 | 141/200 (70%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) |
| aniso2_w12 L=128 | 8 | 1.00 | 128 | 33/200 (16%) | 199/200 (100%) | 197/200 (98%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=192 | 10 | 1.00 | 8 | 158/200 (79%) | 198/200 (99%) | 197/200 (98%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=192 | 10 | 1.00 | 32 | 32/200 (16%) | 189/200 (94%) | 188/200 (94%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=192 | 10 | 1.00 | 128 | 1/200 (0%) | 193/200 (96%) | 188/200 (94%) | 200/200 (100%) | 200/200 (100%) |
| aniso2_w12 L=256 | 12 | 1.00 | 8 | 113/200 (56%) | 193/200 (96%) | 192/200 (96%) | 200/200 (100%) | 200/200 (100%) |
| aniso2_w12 L=256 | 12 | 1.00 | 32 | 8/200 (4%) | 187/200 (94%) | 181/200 (90%) | 200/200 (100%) | 200/200 (100%) |
| aniso2_w12 L=256 | 12 | 1.00 | 128 | 0/200 (0%) | 194/200 (97%) | 192/200 (96%) | 200/200 (100%) | 199/200 (100%) |
| aniso2_w12 L=384 | 14 | 1.00 | 8 | 31/200 (16%) | 198/200 (99%) | 191/200 (96%) | 195/200 (98%) | 199/200 (100%) |
| aniso2_w12 L=384 | 14 | 1.00 | 32 | 8/200 (4%) | 188/199 (94%) | 177/200 (88%) | 197/200 (98%) | 200/200 (100%) |
| aniso2_w12 L=384 | 14 | 1.00 | 128 | 0/200 (0%) | 187/200 (94%) | 162/200 (81%) | 196/200 (98%) | 200/200 (100%) |
| aniso3_w6 L=48 | 3 | 1.15 | 8 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=48 | 3 | 1.15 | 32 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=48 | 3 | 1.15 | 128 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=64 | 5 | 1.62 | 8 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=64 | 5 | 1.62 | 32 | 15/22 (68%) | 22/22 (100%) | 22/22 (100%) | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=64 | 5 | 1.62 | 128 | 1/7 (14%) | 7/7 (100%) | 7/7 (100%) | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=96 | 7 | 1.64 | 8 | 98/124 (79%) | 123/124 (99%) | 124/124 (100%) | 200/200 (100%) | 198/200 (99%) |
| aniso3_w6 L=96 | 7 | 1.64 | 32 | 56/200 (28%) | 198/200 (99%) | 196/200 (98%) | 200/200 (100%) | 199/200 (100%) |
| aniso3_w6 L=96 | 7 | 1.64 | 128 | 2/200 (1%) | 198/200 (99%) | 196/200 (98%) | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=128 | 9 | 1.64 | 8 | 123/196 (63%) | 195/196 (99%) | 195/196 (99%) | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=128 | 9 | 1.64 | 32 | 9/200 (4%) | 195/200 (98%) | 192/200 (96%) | 200/200 (100%) | 198/200 (99%) |
| aniso3_w6 L=128 | 9 | 1.64 | 128 | 0/200 (0%) | 197/200 (98%) | 197/200 (98%) | 200/200 (100%) | 200/200 (100%) |
| aniso3_w6 L=192 | 11 | 1.64 | 8 | 47/200 (24%) | 196/200 (98%) | 194/200 (97%) | 200/200 (100%) | 199/200 (100%) |
| aniso3_w6 L=192 | 11 | 1.64 | 32 | 0/200 (0%) | 194/200 (97%) | 193/200 (96%) | 200/200 (100%) | 197/200 (98%) |
| aniso3_w6 L=192 | 11 | 1.64 | 128 | 0/200 (0%) | 197/200 (98%) | 195/200 (98%) | 200/200 (100%) | 200/200 (100%) |

Control cells: simultaneous-band coverage k/answered (G4: E4c, E4r, E7, E8 at every wpn; E1q for reference):

| cell | exact_J | wpn | E1q | E4c | E4r | E7 | E8 |
|---|---|---|---|---|---|---|---|
| ring L=384 | 12 | 8 | 185/200 (92%) | 199/200 (100%) | 199/200 (100%) | 200/200 (100%) | 199/200 (100%) |
| ring L=384 | 12 | 32 | 189/200 (94%) | 200/200 (100%) | 199/200 (100%) | 200/200 (100%) | 199/200 (100%) |
| ring L=384 | 12 | 128 | 192/200 (96%) | 195/200 (98%) | 196/200 (98%) | 200/200 (100%) | 199/200 (100%) |
| ring L=768 | 16 | 8 | 187/200 (94%) | 198/200 (99%) | 199/200 (100%) | 199/200 (100%) | 198/200 (99%) |
| ring L=768 | 16 | 32 | 183/200 (92%) | 197/200 (98%) | 198/200 (99%) | 196/200 (98%) | 199/200 (100%) |
| ring L=768 | 16 | 128 | 191/200 (96%) | 197/200 (98%) | 197/200 (98%) | 188/200 (94%) | 200/200 (100%) |
| torus2 L=40 | 6 | 8 | 90/95 (95%) | 95/95 (100%) | 95/95 (100%) | 200/200 (100%) | 198/200 (99%) |
| torus2 L=40 | 6 | 32 | 187/198 (94%) | 198/198 (100%) | 198/198 (100%) | 200/200 (100%) | 199/200 (100%) |
| torus2 L=40 | 6 | 128 | 184/200 (92%) | 199/200 (100%) | 199/200 (100%) | 200/200 (100%) | 200/200 (100%) |
| torus2 L=80 | 10 | 8 | 188/200 (94%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) |
| torus2 L=80 | 10 | 32 | 185/200 (92%) | 199/200 (100%) | 197/200 (98%) | 200/200 (100%) | 199/200 (100%) |
| torus2 L=80 | 10 | 128 | 188/200 (94%) | 198/200 (99%) | 198/200 (99%) | 199/200 (100%) | 200/200 (100%) |
| torus3 L=18 | 4 | 8 | n/a | n/a | n/a | 200/200 (100%) | 198/200 (99%) |
| torus3 L=18 | 4 | 32 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| torus3 L=18 | 4 | 128 | n/a | n/a | n/a | 200/200 (100%) | 200/200 (100%) |
| torus3 L=26 | 7 | 8 | 184/195 (94%) | 194/195 (99%) | 194/195 (99%) | 200/200 (100%) | 200/200 (100%) |
| torus3 L=26 | 7 | 32 | 182/200 (91%) | 199/200 (100%) | 198/200 (99%) | 200/200 (100%) | 200/200 (100%) |
| torus3 L=26 | 7 | 128 | 187/200 (94%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) | 200/200 (100%) |

Control false-change declarations (G5 at wpn 32; margin 0.25, any sign), k/answered:

| cell | wpn | E1q | E4c | E4r | E7 | E8 |
|---|---|---|---|---|---|---|
| ring L=384 | 8 | 2/200 (1%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| ring L=384 | 32 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| ring L=384 | 128 | 1/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| ring L=768 | 8 | 2/200 (1%) | 1/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| ring L=768 | 32 | 1/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| ring L=768 | 128 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=40 | 8 | 2/95 (2%) | 0/95 (0%) | 0/95 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=40 | 32 | 1/198 (1%) | 0/198 (0%) | 0/198 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=40 | 128 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=80 | 8 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=80 | 32 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus2 L=80 | 128 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus3 L=18 | 8 | n/a | n/a | n/a | 0/200 (0%) | 0/200 (0%) |
| torus3 L=18 | 32 | n/a | n/a | n/a | 0/200 (0%) | 0/200 (0%) |
| torus3 L=18 | 128 | n/a | n/a | n/a | 0/200 (0%) | 0/200 (0%) |
| torus3 L=26 | 8 | 5/195 (3%) | 0/195 (0%) | 0/195 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus3 L=26 | 32 | 1/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |
| torus3 L=26 | 128 | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) | 0/200 (0%) |


## Amendment 12 results (PRE-REGISTERED: Gašpar–Kukal estimator benchmark)
Bars, arms, cells, seeds and the reading rule are those of PLAN.md Amendment 12, written before any run; analysis is `analyze_amend12.py` (sha256 1439d8ca...), unmodified, run once. Runs on host-D (WSL); manifests `results/manifest_amend12_b{8,32,128}.json` are `complete: true`, 0 errors, 640/640 items, 57,600 rows each, and the raw-CSV sha256 values in the manifests were re-checked against the copies in `results/`; p22lib.py, run_p22.py and p22amend12.py hashes match the manifests. Outputs: `results/amend12_report.txt`, `amend12_table.csv`, `amend12_minL.csv`, `amend12_control_cells.csv` (per-control-cell k/n for every arm and budget). "Coverage" is simultaneous coverage of the exact d(t) curve over non-abstained control datasets; k/n pooled, Clopper–Pearson (CP) 95%. Abstention is about 31% at every budget (J < 6), and identical across arms.

**Verdicts at wpn 32 (gated).** (C) pooled >= 93% and no control cell with CP upper < 93%; (F) FPR <= 5%, no cell with CP lower > 5% (two-segment arms).

| arm | what it is | C | F | pooled coverage (wpn 32) | failing control cells |
|---|---|---|---|---|---|
| GKW | their WSSQ, regression se | **FAIL** | n/a | 1716/1930 = 88.9% (87.4-90.3) | 5 of 10 (9001, 9002, 9006, 9007, 9011) |
| GKL | their LSSQ, regression se | **FAIL** | n/a | 1659/1930 = 86.0% (84.3-87.5) | 9 (9000, 9001, 9002, 9006, 9007, 9008, 9011, 9012, 9013) |
| GKWm | WSSQ, model-based se | PASS | n/a | 1813/1930 = 93.9% (92.8-95.0) | 0 |
| GKWb | WSSQ + our bootstrap | PASS | n/a | 1807/1930 = 93.6% (92.4-94.7) | 0 |
| E0g | reference | PASS | n/a | 1817/1930 = 94.1% (93.0-95.2) | 0 |
| GK2W | two-segment, analytic se | **FAIL** | PASS (32/1930 = 1.7%) | 1255/1930 = 65.0% (62.9-67.2) | 10 |
| GK2L | two-segment, LSSQ | **FAIL** | PASS (36/1930 = 1.9%) | 1221/1930 = 63.3% (61.1-65.4) | 10 |
| GK2Wb | two-segment + bootstrap | **FAIL** | PASS (4/1930 = 0.2%) | 1811/1930 = 93.8% (92.7-94.9) | 1 (torus2 L=96: 176/200 = 88.0%) |
| E0 | reference | PASS | PASS (6/1930 = 0.3%) | 1821/1930 = 94.4% (93.2-95.3) | 0 |

GK2Wb fails C only through the single-cell rule: its pooled coverage (93.8%) is above 93%; the same cell is 180/200 = 90.0% for E0, which passes because its CP upper limit stays above 93%. Reported as found.

**Other budgets (reported, not gated).** GKW fails C at all three budgets: wpn 8 1730/1935 = 89.4%, wpn 32 88.9%, wpn 128 1739/1968 = 88.4%. GKL fails at all three (86.2%, 86.0%, 86.7%). GKWm passes at 8 (93.5%), 32 (93.9%) and fails at 128 on one cell (92.6%, cell 9213). GKWb and E0g pass at 8 and 32 and fail at 128 on cells 9212/9213 (92.5% and 91.9%; E0g does so as well, so this is a property of the budget and not of the weighting). GK2W/GK2L fail C at 8 and 128 (62.6%/61.5% and 66.9%/66.0%) and fail F at wpn 8 (FPR 5.3% and 5.7%, cells 9100, 9106, 9110-9112). GK2Wb: wpn 8 PASS C and F (93.6%, FPR 0.7%), wpn 128 PASS C and F (94.5%, 0.0%). E0 passes C and F at 8 and 128 (93.8%/0.9% and 94.1%/0.0%). Per-cell counts for all arms and budgets are in `amend12_control_cells.csv`.

**Nominal constant exponent on controls (ungated; the target of their Z test).** GKW covers the nominal d (1, 2, 3) in 86.9% (wpn 32), 88.8% (wpn 8), 82.4% (wpn 128) of datasets, i.e. below 95% and worse at the largest budget; GKWb 92.5%, 93.9%, 86.1%.

**M (minimum detectable L, power >= 0.8 at that and all larger L; abstention = miss; no pass threshold set in the amendment).** wpn 32: GK2W, GK2L, GK2Wb and E0 all give aniso2_w16 = 192 and aniso3_w8 = 80. wpn 8: GK2W 256 / 96, GK2L 224 / 96, GK2Wb 224 / 96, E0 224 / 96. wpn 128: all four 160 / 80. M is only informative for arms that pass C; GK2W and GK2L, which do not, match E0 at 32 and 128 while sitting at about 65% coverage.

**Coverage of realised crossovers (X, descriptive, one-scale arms).** GKW 37.0%, GKWb 16.9%, E0g 15.0% at wpn 32: a single global exponent cannot cover a changing d(t), as expected; GKW's higher figure reflects its narrower intervals, not a better target.

**Pre-registered predictions.** (1) GKW points nearly identical to E0g: held (mean bias -0.0018 vs -0.0019; wpn 8 -0.0137 vs -0.0138; wpn 128 0.0013 both). (2) analytic se narrower than the bootstrap, GKW coverage <= GKWb: held at all budgets (88.9% vs 93.6% at wpn 32). (3) GKW/GKL fail C at one or more budgets: held (all three). (4) weighting changes coverage and minimum L little relative to E0g/E0: held for point estimates and minimum L; coverage differences come from the se, not the weights (GKW vs GKWm differ only in the se: 88.9% vs 93.9%).

**Reading (rule applied as written).** The published estimator with its regression standard error is not calibrated on these graphs at any budget tested: 88-89% coverage of the exact curve and 82-89% of the nominal constant, consistent with (not separately tested) the regression se treating time points that share walkers as independent. Replacing only the se (model-based, or our cluster bootstrap) restores calibration at wpn 8 and 32 with the same point estimates. Their estimator was designed for, and evaluated on, one very host-D substrate per case with an oracle window; none of this speaks to that setting. What this benchmark does not measure is their estimator in their setting. Amendment 12 changes no earlier verdict. Stale text in this file ("Gaspar-Kukal not implemented", in the notes above) is superseded by this section.
