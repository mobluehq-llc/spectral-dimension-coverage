# RESULTS_BIAS_NOISE (Amendment 14, 2026-11-24; EXPLORATORY, descriptive, stored data only, no walks simulated)
Code `bias_noise.py` (sha256 fe281abb...), data `archive/data/raw_amend10c.csv` (sha256 fe21111c...), E1q fixed-window rows, 16 cells x 3 budgets (aniso2_w16 L 128-768, aniso3_w8 L 80-256), 200 datasets each, all answered. aniso3_w8 L = 48 and 64 dropped by the registered rule (all datasets abstained). Outputs: `results/bias_noise_{cells,by_h,rising,frozen_partial}.csv`, `bias_noise_rule.json`.

## Registered rule (literal): NOT EXPLAINED
- Clause (i): R1 (|window-mean deterministic bias| / SD of window-mean noise) rises strictly 8 < 32 < 128 in 7/16 cells = 43.8% (bar >= 75%): FAIL.
- Clause (ii): cell-budget pairs with R1 > 2 hold 0 of 7,567 misses (bar >= 80%): FAIL. R1 > 1 holds 7.9%. Max R1 is 1.96.
- Reason: the window-mean bias cancels in sign along the window (R1 numerator median about 0.02, against a per-time bias of up to 0.4-0.6), so the registered functional hides the bias that breaks a simultaneous band. Under the rule, the paper may make no mechanism claim from this amendment.

## Descriptive, not registered (labelled post hoc)
- Per-time bias is large: mean over datasets of max_j |b_j(h_ds)| is 0.15-0.67 (median 0.25 at wpn 32), budget-invariant for fixed h by construction.
- Sampling spread shrinks about 3.3x from wpn 8 to 128 (SD of window-mean noise 0.114 -> 0.065 -> 0.034, medians); mean bootstrap se 0.198 -> 0.106 -> 0.058.
- R2 = max-bias / mean se: median 2.12, 2.90, 3.61 at 8, 32, 128; rises monotonically in 13/16 cells; R2 > 2 in 37/48 cell-budgets holding 88.4% of misses; R2 > 1 in 45/48 holding 98.2%.
- R3 = max-bias / median simultaneous half-width: median 0.88, 1.09, 1.29; R3 > 0.5 holds 97.0% of misses, R3 > 1 holds 75.8%.
- Misses rise 1,908 (wpn 8) -> 2,573 (32) -> 3,086 (128) of 3,200 answered per budget.
- The CV bandwidth shrinks with budget (mean h 2.41 -> 1.90 -> 0.97; share h <= 1: 40% -> 54% -> 90%), which lowers bias and raises variance, yet misses still rise: consistent with bias staying above a shrinking band. Not a test of cause.
- Identity check: stored (bias - bias_pop) versus population-weight bias at the stored h: cell median absolute difference <= 0.0045 (prediction b: < 0.02, held).

## Frozen-bandwidth arm: PARTIAL, as registered
Stored data hold no per-time estimates and no frozen-h run, so no frozen-h spread exists. Only the deterministic bias at the wpn-32 modal h (budget-invariant) is available (`bias_noise_frozen_partial.csv`). A real frozen-h run needs new walks and was not done.

## Reading
Registered descriptive bar not met (R1 flaw: the window-mean functional). The per-time ratio, chosen after seeing that failure and unregistered, is consistent with the bias-versus-shrinking-band reading but would need a pre-registered rerun (and a frozen-h walk run) to carry any claim. Predictions: (a) partly held (CV h shrinks, R1 does not track sqrt(wpn)); (b) held; (c) failed (R1 > 2 in no pair).
