# EXPLORATORY: exploratory-log row 29 (E0g misses at wpn 128, Amendment 12)

**Status: EXPLORATORY. No confirmatory claim.** Written after the Amendment 12 results were seen; not pre-registered; no bar, verdict or earlier result is changed by it. Scripts: `row29_explore.py` (run), `row29_analyze.py` (analysis); full numbers in `row29_analysis.txt`.

## Question
At wpn 128, E0g (and GKWb, GKWm) fell below the 93% simultaneous-coverage bar on controls 9212 and 9213 (torus3 L=28 and L=32; the estimator is unexpectedly worse than at wpn 32). Is that bias, an SE problem, or window selection?

## Method and reproduction (on host-T)
- Pinned venv on host-T (`row29_requirements.lock`: numpy 2.5.5, scipy 1.18.1, pandas 3.1.1, python-flint 0.9.0, python-dateutil 4.0.0, tzdata 2026.6), `python -I`, 1 thread per worker, 8 workers. Code: repo HEAD `e255d37` p22 files (frozen `p22lib.py`/`run_p22.py` sha256 asserted equal to Amendment 11 by `run_amend12.verify()`; `p22amend12.py` sha equals the manifest's).
- Reran cells 9212 and 9213 through the unchanged `run_amend12.run_dataset`, same seeds `SeedSequence([22, 12000+id, ds])`. **Datasets 0-199 reproduce the stored `raw_amend12_b128.csv` exactly: 3600 rows, 75,600 numeric values, all bit-equal, no NaN-pattern difference.** Datasets 200-799 are new (600 new = 3x the original 200; 800 total per cell). Platform of the stored run was not recorded; this check says the result is not platform-fragile for these cells.
- Extra diagnostic on the same walker data (own bootstrap stream, nothing stored is disturbed): an E0g-style global weighted fit with bootstrap SE for the data-driven window and for fixed window lengths J=5..10, recording d_hat, boot SE, exact d at the fit's log-mean time (`d_true`, the scoring truth) and the same fit applied to the exact return curve (`d_pop`, what the estimator targets in expectation).

## Results (pooled 800 datasets per cell; original-200 and new-600 agree, see `row29_analysis.txt`)
| | cell 9212 (torus3 L=28) | cell 9213 (torus3 L=32) |
|---|---|---|
| E0g coverage, original 200 / new 600 / pooled | 88.0 / 90.5 / **89.9%** (CP 87.6-91.9) | 84.5 / 81.8 / **82.5%** (CP 79.7-85.1) |
| E0 (two-segment) pooled | 94.0% | 93.0% |
| mean(d_hat - d_true) | +0.0095 (+16.5 se of the mean) | +0.0118 (+26.5 se) |
| mean(d_hat - d_pop) | -0.0007 (about -1.2 se) | +0.0005 (about +1.1 se) |
| empirical SD(d_hat) vs mean bootstrap SE | 0.01637 vs 0.01646 (ratio 1.005) | 0.01253 vs 0.01256 (ratio 1.002) |
| z = err/SE: mean, SD | +0.58, 1.005 | +0.94, 1.003 |
| coverage of d_true by +-1.96 SE | 90.8% | 82.9% |
| coverage of d_pop by +-1.96 SE | 96.2% | 94.8% |
| window J chosen by the rule | 7 on 800/800 datasets | 8 on 800/800 datasets |

## Diagnosis
1. **The misses are real and not noise**: the 3x larger rerun gives the same or slightly worse coverage (9213: 81.8% on the new 600). The original 200 were not a fluke.
2. **It is not SE miscalibration.** The bootstrap SE matches the empirical SD of d_hat to within 0.5% in both cells, and d_hat covers its own population value d_pop at 96.2% and 94.8% (nominal 95%).
3. **It is a fixed offset between the estimand and the scoring truth, not estimator bias.** d_hat is unbiased for d_pop (differences of order 1 se of the mean), but d_pop differs from the scoring truth d_true by about +0.0102 (9212) and +0.0113 (9213) (the `d_true - d_pop` column is -0.0102 / -0.0113). The exact d(t) on these controls is still falling across the window (d about 3.058 down to 3.003), so the log-log slope over the window is not the local exponent at the log-mean time. That offset is 0.6 SE (9212) and 0.9 SE (9213), which with perfectly calibrated SE gives the observed coverage: normal-theory coverage with the observed bias/SE and SD/SE is 91.2% and 84.6%, against measured 90.8% and 82.9%.
4. **Window selection is not the cause here.** The rule picked the same J on every dataset (7 and 8), so there is no data-dependent selection to bias anything. Window length does matter through the offset: fixed shorter windows shrink it and coverage rises (9213: J=5 gives 92.9%, J=6 90.6%, J=8 82.9%, J=9 81.5%), with SE growing; 9212 is non-monotone (J=5/6/10 about 93%, J=7 90.8%, J=8/9 about 89-90%).
5. **Why wpn 128 and not wpn 32 (inference, not measured here):** the SE scales roughly as 1/sqrt(walkers) while the offset does not, so the same offset is a larger fraction of the SE at wpn 128. I did not run wpn 32 in this follow-up, so this is consistent with the numbers but not tested.
6. Why E0 (two segments) passes: its bands are about 4x wider (hw 0.11 vs 0.03), a different, much less sharp statement; not evidence the offset is absent there.

## What this does not say
- It does not show the E0g/GK arms are "wrong": by this analysis they estimate the window-averaged slope correctly. The failure is against the Amendment's chosen truth (exact d at the log-mean time of the window) on controls whose d(t) is still curving. Whether that truth definition is the right one for a bar is a design choice for the paper; nothing was changed here.
- Two controls only; one estimator family; windows are the pipeline's, not tuned. Not run: other torus3 sizes, other wpn, an alternative truth definition.

## Reproduce
`<remote command omitted>` (resumable, about 2.4 h at 8 workers under load) then `../venv/bin/python -I row29_analyze.py ../out ../stored/raw_amend12_b128.csv`. Raw units: host-T `<local path omitted> (320 pickles).
