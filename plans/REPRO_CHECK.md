# REPRO_CHECK (paper A, round 3) -- 2026-11-23

**What this is.** One small cell of each of Amendments 9b, 10b, 10c and 11, re-run with the p22 code at the current commit in a fresh pinned environment on `host-D`, compared with the stored raw results. No estimator, bar, seed or cell definition was changed; `repro_check.py` only calls the runners' own per-dataset functions (`run_amend1011.py` helpers, `run_p22i/j/k.py`). It is a reproducibility spot check, not a re-run of the amendments.

## Environment (recorded, pinned)
- Host `host-D`: WSL2 Ubuntu (Linux 6.18.58.1-microsoft-standard-WSL2), AMD Ryzen 7 255, Python 3.12.3 (GCC 13.3.0). Venv created with `python3 -m venv --without-pip` + `pip --python`, `python -I` (no user site, no PYTHONPATH, so the older hand-built `<local path omitted>` is NOT used).
- Exact versions (`repro_out/requirements.lock`): numpy 2.5.5, scipy 1.18.1, pandas 3.1.1, python-flint 0.9.0, python-dateutil 4.0.0, tzdata 2026.6. These equal the versions in the old pylib on host-D, which produced the stored 10b/10c/11 results. The stored 9b results were produced on host-T; its versions were not recorded here.
- Code: repo HEAD `89a381b`; the p22 directory last changed in `f23021e` (the Amendment 10b/10c/11 commit). Tarball of the tracked p22 files plus `repro_check.py`; `run_amend1011.verify()` confirmed the six frozen estimator files byte-identical to 1ba9ca6 before every run.
- Stored comparands: `results/raw_amend9_b32.csv` (git, sha256 cf5bfd10...), `results_amend10_11/raw_amend10b.csv` (git, 2adfb2c6...), `raw_amend10c.csv` (fe21111c...) and `raw_amend11.csv` (533ef4ce...), the last two NOT in git, taken from host-D `<local path omitted>` with sha256 asserted equal to the committed manifests before use.
- Each run was a single process, 1 thread (OMP/OPENBLAS/MKL = 1), 4 cores total.

## Comparison rule (fixed before the numbers below were read, with one correction disclosed after)
Timing columns excluded. Strings and verdict-bearing discrete columns (abstain, J, coverage and declaration flags, statuses, counts of LP/gate outcomes) must be EXACT. Floats are counted as exact / within 1e-9 / within 1e-6 (relative+absolute) / larger. Verdicts: MATCH_EXACT, MATCH_WITHIN_TOL, DISCRETE_MATCH_FLOAT_DIFFERS (all verdict-bearing columns identical but solver-level floats differ), MISMATCH.

## Results
| Amendment | Cell re-run | Rows | Verdict |
|---|---|---|---|
| 9b | ring L=128, wpn 32 (control), datasets 0-9, all 7 estimators | 70/70 matched | **DISCRETE_MATCH_FLOAT_DIFFERS** (E7 only) |
| 10b | twostep_8_32 L=320, datasets 32 and 48 (group sens), certified bands | 2/2 | **MATCH_EXACT** (every non-timing value bit-identical) |
| 10c | aniso2_w16 L=128, wpn 8, datasets 0-5, original and fixed window (J=8), E1q/E4c/E4r | 36/36 | **MATCH_EXACT** |
| 11 | ring L=384, wpn 8 (smallest cell), datasets 0-5, E1q/E4c/E4r/E7/E8 | 36/36 | **MATCH_EXACT** |

**9b detail.** Every verdict-bearing discrete value (576 values: abstain, J, coverage, declarations) is identical, and every estimator except E7 is bit-identical or within 1e-9. For E7 on this cell, 51 float values differ by more than 1e-6 (columns c_box, gate1, hw, width_mean, a1p_maxgap, and the diagnostic counts n_boot, n_cp): e.g. dataset 0, c_box 3.5176 now vs 3.5541 stored, n_boot 15 vs 16; dataset 1 matches exactly. Largest absolute difference 1.0 (a diagnostic count). Re-running dataset 0 and 1 twice on host-D gives identical E7 output both times, so the code is deterministic here; the difference is between platforms/libraries (stored 9b came from host-T), most likely the E7 box construction flipping a near-tie threshold count. The earlier Amendment-10 `repro` step (3 cells x 10 datasets, `results/repro_amend9_b32.json`) reported the same pattern: discrete mismatches 0, continuous max diff 0.017. **Reading:** 9b reproduces at the level of every reported coverage/abstention/declaration value on this cell; E7's continuous band widths on 9b are not bit-reproducible across platforms. This was not tested further (single cell, 10 datasets) and no claim in the paper rests on E7 hw at the 1e-6 level.

## Disclosures
- Runs 1 and 2 of this script (first-run, before the comparison rule above) flagged 9b as MISMATCH under a tolerance that counted n_boot/n_cp as discrete and lp_seconds as a value, and 10c used a cell (aniso3_w8 L=48) on which every estimator abstained (a vacuous comparison). I corrected both (timing columns excluded; n_boot/n_cp treated as diagnostic counts; 10c cell changed to aniso2_w16 L=128) before accepting results; the discarded outputs are kept on host-D in `<local path omitted>`, `out_v2_discarded`. The 11 and 10b verdicts were identical under every version.
- One cell per amendment is a spot check: it does not establish that all cells reproduce.
- Not run: Amendment 9a/9c, 10a, the exact-q 10b arm, or the 10c windows/Amendment 11 exact-table steps.

## Reproduce
`ssh host-D 'wsl -e bash -c "<local path omitted> 9b 10b 10c 11"'` (skips finished); or in a checkout with the pinned venv: `python -I repro_check.py {9b|10b|10c|11} --stored <dir with raw_amend9_b32.csv, raw_amend10c.csv, raw_amend11.csv> --out repro_out`. Raw JSON per amendment and logs: `repro_out/`. Wall time: 9b 5 s, 10c 1 s, 11 5 s, 10b 920 s.
