# ROUND5_CHECKS: round-5 referee items (paper A v5b) -- 2026-11-24

All checks below are EXPLORATORY or corrective reporting from stored results; none feeds a gate, bar or verdict and none changes an estimator, seed or cell. Outputs are in `round5_checks/`; every number in the paper traces through `ideas/A_spectral_dimension_coverage/CLAIMS.md` rows 248-307 and `verify.py` (302/302 verified).

## 0. Priority (ChatGPT s8): cross-time dependence of the simulated counts -- `round5_dep_check.py`, `round5_checks/round5_dep_check.txt`
**Shared walkers CONFIRMED.** `p22lib.Graph._direct` / `_steps` advance one position array through all grid times, so a block's counts at different times are sums over the same walkers. Ring C128, lazy, times 10/24/34 (run on host-D, numpy there): exact per-walker cov(ret24, ret34) = q(24)[q(10)-q(34)] = 0.009142; empirical 0.009132 (se 0.000013, 40 x 2e6 walkers, z = -0.8); marginals match (q24 0.114567 vs 0.114584, q34 0.096403 vs 0.096346). Blocks of m = 20 (800,000 blocks): 0.1843 (se ~0.002) vs exact 0.1828. Independent draws would give 0. (Two earlier smaller runs sat about 2 se low; the 40-run pool removes that.)

## 1. Gemini typo list
Checked in source and PDF text: "E0 and E2" is correct in the v5 source; no "IIII" (cross-references render "Sec. III I"); `\ln` renders; no doubled "( (" in text or Table I; SM Table S3 parentheses and Tables S5-S7 columns render cleanly at 80 dpi. Nothing to fix; the referee's extraction artefacts are not in the PDFs. (Only change: none.)

## 2. E2 abstention: 57% vs (2800-1108)/2800 -- `round5_misc.py`, `round5_checks/round5_misc.txt`
The referee read 1108 as the answered count; it is the COVERED count. E2 on the 14 control cells: n = 2800, answered 1208, covered 1108, abstained 1592 = 56.9%. 57% was right; the text now says "(1592 of 2800)".

## 3. SD of studentised stochastic errors -- `round5_sdz.py`, `round5_checks/round5_sdz.txt`
Stored fresh-seed run (ring L=512, 1000 datasets, 13,653 (n_j, sd_j) pairs, `markowitz_fresh_perj.csv` on host-D `<local path omitted> no rerun needed). SD(z_j) = 1.167 (V0), 1.118 (V1, B=1000), **1.178 (V2, independent split)**; median over 27 strata (>=100 datasets) 1.165 / 1.124 / 1.166; |z|>1.96 in 9.2 / 8.1 / 9.7% of pairs. **About 1.17, not 1.00:** the 0.86 ratio is not a marginal-vs-conditional (scale-mixture) artefact; the bootstrap sd is genuinely too small by about that factor, also with the independent covariance. Stated in main Sec. III G and SM S7B.

## 4. L_min-bracketing counts and holdout gate inputs -- `round5_tables.py`
SM Table S11 (E1q, E4c, E4r, E6, Mb, CL: counts and CP at L_min and adjacent tested sizes; the script asserts its recomputed L_min equals Table II for all 12 arm-family entries and exits non-zero otherwise) and SM Table S12 (G3 covered/answered; G4 pooled counts and the control cell with the lowest CP upper limit). E4c 2D reaches no L_min within the sweep (79.0% at 768).

## 5. E7 relaxations A/B (ChatGPT s3)
E7j exists only for relaxation B. Endpoint test at 3D L = 160/192/256: A 0/0/0 of 200, B 35/46/147. Joint (B) 170/186/198. All three joint point estimates meet 0.8; only the L=160 CP interval (79.3-89.6) crosses below.
