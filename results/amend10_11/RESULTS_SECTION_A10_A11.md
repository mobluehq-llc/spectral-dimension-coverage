
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

