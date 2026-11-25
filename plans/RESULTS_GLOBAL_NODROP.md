# RESULTS: global no-drop LP (2026-11-24). Pre-registration PREREG_GLOBAL_NODROP.md (commit 30b2514); code rt_nodrop.py, rt_nodrop_analyze.py

Environment: archive/.venv (numpy 2.5.5, scipy 1.18.1), host-W CPU, 5 workers, nice. Total wall about 20 min (budget 3 h); no dataset count was cut (100 per cell as registered). Raw: results/nodrop_exact.csv, nodrop_data.csv, nodrop_controls.csv (+ logs). EXPLORATORY; the paper is not edited.

## 1. Exact q (sign test) PASSED
Box = exact q (1 +- eps) at every grid time, Jw = exact_J. eps = 1e-4, 1e-3, 1e-2 gave the same outcomes (12 LPs each, no gate disagreement).
| cell | exact max drop in window | outcome (all 3 eps) | joint | E7 endpoint |
|---|---|---|---|---|
| C768 x C16 (Jw 18) | 1.004 | **DROP** | declares | declares |
| C160 x C8^2 (Jw 10) | 2.020 | **DROP** | declares | declares |
| controls ring 1024, torus2 96, torus3 24 | 0, 0.004, 0.053 | WITNESS (all 9 rows) | no | no |
(The secant certificate of 0.492 is a lower bound; the exact in-window maximum drop is 1.004.) Narrower boxes (+-1e-6) are BOXINFEASIBLE on the 800-atom grid in the unit test, hence the ladder.

## 2. Stored-type datasets, 32 walkers per node, first 100 per cell, same box as the joint test
Data regenerated from the frozen seeds (stored 2.4.6 walks are not bit-reproducible under 2.5.5). Regenerated joint decisions agree with rt_monolp_p2.csv B0 on 97/100 (2D; the 3 differences are the 3 datasets where the frozen E7 LP returns numerical status 4, counted infeasible by the joint test) and 100/100 (3D); against the original stored raw_amend9_a.csv: 2D 100/100 (0 declarations both), 3D 96/100 decisions, 87 declarations both.
| cell | global DROP | WITNESS | INCONCLUSIVE | joint-test drops | E7 endpoint drops |
|---|---|---|---|---|---|
| C768 x C16, 32/node | **37** | 58 | 5 | **0** (3 infeasible, status 4) | 0 |
| C160 x C8^2, 32/node | **100** | 0 | 0 | 87 | 27 |
Cross-tab 3D: joint drop 87 -> global DROP 87; joint no-drop 13 -> global DROP 13. 2D: joint 0 -> global DROP 37 / WITNESS 58 / INCONCLUSIVE 5. Monotonicity holds (every joint declaration is also a global DROP; 0 exceptions). Every DROP passed the 1597-atom gate (gate disagreements 0). AMBIGUOUS (no-drop witness AND drop witness both verified): 58/100 in 2D, 0/100 in 3D. Witnesses sit on the boundary of the null (min Delta = -0.25 to 1e-9, median -0.250000), i.e. the feasible no-drop measures are barely no-drop. Cost: 2D 6.7 s and 31 LPs per dataset (max 401 = cap, the 5 INCONCLUSIVE), 3D 1.6 s and 13 LPs.
Caveat: of the 37 2D DROP datasets, one (ds 71) is among the 3 whose frozen E7 gate returns LP status 4; the global solver itself reported no numerical status there (36 excluding it).

## 3. Control run (post-hoc, labelled in the prereg Addendum 1, added after seeing 2.)
Ring L=1024, torus2 L=96, torus3 L=24, first 100 datasets each at 32/node: **WITNESS 300/300, DROP 0/300** (joint 0/300). So the global test made no false declaration on the three no-drop controls (95% upper bound on the false rate 1.0%, rule of three).

## 4. Against the pre-registered bars and predictions
- Exact-q sign test: passed.
- "Identifiability limit" bar (keep wording only if 0 of 100 DROP and WITNESS in >= 95/100): **FAILED**. 37/100 datasets are DROP-certified at 32 walkers per node on C768 x C16 (WITNESS 58, INCONCLUSIVE 5). The paper sentence "the implemented joint test declared no 2D drop through L=768" stays true of the joint test, but a statement that the 32-walker box cannot resolve the 2D drop is false for the global test: it certifies the drop in about 37% (up to 42% if the 5 INCONCLUSIVE were all drops).
- Prediction (0 DROP in 2D) was **wrong**; the 3D prediction (global >= joint) held, with 100 vs 87.
- Reading: the joint test asks that ONE time t_j drop for every measure; the global test lets the dropping time depend on the measure, which removes the slack that the single-j test leaves. For the 2D cell this converts 0 into 37 certified declarations without any false declaration on the controls.

## 5. Not done / limits
No exact rational audit of the bases; DROP is certified up to LP tolerance with margin 0.005 and on a finite atom grid (800, gate 1597). One cell per family, 100 datasets, no coverage-of-truth analysis beyond the controls. Walker budgets 8 and 128 per node were not run (not requested). The 5 2D INCONCLUSIVE datasets hit the 400-LP cap, not the interval floor; a larger cap was not tried.
