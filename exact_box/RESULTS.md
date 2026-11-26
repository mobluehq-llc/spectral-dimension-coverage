# BD RESULTS: exact finite-sample box, near-null controls, power map

Run 2026-11-25 on host D (WSL, `bd_chain.sh 10`, harness then a..f, all modes exit 0, "chain complete"). Runners hashed and committed before launch (PLAN Amendment 1; `src/SHA256SUMS`; host D copies re-checked: rt_cpbox.py e156eb4c..., analyzer 24f56d53...). Analysis by the hashed `rt_cpbox_analyze.py`; full tables `results/BD_tables.md`; raw rows (`bd_*.jsonl`, 8,947 lines, 0 ERROR rows) are in `data/` (this release).
**Status: CONFIRMATORY for bars 1-4 (pre-registered, unchanged). Arm (f) is EXPLORATORY.**

## Bars
1. **Harness: PASS.** CP box contains exact q at all K=28 times in 1917/2000 = 0.9585; one-sided CP lower 0.9504 (needs 0.94).
2. **Level: PASS.** False DROP 0/900 [0, 0.004] pooled over (a) controls and (b) sub-margin graphs. Box coverage over all (a) datasets 862/900 = 0.958 [0.943, 0.970]. Caveat: arm (b) has no power at all (below), so it adds little beyond the controls; the evidence for the exact level is the coverage of the box (bar 1 and the 0.958 above), as section 2 says.
3. **Power GO: PASS.** C768xC16, CP box, global: 1, 85, 100, 100 of 100 at 8, 32, 128, 512 per node. GO.
4. **Power map: PASS.** All 8 cells monotone in walkers (aniso2 L=320: 0,3,96,100; L=768: 1,85,100,100; aniso3 L=128: 8,100,100,100; L=256: 93,100,100,100); every cell reaches 100/100 at 512. Arm (d) (exact q at eps 1e-4, 1e-3, 1e-2): DROP in all 24 cells.

## Other arms
- (a) Same seeds, 32/node: CP box global DROP 85/100 (aniso2) against 37/100 for the stored nominal box; all 37 nominal DROPs also CP DROPs, 48 CP-only. So the exact box is MORE powerful here, not less (hw ratio CP/nominal about 0.85). aniso3: 100 vs 100. Joint, pointwise and box-secant tests are far weaker than the global test at 32/node (5, 3, 0 of 100 in 2D; at 128: 95, 95, 0). Controls (ring, torus2, torus3): 0/100 for all four procedures.
- (b) Near-null: **0 DROP in every graph at 32 and 128 per node, both boxes**, including the three super-margin graphs (exact drop 0.272, 0.306, 0.339). Windows are 2 to 7 grid times (N 287 to 1785), and 30-53 of 100 datasets on C39xC8 and C41xC7 had no window. So (b) measures the level on easy cases only; it says nothing about power near the margin, and it does not test the level just under 0.25 with a usable window.
- (e) All 423 CP-box global DROPs of (a) CERTIFIED by the continuum certificate; ball-arithmetic re-check of the seeded 10: 10/10 verified.
- (f) EXPLORATORY: exact q(1..3) as equalities raises aniso2 L=768 at 32/node from 85 to 100 of 100 (aniso3 100 to 100).

## Predictions (section 7) against outcome
2D at 32/node predicted 8-20, **observed 85 (wrong, in the favourable direction)**; at 128 predicted >= 50, observed 100. 3D >= 95 at 32: 100 (right). Near-null zero up to 0.23: right; "non-zero power only from 0.30": **wrong, zero at 0.30 and 0.34 too** (short windows). Arm (f) at most 3 more: **wrong, +15** (exploratory).

## Verdict
**GO** per section 8: bars 1-3 hold. The declaration level of the global no-drop test can be called exact (Clopper-Pearson box at level 0.05/K over the full grid), within the section-2 conditions (uniform start, independent walkers, lazy walk, connected graph; all BD cells are tori). Not shown: level or power for random graphs (all cells are tori, via `_direct`), or near the 0.25 margin with a usable window. Paper A is not edited by this run; a wording change goes through review first. Public release needs Adam's go.
