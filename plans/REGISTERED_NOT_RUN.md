# Checks registered and not run (moved from SM Sec. S4, arXiv v3)

(i) A finite-sample box (Clopper-Pearson at every grid time, fixed Bonferroni allocation) under the same certificate and for all four procedures.
(ii) Near-null controls: product graphs with exact drops of 0.15 to 0.35, chosen from the spectrum before any walk, to measure the declaration rate as the exact drop crosses 0.25 (PREREG_REDTEAM_FINAL.md, P4).
(iii) An ellipsoidal box from the walker-block covariance, a second-order cone programme, with a level check on controls first (P5).
(iv) A two-stage ensemble MonoLP, bootstrapping graphs and then walkers, for the spectral dimension of the mean return curve.
(v) For SM Sec. S7, robust bias correction with a bandwidth chosen for the derivative on fresh walks recorded at every second step, on C_256 x C_16 (this is the "first" and, with another time grid, the "last" of the untested items named there).
