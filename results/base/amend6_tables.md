groups present: {'big2d': True, 'sep': True, 'fresh': True, 'orig': True}
## Item 1: minimum-size recipe (fit on 3D only; power = declaration rate, abstention counted as a miss)
| est | predicted_minL | measured_minL | index_distance | within_one_step | power_by_L |
|---|---|---|---|---|---|
| E0 | 128 | 192 | 1 | True | 32:0.00 48:0.00 64:0.01 96:0.05 128:0.25 192:0.99 256:1.00 |
| E1q | 192 | 192 | 0 | True | 32:0.00 48:0.00 64:0.00 96:0.01 128:0.17 192:0.83 256:0.94 |
| E6 | 192 | 256 | 1 | True | 32:0.00 48:0.00 64:0.00 96:0.01 128:0.01 192:0.59 256:0.98 |
| Mb | >256 | 256 | 1 | True | 32:0.00 48:0.00 64:0.00 96:0.00 128:0.23 192:0.74 256:0.87 |
| CL | 192 | 128 | 1 | True | 32:0.10 48:0.28 64:0.45 96:0.71 128:0.81 192:0.90 256:0.93 |
| E4c | 768 | >768 | 1 | True | 32:0.00 48:0.00 64:0.00 96:0.00 128:0.01 192:0.10 256:0.23 384:0.41 512:0.55 768:0.79 |
| E4r | 768 | 768 | 0 | True | 32:0.00 48:0.00 64:0.00 96:0.00 128:0.00 192:0.12 256:0.32 384:0.49 512:0.71 768:0.92 |

First five (consistency): 5 of 5 within one grid step. E4c/E4r (out of sample): 2 of 2.

| est | Jmin | n_fit | status | a | delta | Tlo | Thi | T | flag | Jfit |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 | 6 | 3 | ok | 0.369 | 1.864 | 1.0 | 4.881 | 2.209 |  | 7,9,11 |
| E1q | 6 | 3 | ok | 0.405 | 0.848 | 1.0 | 4.327 | 2.08 |  | 7,9,11 |
| E6 | 6 | 3 | ok | 1.277 | 2.307 | 1.0 | 1.417 | 1.19 |  | 7,9,11 |
| Mb | 8 | 2 | ok | 0.311 | -1.176 | 1.0 | 3.346 | 1.829 |  | 9,11 |
| CL | 5 | 4 | ok | 0.451 | 0.592 | 1.0 | 2.758 | 1.661 |  | 5,7,9,11 |
| E4c | 6 | 3 | ok | 1.551 | 1.459 | 1.608 | 2.238 | 1.897 |  | 7,9,11 |
| E4r | 6 | 3 | ok | 1.462 | 1.621 | 1.764 | 2.563 | 2.126 |  | 7,9,11 |

E4 half-width and power in the large 2D cells (hw_pred is the recipe's J-law, ungated):
| est | L | J | hw_pred | hw_meas | nonabs | power |
|---|---|---|---|---|---|---|
| E4c | 192 | 10 | 0.921 | 0.79 | 1.0 | 0.1 |
| E4c | 256 | 12 | 0.706 | 0.654 | 1.0 | 0.235 |
| E4c | 384 | 14 | 0.564 | 0.56 | 1.0 | 0.41 |
| E4c | 512 | 16 | 0.464 | 0.509 | 1.0 | 0.545 |
| E4c | 768 | 18 | 0.391 | 0.429 | 0.98 | 0.79 |
| E4r | 192 | 10 | 0.82 | 0.725 | 1.0 | 0.115 |
| E4r | 256 | 12 | 0.61 | 0.579 | 1.0 | 0.32 |
| E4r | 384 | 14 | 0.475 | 0.514 | 1.0 | 0.49 |
| E4r | 512 | 16 | 0.383 | 0.456 | 1.0 | 0.71 |
| E4r | 768 | 18 | 0.316 | 0.376 | 1.0 | 0.92 |

Coverage of exact d(t) in the large 2D cells (non-abstained):
| est | L | cov | n |
|---|---|---|---|
| E4c | 384 | 0.95 | 200 |
| E4c | 512 | 0.98 | 200 |
| E4c | 768 | 0.954 | 196 |
| E4r | 384 | 0.93 | 200 |
| E4r | 512 | 0.96 | 200 |
| E4r | 768 | 0.92 | 200 |
| E6 | 384 | 0.96 | 200 |
| E6 | 512 | 0.9 | 200 |
| E6 | 768 | 0.645 | 200 |

## Item 2: E6 misfit alarms (original seeds, 37 cells x 200)
E4c abstained where E6 did not: 4 datasets (no alarm).
Reproduction of E6 vs raw_amend5: 7400 datasets, abstain identical 1.0000, cov_simult identical 1.0000, max |hw diff| 0.00e+00
| bar | A1 | A2 | LoF |
|---|---|---|---|
| alarm | A1 E6-vs-E4c band disagreement | A2 one-vs-two-step (fallback) | LoF (Amendment 5, for reference) |
| a_pooled_P_alarm_given_miss | 0.7 (0.1-2.0)  [misses 441] | 36.5 (32.0-41.2)  [misses 441] | 18.8 (15.3-22.8)  [misses 441] |
| a_pass | False | False | False |
| b_twostep320_P_alarm_given_miss | 1.6 (0.3-4.5)  [misses 193] | 56.0 (48.7-63.1)  [misses 193] | 30.6 (24.2-37.6)  [misses 193] |
| b_pass | False | False | False |
| c_control_alarm | 0.0 (0.0-0.2) | 3.3 (2.6-4.2) | 3.7 (2.9-4.7) |
| c_pass | True | True | True |
| d_xover_power_E6 | 0.569 | 0.569 | 0.569 |
| d_power_E6g | 0.569 | 0.555 | 0.552 |
| d_loss_points | 0.0 | 1.4 | 1.6 |
| d_pass | True | True | True |
| ALL_PASS | False | False | False |

Per cell (non-abstained E6; a1/a2/lof = alarm rates):
| cell | group | n | cov | misses | a1 | a2 | lof | a1_given_miss | a2_given_miss |
|---|---|---|---|---|---|---|---|---|---|
| aniso2_w16_128 | xover_real | 200 | 1.0 | 0 | 0.0 | 0.02 | 0.05 |  |  |
| aniso2_w16_192 | xover_real | 200 | 0.995 | 1 | 0.0 | 0.02 | 0.025 | 0.0 | 0.0 |
| aniso2_w16_256 | xover_real | 200 | 0.99 | 2 | 0.0 | 0.01 | 0.045 | 0.0 | 0.0 |
| aniso2_w16_64 | xover_notreal | 5 | 1.0 | 0 | 0.0 | 0.0 | 0.0 |  |  |
| aniso2_w16_96 | xover_notreal | 198 | 1.0 | 0 | 0.0 | 0.03 | 0.02 |  |  |
| aniso3_w8_128 | xover_real | 200 | 0.985 | 3 | 0.0 | 0.02 | 0.03 | 0.0 | 0.0 |
| aniso3_w8_192 | xover_real | 200 | 0.98 | 4 | 0.0 | 0.05 | 0.04 | 0.0 | 0.0 |
| aniso3_w8_64 | xover_real | 32 | 1.0 | 0 | 0.0 | 0.0 | 0.0 |  |  |
| aniso3_w8_96 | xover_real | 199 | 1.0 | 0 | 0.0 | 0.015 | 0.005 |  |  |
| brush_f16_12 | brush | 200 | 0.965 | 7 | 0.0 | 0.07 | 0.035 | 0.0 | 0.0 |
| brush_f16_16 | brush | 200 | 0.895 | 21 | 0.0 | 0.045 | 0.055 | 0.0 | 0.048 |
| brush_f24_8 | brush | 200 | 0.94 | 12 | 0.0 | 0.07 | 0.06 | 0.0 | 0.083 |
| gasket_5 | gasket | 158 | 0.987 | 2 | 0.0 | 0.019 | 0.063 | 0.0 | 0.0 |
| gasket_6 | gasket | 200 | 0.98 | 4 | 0.0 | 0.065 | 0.05 | 0.0 | 0.0 |
| gasket_7 | gasket | 200 | 0.855 | 29 | 0.0 | 0.065 | 0.08 | 0.0 | 0.034 |
| ring_1024 | const | 200 | 0.995 | 1 | 0.0 | 0.075 | 0.075 | 0.0 | 0.0 |
| ring_128 | const | 127 | 0.992 | 1 | 0.0 | 0.008 | 0.008 | 0.0 | 0.0 |
| ring_256 | const | 200 | 0.99 | 2 | 0.0 | 0.04 | 0.025 | 0.0 | 0.0 |
| ring_512 | const | 200 | 0.995 | 1 | 0.0 | 0.06 | 0.065 | 0.0 | 0.0 |
| torus2_32 | const | 3 | 1.0 | 0 | 0.0 | 0.0 | 0.0 |  |  |
| torus2_48 | const | 200 | 0.995 | 1 | 0.0 | 0.02 | 0.025 | 0.0 | 0.0 |
| torus2_64 | const | 200 | 1.0 | 0 | 0.0 | 0.025 | 0.025 |  |  |
| torus2_96 | const | 200 | 1.0 | 0 | 0.0 | 0.045 | 0.045 |  |  |
| torus3_24 | const | 200 | 0.995 | 1 | 0.0 | 0.015 | 0.035 | 0.0 | 0.0 |
| torus3_28 | const | 200 | 0.995 | 1 | 0.0 | 0.02 | 0.02 | 0.0 | 0.0 |
| torus3_32 | const | 200 | 1.0 | 0 | 0.0 | 0.015 | 0.04 |  |  |
| twostep_8_32_128 | twostep | 200 | 0.995 | 1 | 0.0 | 0.03 | 0.03 | 0.0 | 0.0 |
| twostep_8_32_192 | twostep | 200 | 0.83 | 34 | 0.0 | 0.08 | 0.07 | 0.0 | 0.147 |
| twostep_8_32_256 | twostep | 200 | 0.31 | 138 | 0.0 | 0.285 | 0.075 | 0.0 | 0.326 |
| twostep_8_32_320 | twostep | 200 | 0.035 | 193 | 0.015 | 0.555 | 0.3 | 0.016 | 0.56 |

E4c (re-run) simultaneous coverage by group: brush: 97.7 (96.1-98.7); const: 99.3 (98.8-99.6); gasket: 98.7 (97.4-99.5); twostep: 95.2 (93.5-96.6); xover_notreal: 100.0 (98.2-100.0); xover_real: 97.3 (96.3-98.1)

## Item 3: fresh-seed replication (seeds [22, cell, 1000+ds])
| est | check | original | replicate | inside |
|---|---|---|---|---|
| E6 | coverage, realised crossover pooled | 99.2 (98.5-99.6) | 99.0 (98.3-99.5) | True |
| E6 | coverage, twostep_8_32 L=128 | 99.5 (97.2-100.0) | 99.0 (96.4-99.9) | True |
| E6 | coverage, twostep_8_32 L=192 | 83.0 (77.1-87.9) | 85.0 (79.3-89.6) | True |
| E6 | coverage, twostep_8_32 L=256 | 31.0 (24.7-37.9) | 26.0 (20.1-32.7) | True |
| E6 | coverage, twostep_8_32 L=320 | 3.5 (1.4-7.1) | 1.5 (0.3-4.3) | True |
| E6 | coverage, twostep pooled | 54.2 (50.7-57.7) | 52.9 (49.3-56.4) | True |
| E6 | coverage, brush_f24 Lb=8 | 94.0 (89.8-96.9) | 95.0 (91.0-97.6) | True |
| E6 | coverage, brush_f16 Lb=12 | 96.5 (92.9-98.6) | 95.5 (91.6-97.9) | True |
| E6 | coverage, brush_f16 Lb=16 | 89.5 (84.4-93.4) | 92.5 (87.9-95.7) | True |
| E6 | coverage, brush pooled | 93.3 (91.0-95.2) | 94.3 (92.2-96.0) | True |
| E6 | correct-sign power, realised crossover (abstain=miss) | 56.9 (54.4-59.3) | 55.6 (53.1-58.0) | True |
| E4c | coverage, realised crossover pooled | 96.4 (95.2-97.4) | 96.9 (95.8-97.8) | True |
| E4c | coverage, twostep_8_32 L=128 | 97.0 (93.6-98.9) | 96.0 (92.3-98.3) | True |
| E4c | coverage, twostep_8_32 L=192 | 95.5 (91.6-97.9) | 93.0 (88.5-96.1) | True |
| E4c | coverage, twostep_8_32 L=256 | 93.5 (89.1-96.5) | 94.5 (90.3-97.2) | True |
| E4c | coverage, twostep_8_32 L=320 | 95.0 (91.0-97.6) | 96.0 (92.3-98.3) | True |
| E4c | coverage, twostep pooled | 95.2 (93.5-96.6) | 94.9 (93.1-96.3) | True |
| E4c | coverage, brush_f24 Lb=8 | 98.5 (95.7-99.7) | 96.5 (92.9-98.6) | True |
| E4c | coverage, brush_f16 Lb=12 | 98.0 (95.0-99.5) | 99.0 (96.4-99.9) | True |
| E4c | coverage, brush_f16 Lb=16 | 96.5 (92.9-98.6) | 97.5 (94.3-99.2) | True |
| E4c | coverage, brush pooled | 97.7 (96.1-98.7) | 97.7 (96.1-98.7) | True |
| E4c | correct-sign power, realised crossover (abstain=miss) | 25.0 (22.9-27.2) | 23.8 (21.7-25.9) | True |

E6 gated checks inside the original CP interval: 11 of 11.  E4c (ungated): 11 of 11.

Fresh seeds, descriptive alarm check (twostep+brush only, no gasket/controls): misses 411; P(A1|miss) = 0.0 (0.0-0.9); P(A2|miss) = 39.2 (34.4-44.1)
  twostep_8_32 L=320: misses 197, P(A1|miss) = 0.0 (0.0-1.9), P(A2|miss) = 58.9 (51.7-65.8)

## Step-separation scan (twostep C_256 x C_8 x C_b; E6 simultaneous coverage of exact d(t), non-abstained; b=32 is original+fresh pooled)
| b | ratio | est | n | abstain | cov | alarm_A1 | alarm_A2 | hw |
|---|---|---|---|---|---|---|---|---|
| 12 | 1.5 | E4c | 200 | 0.0 | 94.0 (89.8-96.9) |  |  | 0.623 |
| 16 | 2.0 | E4c | 200 | 0.0 | 97.0 (93.6-98.9) |  |  | 0.519 |
| 20 | 2.5 | E4c | 200 | 0.0 | 89.0 (83.8-93.0) |  |  | 0.483 |
| 24 | 3.0 | E4c | 200 | 0.0 | 89.5 (84.4-93.4) |  |  | 0.498 |
| 32 | 4.0 | E4c | 398 | 0.005 | 94.0 (91.2-96.1) |  |  | 0.518 |
| 12 | 1.5 | E6 | 200 | 0.0 | 95.5 (91.6-97.9) | 0.0 | 0.015 | 0.341 |
| 16 | 2.0 | E6 | 200 | 0.0 | 92.5 (87.9-95.7) | 0.0 | 0.03 | 0.332 |
| 20 | 2.5 | E6 | 200 | 0.0 | 42.0 (35.1-49.2) | 0.0 | 0.11 | 0.304 |
| 24 | 3.0 | E6 | 200 | 0.0 | 12.0 (7.8-17.3) | 0.0 | 0.215 | 0.274 |
| 32 | 4.0 | E6 | 400 | 0.0 | 28.5 (24.1-33.2) | 0.0 | 0.275 | 0.238 |

E6: b with CP-lower >= 93%: []; b with point estimate < 93%: [16, 20, 24, 32]; b with CP-upper < 93%: [20, 24, 32]