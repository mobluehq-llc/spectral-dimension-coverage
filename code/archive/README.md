# Regeneration of every table and gate verdict (paper A)

Run from the bundle root:

```
code/archive/setup_env.sh                                 # once: Python 3.12 + exact versions in requirements.lock (needs uv)
code/archive/.venv/bin/python -s code/archive/regenerate_all.py
```

It (1) checks the interpreter against `requirements.lock`, (2) checks the sha256 of all 72 per-dataset inputs against `data_manifest.json`
(54 under `results/base` and `results/amend10_11`, 18 under `results/external`), (3) copies code and inputs to a work directory in the original
layout so the stored tables are never touched, (4) runs 19 analysis steps (original replication, Amendments 1, 3-9, 10a, 10b/10c/11, 12, and the
round-3 tables), (5) compares each regenerated table with the stored one (IDENTICAL / EQUAL-1e-9 / DIFFERS / NOT-WRITTEN; stale copies are
deleted before each step), and (6) writes `code/archive/out/` (`regeneration_report.csv`, `gate_verdicts.txt`, `gate_verdicts_10b_10c_11.csv`, one log
per step). Exit status 0 only if all hashes match, all steps exit 0 and no table differs. Run time about 30 s on one core. A copy of the report from the
release run is in `results/reproduction/`.

## What it does NOT do
* It regenerates tables and verdicts from per-dataset results. It does not re-run any simulation. Per-dataset files are platform-dependent
  (NumPy 2.4.6 and 2.5.5 draw different walks from the same seed), so re-simulating on another machine gives a different realisation of the same seeds.
  The files in the manifest are the record.
* Heavy compute outputs of round 3 are inputs: `results/external/e4diag/*` (rerun with `code/e4_oracle_diag.py`) and `results/round3/lp_verify_results*.json`
  (exact LP verification, rerun with `code/round3_lp_verify.py run`, about 35 min on 4-5 cores).
* The causal-set stress test and the track audit (`code/analyze_causets.py`, `code/audit/`) are not part of this regeneration.

## Pinned environment
See `requirements.lock`. NumPy 2.5.5, SciPy 1.18.1, pandas 3.1.1, python-flint 0.9.0, highspy 1.16.0, matplotlib 3.11.2 (figures only), Python 3.12.
The stored Amendment 3-9 results were simulated under NumPy 2.4.6 and SciPy 1.17.1; their analysis here is pure pandas/scipy on stored rows and reproduces them byte for byte.
