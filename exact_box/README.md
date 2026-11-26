# Exact Clopper-Pearson box run (paper A, Sec. III.E, Table III and Supplement Sec. S4)

New in bundle v1.2. v1.1 regenerates every table of the v1 analyses (`code/archive/regenerate_all.py`); it did not contain the
exact-box run, which was registered and run after v1.1. This folder holds it, with the later ball-arithmetic and window checks.

* `PLAN.md` - the registration (committed, and the runners hashed, before launch).
* `src/` - `rt_cpbox.py` (runner, modes harness, a..f), `rt_cpbox_analyze.py` (hashed analyzer), `bd_chain.sh`; `SHA256SUMS`
  holds the pre-launch hashes of these three, which still match (their comments keep the original host and path names).
* `BD_tables.md`, `bd_b_graphs.json`, `bd_harness_summary.json` - the result tables behind Table III and Table S-cpmap.
* `RESULTS.md` - bars, outcomes and failed predictions.
* `data/` - the 7 raw jsonl files of the run, 8,947 rows, 0 ERROR rows; `data/SHA256SUMS` holds the sha256 computed on the run host.
* `src/rt_cpbox_arball.py` -> `results/bd_arball.{jsonl,log}`: full-chain Arb ball recheck (the P9 check of
  `code/rt_certify_arbfull.py`) of every DROP outcome of the run, 2,330 of them (the 423 of Table III and all DROPs of the map).
* `src/rt_cpbox_arball_rig.py` -> `results/bd_arball_rig.{jsonl,log}`: the same recheck with the Clopper-Pearson endpoints widened
  outward and the widening proven sufficient in Arb ball arithmetic; 2,330/2,330 CERTIFIED.
* `src/rt_null_windows.py` -> `results/null_windows.jsonl`: the exact drop over the window the global test actually used, for each
  of the 900 level datasets (the 0/900 false DROP count of bar 2; max exact drop 0.229, 71 datasets with no window). The file
  in v1.2 was regenerated for this release (on host W, and again on host T with the same windows); its window lengths `Jw` agree
  with those recorded by the run itself (`data/bd_a.jsonl`, `bd_b.jsonl`) in 900/900 datasets. An earlier copy disagreed in 63
  control datasets (rings, cells 3 and 8; 14 tested drops); the counts above are the same in both.

Host names "host D", "host T" and "host W" stand for compute machines, as elsewhere in the bundle. The text of `PLAN.md` and
`RESULTS.md` differs from the committed originals only in those names.

## Regenerate the tables
With `data/*`, `bd_harness_summary.json`, `bd_b_graphs.json` and `../results/v10/nodrop_{data,controls}.csv` in one `results/` folder:

```
uv run --no-project --python 3.12 --with numpy==2.5.5,pandas==3.1.1,scipy==1.18.1,tabulate==0.10.0 python src/rt_cpbox_analyze.py results
```

writes a `BD_tables.md` byte-identical to the stored one (seconds; the analyzer needs `tabulate`, not in v1.1's requirements.lock).
`python3 -I verify.py CLAIMS.md` re-checks rows 489 onward against the stored files here.

## Re-run the checks
The three check scripts import the frozen `code/` files and `src/rt_cpbox.py`; run them from a directory holding `code/*.py`,
`src/*.py` and `results/` with the raw `bd_*.jsonl` (NumPy 2.5.5 / SciPy 1.18.1 / python-flint 0.9.0, as pinned). The datasets are
regenerated from the frozen seeds, so other NumPy versions draw other walks.
