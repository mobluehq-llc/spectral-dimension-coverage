# Do error bands for the scale-dependent spectral dimension cover the exact curve?

Reproducibility bundle for paper A (version 6, 24 November 2026): *A benchmark on finite graphs with a frozen replication and a holdout.*

**Author:** Adam Faust, MOBLUEHQ, Inc. **Contact:** adam@mobluehq.com
**Citation:** DOI: [10.5281/zenodo.24488979](https://doi.org/10.5281/zenodo.24488979). See `CITATION.cff`.

## What the paper is
Claims that the spectral dimension of a network or quantum-gravity model flows with scale rest on the log-derivative of a noisy return probability on a
finite graph. The paper tests error bands against a known truth, on tori and products of cycles whose exact d_s(t) follows from the spectrum, with bars fixed
before each run, using a frozen-code replication, a holdout, and corrective reanalyses. Gated and exploratory results are labelled as such. The pre-registered
plan with all twelve amendments is `plans/PLAN.md`; every number in the paper is mapped to a stored file in `CLAIMS.md`.

## Reproduce (one command)
Requires `uv` and Python 3.12 (uv fetches it). From this directory:

```
code/archive/setup_env.sh && code/archive/.venv/bin/python -s code/archive/regenerate_all.py
python3 verify.py CLAIMS.md          # re-checks every numbered claim against the stored evidence (stdlib only)
```

`regenerate_all.py` verifies the sha256 of the 72 stored per-dataset inputs, re-runs the 19 analysis steps, and requires every regenerated table to match
the stored one; it also writes every gate verdict to `code/archive/out/`. About 30 s on one core. The release run reproduced all 61 stored tables
IDENTICAL and 44 verdict lines; see `results/reproduction/`. This reproduces tables and verdicts from stored per-dataset results; it does not repeat the simulations.

## Environment pinning and the NumPy-version note
Exact versions are in `requirements.txt` and `code/archive/requirements.lock` (NumPy 2.5.5, SciPy 1.18.1, pandas 3.1.1, python-flint 0.9.0, highspy 1.16.0,
matplotlib 3.11.2, Python 3.12). **NumPy note:** the stored Amendment 3-9 simulations were run under NumPy 2.4.6 / SciPy 1.17.1 (host T), and the Amendment 10b/10c/11/12
simulations under NumPy 2.5.5 / SciPy 1.18.1 (host D). The same seed draws different random walks under the two versions, so a fresh re-simulation is a new
realisation, not a bit-copy; this is why the stored per-dataset files, not the seeds, are the record. `plans/REPRO_CHECK.md` and `plans/ROUND4_CHECKS.md`
document the spot checks of this. In the plans and logs, "host T" and "host D" stand for the two compute machines used (names replaced for the release);
paths to those machines are not reproducible by readers and are not needed.

## File map
- `paper/` main.tex, supplement.tex, refs.bib, figures (`figs/`), table includes (`tab_*.tex`, `round3_tables_numbers.txt`), and the v6 PDFs `A_spectral_dimension_v6.pdf`, `A_supplement_v6.pdf`.
- `code/` estimators and simulation: `p22lib.py` (graphs, walks, exact d_s(t)), `p22*.py` (estimators and amendment variants), `run_p22*.py` / `run_amend*.py` (drivers), `analyze*.py` and `verdict_amend10_11.py` (analysis and verdicts), `round3_*`/`round4_*`/`round5_*`/`round6_*` (referee-response checks), `paper_scripts/` (figure/table builders), `archive/` (`regenerate_all.py`, pinned lock, input manifest), `audit/`, `lean/`, `verify_sage_*.sage` (independent checks).
- `plans/` PLAN.md (all amendments, each written before its run), IMPLEMENTATION_NOTES_A10_A11.md, REPRO_CHECK.md, ROUND4/5_CHECKS*.md, RESULTS.md, VERIFY_SAGE.md.
- `results/` stored tables and per-dataset files: `base/` (original replication and Amendments 1-9, 12), `amend10_11/`, `round3/`, `external/` (inputs not otherwise tracked: raw_amend10c/11 and E4 diagnosis), `reproduction/` (release-run report).
- `CLAIMS.md` every claim to its evidence file and check; `verify.py` the checker.

Scripts in `code/` that were run on the original machines (drivers, `repro_check.py`) refer to the original directory layout; `regenerate_all.py` rebuilds it in a work directory.

## Licences
Code (`code/`, `verify.py`): MIT, see `LICENSE`. Paper text, figures, plans and data: CC BY 4.0, see `LICENSE-DATA-TEXT`.
