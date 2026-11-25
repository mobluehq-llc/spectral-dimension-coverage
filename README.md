# Do error bands for the scale-dependent spectral dimension cover the exact curve?

Reproducibility bundle for paper A, **version 1.1 (24 November 2026)**, matching the arXiv version 8 text. v1.0 (the version-6 text) is Zenodo DOI 10.5281/zenodo.24488979; the DOI of v1.1 is on the Zenodo record linked from this repository's release page. Paper: *More walkers, worse error bars: bias-driven undercoverage and structural bounds for the scale-dependent spectral dimension* (`paper/A_arxiv_v8.pdf`, supplement `paper/A_arxiv_supplement_v8.pdf`).

**Author:** Adam Faust, MOBLUEHQ, Inc. **Contact:** adam@mobluehq.com
**Citation:** concept DOI (all versions): [10.5281/zenodo.24488979](https://doi.org/10.5281/zenodo.24488979) resolves to the latest version. See `CITATION.cff`.

## What the paper is
Claims that the spectral dimension of a network or quantum-gravity model flows with scale rest on the log-derivative of a noisy return probability on a
finite graph. The paper tests error bands against a known truth, on tori and products of cycles whose exact d_s(t) follows from the spectrum, with bars fixed
before each run, using a frozen-code replication, a holdout, and corrective reanalyses. Gated and exploratory results are labelled as such. The pre-registered
plan with all twelve amendments is `plans/PLAN.md`; every number in the paper is mapped to a stored file in `CLAIMS.md`.

## What v1.1 adds to v1
v1 predated the split-covariance and random-graph runs and several later checks. v1.1 adds their code, plans, outputs and ledger rows:
- **BiasAware-Split (Amendment 13a)** and **random graphs (13b)**, and the step/regression arm (13c): `code/{p22amend13,run_amend13a,run_amend13b,analyze_amend13,a13c_scaling}.py`, `results/amend13/` (per-dataset raw files, manifests, calibration and coverage tables, reports), `plans/RESULTS_AMEND13.md`.
- **Mechanism of the budget effect (paper Sec. III.B)**: `code/bias_noise.py`, `code/rt_rbc.py`, `code/rt_test_rbc.py`, `code/rt_monolp.py`, `code/rt_analyze.py`, outputs `results/v10/bias_noise_*`, `rt_rbc*`, `rt_monolp*`, `plans/RESULTS_BIAS_NOISE.md`, `plans/RESULTS_REDTEAM_GEMINI.md`.
- **Global no-drop test, minimum compatible drop, continuum certificates (Sec. III.E)**: `code/rt_nodrop.py`, `rt_nodrop_analyze.py`, `rt_mindrop.py`, `rt_mindrop_analyze.py`, `rt_certify.py`, `rt_certify_falsifier.py`; outputs `results/v10/nodrop_*`, `mindrop_*`, `certify_*` (including `certify_2d_mpall.csv`, the Arb ball recheck of all 37 certificates, pre-registered as P8 in `plans/PREREG_REDTEAM_FINAL.md`), summaries `results/v10/v9_nodrop_numbers.txt`, `v10_mindrop_cert_numbers.txt`, `v8_mechanism_numbers.txt`; plans `plans/PREREG_GLOBAL_NODROP.md`, `PREREG_REDTEAM_*.md`, `RESULTS_GLOBAL_NODROP.md`, `RESULTS_MINDROP_CERT.md`.
- Paper sources of the arXiv version (`paper/main.tex`, `paper/supplement.tex`, tables, figures, and the built PDFs) and the scripts that build its numbers and figures (`code/paper_scripts/`); the complete prior-art note `plans/PRIORART_MONOLP.md`.
- `CLAIMS.md` now has 483 checked rows (v1: 349); `python3 -I verify.py CLAIMS.md` (or `make verify`) re-checks all of them against the stored evidence. Rows 462-488 cover the items above, including the P9 full ball-arithmetic recheck (`code/rt_certify_arbfull.py`, `results/v10/certify_2d_arbfull.*`, `certify_arbfull_falsifier.log`).

**Scope of the one-command regeneration.** `code/archive/regenerate_all.py` still runs the 19 analysis steps of v1 (it was not extended). The v1.1 outputs are verified by `verify.py` against the stored files and can be re-run with the scripts named above in the pinned environment; the per-dataset files are the record, because the same seed draws different walks under different NumPy versions (see below). Host names: "host T", "host D" and "host W" stand for three compute machines; their names and local paths are omitted. The LLM critiques of earlier drafts are not included.

## Reproduce (one command)
Requires `uv` and Python 3.12 (uv fetches it). From this directory:

```
code/archive/setup_env.sh && code/archive/.venv/bin/python -s code/archive/regenerate_all.py
python3 -I verify.py CLAIMS.md       # re-checks every numbered claim against the stored evidence (stdlib only)
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
- `paper/` main.tex, supplement.tex, refs.bib, figures (`figs/`), table includes (`tab_*.tex`, `v7tabs/`), the cross-reference files `xref_M.tex` and `xref_S.tex`, and the arXiv-version PDFs `A_arxiv_v8.pdf`, `A_arxiv_supplement_v8.pdf`. The v6 PDFs of v1.0 are in the v1.0 release.
- `code/` estimators and simulation: `p22lib.py` (graphs, walks, exact d_s(t)), `p22*.py` (estimators and amendment variants), `run_p22*.py` / `run_amend*.py` (drivers), `analyze*.py` and `verdict_amend10_11.py` (analysis and verdicts), `round3_*`/`round4_*`/`round5_*`/`round6_*` (referee-response checks), `paper_scripts/` (figure/table builders), `archive/` (`regenerate_all.py`, pinned lock, input manifest), `audit/`, `lean/`, `verify_sage_*.sage` (independent checks).
- `plans/` PLAN.md (all amendments, each written before its run), PREREG_*.md, RESULTS*.md, REGISTERED_NOT_RUN.md, IMPLEMENTATION_NOTES_A10_A11.md, REPRO_CHECK.md, ROUND4/5_CHECKS*.md, RESULTS.md, VERIFY_SAGE.md.
- `results/` stored tables and per-dataset files: `v10/` and `amend13/` (new in v1.1), `base/` (original replication and Amendments 1-9, 12), `amend10_11/`, `round3/`, `external/` (inputs not otherwise tracked: raw_amend10c/11 and E4 diagnosis), `reproduction/` (release-run report).
- `CLAIMS.md` every claim to its evidence file and check; `verify.py` the checker.

Scripts in `code/` that were run on the original machines (drivers, `repro_check.py`) refer to the original directory layout; `regenerate_all.py` rebuilds it in a work directory.

## Licences
Code (`code/`, `verify.py`): MIT, see `LICENSE`. Paper text, figures, plans and data: CC BY 4.0, see `LICENSE-DATA-TEXT`.
