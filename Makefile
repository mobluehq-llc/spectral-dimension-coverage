# Copy this directory's Makefile and verify.py into an idea/paper repo root (next to CLAIMS.md).
# `make verify` re-runs every CLAIMS.md row's check from its evidence file; non-zero exit on any mismatch.
PY ?= python3
.PHONY: verify
verify:
	$(PY) -I verify.py CLAIMS.md
