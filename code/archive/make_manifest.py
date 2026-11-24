#!/usr/bin/env python3
"""Writes archive/data_manifest.json: sha256 of every per-dataset input (in git and external). Run once after the inputs are final."""
import os, json, hashlib, glob
HERE = os.path.dirname(os.path.abspath(__file__)); P22 = os.path.dirname(HERE)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
ing = {}
for pat in ("results/raw*.csv", "results/cells_meta*.csv", "results_amend10_11/raw_amend10b*.csv", "results_amend10_11/*.json", "results_amend10_11/cells_meta*.csv", "results_amend10_11/agg_*.csv"):
    for p in sorted(glob.glob(os.path.join(P22, pat))):
        rel = os.path.relpath(p, P22)
        if "smoke" in rel or "partial" in rel: continue
        ing[rel] = sha(p)
ext = {os.path.relpath(p, os.path.join(HERE, "data")): sha(p) for p in sorted(glob.glob(os.path.join(HERE, "data", "**", "*.csv"), recursive=True))}
json.dump(dict(in_git=ing, external=ext, note="in_git: committed per-dataset files (hash checked in the checkout); external: not in git, kept in archive/data (hash checked)"), open(os.path.join(HERE, "data_manifest.json"), "w"), indent=1, sort_keys=True)
print(len(ing), "in-git inputs,", len(ext), "external inputs")
