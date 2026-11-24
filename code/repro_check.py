#!/usr/bin/env python3
"""Round-3 reproducibility check for paper A (p22): re-run ONE small cell of each of Amendments 9b, 10b, 10c, 11 with the code at the
current commit, in a fresh pinned venv, and compare to the stored raw results. Adds NO estimator, bar, cell definition or seed; it only
calls the runners' own per-dataset functions (run_amend1011 helpers) on a handful of datasets.
usage (cwd = checkout root, PYTHONPATH unset, venv python):  repro_check.py {9b|10b|10c|11} [--stored DIR] [--out DIR]
Stored raws: 9b results/raw_amend9_b32.csv ; 10b results_amend10_11/raw_amend10b.csv ; 10c/11 raw_amend10c.csv / raw_amend11.csv in --stored
(10c/11 raw CSVs are not in git; their sha256 is checked against the committed manifest before use).
Tolerance (documented): non-timing columns compared after CSV round trip. Strings and discrete columns (abstain, J, coverage flags, declarations, status, counts)
must be EXACT. Floats: EXACT counted separately; otherwise |a-b| <= 1e-9 + 1e-9*|b| counts as WITHIN_TOL (LP solver / BLAS ties); anything else is a MISMATCH.
Timing columns (seconds, t_*) are excluded."""
import os, sys, json, time, hashlib
ARGV = list(sys.argv); WHICH = ARGV[1]
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
STORED = ARGV[ARGV.index("--stored") + 1] if "--stored" in ARGV else os.path.join(HERE, "results")
OUT = ARGV[ARGV.index("--out") + 1] if "--out" in ARGV else os.path.join(HERE, "repro_out")
os.makedirs(OUT, exist_ok=True)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(v, "1")
sys.argv = [ARGV[0], "verify"]
import numpy as np, pandas as pd
import run_amend1011 as M
sys.argv = ARGV
R, RI, RJ, RK, A, A6, H = M.R, M.RI, M.RJ, M.RK, M.A, M.A6, M.H

import re
DISC_RE = re.compile(r"(^abstain$|^J$|^cov|^decl|status$|_neg$|_pos$|_amb$|_certnon$|nlp$|nbad$|^gate_fail$|^nowin$|^infeasible$|natoms$|^error$|^why$|^group$|^fam$|^kind$|^mode$|^J_fixed$|^wpn$|contains|^m_viol$|nopen$)")
def is_timing(c): return c == "seconds" or "seconds" in c or c in ("t_orig", "t_B", "t_A", "t_Bj")

def compare(new, old, keys, name):
    """Tiers: float columns EXACT / <=1e-9 / <=1e-6 / >1e-6 (relative+absolute); discrete (verdict-bearing) columns and strings must be EXACT."""
    new = pd.DataFrame(new)
    m = new.merge(old, on=keys, suffixes=("_n", "_o"), how="outer", indicator=True)
    mm = m[m["_merge"] == "both"]
    cols = [c for c in new.columns if c not in keys and c + "_o" in mm.columns and not is_timing(c)]
    T = dict(exact=0, tol1e9=0, tol1e6=0, big=0, both_nan=0, disc_exact=0, disc_bad=0, str_exact=0, str_bad=0)
    bigcols = {}; badcols = {}; worst = 0.0; worst_rel = 0.0
    for c in cols:
        disc = bool(DISC_RE.search(c))
        for x, y, e in zip(mm[c + "_n"].tolist(), mm[c + "_o"].tolist(), mm["est"].tolist() if "est" in mm else [""] * len(mm)):
            xn = x is None or (isinstance(x, float) and np.isnan(x)); yn = y is None or (isinstance(y, float) and np.isnan(y))
            if xn and yn: T["both_nan"] += 1; continue
            if xn != yn:
                (T.__setitem__("disc_bad", T["disc_bad"] + 1) if disc else T.__setitem__("big", T["big"] + 1)); (badcols if disc else bigcols).setdefault(f"{c}[{e}]", 0); (badcols if disc else bigcols)[f"{c}[{e}]"] += 1; continue
            if isinstance(x, str) or isinstance(y, str):
                if str(x).strip() == str(y).strip(): T["str_exact"] += 1; continue
                try:
                    xs = np.array([float(t) for t in str(x).split(";")]); ys = np.array([float(t) for t in str(y).split(";")]); d = np.abs(xs - ys)
                    if len(xs) == len(ys) and (d <= 1e-9 + 1e-9 * np.abs(ys)).all(): T["tol1e9"] += 1; worst = max(worst, float(d.max())); continue
                except Exception: pass
                T["str_bad"] += 1; badcols[f"{c}[{e}]"] = badcols.get(f"{c}[{e}]", 0) + 1; continue
            x = float(x); y = float(y)
            if disc:
                if x == y: T["disc_exact"] += 1
                else: T["disc_bad"] += 1; badcols[f"{c}[{e}]"] = badcols.get(f"{c}[{e}]", 0) + 1
                continue
            d = abs(x - y); rel = d / max(abs(y), 1e-300)
            if d == 0: T["exact"] += 1
            elif d <= 1e-9 + 1e-9 * abs(y): T["tol1e9"] += 1; worst = max(worst, d)
            elif d <= 1e-6 + 1e-6 * abs(y): T["tol1e6"] += 1; worst = max(worst, d)
            else: T["big"] += 1; bigcols[f"{c}[{e}]"] = bigcols.get(f"{c}[{e}]", 0) + 1; worst = max(worst, d); worst_rel = max(worst_rel, rel if abs(y) > 1e-6 else 0)
    rows_new = len(new); matched = len(mm); keys_ok = matched == rows_new and rows_new > 0
    disc_ok = T["disc_bad"] == 0 and T["str_bad"] == 0
    if keys_ok and disc_ok and T["tol1e6"] == 0 and T["big"] == 0: verdict = "MATCH_EXACT" if T["tol1e9"] == 0 else "MATCH_WITHIN_TOL_1e-9"
    elif keys_ok and disc_ok and T["big"] == 0: verdict = "MATCH_WITHIN_TOL_1e-6"
    elif keys_ok and disc_ok: verdict = "DISCRETE_MATCH_FLOAT_DIFFERS (verdict-bearing columns identical; solver-level floats differ > 1e-6)"
    else: verdict = "MISMATCH"
    return dict(amendment=name, verdict=verdict, rows_rerun=rows_new, rows_matched_to_stored=matched, counts=T, float_diffs_gt_1e6_by_column_est=bigcols, discrete_mismatch_by_column_est=badcols,
                worst_abs_diff_nondiscrete=worst, worst_rel_diff_gt1e6=worst_rel, columns_compared=len(cols), excluded_timing_columns=True)

def sha(p): return M.sha(p)

def r9b():
    RI.GROUP = "b32"; RI.SEED_OFF = 2000
    sp = RI.select(); s = [s for s in sp if s["kind"] == "control" and s["fam"] == "ring" and s["L"] == 128][0]
    R.prep(s, np.random.default_rng(2200)); RI.CELLS = {s["id"]: s}
    nd = 10; rows = []
    for ds in range(nd): rows += RI.run_dataset(s, ds)
    p = os.path.join(STORED, "raw_amend9_b32.csv"); old = pd.read_csv(p, low_memory=False); old = old[(old.cell == s["id"]) & (old.ds < nd)]
    r = compare(rows, old, ["cell", "ds", "est"], "9b"); r.update(cell=f"{s['fam']} L={s['L']} id={s['id']} wpn=32 ({s['kind']})", datasets=nd, stored_file=os.path.basename(p), stored_sha256=sha(p)); return r

def r10b():
    sel = json.load(open(os.path.join(HERE, "results_amend10_11", "amend10b_selection.json")))
    tasks = [(103, 32, "sens"), (103, 48, "sens")]; assert all([c, d] in sel["sens"] for c, d, _ in tasks)
    rows = [RJ.run_one(t) for t in tasks]
    p = os.path.join(HERE, "results_amend10_11", "raw_amend10b.csv"); old = pd.read_csv(p, low_memory=False)
    old = old.merge(pd.DataFrame([dict(cell=c, ds=d) for c, d, _ in tasks]), on=["cell", "ds"])
    new = pd.DataFrame(rows)
    r = compare(new, old, ["cell", "ds"], "10b"); c = RJ.get_cell(103)
    r.update(cell=f"{c['fam']} L={c['L']} id=103, datasets 32 and 48 (group sens)", datasets=len(tasks), stored_file="results_amend10_11/raw_amend10b.csv", stored_sha256=sha(p)); return r

def _raw_checked(fname, mname):
    p = os.path.join(STORED, fname); want = json.load(open(os.path.join(HERE, "results_amend10_11", mname)))["csv_sha256"]; got = sha(p)
    assert got == want, f"stored {fname} sha {got} != committed manifest {want}"
    return p, got

def r10c():
    p, h = _raw_checked("raw_amend10c.csv", "manifest_amend10c.json")
    wins = json.load(open(os.path.join(HERE, "results_amend10_11", "amend10c_windows.json")))["cells"]
    cc = M.cells_10c(); C = {}
    for s in cc[8]:
        if s["fam"] == "aniso2_w16" and s["L"] == 128:
            e = wins[f"{s['fam']}_{s['L']}"]; assert e["ids"]["8"] == s["id"]; s["J_fixed"] = int(e["exact_J_wpn32"]); C[s["id"]] = s
    assert len(C) == 1; M.C10 = C; s = list(C.values())[0]; M.g3cm(s)
    nd = 6; rows = M.run_10c((s["id"], 0, nd))
    old = pd.read_csv(p, low_memory=False); old = old[(old.cell == s["id"]) & (old.ds < nd)]
    r = compare(rows, old, ["cell", "ds", "est", "mode"], "10c")
    r.update(cell=f"{s['fam']} L={s['L']} id={s['id']} wpn=8, orig + fixed window J={s['J_fixed']}", datasets=nd, stored_file="raw_amend10c.csv (host-D, sha matches committed manifest)", stored_sha256=h); return r

def r11():
    p, h = _raw_checked("raw_amend11.csv", "manifest_amend11.json")
    sp = RK.specs(); rr = np.random.default_rng(2200)
    cand = [x for x in sp if x["wpn"] == 8]
    for x in cand: R.prep(x, rr)
    s = min(cand, key=lambda x: x["N"] * x["wpn"] * len(x["ts_w"]))
    RK.CELLS = {s["id"]: s}; nd = 6; rows = []
    for ds in range(nd): rows += RK.run_dataset(s, ds)
    old = pd.read_csv(p, low_memory=False); old = old[(old.cell == s["id"]) & (old.ds < nd)]
    r = compare(rows, old, ["cell", "ds", "est"], "11"); r.update(cell=f"{s['fam']} L={s['L']} id={s['id']} wpn={s['wpn']}", datasets=nd, stored_file="raw_amend11.csv (host-D, sha matches committed manifest)", stored_sha256=h); return r

if __name__ == "__main__":
    t0 = time.time(); out = os.path.join(OUT, f"repro_{WHICH}.json")
    if os.path.exists(out): print("already done", out); sys.exit(0)
    r = {"9b": r9b, "10b": r10b, "10c": r10c, "11": r11}[WHICH]()
    r["seconds"] = round(time.time() - t0, 1); r["finished"] = time.strftime("%F %T")
    json.dump(r, open(out, "w"), indent=1, sort_keys=True); print(json.dumps(r, indent=1, sort_keys=True))
