#!/usr/bin/env python3
"""ONE COMMAND that regenerates every table and gate verdict of paper A (p22) from the per-dataset result files.

    code/archive/.venv/bin/python -s code/archive/regenerate_all.py [--data DIR] [--out DIR] [--steps a,b,...] [--keep-work]

Nothing in the committed tree is modified: the analysis code and the per-dataset inputs are copied to a work directory, every analysis script is run
there, and each file it writes is compared with the committed copy (byte-identical / numerically equal at 1e-9 / DIFFERS / NEW). Gate verdicts are
collected from the scripts' own output. No estimator is re-run and no simulation is repeated: this reproduces TABLES and VERDICTS from per-dataset results.
Re-running the simulations themselves is a separate step (see archive/README.md, 'Regenerating the per-dataset files').

--data DIR  directory with the per-dataset files that are NOT in git (default results/external): raw_amend10c.csv, raw_amend11.csv, and the round-3 compute
            outputs e4diag_*.csv.  Their sha256 values are checked against archive/data_manifest.json before anything runs (mismatch = exit 2).
Exit status: 0 only if every input hash matches, every step exits 0, and no regenerated table DIFFERS from the committed one.
"""
import os, sys, json, hashlib, shutil, subprocess, time, re, csv, math
HERE = os.path.dirname(os.path.abspath(__file__)); P22 = os.path.dirname(HERE); BUNDLE = os.path.dirname(P22)
# Bundle layout: code in code/, stored results in results/{base,amend10_11,round3,external}. The scripts expect the original layout
# (results/, results_amend10_11/, results_round3/ beside the code), which is rebuilt in the work directory.
RESMAP = {"results/": "base/", "results_amend10_11/": "amend10_11/", "results_round3/": "round3/"}
def committed(rel):
    for k, v in RESMAP.items():
        if rel.startswith(k): return os.path.join(BUNDLE, "results", v + rel[len(k):])
    return os.path.join(P22, rel)

args = sys.argv[1:]
def opt(n, d=None): return args[args.index(n) + 1] if n in args else d
DATA = os.path.abspath(opt("--data", os.path.join(BUNDLE, "results", "external"))); OUT = os.path.abspath(opt("--out", os.path.join(HERE, "out")))
ONLY = opt("--steps"); KEEP = "--keep-work" in args
WORK = os.path.join(OUT, "work")
PY = sys.executable


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def env_check():
    lock = {l.split("==")[0].lower().replace("_", "-"): l.strip().split("==")[1] for l in open(os.path.join(HERE, "requirements.lock")) if "==" in l}
    import importlib.metadata as md
    bad = []
    for k, v in lock.items():
        try: got = md.version(k)
        except Exception: got = None
        if got != v: bad.append(f"{k}: want {v} got {got}")
    print("environment:", "matches requirements.lock" if not bad else "MISMATCH " + "; ".join(bad), "| python", sys.version.split()[0])
    return not bad


# (name, script (relative to p22 dir), argv, file globs it writes (relative to work dir), needs)
STEPS = [
    ("replication", "analyze.py", [], ["results/cell_summary.csv", "results/power_table.csv", "results/verdicts.csv", "results/ablation_parity.csv", "results/convergence_halves.csv"]),
    ("amend1", "analyze_amend.py", [], ["results/amend_final_table.csv", "results/amend_cell_summary.csv"]),
    ("amend3", "analyze_amend3.py", [], ["results/amend3_final_table.csv", "results/amend3_cell_summary.csv"]),
    ("amend4", "analyze_amend4.py", [], ["results/amend4_final_table.csv", "results/amend4_cell_summary.csv"]),
    ("amend5", "analyze_amend5.py", [], ["results/amend5_e6_cells.csv", "results/amend5_passbar.csv", "results/amend5_groups_all_estimators.csv", "results/amend5_tables.md"]),
    ("amend6", "analyze_amend6.py", [], ["results/amend6_alarm_bars.csv", "results/amend6_alarm_cells.csv", "results/amend6_replication.csv", "results/amend6_separation_scan.csv", "results/amend6_tables.md"]),
    ("amend7", "analyze_amend7.py", [], ["results/amend7_calib.csv", "results/amend7_minL.csv", "results/amend7_power.csv", "results/amend7_scores.csv"]),
    ("amend8", "analyze_amend8.py", [], ["results/amend8_cell_table.csv", "results/amend8_report.txt"]),
    ("amend8_table_iv", "paper_scripts/make_e7_outcomes.py", [], ["results/amend8_e7_outcomes.csv"]),
    ("amend9", "analyze_amend9.py", [], ["results/amend9_b_table.csv", "results/amend9_b_minL.csv", "results/amend9_c_table.csv", "results/amend9a_power.csv", "results/amend9_report.txt"]),
    ("amend10a", "analyze_amend10.py", ["a"], ["results/amend10a_table.csv", "results/amend10a_report.txt"]),
    ("amend10b_10c_11", "verdict_amend10_11.py", ["--raw", "{WORK}"], ["results_amend10_11/verdicts.csv", "results_amend10_11/verdict_report.txt", "results_amend10_11/tables_10b_cell_outcomes.csv",
                                                                     "results_amend10_11/tables_10b_3d_power.csv", "results_amend10_11/tables_10c_pooled.csv", "results_amend10_11/tables_10c_cells_fixed_vs_orig.csv",
                                                                     "results_amend10_11/tables_10c_e4_diagnosis.csv", "results_amend10_11/tables_11_cells.csv"]),
    ("amend12", "analyze_amend12.py", [], ["results/amend12_table.csv", "results/amend12_minL.csv", "results/amend12_report.txt"]),
    ("round3_e7_crosswalk", "round3_e7_crosswalk.py", ["--out", "results_round3"], ["results_round3/r3_e7_crosswalk_per_dataset.csv", "results_round3/r3_e7_crosswalk_summary.csv"]),
    ("round3_lp_tolerance", "round3_lp_verify.py", ["tolerance", "--out", "results_round3"], ["results_round3/lp_tolerance_margins.csv", "results_round3/lp_tolerance_counters.json"]),
    ("round3_lp_summary", "round3_lp_summary.py", ["--out", "results_round3"], ["results_round3/lp_verify_summary.csv", "results_round3/lp_verify_totals.json"]),
    ("round3_10c_11_tables", "round3_10c_tables.py", ["--raw", "{WORK}", "--out", "results_round3"], ["results_round3/r3_10c_cell_arm_budget_window.csv", "results_round3/r3_10c_pooled.csv",
                                                                                                   "results_round3/r3_11_cell_arm_budget.csv", "results_round3/r3_11_G5_denominators.csv", "results_round3/r3_11_E7_solver_accounting.csv"]),
    ("round3_native_targets", "round3_native_targets.py", ["--out", "results_round3"], ["results_round3/r3_native_E0_E2_all_cells.csv", "results_round3/r3_native_E0_E2_threshold_brackets.csv"]),
    ("round3_e4_oracle_diag", "e4_oracle_diag_report.py", ["--in", "{WORK}/e4diag", "--out", "results_round3"], ["results_round3/e4diag_summary.csv", "results_round3/e4diag_report.txt", "results_round3/e4diag_replay_each_miss.csv"]),
]
NEEDS = {"amend10b_10c_11": ["raw_amend10c.csv", "raw_amend11.csv"], "round3_10c_11_tables": ["raw_amend10c.csv", "raw_amend11.csv"], "round3_e4_oracle_diag": ["e4diag/e4diag_replay_cell0.csv"]}


def numeq(a, b, tol=1e-9):
    try:
        ra = list(csv.reader(open(a, newline=""))); rb = list(csv.reader(open(b, newline="")))
    except Exception: return None
    if len(ra) != len(rb): return False
    for x, y in zip(ra, rb):
        if len(x) != len(y): return False
        for u, v in zip(x, y):
            if u == v: continue
            try:
                fu, fv = float(u), float(v)
                if not (math.isclose(fu, fv, rel_tol=tol, abs_tol=tol) or (math.isnan(fu) and math.isnan(fv))): return False
            except ValueError: return False
    return True


def compare(rel, work):
    a, b = committed(rel), os.path.join(work, rel)
    if not os.path.exists(b): return "NOT-WRITTEN"
    if not os.path.exists(a): return "NEW"
    if sha(a) == sha(b): return "IDENTICAL"
    if rel.endswith(".csv"):
        r = numeq(a, b)
        if r: return "EQUAL-1e-9"
    return "DIFFERS"


def stage_work():
    if os.path.exists(WORK): shutil.rmtree(WORK)
    ign = shutil.ignore_patterns("archive", "__pycache__", "lean", "repro_out", "causet*", "*.pdf", "figs", ".git", "audit", "round4_checks", "round5_checks", "round6_checks")
    shutil.copytree(P22, WORK, ignore=ign, symlinks=False)
    for k, v in RESMAP.items():
        shutil.copytree(os.path.join(BUNDLE, "results", v), os.path.join(WORK, k.rstrip("/")), ignore=shutil.ignore_patterns("__pycache__"))
    for name in ("raw_amend10c.csv", "raw_amend11.csv"):
        shutil.copy(os.path.join(DATA, name), os.path.join(WORK, name))
    for f in os.listdir(os.path.join(BUNDLE, "results", "amend10_11")):
        if f.startswith("raw_amend10b"): shutil.copy(os.path.join(BUNDLE, "results", "amend10_11", f), os.path.join(WORK, f))
    if os.path.isdir(os.path.join(DATA, "e4diag")): shutil.copytree(os.path.join(DATA, "e4diag"), os.path.join(WORK, "e4diag"))
    os.makedirs(os.path.join(WORK, "results_round3"), exist_ok=True)


def main():
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    env_ok = env_check()
    # 1. inputs: hash check (data outside git) and committed per-dataset files
    man = json.load(open(os.path.join(HERE, "data_manifest.json")))
    bad = []
    for rel, h in man["external"].items():
        p = os.path.join(DATA, rel)
        if not os.path.exists(p): bad.append(f"MISSING {rel}"); continue
        if sha(p) != h: bad.append(f"MISMATCH {rel} want {h[:12]} got {sha(p)[:12]}")
    for rel, h in man["in_git"].items():
        p = committed(rel)
        if not os.path.exists(p): bad.append(f"MISSING {rel}"); continue
        if sha(p) != h: bad.append(f"MISMATCH {rel} want {h[:12]} got {sha(p)[:12]}")
    if bad:
        print("INPUT CHECK FAILED:\n  " + "\n  ".join(bad)); sys.exit(2)
    print(f"inputs: {len(man['external'])} external + {len(man['in_git'])} in-git per-dataset files, all sha256 match")
    stage_work()
    rows = []; gates = []; fail = False
    for name, script, argv, outs in STEPS:
        if ONLY and name not in ONLY.split(","): continue
        miss = [n for n in NEEDS.get(name, []) if not os.path.exists(os.path.join(WORK, n))]
        if miss: rows.append((name, "SKIPPED (missing " + ",".join(miss) + ")", 0, [])); fail = True; continue
        for o in outs:                                   # delete stale copies so a file that is not rewritten shows as NOT-WRITTEN
            try: os.remove(os.path.join(WORK, o))
            except FileNotFoundError: pass
        cmd = [PY, "-E", "-s", os.path.join(WORK, script)] + [a.replace("{WORK}", WORK) for a in argv]
        t1 = time.time(); r = subprocess.run(cmd, cwd=WORK, capture_output=True, text=True)
        log = os.path.join(OUT, f"log_{name}.txt"); open(log, "w").write(r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr.strip() else ""))
        cmp_ = {o: compare(o, WORK) for o in outs}
        for line in r.stdout.splitlines():
            if re.search(r"\[(PASS|FAIL|REPORTED)\]|VERDICT|PASS_CF|bar (met|NOT met)", line): gates.append((name, line.strip()))
        ok = r.returncode == 0 and all(v in ("IDENTICAL", "EQUAL-1e-9", "NEW") for v in cmp_.values())
        fail |= not ok
        rows.append((name, "ok" if ok else f"exit={r.returncode} " + ",".join(f"{os.path.basename(k)}:{v}" for k, v in cmp_.items() if v not in ("IDENTICAL", "EQUAL-1e-9", "NEW")), time.time() - t1, cmp_))
        print(f"[{'ok' if ok else 'FAIL'}] {name:24s} {time.time()-t1:5.1f}s  " + " ".join(f"{os.path.basename(k)}={v}" for k, v in cmp_.items()))
    # gate verdicts file (everything the scripts print as a verdict, with its step)
    with open(os.path.join(OUT, "gate_verdicts.txt"), "w") as f:
        for n, l in gates: f.write(f"{n}\t{l}\n")
    vp = os.path.join(WORK, "results_amend10_11", "verdicts.csv")
    if os.path.exists(vp): shutil.copy(vp, os.path.join(OUT, "gate_verdicts_10b_10c_11.csv"))
    with open(os.path.join(OUT, "regeneration_report.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step", "status", "seconds", "file", "comparison"])
        for n, s, t, c in rows:
            for k, v in (c.items() if isinstance(c, dict) and c else [("", "")]): w.writerow([n, s, f"{t:.1f}", k, v])
    print(f"{len(gates)} verdict lines collected -> {os.path.join(OUT, 'gate_verdicts.txt')}; total {time.time()-t0:.0f}s; environment {'pinned' if env_ok else 'NOT PINNED'}")
    if not KEEP: shutil.rmtree(WORK, ignore_errors=True)
    print("REGENERATION", "COMPLETE: all steps ran and every table matches the committed one" if not fail else "INCOMPLETE OR DIFFERS: see regeneration_report.csv")
    sys.exit(0 if not fail else 1)


if __name__ == "__main__": main()
