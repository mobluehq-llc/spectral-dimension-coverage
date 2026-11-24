#!/usr/bin/env python3
"""Paper A v4: write the LaTeX bodies of the round-3 tables straight from the committed CSVs in ../results_round3/,
so that no number in those tables is typed by hand.  No analysis is done here beyond pooling counts.
Usage: python3 make_round3_tables.py   (writes tab_r3_*.tex and round3_tables_numbers.txt next to this file)"""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
R3 = os.path.join(HERE, "..", "results_round3")


def rows(name):
    return list(csv.DictReader(open(os.path.join(R3, name))))


def pct(x, nd=1):
    return f"{100 * float(x):.{nd}f}"


def write(name, lines):
    with open(os.path.join(HERE, name), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", name, len(lines), "rows")


# ---------------------------------------------------------------- 10c pooled (both families)
p = [x for x in rows("r3_10c_pooled.csv") if x["fam"] == "both"]
out = []
for est in ("E1q", "E4c", "E4r"):
    for mode, lab in (("fixed", "fixed"), ("orig", "own")):
        for w in ("8", "32", "128"):
            x = [y for y in p if y["est"] == est and y["mode"] == mode and y["wpn"] == w]
            assert len(x) == 1; x = x[0]
            assert int(x["n"]) == 3600 and int(x["answered"]) == 3600 - int(x["abstain"])
            out.append(f"{est} & {lab} & {w} & {x['abstain']} & {x['answered']} & {x['k_cov']} & "
                       f"{pct(x['cov_cond'])} ({pct(x['cov_cp_lo'])}--{pct(x['cov_cp_hi'])}) & {pct(x['cov_strict'])} & "
                       f"{x['k_drop']} & {pct(x['power_ans'])} ({pct(x['power_cp_lo'])}--{pct(x['power_cp_hi'])}) & {pct(x['power_all'])} \\\\")
write("tab_r3_10c_pooled.tex", out)

# ---------------------------------------------------------------- 10c per cell, one table per arm
c = rows("r3_10c_cell_arm_budget_window.csv")
fams = (("aniso2_w16", r"$C_{%s}\times C_{16}$"), ("aniso3_w8", r"$C_{%s}\times C_8^2$"))
for est in ("E1q", "E4c", "E4r"):
    out = []
    for fam, fmt in fams:
        Ls = sorted({int(x["L"]) for x in c if x["fam"] == fam})
        for L in Ls:
            cells = []
            for mode in ("fixed", "orig"):
                for w in ("8", "32", "128"):
                    x = [y for y in c if y["fam"] == fam and int(y["L"]) == L and y["est"] == est and y["mode"] == mode and y["wpn"] == w]
                    assert len(x) == 1; x = x[0]
                    assert int(x["n"]) == 200
                    cells.append(f"{x['k_cov']}/{x['answered']} & {x['k_drop']}")
                    J = int(float(x["J_fixed"]))
            out.append((fmt % L) + f" & {J} & " + " & ".join(cells) + r" \\")
    write(f"tab_r3_10c_cells_{est}.tex", out)

# ---------------------------------------------------------------- native targets for E0 and E2 (threshold brackets)
nb = rows("r3_native_E0_E2_threshold_brackets.csv")
out = []
for x in nb:
    fam = r"2D" if x["fam"] == "aniso2_w16" else r"3D"
    a = int(x["answered"])
    if a == 0:
        cp_, cn = "--", "--"
    else:
        cp_ = f"{x['k_cov_pointwise']}/{a} = {pct(x['cov_pointwise'])}"
        cn = f"{x['k_cov_native']}/{a} = {pct(x['cov_native'])} ({pct(x['cov_native_lo'])}--{pct(x['cov_native_hi'])})"
    out.append(f"{x['est']} & {fam} & {x['L']} & {a} & {cp_} & {cn} & "
               f"{x['k_drop']}/200 = {pct(x['power_all'])} ({pct(x['power_all_lo'])}--{pct(x['power_all_hi'])}) \\\\")
write("tab_r3_native.tex", out)

# pooled native-target coverage over all main-benchmark cells (counts summed, nothing else)
na = rows("r3_native_E0_E2_all_cells.csv")
notes = []
for est in ("E0", "E2"):
    for kind, real, lab in (("control", "False", "controls"), ("crossover", "True", "realised crossovers")):
        s = [x for x in na if x["est"] == est and x["kind"] == kind and x["real"] == real]
        a = sum(int(x["answered"]) for x in s); kp = sum(int(x["k_cov_pointwise"]) for x in s); kn = sum(int(x["k_cov_native"]) for x in s)
        notes.append(f"{est} {lab}: pointwise {kp}/{a} = {100 * kp / a:.1f}%   native {kn}/{a} = {100 * kn / a:.1f}%   (from r3_native_E0_E2_all_cells.csv)")

# ---------------------------------------------------------------- E4 oracle diagnosis
e = rows("e4diag_summary.csv")
out = []
tot_ans = tot_bound = 0
for arm, lab in (("replay", "original"), ("fresh", "fresh")):
    for L in ("128", "256", "512", "1024"):
        x = [y for y in e if y["arm"] == arm and y["cell"] == f"ring_{L}"]
        assert len(x) == 1; x = x[0]
        tot_ans += int(x["answered"]); tot_bound += int(x["miss_bound"])
        assert int(x["misses"]) == int(x["miss_bound"]) + int(x["miss_calibration"])
        out.append(f"{L} & {lab} & {x['answered']} & {pct(x['cov_cond'])} ({pct(x['cov_lo'])}--{pct(x['cov_hi'])}) & {x['misses']} & {x['miss_bound']} & "
                   f"{x['miss_calibration']} & {float(x['max_b_over_B']):.2f} & {pct(x['noise_cover'])} & {float(x['median_sd_boot_over_sd_emp']):.2f} \\\\")
write("tab_r3_e4o.tex", out)
notes.append(f"E4o diagnosis: answered datasets in all eight rows {tot_ans}; bound misses {tot_bound}   (from e4diag_summary.csv)")
notes.append("E4o diagnosis: min over rows of min_sd_ratio = %.2f" % min(float(x["min_sd_ratio"]) for x in e))
with open(os.path.join(HERE, "round3_tables_numbers.txt"), "w") as f:
    f.write("\n".join(notes) + "\n")
print("\n".join(notes))
