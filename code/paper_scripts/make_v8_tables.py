"""v8 tables for the supplement (presentation only: parsed from stored report files; nothing is simulated or recomputed).
Run from p22_spectral_dimension_flow:  python3 paper/make_v8_tables.py
Reads results_amend13/amend13a_report.txt, results_amend13/amend13b_report.txt, results/rt_rbc_analysis.txt
Writes paper/tab_v8_a13a.tex, paper/tab_v8_a13b.tex, paper/tab_v8_rbc_cells.tex"""
import re
NAME = {"E4o": "BiasAware-Oracle", "E4c": "BiasAware-CV", "E4o_V2": "BiasAware-Split (exact $M$)", "E4c_V2": "BiasAware-Split (CV $M$)",
        "E1q": "LocQuad", "E7": "MonoLP"}
# ---- 13a
rx = re.compile(r"\[(rep|hold) wpn\s+(\d+)\] (\S+)\s+C (\d+)/(\d+) = ([\d.]+)% \(CP ([\d.]+)-([\d.]+)\) -> (PASS|FAIL) \(cells CP-upper<93: \[(.*?)\]\) \| "
                r"F (\d+)/(\d+) = [\d.]+% .*?-> (PASS|FAIL) \| X (\d+)/(\d+) = ([\d.]+)% \(CP ([\d.]+)-([\d.]+)\).*?minL \S+ (\S+), \S+ (\S+)")
out = []; n = 0
for l in open("results_amend13/amend13a_report.txt"):
    m = rx.search(l)
    if not m: continue
    n += 1
    des, w, arm, ck, cn, cp, clo, chi, cv, cells, fk, fn, fv, xk, xn, xp, xlo, xhi, m2, m3 = m.groups()
    ncell = len([c for c in cells.split(",") if c.strip()])
    mm = lambda v: "--" if v == "inf" else v
    out.append(f"{'Replication' if des == 'rep' else 'Holdout'} & {w} & \\textsf{{{NAME[arm]}}} & {ck}/{cn} & {cp} ({clo}--{chi}) & "
               f"{'met' if cv == 'PASS' else 'not met'} ({ncell}) & {fk}/{fn} & {xk}/{xn} & {xp} ({xlo}--{xhi}) & {mm(m2)}, {mm(m3)} \\\\")
assert n == 24, n
open("paper/tab_v8_a13a.tex", "w").write("\n".join(out) + "\n")
# ---- 13b
rx = re.compile(r"^\s+(E1q|E4c|E4c_V2|E7)\s+wpn\s+(\d+) (crossover_like|flat_like)\s+(\d+)/(\d+) = ([\d.]+)% \(CP ([\d.]+)-([\d.]+)\)")
tgt = None; D = {}
for l in open("results_amend13/amend13b_report.txt"):
    if l.startswith("--- target:"): tgt = "S" if "single" in l else "E"; continue
    m = rx.match(l)
    if m: D[(tgt, m[1], m[3], int(m[2]))] = (m[4], m[5], m[6], m[7], m[8])
assert len(D) == 48, len(D)
out = []
for arm in ("E1q", "E4c", "E4c_V2", "E7"):
    for g in ("crossover_like", "flat_like"):
        for w in (8, 32, 128):
            a = D[("S", arm, g, w)]; b = D[("E", arm, g, w)]
            nm = NAME[arm].replace(" (CV $M$)", "")
            out.append(f"\\textsf{{{nm}}} & {g.split('_')[0]} & {w} & {a[0]}/{a[1]} & {a[2]} ({a[3]}--{a[4]}) & {b[0]}/{b[1]} & {b[2]} ({b[3]}--{b[4]}) \\\\")
open("paper/tab_v8_a13b.tex", "w").write("\n".join(out) + "\n")
# ---- RBC by cell
rows = {}; fam = L = None
for l in open("results/rt_rbc_analysis.txt"):
    m = re.match(r"^(aniso2_w16|aniso3_w8)?\s*(\d+)?\s+(8|32|128)\s+(\d\.\d\d)\s+(\d\.\d\d)\s+(\d\.\d\d)\s*$", l)
    if not m: continue
    if m[1]: fam = m[1]
    if m[2]: L = int(m[2])
    rows[(fam, L, int(m[3]))] = (m[4], m[5], m[6])
assert len(rows) == 18, len(rows)
out = []
for (f, L) in sorted({k[:2] for k in rows}):
    c = lambda i: ", ".join(str(int(round(100 * float(rows[(f, L, w)][i])))) for w in (8, 32, 128))
    out.append(f"{'$C_L\\times C_{16}$' if f == 'aniso2_w16' else '$C_L\\times C_8^2$'} & {L} & {c(0)} & {c(1)} & {c(2)} \\\\")
open("paper/tab_v8_rbc_cells.tex", "w").write("\n".join(out) + "\n")
print("wrote 3 tables")
# ---- 13b cells
import csv
out = []
for r in csv.DictReader(open("results_amend13/amend13b_cells.csv")):
    ws = r["fam"] == "ws"; a, b_ = r["name"].split("_")[1:]
    out.append((0 if ws else 1, f"{'Small-world' if ws else 'Percolation'} & ${'N' if ws else 'L'}={a[1:]}$, $p={b_[1:]}$ & {float(r['N_mean']):.0f} ({r['N_min']}--{r['N_max']}) & {r['exact_J']} & "
               f"{float(r['max_abs_delta']):.2f} & {float(r['sd_d_between_graphs_mean_over_window']):.3f} & {'crossover-like' if r['crossover_like'] == 'True' else 'flat-like'} \\\\"))
open("paper/tab_v8_a13b_cells.tex", "w").write("\n".join(x[1] for x in sorted(out)) + "\n")
