"""v7 presentation layer: map the archive's estimator codes (E0, E1q, E4c, E7, GKW, ...) to the self-documenting names used in
the paper. The archive, the plans and the result files keep the codes; SM Table S1 is the crosswalk.
Usage: python3 apply_names.py            -> writes v7tabs/<name>.tex for every generated tab_*.tex (inputs are not modified)
       python3 apply_names.py FILE.tex   -> prints FILE with codes replaced (used once to convert supplement.tex for v7)"""
import re, sys, glob, os
MAP = [("GK2Wb", "WLS-2seg-Boot"), ("GK2W", "WLS-2seg"), ("GK2L", "OLS-2seg"), ("GKWm", "WLS-loglog-Model"), ("GKWb", "WLS-loglog-Boot"),
       ("GKW", "WLS-loglog"), ("GKL", "OLS-loglog"), ("E0g", "OneLine"), ("E0", "TwoLine"), ("E1q", "LocQuad"), ("E1x", "LocLin-NoSplit"),
       ("E1u", "LocLin-Under"), ("E1", "LocLin"), ("E2", "TraceSecant"), ("E3s", "Collapse1"), ("E3", "Collapse3"), ("E4o", "BiasAware-Oracle"),
       ("E4c", "BiasAware-CV"), ("E4r", "BiasAware-Rule"), ("E4", "BiasAware"), ("E6", "StepParam"), ("E7j", "MonoLP-Joint"), ("E7", "MonoLP"),
       ("E8", "ConvexSecant"), ("CL", "Collapse-Opt"), ("Mb", "WML-Boot")]
RX = re.compile(r"(?<![A-Za-z0-9_\\{:])(" + "|".join(k for k, _ in MAP) + r")(?![A-Za-z0-9_])")
D = dict(MAP)
def conv(s, wrap=True):
    s = RX.sub(lambda m: ("\\textsf{%s}" % D[m.group(1)]) if wrap else D[m.group(1)], s)
    return s
def fix_brk(line):
    if line.count("&") != 6: return line
    c = [x.strip() for x in line.rstrip().rstrip("\\").split("&")]
    if "$^\\ast$" in c[4]:
        assert c[6] == "", line
        c[4:7] = ["", c[4], c[5]]
    c[4:7] = [x if x else "--" for x in c[4:7]]
    return " & ".join(c) + " \\\\"
if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.stdout.write(conv(open(sys.argv[1]).read()))
    else:
        os.makedirs("v7tabs", exist_ok=True)
        for f in sorted(glob.glob("tab_*.tex")):
            if f == "tab_v7_joint.tex": continue
            s = conv(open(f).read())
            s = re.sub(r"(?m)^M &", r"\\textsf{WML} &", s)      # the windowed-ML arm, archive code M (not the bound $M$)
            if f.startswith("tab_r6_brk"):                      # v8: when L_min is the smallest tested size its entry goes in the
                s = "\n".join(fix_brk(l) for l in s.split("\n"))  # "At L_min" column (it was printed under "Below"); empty cells -> "--"
            open(os.path.join("v7tabs", f), "w").write(s)
        print("wrote", len(glob.glob("v7tabs/*.tex")), "tables")
