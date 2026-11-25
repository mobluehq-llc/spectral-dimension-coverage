"""Build paper.md from main.tex via pandoc (no TeX engine needed). Usage: python3 make_md.py  (run in paper/)."""
import re, subprocess
s = open("main.tex").read()
s = s.replace("\\documentclass[aps,pre,twocolumn,superscriptaddress,nofootinbib,10pt]{revtex4-2}", "\\documentclass{article}\n\\usepackage{natbib}")
s = s.replace("\\author{[Author names to be added]}\n\\affiliation{[Affiliation to be added]}", "\\author{[Author names and affiliation to be added]}")
s = s.replace("\\maketitle\n", "").replace(" \\nocite{andi2021}", "")
s = s.replace("table*", "table").replace("figure*", "figure")
s = re.sub(r"(figs/fig[0-9a-z_]+)\}", r"\1.png}", s)
refs = {"tab:main": "1", "fig:ex": "1", "fig:front": "2", "fig:sep": "3", "fig:cs": "4", "sec:bench": "2", "sec:res": "3", "sec:audit": "4", "sec:rec": "5", "sec:limits": "Limitations"}
s = re.sub(r"\\ref\{([^}]*)\}", lambda m: refs[m.group(1)], s)
s = s.replace("Sec.~Limitations", "the Limitations section")
open("/tmp/p22_md_src4.tex", "w").write(s)
subprocess.run(["pandoc", "/tmp/p22_md_src4.tex", "-f", "latex", "-t", "markdown-raw_html-native_divs+implicit_figures", "--citeproc", "--bibliography=refs.bib", "--wrap=none", "-s", "-o", "paper.md"], check=True)
t = open("paper.md").read()
t = t.replace("<local path omitted> "refs.bib").replace("title:", "nocite: '@andi2021'\ntitle:", 1)
pat = re.compile(r':{3,} \{#fig:(\w+) \.figure[^\n]*\}\n!\[\]\((figs/[^)]+)\)\{[^}]*\}\n\n::: caption\n(.*?)\n:::\n:{3,}\n', re.S)
t = pat.sub(lambda m: f"![**Figure ?.** {m.group(3).strip()}]({m.group(2)}){{width=85%}}\n", t)
t = re.sub(r"\{width=\"[^\"]*\"\}", "{width=85%}", t)
t = t.replace("![Coverage--power", "![**Figure ?.** Coverage--power")
k = [0]
def num(m): k[0] += 1; return f"![**Figure {k[0]}.**"
t = re.sub(r"!\[\*\*Figure \?\.\*\*", num, t)
t = t.replace("  : Estimator summary.", "  : **Table 1.** Estimator summary.")
open("paper.md", "w").write(t)
