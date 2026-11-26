#!/usr/bin/env python3
"""ChatGPT R4 check (Paper A v16): exact drop of the generating graph over the window the global test ACTUALLY used,
for every one of the 900 level datasets of the exact-CP-box run (BD bar 2): 300 controls of arm (a) at 32/node and
600 datasets of arm (b) (the three graphs with reference-window exact drop < 0.245, at 32 and 128/node; the window
does not depend on the box). Regenerates each dataset from the frozen seeds with rt_cpbox.gen (no LP is solved) and
records Jw, the exact drop max_{j<Jw}[d(t_0)-d(t_j)], and the exact drop over every window the rule could select
(prefixes of length 2..K; the largest is the full-grid prefix). Rows: results/null_windows.jsonl."""
import os, sys, json
sys.argv = [sys.argv[0], "a"] + sys.argv[1:]
import numpy as np
import rt_cpbox as B

def drop_over(d, J): return float((d[0] - d[:J]).max()) if J >= 2 else float("nan")

def main():
    out = []
    a = B.core_cells(B.A_CTRL)
    jobs = [(c, 32, "a") for c in a] + [(c, w, "b") for c in B.load_b() if c["max_drop"] < 0.245 for w in B.B_WPN]
    for c, w, arm in jobs:
        ts = c["ts_w"]; d = c["graph"].exact(ts)[1]; K = len(ts)
        allw = [drop_over(d, J) for J in range(2, K + 1)]
        for ds in range(100):
            _, U, m, Jw = B.gen(c, w, ds)
            out.append(dict(arm=arm, cell=int(c["id"]), fam=c["fam"], L=int(c["L"]), N=int(c["N"]), wpn=w, ds=ds, K=K, Jw=int(Jw),
                            ref_J=int(c["exact_J"]), ref_drop=float(c["max_drop"]), tested_drop=drop_over(d, Jw), max_any_window=max(allw)))
        r = [x for x in out if x["cell"] == c["id"] and x["wpn"] == w]
        print(c["fam"], c["L"], w, "ref", round(c["max_drop"], 4), "Jw range", min(x["Jw"] for x in r), max(x["Jw"] for x in r),
              "max tested", np.nanmax([x["tested_drop"] for x in r] + [-1]), "max any", r[0]["max_any_window"], flush=True)
    os.makedirs("results", exist_ok=True)
    with open("results/null_windows.jsonl", "w") as f:
        for x in out: f.write(json.dumps(x) + "\n")
    n = len(out); tested = [x for x in out if x["Jw"] >= 2]
    print("datasets", n, "with window", len(tested), "tested drop > 0.25:", sum(x["tested_drop"] > 0.25 for x in tested),
          "> 0.245:", sum(x["tested_drop"] > 0.245 for x in tested), "max", max(x["tested_drop"] for x in tested))

if __name__ == "__main__": main()
