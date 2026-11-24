# ROUND 6: E7 control datasets with no band (infeasible) per budget in the holdout (Amendment 11). Reads archive/data/raw_amend11.csv.
import os, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(H, "archive/data/raw_amend11.csv"), low_memory=False)
m = pd.read_csv(os.path.join(H, "results_amend10_11/cells_meta_amend11.csv"))
meta = pd.read_csv(os.path.join(H, "results_amend10_11/tables_11_cells.csv"))[["cell", "kind"]].drop_duplicates()
e = d[d.est == "E7"].merge(meta, on="cell"); c = e[e.kind == "control"]
for w, g in c.groupby("wpn"): print(f"E7 holdout controls wpn={w}: n={len(g)} no_band={int(g.infeasible.sum())}")
