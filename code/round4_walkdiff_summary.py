"""ROUND 4 EXPLORATORY: summarise e7_walk_diff_by_dataset.csv by family and by whether the first differing time index lies inside the exact window (exact_J points) -- approximate reach of window-based estimators."""
import sys, pandas as pd
d = pd.read_csv(sys.argv[1]); out = []
for g in ("b8", "b32"):
    m = pd.read_csv(sys.argv[2] % g).rename(columns={"id": "cell"})[["cell", "fam", "kind", "exact_J"]]
    x = d[d.group == g].merge(m, on="cell")
    x["in_window"] = (x.differs == 1) & (x.first_diff_t < x.exact_J)
    t = x.groupby("fam").agg(cells=("cell", "nunique"), datasets=("ds", "size"), walks_differ=("differs", "sum"), first_diff_inside_exact_window=("in_window", "sum"), min_first_diff_from_end=("first_diff_t", lambda s: 0)).drop(columns="min_first_diff_from_end")
    t["group"] = g; out.append(t.reset_index())
    print(g); print(t.to_string())
    print("cells with any difference:", sorted(x[x.differs == 1].cell.unique().tolist())); print("cells with no difference:", sorted(set(x.cell) - set(x[x.differs == 1].cell)))
    print("TOTAL", len(x), "differ", int(x.differs.sum()), "first diff inside exact window", int(x.in_window.sum()))
pd.concat(out).to_csv(sys.argv[3], index=False)
