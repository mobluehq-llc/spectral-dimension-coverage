"""v9 evidence lines for the global no-drop test (no simulation: stored per-dataset rows only).
Run from p22_spectral_dimension_flow:  archive/.venv/bin/python paper/make_v9_numbers.py  -> results/v9_nodrop_numbers.txt
Inputs: results/nodrop_exact.csv, nodrop_data.csv, nodrop_controls.csv (rt_nodrop.py; PREREG_GLOBAL_NODROP.md, commit 30b2514)."""
import pandas as pd
from scipy.stats import beta
def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1); hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return 100 * lo, 100 * hi
out = []
e = pd.read_csv("results/nodrop_exact.csv")
for (f, L), g in e.groupby(["fam", "L"], sort=False):
    oc = sorted(set(g.g_outcome)); assert len(oc) == 1, (f, L, oc)
    out.append(f"exact {f} L={L}: eps {sorted(g.eps)} | outcome {oc[0]} at all {len(g)} eps | exact max drop {g.max_drop_exact.iloc[0]:.3f} | "
               f"joint drop {int(g.j_drop.min())} | endpoint drop {int(g.e7_neg.min())} | gate disagreements {int(g.g_gate.sum())}")
d = pd.read_csv("results/nodrop_data.csv")
for (f, L), g in d.groupby(["fam", "L"], sort=False):
    n = len(g); c = g.g_outcome.value_counts().to_dict(); k = c.get("DROP", 0); lo, hi = cp(k, n); kj = int(g.j_drop.sum()); lj, hj = cp(kj, n)
    amb = int(((g.g_outcome == "WITNESS") & (g.dw_found == 1)).sum())
    out.append(f"data {f} L={L} wpn32: n {n} | DROP {k} (CP {lo:.1f}-{hi:.1f}) | WITNESS {c.get('WITNESS', 0)} | INCONCLUSIVE {c.get('INCONCLUSIVE', 0)} | "
               f"other {n - k - c.get('WITNESS', 0) - c.get('INCONCLUSIVE', 0)} | joint drop {kj} (CP {lj:.1f}-{hj:.1f}) | joint infeasible {int(g.j_infeasible.sum())} | "
               f"endpoint drop {int(g.e7_neg.sum())} | joint drop and global not DROP {int(((g.j_drop == 1) & (g.g_outcome != 'DROP')).sum())} | "
               f"both witnesses {amb} | gate disagreements {int(g.g_gate.sum())} | LPs mean {g.g_nlp.mean():.0f} max {int(g.g_nlp.max())} | seconds mean {g.g_seconds.mean():.1f}")
    w = g[g.g_outcome == 'WITNESS']
    if len(w): out.append(f"data {f} L={L} wpn32: no-drop witnesses, smallest change min {w.g_mindelta.min():.3f} median {w.g_mindelta.median():.3f} max {w.g_mindelta.max():.3f}")
    out.append(f"data {f} L={L} wpn32: DROP if every INCONCLUSIVE were a drop {k + c.get('INCONCLUSIVE', 0)} | DROP among joint no-drop {int(((g.j_drop == 0) & (g.g_outcome == 'DROP')).sum())} of {int((g.j_drop == 0).sum())}")
c = pd.read_csv("results/nodrop_controls.csv"); n = len(c); k = int((c.g_outcome == "DROP").sum())
for (f, L), g in c.groupby(["fam", "L"], sort=False):
    out.append(f"control {f} L={L} wpn32: n {len(g)} | DROP {int((g.g_outcome == 'DROP').sum())} | WITNESS {int((g.g_outcome == 'WITNESS').sum())} | joint drop {int(g.j_drop.sum())}")
out.append(f"controls pooled: DROP {k}/{n} | WITNESS {int((c.g_outcome == 'WITNESS').sum())}/{n} | joint drop {int(c.j_drop.sum())}/{n} | CP 95% upper limit on the false-drop rate {cp(k, n)[1]:.1f}%")
open("results/v9_nodrop_numbers.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
