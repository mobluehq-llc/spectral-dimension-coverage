"""POST-HOC diagnostic (added after E1q and E2 disagreed in sign on the causet graphs; not pre-registered): numerical reference d_ref(t) = -2 dln qbar/dln t
by central differences of the sprinkling-ensemble-mean q on the stored time union (noise-free given the ensemble; NOT an exact truth: it has secant
bias and is an ensemble mean), restricted to the p22 window N*qbar >= 1.3. Reports peak value / location, and the value at the last window point."""
import pickle, glob, os, numpy as np, pandas as pd
import causet_lib as cl
rows = [pickle.load(open(f, "rb")) for f in glob.glob("results/causet/r_*.pkl") if os.path.getsize(f) > 0]
by = {}
for r in rows: by.setdefault((r["dim"], r["N"]), []).append(r)
g = cl.GRID.astype(float); out = []
for (dim, N), rs in sorted(by.items()):
    for con in ("k0", "k3", "k15", "rgg0", "rgg3", "rgg15"):
        q = np.mean([r[f"q_{con}"] for r in rs], 0)
        ok = np.nonzero(N * q >= 1.3)[0]; ok = ok[(ok >= 1) & (ok < len(g) - 1)]
        if len(ok) < 3: continue
        d = np.array([-2 * (np.log(q[j + 1]) - np.log(q[j - 1])) / (np.log(g[j + 1]) - np.log(g[j - 1])) for j in ok])
        k = int(np.argmax(d))
        out.append(dict(dim=dim, N=N, con=con, R=len(rs), d_peak=round(d[k], 2), t_peak=g[ok[k]], d_last=round(d[-1], 2), t_last=g[ok[-1]], n_window_pts=len(ok)))
o = pd.DataFrame(out); o.to_csv("results/causet_reference_peak.csv", index=False)
pv = o.pivot_table(index=["dim", "con"], columns="N", values="d_peak"); print(pv.round(2).to_string())
for (dim, con), x in o.groupby(["dim", "con"]):
    x = x.sort_values("N"); print(dim, con, "peak-d exponent", round(np.polyfit(np.log(x.N), np.log(x.d_peak), 1)[0], 3), "t_peak", x.t_peak.tolist())
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), sharey=False)
for ax, dim in zip(axs, (1, 3)):
    for con, ls in (("k0", "-"), ("rgg0", "--")):
        for N in (500, 1000, 2000, 4000, 8000):
            rs = by[(dim, N)]; q = np.mean([r[f"q_{con}"] for r in rs], 0)
            ok = np.nonzero(N * q >= 1.3)[0]; ok = ok[(ok >= 1) & (ok < len(g) - 1)]
            d = [-2 * (np.log(q[j + 1]) - np.log(q[j - 1])) / (np.log(g[j + 1]) - np.log(g[j - 1])) for j in ok]
            ax.plot(g[ok], d, ls, marker="o", ms=2.5, color=plt.cm.viridis(np.log2(N / 500) / 4), label=f"{con} N={N}" if con == "k0" else None)
    ax.set_xscale("log"); ax.set_xlabel("lazy-walk steps t"); ax.set_ylabel("d_ref(t) (central difference of ensemble-mean q)"); ax.grid(alpha=.3)
    ax.set_title(f"{'1+1' if dim == 1 else '3+1'}: Hasse/link graph (solid) vs matched-degree RGG (dashed)", fontsize=9)
axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig("results/fig7_causet_reference_dt.png", dpi=130)
