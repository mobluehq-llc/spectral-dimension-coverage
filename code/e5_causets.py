"""Amendment 4: E5qu monotone-band diagnostic on the causal-set datasets (no truth). Both directions reported (two looks). Writes results/causet_E5.csv"""
import glob, pickle, os
import numpy as np, pandas as pd
import p22lib as L, p22amend as A, p22shape as S, causet_lib as cl
CONS = ["k0", "k3", "k15", "rgg0", "rgg3", "rgg15"]
rows = [pickle.load(open(f, "rb")) for f in glob.glob("results/causet/r_*.pkl") if os.path.getsize(f) > 0]
by = {}
for r in rows: by.setdefault((r["dim"], r["N"]), []).append(r)
idx = {int(t): i for i, t in enumerate(cl.GRID)}
out = []
for (dim, N), rs in sorted(by.items()):
    rs = sorted(rs, key=lambda x: x["rep"])
    for con in CONS:
        for thr in (6.0, 1.3):
            grid = cl.GRID4; ci = [idx[int(t)] for t in grid]
            U = np.array([r[f"q_{con}"][ci] for r in rs]); n = len(rs); mw = 32.0 * N / n
            rec = dict(dim=dim, N=N, con=con, thr=thr, R=n)
            p = U[0::2].mean(0); W = U[0::2].shape[0] * mw; cnt = p * W
            J = L.window_len(cnt, p, N, thr)
            if J < 6: rec["abstain"] = 1; out.append(rec); continue
            x = np.log(grid[:J].astype(float)); yA = np.log(p[:J])
            best = S.slope_ops(yA, cnt[:J], x)
            if best is None: rec["abstain"] = 1; out.append(rec); continue
            rng = np.random.default_rng(np.random.SeedSequence([22, 4, dim, N, CONS.index(con), int(thr * 10)]))
            UB = U[1::2]; floor = 0.5 / (UB.shape[0] * mw)
            yB = np.log(np.maximum(UB.mean(0), floor)[:J])
            Ys = np.log(np.maximum(L.boot_means([UB], 200, rng)[:, :J], floor))
            _, sl = A.poly_rows(x, cnt[:J], 0.7 * best[1], 2)
            dh = -2 * sl @ yB; Ds = -2 * Ys @ sl.T
            rec.update(abstain=0, J=J, d_first=float(dh[0]), d_last=float(dh[-1]), d_max=float(dh.max()))
            for nm, dec in (("drop", True), ("rise", False)):
                r = S.e5_core(dh, Ds, None, "E5qu", dec)
                rec[f"{nm}_declared"] = r["decl_neg"] if dec else r["decl_pos"]
                rec[f"{nm}_band_gap"] = r["band_drop_lo"]
            out.append(rec)
df = pd.DataFrame(out); df.to_csv("results/causet_E5.csv", index=False)
d = df[df.abstain == 0]
print(d.groupby(["dim", "con", "thr"])[["drop_declared", "rise_declared"]].mean())
print(len(d), "non-abstained of", len(df))
