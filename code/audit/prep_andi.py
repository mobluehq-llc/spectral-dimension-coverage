"""Build a compact npz of AnDi-2 track_1 tracks (exp_1, exp_7, exp_8; length >= 100) from the extracted csv. Reads with python -I."""
import sys, glob, os, numpy as np, pandas as pd
R = sys.argv[1]; OUT = sys.argv[2]
rows = []
for e in (1, 7, 8):
    for f in sorted(glob.glob(f"{R}/track_1/exp_{e}/trajs_fov_*.csv"), key=lambda s: int(s.split("_")[-1][:-4])):
        fov = int(f.split("_")[-1][:-4]); d = pd.read_csv(f)
        for k, t in d.groupby("traj_idx"):
            if len(t) < 100: continue
            t = t.sort_values("frame"); assert (np.diff(t.frame.values) == 1).all()
            n = min(len(t), 200)
            xy = np.full((200, 2), np.nan); al = np.full(200, np.nan); D = np.full(200, np.nan)
            xy[:n] = t[["x", "y"]].values[:n]; al[:n] = t.alpha.values[:n]; D[:n] = t.D.values[:n]
            rows.append((e, fov, int(k), n, xy, al, D))
np.savez_compressed(OUT, exp=np.array([r[0] for r in rows]), fov=np.array([r[1] for r in rows]), idx=np.array([r[2] for r in rows]), n=np.array([r[3] for r in rows]),
                    xy=np.stack([r[4] for r in rows]), alpha=np.stack([r[5] for r in rows]), D=np.stack([r[6] for r in rows]))
print(len(rows), "tracks ->", OUT)
