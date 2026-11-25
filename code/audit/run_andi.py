"""Part A runner: AnDi-2 tracks through E0/E1q/E4c/E6.  usage: run_andi.py tracks.npz out.csv [NCPU] [LIMIT]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, multiprocessing as mp, time
import audit_lib as AL, andi_lib as AN

B = 200; TS = (100, 150, 200); CAP = 400
SEG_MIN = {100: 10, 150: 15, 200: 20}


def select(Z):
    sel = []
    for i in range(len(Z["n"])):
        e = int(Z["exp"][i])
        for T in TS:
            if Z["n"][i] < T: continue
            a = Z["alpha"][i, :T]; D = Z["D"][i, :T]
            ch = np.nonzero((np.abs(np.diff(a)) > 1e-9) | (np.abs(np.diff(D)) > 1e-9))[0] + 1
            if len(ch) == 0 and e in (1, 7, 8):
                sel.append((i, T, "SS", np.nan))
            elif e == 7 and 1 <= len(ch) <= 5:
                seg = np.diff(np.concatenate([[0], ch, [T]]))
                if seg.min() < SEG_MIN[T]: continue
                lags = AN.lag_grid(T)
                if len(lags) < 6: continue
                dl, _ = AN.label_expected_exponent(a, D, lags)
                S = float(np.abs(dl - dl[0]).max())
                if S >= 0.30: sel.append((i, T, "SW", S))
    out = []
    rng = np.random.default_rng(2201)
    for T in TS:
        for cls in ("SS", "SW"):
            s = [x for x in sel if x[1] == T and x[2] == cls]
            if len(s) > CAP: s = [s[j] for j in sorted(rng.choice(len(s), CAP, replace=False))]
            out += s
    return out


def work(args):
    i, T, cls, Slab, Zi = args
    xy, al, D = Zi
    uid = int(Zi[3]) if False else None
    rng = np.random.default_rng(np.random.SeedSequence([2201, int(i), T]))
    lags, inp = AN.track_inputs(xy[:T], T, B, rng)
    row = dict(track=i, T=T, cls=cls, J=len(lags), S_lab=Slab)
    if inp is None:
        row["abstain_all"] = 1; return [row]
    x = np.log(lags.astype(float))
    dl, nseg = AN.label_expected_exponent(al[:T], D[:T], lags)
    truth = dl
    row.update(nseg=nseg, alpha_true=float(al[0]) if cls == "SS" else np.nan, sign_lab=float(np.sign(dl[np.argmax(np.abs(dl - dl[0]))] - dl[0])))
    A_, B_, All_ = inp
    res = AL.run_all(x, A_, B_, All_, truth_d=truth, delta=0.25)
    rows = []
    for tag, (r, cur) in res.items():
        rr = dict(row); rr.update(AL.summarize(r)); rr["est"] = tag
        if cur is not None:
            rr["d"] = ";".join("%.4f" % v for v in cur["d"]); rr["lo"] = ";".join("%.4f" % v for v in cur["lo"]); rr["hi"] = ";".join("%.4f" % v for v in cur["hi"])
        rr["truth"] = ";".join("%.4f" % v for v in truth); rr["lags"] = ";".join(map(str, lags))
        rows.append(rr)
    return rows


if __name__ == "__main__":
    Z = dict(np.load(sys.argv[1])); out = sys.argv[2]
    ncpu = min(int(sys.argv[3]) if len(sys.argv) > 3 else 8, 16); lim = int(sys.argv[4]) if len(sys.argv) > 4 else None
    sel = select(Z)
    if lim: sel = sel[:lim]
    print("tracks selected:", pd.Series([(t, c) for _, t, c, _ in sel]).value_counts().sort_index().to_dict(), flush=True)
    tasks = [(i, T, c, S, (Z["xy"][i], Z["alpha"][i], Z["D"][i])) for i, T, c, S in sel]
    t0 = time.time(); allrows = []
    with mp.Pool(ncpu) as p:
        for k, rr in enumerate(p.imap_unordered(work, tasks, chunksize=2)):
            allrows += rr
            if k % 100 == 0: print(k, len(tasks), round(time.time() - t0), flush=True)
    pd.DataFrame(allrows).to_csv(out, index=False); print("done", len(allrows), round(time.time() - t0), "s")
