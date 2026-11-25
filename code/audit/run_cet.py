"""Part B runner: HadCET DFA-2 audit and synthetic pipeline calibration.
usage: run_cet.py real HadCET.txt out_prefix | run_cet.py synth out_prefix [NCPU] [REPS] | run_cet.py pop out_prefix [NCPU]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, multiprocessing as mp, time, json
import audit_lib as AL, dfa_lib as D

B = 300; DELTAS = (0.10, 0.25)
TRAIN = (1772, 1850); P1 = (1851, 1950); P2 = (1951, 2025)
FINE = np.array(sorted({int(round(7 * 2 ** (k / 8))) for k in range(0, 57)}))
FINE = FINE[FINE <= D.SCALES[-1]]


def load_cet(path):
    rows = []
    for ln in open(path):
        p = ln.split()
        if len(p) == 2 and p[0][:2] in ("17", "18", "19", "20") and p[0][4] == "-":
            rows.append((p[0], float(p[1])))
    df = pd.DataFrame(rows, columns=["date", "v"]); df["date"] = pd.to_datetime(df.date)
    return df


def harmonics(doy, K):
    w = 2 * np.pi * doy / 365.25
    cols = [np.ones_like(w)]
    for k in range(1, K + 1): cols += [np.cos(k * w), np.sin(k * w)]
    return np.column_stack(cols)


def deseason(df):
    tr = df[(df.date.dt.year >= TRAIN[0]) & (df.date.dt.year <= TRAIN[1])]
    doy = df.date.dt.dayofyear.values.astype(float); dtr = tr.date.dt.dayofyear.values.astype(float)
    c1, *_ = np.linalg.lstsq(harmonics(dtr, 4), tr.v.values, rcond=None)
    res = tr.v.values - harmonics(dtr, 4) @ c1
    c2, *_ = np.linalg.lstsq(harmonics(dtr, 2), res ** 2, rcond=None)
    mu = harmonics(doy, 4) @ c1; sd = np.sqrt(np.maximum(harmonics(doy, 2) @ c2, 1e-6))
    return (df.v.values - mu) / sd, dict(train_resid_sd=float(res.std()), sd_min=float(sd.min()), sd_max=float(sd.max()))


def inputs_from_series(z, years, rng):
    tab = D.cluster_tables(z, years)
    tA, tB, tAll = D.split_halves(tab)
    return [D.curve_and_boot(t, B, rng) for t in (tA, tB, tAll)], tab


def run_one(z, years, rng, truth_d=None, deltas=DELTAS, tags=("E0", "E1q", "E4c", "E6")):
    inp, tab = inputs_from_series(z, years, rng)
    x = np.log(D.SCALES.astype(float)); out = {}
    for dl in deltas:
        out[dl] = AL.run_all(x, inp[0], inp[1], inp[2], truth_d=truth_d, delta=dl, tags=tags)
    return x, inp, tab, out


def flat(r):
    return {k: v for k, v in r.items() if isinstance(v, (int, float, np.integer, np.floating, str)) or v is None}


def cmd_real(path, prefix):
    df = load_cet(path); assert len(df) > 90000 and df.date.is_monotonic_increasing and (df.date.diff().dropna() == pd.Timedelta(days=1)).all(), "missing days"
    assert (df.v > -50).all()
    z, info = deseason(df); yrs = df.date.dt.year.values
    meta = dict(n_days=len(df), first=str(df.date.iloc[0].date()), last=str(df.date.iloc[-1].date()), sha256=os.popen(f"shasum -a 256 {path}").read().split()[0], **info)
    json.dump(meta, open(prefix + "_meta.json", "w"), indent=1)
    rows = []; curves = {}
    for name, (a, b) in (("P1", P1), ("P2", P2)):
        m = (yrs >= a) & (yrs <= b)
        rng = np.random.default_rng(np.random.SeedSequence([2203, a]))
        x, inp, tab, out = run_one(z[m], yrs[m], rng)
        y = inp[2][0]; Fs = np.exp(-2 * y)
        curves[name] = dict(scales=D.SCALES.tolist(), F=Fs.tolist(), logF_sd=(2 * inp[2][1].std(0)).tolist(), boxes=tab["cnt"].sum(0).tolist())
        for dl, res in out.items():
            for tag, (r, cur) in res.items():
                rr = flat(r); rr.update(period=name, delta=dl, est=tag, years=f"{a}-{b}")
                if cur is not None: rr.update(d=";".join("%.5f" % v for v in cur["d"]), lo=";".join("%.5f" % v for v in cur["lo"]), hi=";".join("%.5f" % v for v in cur["hi"]), xs=";".join("%.5f" % v for v in cur["x"]))
                rows.append(rr)
    pd.DataFrame(rows).to_csv(prefix + "_real.csv", index=False); json.dump(curves, open(prefix + "_curves.json", "w"))
    print("real done", len(rows))


def gen(kind, n, rng):
    if kind == "N0": return D.fgn(n, 0.65, rng)
    if kind == "N1":
        e = rng.standard_normal(n); a = np.empty(n); a[0] = e[0]
        phi = 0.85; s = np.sqrt(1 - phi ** 2)
        for i in range(1, n): a[i] = phi * a[i - 1] + s * e[i]
        return a + 0.5 * D.fgn(n, 0.65, rng)


def pop_one(args):
    kind, r = args
    rng = np.random.default_rng(np.random.SeedSequence([2204, r, 0 if kind == "N0" else 1]))
    n = 36500; z = gen(kind, n, rng); yrs = 1851 + np.arange(n) // 365
    tab = D.cluster_tables(z, yrs, FINE)
    return tab["sumF"].sum(0), tab["cnt"].sum(0)


def cmd_pop(prefix, ncpu):
    out = {}
    with mp.Pool(ncpu) as p:
        for kind in ("N0", "N1"):
            res = p.map(pop_one, [(kind, r) for r in range(200)])
            f2 = sum(r[0] for r in res) / sum(r[1] for r in res)
            lx = np.log(FINE.astype(float)); ly = 0.5 * np.log(f2)           # log F
            d = np.gradient(ly, lx)                                          # central differences on the fine grid
            out[kind] = dict(fine=FINE.tolist(), logF=ly.tolist(), d_coarse=np.interp(np.log(D.SCALES.astype(float)), lx, d).tolist())
    json.dump(out, open(prefix + "_pop.json", "w"))
    print({k: np.round(v["d_coarse"], 3).tolist() for k, v in out.items()})


def synth_one(args):
    kind, r, truth = args
    rng = np.random.default_rng(np.random.SeedSequence([2202, r, 0 if kind == "N0" else 1]))
    n = 36500; z = gen(kind, n, rng); yrs = 1851 + np.arange(n) // 365
    x, inp, tab, out = run_one(z, yrs, rng, truth_d=np.array(truth))
    rows = []
    for dl, res in out.items():
        for tag, (rr, cur) in res.items():
            q = flat(rr); q.update(kind=kind, rep=r, delta=dl, est=tag); rows.append(q)
    return rows


def cmd_synth(prefix, ncpu, reps):
    pop = json.load(open(prefix + "_pop.json")); tasks = []
    for kind in ("N0", "N1"):
        for r in range(reps): tasks.append((kind, r, pop[kind]["d_coarse"]))
    t0 = time.time(); rows = []
    with mp.Pool(ncpu) as p:
        for k, rr in enumerate(p.imap_unordered(synth_one, tasks)):
            rows += rr
            if k % 20 == 0: print(k, len(tasks), round(time.time() - t0), flush=True)
    pd.DataFrame(rows).to_csv(prefix + "_synth.csv", index=False); print("synth done")


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "real": cmd_real(sys.argv[2], sys.argv[3])
    elif c == "pop": cmd_pop(sys.argv[2], min(int(sys.argv[3]), 16))
    elif c == "synth": cmd_synth(sys.argv[2], min(int(sys.argv[3]), 16), int(sys.argv[4]) if len(sys.argv) > 4 else 100)
