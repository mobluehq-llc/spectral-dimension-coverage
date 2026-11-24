"""Run the causal-set extension (PLAN.md AMENDMENT 2). Usage: run_causets.py [--smoke]. Resumable: one pickle per (dim,N,rep)."""
import sys, os, time, pickle
import multiprocessing as mp
import causet_lib as cl

PLAN = {500: 100, 1000: 100, 2000: 100, 4000: 40, 8000: 12}
if "--smoke" in sys.argv:
    PLAN = {500: 8, 1000: 8, 2000: 8}
OUT = "results/causet" + ("_smoke" if "--smoke" in sys.argv else "")
os.makedirs(OUT, exist_ok=True)

def run(a):
    f = f"{OUT}/r_{a[0]}_{a[1]}_{a[2]}.pkl"
    if os.path.exists(f):
        return f, 0.0
    t = time.time(); o = cl.job(a)
    pickle.dump(o, open(f + ".tmp", "wb")); os.replace(f + ".tmp", f)
    return f, time.time() - t

if __name__ == "__main__":
    tasks = sorted([(d, N, r) for N, R in PLAN.items() for d in (1, 3) for r in range(R)], key=lambda x: (-x[1], x[0], x[2]))
    t0 = time.time()
    with mp.Pool(2 if "--smoke" in sys.argv else 9) as p:
        for i, (f, dt) in enumerate(p.imap_unordered(run, tasks, chunksize=1)):
            if i % 20 == 0:
                print(i, len(tasks), f, round(dt, 1), round(time.time() - t0), flush=True)
    print("DONE", round(time.time() - t0), flush=True)
