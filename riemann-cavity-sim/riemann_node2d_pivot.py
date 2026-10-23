"""Iteration-5 pivot: bracket the MID (left-edge) and GP (right-edge)
node minima. Pre-stated in ITERATION_LOG 2026-10-22 18:26 EDT."""
import csv, time
import numpy as np
import riemann_node2d_sim as N5
import riemann_tp_cylinder_sim as S

EXTRA = {}
for d in (0.02, 0.05, 0.08, 0.10):
    EXTRA[f"MID_d{d:.2f}"] = N5._acell([N5.TONE_MID], round(0.79 - d, 4), 0.0)
for d in (0.30, 0.34, 0.40):
    EXTRA[f"GP_d{d:.2f}"] = N5._acell([N5.TONE_GP], round(0.79 - d, 4), 0.0)

rows = []
t0 = time.time()
for name, cell in EXTRA.items():
    for seed in N5.SEEDS:
        res = N5.run_node(S.ARMS["LADDER"], seed, cell)
        ex = res["exch"]
        rows.append({"cell": name, "seed": seed, "E_inj": ex["E_inj"],
                     "E_ind": ex["E_ind"],
                     "ind_frac": ex["E_ind"] / ex["E_inj"],
                     "closure_frac": ex["closure_frac"],
                     "sat_frac": ex["sat_frac"],
                     "clamp_count": ex["clamp_count"]})
        print(name, seed, f"ind_frac={rows[-1]['ind_frac']:.5f}",
              f"elapsed={time.time()-t0:.0f}s", flush=True)
with open(f"{N5.HD}/riemann_node2d_pivot_runs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
for pfx in ("MID", "GP"):
    g = {}
    for r in rows:
        if r["cell"].startswith(pfx):
            g.setdefault(r["cell"], []).append(r["ind_frac"])
    print(pfx, {c: round(float(np.median(v)), 5) for c, v in sorted(g.items())})
