"""
riemann_tp_cylinder_audit2.py — second-delta items 2 + 3 (log 14:56 EDT):
(a) AUDIT EXTENSION: flow-coupling combinations (fi +- k*s_i, fi+-fj +- s,
    split-pair midpoints; s_i = FP * i) applied to peaks that the BASE
    audit leaves UNEXPLAINED at tol 0.15, in the SETTLED window (main
    grid) and the fine-pass ring-down windows. Migration counts +
    flow-labelled cluster excess test (LADDER vs controls).
(b) MODE-SPLITTING: per component, amp^2-weighted spectral width within
    +-2.0 of the component at A interior pickups (NONE, DRIVE_END and
    SETTLED), FP in {0, 0.025, 0.15}; band medians per FP per arm.
Merges results into riemann_tp_cylinder_results.json.
"""

import itertools
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import riemann_tp_cylinder_sim as S

HD = S.HD
FP = 0.025
ARMS = S.ARMS
TOL = 0.15

def base_label(fc, F, mod, rate):
    return S.audit_label(fc, F, mod, rate, tol=TOL)

def flow_label(fc, F):
    s = FP * np.arange(1, len(F) + 1)
    def hit(t):
        return t > 0.5 and abs(fc - t) <= TOL
    for i, fi in enumerate(F):
        for k in (1, 2):
            if hit(fi + k * s[i]) or hit(fi - k * s[i]):
                return "flow_adv"
    for (i, fi), (j, fj) in itertools.combinations(enumerate(F), 2):
        if hit((fi + s[i] + fj + s[j]) / 2.0):
            return "flow_mid"
        for sv in (s[i], s[j]):
            if hit(fi + fj + sv) or hit(fi + fj - sv) \
                    or hit(abs(fi - fj) + sv) or hit(abs(fi - fj) - sv):
                return "flow_adv"
    return None

def migration(df, label):
    out = {}
    cand = df[~df["is_drive"]]
    for (arm, mod, rate), cell in cand.groupby(["arm", "modulation", "mod_rate"]):
        F = ARMS[arm]
        mig = Counter()
        flow_rows = []
        for _, r in cell.iterrows():
            b = base_label(r["peak_freq"], F, mod, rate)
            if b != "UNEXPLAINED":
                mig["base_explained"] += 1
                continue
            fl = flow_label(r["peak_freq"], F)
            if fl:
                mig[fl] += 1
                flow_rows.append((arm, mod, rate, r["run_id"], r["seed"],
                                  r["peak_freq"]))
            else:
                mig["still_unexplained"] += 1
        out[f"{arm}|{mod}|{rate:.3f}"] = dict(mig)
        out.setdefault("_flow_rows", []).extend(flow_rows)
    flow_rows = out.pop("_flow_rows")
    # cluster the flow-labelled peaks (>=3/5 seeds, tol 0.5) per arm
    cl = defaultdict(int)
    byarm = defaultdict(list)
    for row in flow_rows:
        byarm[row[0]].append(row)
    for arm, rows_ in byarm.items():
        rows_ = sorted(rows_, key=lambda z: z[5])
        used = [False] * len(rows_)
        n_cl = 0
        for i in range(len(rows_)):
            if used[i]:
                continue
            members = [rows_[i]]
            used[i] = True
            for j in range(i + 1, len(rows_)):
                if not used[j] and abs(rows_[j][5] - rows_[i][5]) <= 0.5 \
                        and rows_[j][1] == rows_[i][1]:
                    members.append(rows_[j])
                    used[j] = True
            if len({m[4] for m in members}) >= 3:
                n_cl += 1
        cl[arm] = n_cl
    out["_flow_clusters_by_arm"] = dict(cl)
    return out

pk = pd.read_csv(f"{HD}/riemann_tp_cylinder_pickups.csv")
settled = pk[pk["window"] == "SETTLED"]
mig_settled = migration(settled, "settled")
fine = pd.read_csv(f"{HD}/riemann_fine_pickups.csv")
ring = fine[fine["window"].isin(["RING1", "RING2", "SETTLED"])]
mig_ring = migration(ring, "ring")

# ------------------------------------------------------------- splitting --
def widths(path, fp):
    d = pd.read_csv(path)
    d = d[(d["modulation"] == "NONE") & (d["cavity"] == "A")
          & (d["pickup_kind"] == "in")]
    res = {}
    for wname in ("DRIVE_END", "SETTLED"):
        dw = d[d["window"] == wname]
        for arm in dw["arm"].unique():
            F = ARMS[arm]
            per_comp = []
            for fd in F:
                ws = []
                for rid, g in dw[dw["arm"] == arm].groupby("run_id"):
                    near = g[np.abs(g["peak_freq"] - fd) <= 2.0]
                    if len(near) < 1:
                        continue
                    a2 = near["peak_amp"] ** 2
                    cen = float(np.sum(a2 * near["peak_freq"]) / np.sum(a2))
                    ws.append(float(np.sqrt(np.sum(a2 * (near["peak_freq"] - cen) ** 2)
                                            / np.sum(a2))))
                per_comp.append(float(np.median(ws)) if ws else float("nan"))
            bands = {}
            for bi, (lo, hi) in enumerate(((0, 10), (10, 20), (20, 30))):
                bands[f"n{lo+1}_{hi}"] = float(np.nanmedian(per_comp[lo:hi]))
            res[f"{arm}|{wname}"] = {"per_component": [round(x, 3) if x == x else None
                                                      for x in per_comp],
                                     "bands": bands}
    return res

split = {
    "FP0": widths(f"{HD}/riemann_fp0_pickups.csv", 0.0),
    "FP025": widths(f"{HD}/riemann_tp_cylinder_pickups.csv", 0.025),
    "FP015": widths(f"{HD}/riemann_fp015_pickups.csv", 0.15),
}

R = json.load(open(f"{HD}/riemann_tp_cylinder_results.json"))
R["audit_extension"] = {"settled_window": mig_settled,
                        "fine_ringdown": mig_ring, "tol": TOL}
R["mode_splitting_widths"] = split
json.dump(R, open(f"{HD}/riemann_tp_cylinder_results.json", "w"), indent=1)
print("SETTLED migration:", json.dumps(mig_settled, indent=0)[:1500])
print("RING migration:", json.dumps(mig_ring, indent=0)[:1500])
for fp, d in split.items():
    print(fp, {k: {kk: round(vv, 3) for kk, vv in v["bands"].items()}
               for k, v in d.items()})
