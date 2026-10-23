"""
riemann_dual_cavity_pass2.py — PASS 2 (pre-stated ITERATION_LOG.md
2026-10-22 13:47 EDT, the ONE authorized follow-up pass):
(a) tight re-audit of the main-grid candidates at tol 0.30 / 0.15;
(b) analysis of the NONE-cell high-resolution runs
    (riemann_dual_cavity_pickups_none_hi.csv) audited at tol 0.15.
Same G4 clustering / eddy rules as the main analysis. Merges results
into riemann_dual_cavity_results.json under "pass2".
"""

import importlib.util
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

HD = "/home/hatch/workspace/harness-dev"
spec = importlib.util.spec_from_file_location("rdc", f"{HD}/riemann_dual_cavity_sim.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def reaudit(df, tol):
    m.PEAK_TOL = tol
    out = []
    for _, r in df.iterrows():
        if r["is_drive"]:
            out.append("drive")
        else:
            out.append(m.audit_label(float(r["peak_freq"]), m.ARMS[r["arm"]],
                                     r["modulation"], float(r["mod_rate"])))
    return out


def cluster_ue(df, label_col):
    """G4 clustering of UNEXPLAINED-labelled peaks; returns cluster list."""
    cand = df[df[label_col] == "UNEXPLAINED"]
    clusters = []
    for (arm, mod, rate), cell in cand.groupby(["arm", "modulation", "mod_rate"]):
        cell = cell.sort_values("peak_freq")
        freqs = cell["peak_freq"].to_numpy()
        used = np.zeros(len(cell), dtype=bool)
        for i in range(len(cell)):
            if used[i]:
                continue
            members = [i]
            used[i] = True
            for j in range(i + 1, len(cell)):
                if not used[j] and abs(freqs[j] - freqs[i]) <= 0.5:
                    members.append(j)
                    used[j] = True
            sub = cell.iloc[members]
            seeds = sorted(sub["seed"].unique())
            if len(seeds) < 3:
                continue
            lit = []
            for (cav, pik), g in sub.groupby(["cavity", "pickup_idx"]):
                if g["seed"].nunique() / len(seeds) >= 0.5:
                    lit.append(f"{cav}{int(pik)}")
            clusters.append({
                "arm": arm, "modulation": mod, "mod_rate": float(rate),
                "freq_median": round(float(sub["peak_freq"].median()), 3),
                "amp_median": round(float(sub["peak_amp"].median()), 5),
                "n_peaks": int(len(sub)), "seeds_present": [int(s) for s in seeds],
                "n_lit_pickups": len(lit), "lit_pickups": lit,
                "eddy": bool(2 <= len(lit) <= 7),
            })
    return clusters


def arm_breakdown(df, label_col):
    cand = df[~df["is_drive"]]
    ue = cand[cand[label_col] == "UNEXPLAINED"]
    return {
        "candidates": int(len(cand)),
        "unexplained_peaks": int(len(ue)),
        "unexplained_frac": round(float(len(ue) / max(len(cand), 1)), 4),
        "unexplained_by_arm": {a: int((ue["arm"] == a).sum())
                               for a in ("LADDER", "DETUNE", "SHIFT")},
        "label_counts": dict(Counter(cand[label_col])),
    }


main_df = pd.read_csv(f"{HD}/riemann_dual_cavity_pickups.csv")
pass2 = {}
for tol in (0.30, 0.15):
    col = f"audit_{tol}"
    main_df[col] = reaudit(main_df, tol)
    cl = cluster_ue(main_df, col)
    pass2[f"main_reaudit_tol{tol}"] = {
        **arm_breakdown(main_df, col),
        "ue_clusters": len(cl),
        "ue_eddy_clusters": int(sum(1 for c in cl if c["eddy"])),
        "ue_eddy_by_arm": {a: int(sum(1 for c in cl if c["eddy"] and c["arm"] == a))
                           for a in ("LADDER", "DETUNE", "SHIFT")},
        "ue_clusters_list": cl,
    }

hi_df = pd.read_csv(f"{HD}/riemann_dual_cavity_pickups_none_hi.csv")
hi_df["audit_0.15"] = reaudit(hi_df, 0.15)
cl_hi = cluster_ue(hi_df, "audit_0.15")
pass2["none_hi_tol0.15"] = {
    **arm_breakdown(hi_df, "audit_0.15"),
    "ue_clusters": len(cl_hi),
    "ue_eddy_clusters": int(sum(1 for c in cl_hi if c["eddy"])),
    "ue_eddy_by_arm": {a: int(sum(1 for c in cl_hi if c["eddy"] and c["arm"] == a))
                       for a in ("LADDER", "DETUNE", "SHIFT")},
    "ue_clusters_list": cl_hi,
}

res_path = f"{HD}/riemann_dual_cavity_results.json"
res = json.load(open(res_path))
res["pass2"] = pass2
json.dump(res, open(res_path, "w"), indent=1)
for k, v in pass2.items():
    print(k, {kk: vv for kk, vv in v.items() if kk != "ue_clusters_list"})
