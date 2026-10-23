"""
riemann_dual_cavity_analyze.py — verdict analysis for the Riemann
dual-cavity sim. Consumes riemann_dual_cavity_pickups.csv (written by
riemann_dual_cavity_sim.py) and applies the pre-registered gates
(ITERATION_LOG.md 2026-10-22 13:26 EDT):

- Candidate cluster: candidate peaks (non-drive) clustered across the 5
  seeds of a cell by frequency (+-0.5); counts only if present in >= 3/5
  seeds (G4).
- EDDY: lit pickups (presence >= 0.5 across the cluster's seeds) number
  2..7 of 12 (G4).
- A cluster is UNEXPLAINED only if every member peak audited UNEXPLAINED.
- G3 VERDICT: per-run counts of unexplained+eddy peaks; ladder arm mean
  must exceed EACH control arm mean by > 2 pooled SD of the difference,
  else NULL.
Writes riemann_dual_cavity_results.json; prints the verdict table.
"""

import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

HD = "/home/hatch/workspace/harness-dev"
TOL_CLUSTER = 0.5

df = pd.read_csv(f"{HD}/riemann_dual_cavity_pickups.csv")
cand = df[~df["is_drive"]].copy()

clusters = []
for (arm, mod, rate), cell in cand.groupby(["arm", "modulation", "mod_rate"]):
    cell = cell.sort_values("peak_freq")
    used = np.zeros(len(cell), dtype=bool)
    freqs = cell["peak_freq"].to_numpy()
    for i in range(len(cell)):
        if used[i]:
            continue
        members = [i]
        used[i] = True
        for j in range(i + 1, len(cell)):
            if not used[j] and abs(freqs[j] - freqs[i]) <= TOL_CLUSTER:
                members.append(j)
                used[j] = True
        sub = cell.iloc[members]
        seeds_present = sorted(sub["seed"].unique())
        if len(seeds_present) < 3:
            continue  # G4: reproducibility floor
        # pickup presence: fraction of cluster seeds detecting at pickup
        lit = []
        for (cav, pik), g in sub.groupby(["cavity", "pickup_idx"]):
            frac = g["seed"].nunique() / len(seeds_present)
            if frac >= 0.5:
                lit.append(f"{cav}{int(pik)}")
        labels = Counter(sub["audit_label"])
        unexplained = set(sub["audit_label"]) == {"UNEXPLAINED"}
        clusters.append({
            "arm": arm, "modulation": mod, "mod_rate": float(rate),
            "freq_median": round(float(sub["peak_freq"].median()), 3),
            "amp_median": round(float(sub["peak_amp"].median()), 5),
            "amp_max": round(float(sub["peak_amp"].max()), 5),
            "n_peaks": int(len(sub)), "seeds_present": [int(s) for s in seeds_present],
            "n_lit_pickups": len(lit), "lit_pickups": lit,
            "eddy": bool(2 <= len(lit) <= 7),
            "unexplained": bool(unexplained),
            "audit_labels": dict(labels),
        })

# per-run unexplained+eddy peak counts (G3 statistic)
ue_keys = {(c["arm"], c["modulation"], c["mod_rate"], round(c["freq_median"], 1))
           for c in clusters if c["unexplained"] and c["eddy"]}
run_counts = defaultdict(int)
for _, r in cand.iterrows():
    key = (r["arm"], r["modulation"], r["mod_rate"], round(r["peak_freq"], 1))
    # match peak to a UE cluster by cell + nearest cluster freq within tol
    for c in clusters:
        if (c["unexplained"] and c["eddy"] and c["arm"] == r["arm"]
                and c["modulation"] == r["modulation"]
                and abs(c["mod_rate"] - r["mod_rate"]) < 1e-9
                and abs(c["freq_median"] - r["peak_freq"]) <= TOL_CLUSTER):
            run_counts[(r["arm"], int(r["run_id"]))] += 1
            break

arm_stats = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    runs = sorted({int(x) for x in df[df["arm"] == arm]["run_id"].unique()})
    counts = np.array([run_counts.get((arm, rid), 0) for rid in runs], dtype=float)
    arm_stats[arm] = {
        "n_runs": len(runs),
        "ue_eddy_peaks_per_run_mean": round(float(counts.mean()), 3),
        "ue_eddy_peaks_per_run_sd": round(float(counts.std(ddof=1)), 3),
        "ue_eddy_peaks_total": int(counts.sum()),
    }

def g3(arm_a, arm_b):
    a, b = arm_stats[arm_a], arm_stats[arm_b]
    diff = a["ue_eddy_peaks_per_run_mean"] - b["ue_eddy_peaks_per_run_mean"]
    pooled = np.sqrt(a["ue_eddy_peaks_per_run_sd"] ** 2 / a["n_runs"]
                     + b["ue_eddy_peaks_per_run_sd"] ** 2 / b["n_runs"])
    return diff, pooled, bool(diff > 2.0 * pooled and diff > 0)

verdict = {}
for ctrl in ("DETUNE", "SHIFT"):
    d, p, ok = g3("LADDER", ctrl)
    verdict[f"LADDER_vs_{ctrl}"] = {"diff": round(d, 3), "pooled_sd": round(p, 3),
                                    "exceeds_2sd": ok}
verdict["ENRICHED"] = all(v["exceeds_2sd"] for v in verdict.values()
                          if isinstance(v, dict))

summary = {
    "total_peaks_logged": int(len(df)),
    "drive_peaks": int(df["is_drive"].sum()),
    "candidate_peaks": int(len(cand)),
    "candidate_peaks_unexplained": int((cand["audit_label"] == "UNEXPLAINED").sum()),
    "audit_label_counts_candidates": dict(Counter(cand["audit_label"])),
    "clusters_total": len(clusters),
    "clusters_explained": int(sum(1 for c in clusters if not c["unexplained"])),
    "clusters_unexplained": int(sum(1 for c in clusters if c["unexplained"])),
    "clusters_unexplained_eddy": int(sum(1 for c in clusters
                                         if c["unexplained"] and c["eddy"])),
    "clusters_by_arm": {
        arm: {
            "total": int(sum(1 for c in clusters if c["arm"] == arm)),
            "unexplained": int(sum(1 for c in clusters
                                    if c["arm"] == arm and c["unexplained"])),
            "unexplained_eddy": int(sum(1 for c in clusters
                                         if c["arm"] == arm and c["unexplained"]
                                         and c["eddy"])),
        } for arm in ("LADDER", "DETUNE", "SHIFT")},
    "arm_stats": arm_stats,
    "g3_verdict": verdict,
}
ue = [c for c in clusters if c["unexplained"]]
ue.sort(key=lambda c: -c["amp_median"])
summary["strongest_unexplained_clusters"] = ue[:10]
summary["all_clusters"] = clusters

json.dump(summary, open(f"{HD}/riemann_dual_cavity_results.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items()
                  if k not in ("all_clusters", "strongest_unexplained_clusters")},
                 indent=1))
print("strongest unexplained clusters:")
for c in ue[:5]:
    print(f"  {c['arm']} {c['modulation']}@{c['mod_rate']:.3f} f={c['freq_median']:.2f} "
          f"amp={c['amp_median']:.4f} seeds={c['seeds_present']} lit={c['n_lit_pickups']} eddy={c['eddy']}")
