"""
riemann_tp_cylinder_analyze.py — verdict analysis for the TP + cylinder
Riemann sim (pre-registration ITERATION_LOG 2026-10-22 14:38; adjustments
14:47). Consumes the sim outputs:
  main grid      riemann_tp_cylinder_{pickups,transients,energy}.csv + raw_meta
  flow sweep     riemann_fp0_*, riemann_fp015_* (FP = 0 / 0.15)
  fine pass      riemann_fine_* (window grid W2, coupling-power series)
Writes riemann_tp_cylinder_results.json and prints section summaries.
"""

import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HD = "/home/hatch/workspace/harness-dev"
import sys
sys.path.insert(0, HD)
from riemann_tp_cylinder_sim import audit_label, ARMS, MEAN_GAP  # noqa: E402

TOL_CLUSTER = 0.5
R = {}

meta = json.load(open(f"{HD}/riemann_tp_cylinder_raw_meta.json"))
F_EIG_A = np.array(meta["F_EIG_A"])
F_EIG_B = np.array(meta["F_EIG_B"])

# ------------------------------------------------------------- G1 / G2 --
R["g1"] = {"hits": meta["g1_drive_hits"], "total": meta["g1_drive_total"],
           "frac": meta["g1_drive_hits"] / max(meta["g1_drive_total"], 1)}
runs_meta = meta["runs"]
cf = np.array([r["closure_frac"] for r in runs_meta])
R["g2"] = {"closure_frac_median": float(np.median(cf)),
           "closure_frac_max": float(cf.max()),
           "runs_within_10pct": int((cf <= 0.10).sum()),
           "n_runs": len(cf)}

# ------------------------------------------------------------- exchange --
exch = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    sel = [r for r in runs_meta if r["arm"] == arm and r["modulation"] == "NONE"]
    exch[arm] = {
        k: float(np.median([r[k] for r in sel])) for k in
        ("E_inj", "E_diss_A", "E_diss_B", "E_rad_A", "E_rad_B",
         "E_stored_end_A", "E_stored_end_B", "closure_frac")}
    tot_out = exch[arm]["E_diss_A"] + exch[arm]["E_diss_B"] + \
        exch[arm]["E_rad_A"] + exch[arm]["E_rad_B"]
    exch[arm]["rad_share_of_output"] = float(
        (exch[arm]["E_rad_A"] + exch[arm]["E_rad_B"]) / tot_out)
    exch[arm]["B_share_of_dissipation"] = float(
        (exch[arm]["E_diss_B"] + exch[arm]["E_rad_B"]) / tot_out)
R["exchange_none_medians"] = exch

tr = pd.read_csv(f"{HD}/riemann_tp_cylinder_transients.csv")
trn = tr[tr["modulation"] == "NONE"]
cross = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    sub = trn[trn["arm"] == arm]
    med = sub.groupby("mode_idx")["transfer_ratio_driveend"].median()
    cross[arm] = {
        "transfer_ratio_by_mode": [round(float(x), 5) for x in med],
        "n_cross_gt_0.1": int((med > 0.1).sum()),
        "n_partial_0.01_0.1": int(((med > 0.01) & (med <= 0.1)).sum()),
        "n_home_lt_0.01": int((med <= 0.01).sum()),
        "modes_cross_gt_0.1": [int(i) for i in med[med > 0.1].index],
    }
lad = trn[trn["arm"] == "LADDER"].groupby("mode_idx")[
    "transfer_ratio_driveend"].median()
f_drive_lad = trn[trn["arm"] == "LADDER"].groupby("mode_idx")["f_drive"].first()
detune_nearest_B = np.array([np.min(np.abs(F_EIG_B - f)) for f in f_drive_lad])
rho, rho_p = spearmanr(lad.values, -detune_nearest_B)
cross["G5_spearman_transfer_vs_neg_detune_to_B"] = {
    "rho": float(rho), "p": float(rho_p), "gate_gt_0.5": bool(rho > 0.5)}
R["crossing"] = cross

# ----------------------------------------------------------- subspectra --
pk = pd.read_csv(f"{HD}/riemann_tp_cylinder_pickups.csv")
pk_de = pk[pk["window"] == "DRIVE_END"].copy()
# tolerance ladder re-audit (DRIVE_END candidates only)
for tol in (0.30, 0.15):
    lab = []
    for (arm, mod, rate), cell in pk_de.groupby(["arm", "modulation", "mod_rate"]):
        F = ARMS[arm]
        for _, r in cell.iterrows():
            lab.append("drive" if r["is_drive"] else
                       audit_label(r["peak_freq"], F, mod, rate, tol=tol))
    pk_de[f"label_{tol}"] = lab

def clusterize(cand, chan_col="pickup"):
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
                if not used[j] and abs(freqs[j] - freqs[i]) <= TOL_CLUSTER:
                    members.append(j)
                    used[j] = True
            sub = cell.iloc[members]
            seeds_present = sorted(sub["seed"].unique())
            if len(seeds_present) < 3:
                continue
            lit = []
            for ch, g in sub.groupby(chan_col):
                if g["seed"].nunique() / len(seeds_present) >= 0.5:
                    lit.append(ch)
            clusters.append({
                "arm": arm, "modulation": mod, "mod_rate": float(rate),
                "freq_median": round(float(sub["peak_freq"].median()), 3),
                "amp_median": round(float(sub["peak_amp"].median()), 5),
                "n_peaks": int(len(sub)),
                "seeds_present": [int(s) for s in seeds_present],
                "n_lit": len(lit), "lit": lit,
                "eddy": bool(2 <= len(lit) <= 7),
                "labels": dict(Counter(sub["audit_label"])),
                "unexplained_045": bool(set(sub["audit_label"]) == {"UNEXPLAINED"}),
                "unexplained_030": bool(set(sub["label_0.3"]) == {"UNEXPLAINED"}),
                "unexplained_015": bool(set(sub["label_0.15"]) == {"UNEXPLAINED"}),
            })
    return clusters

cand_de = pk_de[~pk_de["is_drive"]]
clusters = clusterize(cand_de)
sub_summary = {
    "total_peaks_driveend": int(len(pk_de)),
    "candidates": int(len(cand_de)),
    "unexplained_peaks_by_tol": {
        "0.45": int((cand_de["audit_label"] == "UNEXPLAINED").sum()),
        "0.30": int((cand_de["label_0.3"] == "UNEXPLAINED").sum()),
        "0.15": int((cand_de["label_0.15"] == "UNEXPLAINED").sum())},
    "audit_label_counts": dict(Counter(cand_de["audit_label"])),
    "clusters_total": len(clusters),
    "clusters_by_arm": {},
}
run_counts = defaultdict(int)
for c in clusters:
    if c["unexplained_045"] and c["eddy"]:
        for (arm, mod, rate), cell in cand_de.groupby(["arm", "modulation", "mod_rate"]):
            if arm == c["arm"] and mod == c["modulation"]:
                m = cell[np.abs(cell["peak_freq"] - c["freq_median"]) <= TOL_CLUSTER]
                for rid in m["run_id"].unique():
                    run_counts[(arm, int(rid))] += int((m["run_id"] == rid).sum())
arm_stats = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    rids = sorted(pk_de[pk_de["arm"] == arm]["run_id"].unique())
    counts = np.array([run_counts.get((arm, int(r)), 0) for r in rids], float)
    arm_stats[arm] = {"n_runs": len(rids),
                      "ue_eddy_per_run_mean": float(counts.mean()),
                      "ue_eddy_per_run_sd": float(counts.std(ddof=1))}
    sub_summary["clusters_by_arm"][arm] = {
        "total": int(sum(1 for c in clusters if c["arm"] == arm)),
        "unexplained_045": int(sum(1 for c in clusters
                                   if c["arm"] == arm and c["unexplained_045"])),
        "unexplained_eddy_045": int(sum(1 for c in clusters
                                        if c["arm"] == arm and c["unexplained_045"]
                                        and c["eddy"])),
        "eddy_fraction": float(np.mean([c["eddy"] for c in clusters
                                        if c["arm"] == arm])) if any(
            c["arm"] == arm for c in clusters) else None,
    }
def g3(a, b):
    A, B = arm_stats[a], arm_stats[b]
    diff = A["ue_eddy_per_run_mean"] - B["ue_eddy_per_run_mean"]
    pooled = np.sqrt(A["ue_eddy_per_run_sd"] ** 2 / A["n_runs"]
                     + B["ue_eddy_per_run_sd"] ** 2 / B["n_runs"])
    return {"diff": float(diff), "pooled_sd": float(pooled),
            "exceeds_2sd": bool(diff > 2 * pooled and diff > 0)}
sub_summary["g3"] = {c: g3("LADDER", c) for c in ("DETUNE", "SHIFT")}
sub_summary["g3"]["ENRICHED"] = all(
    v["exceeds_2sd"] for v in sub_summary["g3"].values() if isinstance(v, dict))
sub_summary["arm_stats"] = arm_stats
sub_summary["strongest_unexplained"] = sorted(
    [c for c in clusters if c["unexplained_045"]],
    key=lambda c: -c["amp_median"])[:10]
R["subspectra_driveend"] = sub_summary

# ------------------------------------------------------------ radiatives --
def dirtiness(df):
    tot = (df["peak_amp"] ** 2).sum()
    candp = (df[~df["is_drive"]]["peak_amp"] ** 2).sum()
    return float(candp / tot) if tot > 0 else float("nan")

rad = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    rad[arm] = {}
    for grp, seldf in (
            ("A_in", pk_de[(pk_de["arm"] == arm) & (pk_de["cavity"] == "A")
                           & (pk_de["pickup_kind"] == "in")]),
            ("A_out", pk_de[(pk_de["arm"] == arm) & (pk_de["pickup"] == "A_out")]),
            ("B_in", pk_de[(pk_de["arm"] == arm) & (pk_de["cavity"] == "B")
                           & (pk_de["pickup_kind"] == "in")]),
            ("B_out", pk_de[(pk_de["arm"] == arm) & (pk_de["pickup"] == "B_out")])):
        per_run = [dirtiness(g) for _, g in seldf.groupby("run_id")]
        rad[arm][grp] = float(np.median(per_run))
    # exterior drive recovery (NONE cells): drive freqs seen at out channel
    rec_out = {}
    for cav, out_ch in (("A", "A_out"), ("B", "B_out")):
        hits = tot = 0
        for (rid, seed), g in pk_de[(pk_de["arm"] == arm)
                                    & (pk_de["modulation"] == "NONE")].groupby(
                                        ["run_id", "seed"]):
            F = ARMS[arm]
            seen_out = g[g["pickup"] == out_ch]["peak_freq"].to_numpy()
            seen_in = g[(g["cavity"] == cav) & (g["pickup_kind"] == "in")][
                "peak_freq"].to_numpy()
            for fd in F:
                tot += 1
                if len(seen_out) and np.min(np.abs(seen_out - fd)) <= 0.45:
                    hits += 1
            rec_out.setdefault(cav, 0)
            rec_out[cav] += 0  # placeholder, totals tracked below
        rec_out[cav] = hits / max(tot, 1)
    rad[arm]["exterior_drive_recovery"] = rec_out
rad["G6_ladder_exterior_cleaner"] = {
    "A": bool(rad["LADDER"]["A_out"] < rad["LADDER"]["A_in"]),
    "B": bool(rad["LADDER"]["B_out"] < rad["LADDER"]["B_in"])}
R["radiatives_dirtiness"] = rad

# leakage ratio per component: out amp / mean interior amp at the drive peak.
# UNITS FIX (analysis note): the port series is velocity-weighted
# (sqrt(2*gamma_rad) * v * phi), so its spectral amps carry a factor ~2*pi*f
# vs the displacement interior pickups. Exterior amps are divided by
# 2*pi*f_peak (velocity -> displacement equivalent) before the ratio;
# without this the "ratio" is a units artifact growing linearly in f.
leak = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    F = ARMS[arm]
    sub = pk_de[(pk_de["arm"] == arm)].copy()
    out_mask = sub["pickup_kind"] == "out"
    sub.loc[out_mask, "peak_amp"] = (sub.loc[out_mask, "peak_amp"]
                                     / (2.0 * np.pi * sub.loc[out_mask, "peak_freq"]))
    table = {"A": [], "B": []}
    for cav, out_ch in (("A", "A_out"), ("B", "B_out")):
        for fd in F:
            ratios = []
            for rid, g in sub.groupby("run_id"):
                near = g[np.abs(g["peak_freq"] - fd) <= 0.45]
                o = near[near["pickup"] == out_ch]["peak_amp"]
                ii = near[(near["cavity"] == cav)
                          & (near["pickup_kind"] == "in")]["peak_amp"]
                if len(ii):
                    ratios.append(float(o.max() if len(o) else 0.0)
                                  / float(ii.mean()))
            table[cav].append(round(float(np.median(ratios)), 4) if ratios
                              else None)
    leak[arm] = table
R["leakage_ratio_by_component"] = leak
esc = trn[trn["arm"] == "LADDER"].groupby("mode_idx")["esc_frac_A_mode"].median()
R["modal_escape_fraction_A_ladder"] = [round(float(x), 4) for x in esc]
R["modal_escape_fraction_A_ladder_median"] = float(np.median(esc))

# ------------------------------------------------------------ transients --
tsum = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    sel = [r for r in runs_meta if r["arm"] == arm and r["modulation"] == "NONE"]
    tsum[arm] = {k: float(np.median([r[k] for r in sel])) for k in
                 ("t50_A", "t50_B", "tau_A", "tau_B")}
    tsum[arm]["t50_delay_cavity"] = tsum[arm]["t50_B"] - tsum[arm]["t50_A"]
delays_mode = trn.groupby("arm")["delay_s"].median()
tsum["per_mode_delay_median_s"] = {a: float(delays_mode[a]) for a in delays_mode.index}
R["transients_main"] = tsum
win_counts = {}
for wname in ("EARLY", "DRIVE_END", "SETTLED"):
    sub = pk[(pk["window"] == wname) & (~pk["is_drive"])]
    win_counts[wname] = {
        arm: {"candidate_peaks": int((sub["arm"] == arm).sum()),
              "unexplained_045": int(((sub["arm"] == arm)
                                      & (sub["audit_label"] == "UNEXPLAINED")).sum())}
        for arm in ("LADDER", "DETUNE", "SHIFT")}
R["window_candidate_counts"] = win_counts

# -------------------------------------------------- fine pass: attribution
fine = pd.read_csv(f"{HD}/riemann_fine_pickups.csv")
RING_W = ["RING1", "RING2", "SETTLED"]
ring = fine[fine["window"].isin(RING_W) & (~fine["is_drive"])]
def attribute(row):
    if row["audit_label"] != "UNEXPLAINED":
        return "drive_set_or_product"
    f = row["peak_freq"]
    if np.min(np.abs(F_EIG_A - f)) <= 0.45 or np.min(np.abs(F_EIG_B - f)) <= 0.45:
        return "eigenmode"
    return "unexplained"
ring = ring.assign(attr=[attribute(r) for _, r in ring.iterrows()])
attr_frac = {}
for kind in ("in", "out"):
    s = ring[ring["pickup_kind"] == kind]
    attr_frac[kind] = {k: float(v) for k, v in
                       (s["attr"].value_counts(normalize=True).items())}
    attr_frac[kind]["n"] = int(len(s))
R["ringdown_attribution"] = attr_frac
R["ringdown_prediction_exterior_more_eigen"] = bool(
    attr_frac["out"].get("eigenmode", 0) > attr_frac["in"].get("eigenmode", 0))
fine_win = {}
for wname in ["T1", "T2", "DRIVE_END", "RING1", "RING2", "SETTLED"]:
    s = fine[fine["window"] == wname]
    fine_win[wname] = {
        "peaks": int(len(s)),
        "candidates": int((~s["is_drive"]).sum()),
        "unexplained_045": int(((~s["is_drive"])
                                & (s["audit_label"] == "UNEXPLAINED")).sum()),
        "cand_out_share": float((~s["is_drive"] & (s["pickup_kind"] == "out")).sum()
                                / max((~s["is_drive"]).sum(), 1))}
R["fine_window_counts"] = fine_win

# ------------------------------------------------------- asymmetry (fine) --
fe = pd.read_csv(f"{HD}/riemann_fine_energy.csv")
fen = fe[fe["modulation"] == "NONE"]
asym = {}
for arm in ("LADDER", "SHIFT"):
    rows = []
    for rid, g in fen[fen["arm"] == arm].groupby("run_id"):
        atk = g[(g["t"] >= 0.5) & (g["t"] < 4.0)]
        rng_ = g[(g["t"] >= 4.0) & (g["t"] < 7.6)]
        i_pk = rng_["P_BA"].idxmax()
        rows.append({
            "attack_P_AB": float(atk["P_AB"].mean()),
            "attack_P_BA": float(atk["P_BA"].mean()),
            "ring_P_AB": float(rng_["P_AB"].mean()),
            "ring_P_BA": float(rng_["P_BA"].mean()),
            "backflow_peak_lag_s": float(g.loc[i_pk, "t"] - 4.0)})
    asym[arm] = {k: float(np.median([r[k] for r in rows])) for k in rows[0]}
R["asymmetry_coupling_power"] = asym

# --------------------------------------------------------------- pulling --
def settled_offsets(path, fp):
    d = pd.read_csv(path)
    d = d[(d["window"] == "SETTLED") & (d["arm"] == "LADDER")
          & (d["modulation"] == "NONE") & (d["cavity"] == "A")]
    F = ARMS["LADDER"]
    offs = []
    for n, fd in enumerate(F, start=1):
        near = d[np.abs(d["peak_freq"] - fd) <= 1.5]
        vals = []
        for rid, g in near.groupby("run_id"):
            vals.append(float(g.loc[g["peak_amp"].idxmax(), "peak_freq"]) - fd)
        offs.append(float(np.median(vals)) if vals else float("nan"))
    return np.array(offs)

off0 = settled_offsets(f"{HD}/riemann_fp0_pickups.csv", 0.0)
off25 = settled_offsets(f"{HD}/riemann_tp_cylinder_pickups.csv", 0.025)
off150 = settled_offsets(f"{HD}/riemann_fp015_pickups.csv", 0.15)
n_idx = np.arange(1, 31)
slopes = []
for i in range(30):
    x = np.array([0.0, 0.025 * n_idx[i], 0.15 * n_idx[i]])
    y = np.array([off0[i], off25[i], off150[i]])
    if np.isfinite(y).all():
        slopes.append(float(np.polyfit(x, y, 1)[0]))
    else:
        slopes.append(float("nan"))
slopes = np.array(slopes)
R["pulling"] = {
    "per_component_slope": [round(float(s), 3) if s == s else None for s in slopes],
    "median_slope": float(np.nanmedian(slopes)),
    "frac_components_nearer_eigen_at_fp015": float(np.mean(
        [abs(off150[i] - 0.15 * n_idx[i]) < abs(off150[i])
         for i in range(30) if np.isfinite(off150[i])])),
    "pulling_claim_gate_median_slope_gt_0.5": bool(np.nanmedian(slopes) > 0.5),
    "median_offset_fp0": float(np.nanmedian(np.abs(off0))),
    "median_offset_fp015": float(np.nanmedian(np.abs(off150))),
}
R["all_clusters"] = clusters

json.dump(R, open(f"{HD}/riemann_tp_cylinder_results.json", "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k != "all_clusters"}, indent=1)[:6000])
print("clusters:", len(clusters))
