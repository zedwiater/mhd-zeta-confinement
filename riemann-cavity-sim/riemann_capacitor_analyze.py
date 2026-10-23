"""
riemann_capacitor_analyze.py — iteration-6 verdicts.
Reads riemann_capacitor_{peaks_raw,series,runs}.csv + raw meta +
calib + criterion JSONs.
Pre-registration: ITERATION_LOG.md 2026-10-22 19:3x EDT (gates
GC0-GC6, criterion, audit6 — stated BEFORE the grid; this script
only evaluates them).

Writes:
- riemann_capacitor_scatter.csv (every peak + audit6 + distances)
- riemann_capacitor_grid.csv    (placement search: coordinates +
  criterion values + ranking, incl. refine cells)
- riemann_capacitor_log.csv     (per run: ledger/stability/coords)
- riemann_capacitor_results.json (all verdict numbers)
Audit6: audit5 (base 0.45/0.30/0.15 + flow 0.15 + inducer 0.15)
plus capacitor targets at 0.15 built from BOTH f_bare and the
calibration f_ring: f_cap, 2f_cap, f_cap/2, tone+-k*f_cap,
f_i+-k*f_cap (k=1,2). Engine-expected != unique.
"""

import itertools
import json
from collections import defaultdict

import numpy as np
import pandas as pd

import riemann_tp_cylinder_sim as S

HD = S.HD
FP = 0.025
LAD = S.T_LADDER
EIGA = S.F_EIG_A
EIGB = S.F_EIG_B
BASE12 = list(S.CH_NAMES)

pk = pd.read_csv(f"{HD}/riemann_capacitor_peaks_raw.csv")
runs = pd.read_csv(f"{HD}/riemann_capacitor_runs.csv")
ser = pd.read_csv(f"{HD}/riemann_capacitor_series.csv")
meta = json.load(open(f"{HD}/riemann_capacitor_raw_meta.json"))
crit = json.load(open(f"{HD}/riemann_capacitor_criterion.json"))
calib = json.load(open(f"{HD}/riemann_capacitor_calib.json"))
CELLS = meta["cells"]
TONES = {c: (v["tones"] if v else []) for c, v in CELLS.items()}
F_CAPS = sorted({calib["f_bare"], calib["ring"]["f_ring"]})
WINNER = crit["winner"]


# ---------------------------------------------------------------- audit --
def flow_label(fc, F):
    s = FP * np.arange(1, len(F) + 1)
    def hit(t):
        return t > 0.5 and abs(fc - t) <= 0.15
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


def ind_label5(fc, F, tones):
    if not tones:
        return None
    s = FP * np.arange(1, len(F) + 1)
    def hit(t):
        return t > 0.5 and abs(fc - t) <= 0.15
    for tone in tones:
        if hit(2 * tone) or hit(tone / 2):
            return "ind_harm"
        for fi in F:
            if hit(tone + fi) or hit(abs(tone - fi)) \
                    or hit(2 * tone + fi) or hit(abs(2 * tone - fi)):
                return "ind_mix"
        for sv in s:
            for k in (1, 2):
                if hit(tone + k * sv) or hit(tone - k * sv):
                    return "ind_flow"
    for a, b in itertools.combinations(tones, 2):
        if hit(a + b) or hit(abs(a - b)):
            return "ind_comb"
    return None


def cap_label(fc, F, tones):
    def hit(t):
        return t > 0.5 and abs(fc - t) <= 0.15
    for fc_ in F_CAPS:
        if hit(fc_):
            return "cap_res"
        if hit(2 * fc_) or hit(fc_ / 2):
            return "cap_harm"
        for tone in tones:
            for k in (1, 2):
                if hit(tone + k * fc_) or hit(abs(tone - k * fc_)):
                    return "cap_sb"
        for fi in F:
            for k in (1, 2):
                if hit(fi + k * fc_) or hit(abs(fi - k * fc_)):
                    return "cap_mix"
    return None


def audit6(fc, F, tones):
    """(u045, u030, u015_extended, label6)."""
    u45 = S.audit_label(fc, F, "NONE", 0.0, tol=0.45) == "UNEXPLAINED"
    u30 = S.audit_label(fc, F, "NONE", 0.0, tol=0.30) == "UNEXPLAINED"
    u15 = S.audit_label(fc, F, "NONE", 0.0, tol=0.15) == "UNEXPLAINED"
    lab = "UNEXPLAINED"
    if not u15:
        lab = S.audit_label(fc, F, "NONE", 0.0, tol=0.15)
    if u15 and flow_label(fc, F) is not None:
        u15 = False
        lab = flow_label(fc, F)
    if u15:
        il = ind_label5(fc, F, tones)
        if il is not None:
            u15 = False
            lab = il
    if u15:
        cl = cap_label(fc, F, tones)
        if cl is not None:
            u15 = False
            lab = cl
    return u45, u30, u15, lab


# ------------------------------------------------------- scatter build --
def add_distances(df):
    f = df["peak_freq"].to_numpy()
    df["d_line"] = np.min(np.abs(f[:, None] - LAD[None, :]), axis=1)
    df["d_off"] = np.min(np.abs(f[:, None] - EIGA[None, :]), axis=1)
    df["d_offB"] = np.min(np.abs(f[:, None] - EIGB[None, :]), axis=1)
    return df


cand_mask = (~pk["is_drive"]) & (~pk["is_inducer"])
au = [audit6(r.peak_freq, S.ARMS[r.arm], TONES[r.cell])
      for r in pk[cand_mask].itertuples()]
pk.loc[cand_mask, ["u045", "u030", "u015x"]] = np.array(
    [(a, b, c) for a, b, c, _l in au], dtype=bool)
pk.loc[cand_mask, "audit6_label"] = np.array(
    [l for _a, _b, _c, l in au], dtype=object)
pk.loc[~cand_mask, "audit6_label"] = np.where(
    pk.loc[~cand_mask, "is_inducer"], "inducer_tone", "drive")
pk = add_distances(pk)
pk.to_csv(f"{HD}/riemann_capacitor_scatter.csv", index=False)
cand = pk[cand_mask].copy()

R = {"calibration": calib, "criterion": crit,
     "f_caps_audit": F_CAPS, "winner": WINNER}

# ------------------------------------------------------------- clusters --
def chain_clusters(df, tol=0.45):
    out = []
    if df.empty:
        return out
    d = df.sort_values("peak_freq")
    freqs = d["peak_freq"].to_numpy()
    groups, cur = [], [0]
    for i in range(1, len(freqs)):
        if freqs[i] - freqs[i - 1] <= tol:
            cur.append(i)
        else:
            groups.append(cur); cur = [i]
    groups.append(cur)
    for g in groups:
        sub = d.iloc[g]
        seeds = sorted(sub["seed"].unique())
        lit = sorted(sub["pickup"].unique())
        out.append({
            "freq_median": float(sub["peak_freq"].median()),
            "amp_median": float(sub["peak_amp"].median()),
            "n_peaks": int(len(sub)), "seeds": [int(x) for x in seeds],
            "n_seeds": len(seeds), "lit": lit, "n_lit": len(lit),
            "eddy": bool(2 <= len(lit) <= 7),
            "u015x": bool(sub["u015x"].mean() > 0.5)})
    return out


clusters, unexp_clusters = {}, {}
for (cell, arm, win), cell_df in cand.groupby(["cell", "arm", "window"]):
    cl = [c for c in chain_clusters(cell_df) if c["n_seeds"] >= 2]
    clusters[f"{cell}|{arm}|{win}"] = {
        "n": len(cl), "eddy": int(sum(c["eddy"] for c in cl)),
        "list": cl}
    unexp_clusters[f"{cell}|{arm}|{win}"] = [
        c for c in cl if c["u015x"]]
R["clusters_reproducible"] = clusters
R["unexplained_clusters"] = {k: v for k, v in unexp_clusters.items() if v}
R["gc5_unexplained_by_cell"] = {
    k: len(v) for k, v in unexp_clusters.items() if v}

# ------------------------------------------------------- gates GC0-GC3 --
off = runs[(runs["cell"] == "OFFCAP") & (runs["arm"] == "LADDER")]
R["gc0_offcap_ind_frac_median"] = float((off["E_ind"] / off["E_inj"]).median())
R["gc0_iter5_target"] = 0.016569
R["gc0_pass"] = bool(abs(R["gc0_offcap_ind_frac_median"] / 0.016569 - 1) <= 0.15)
R["gc1_drive_recovery"] = {k: v for k, v in meta["g1"].items()}
R["gc1_min_recovery"] = float(min(v["hits"] / v["total"]
                                  for v in meta["g1"].values()))
R["gc2_max_closure_frac"] = float(runs["closure_frac"].max())
R["gc2_pass"] = bool(runs["closure_frac"].max() <= 0.001)
R["gc3_max_slope_settled"] = float(runs["slope_settled"].max())
R["gc3_max_sat_frac_G"] = float(runs["sat_frac"].max())
R["gc3_clamp_count_total"] = int(runs["clamp_count"].sum())
R["gc3_pass"] = bool(runs["slope_settled"].max() <= 0.005
                     and runs["sat_frac"].max() <= 0.30
                     and runs["clamp_count"].sum() == 0)
cap_runs = runs[runs["cap_x"].notna()]
R["E_cap_flow_resid_max_abs"] = float(cap_runs["E_cap_flow_resid"].abs().max()) \
    if len(cap_runs) else None

# ------------------------------------------------------------ grid CSV --
def persistence_of(df):
    de = df["demod_t0_A_in5_DRIVE_END_amp"]
    st = df["demod_t0_A_in5_SETTLED_amp"]
    return st / de.where(de > 1e-9)


grid_rows = []
for cell, v in CELLS.items():
    if v is None or v.get("cap") is None:
        continue
    g = runs[(runs["cell"] == cell) & (runs["arm"] == "LADDER")]
    if g.empty:
        continue
    ps = persistence_of(g)
    t = crit["table"].get(cell, {})
    grid_rows.append({
        "cell": cell, "stage": ("refine" if cell.startswith("RFC") else "coarse"),
        "cap_x": v["cap"]["x"], "cap_theta_deg": v["cap"]["theta_deg"],
        "ind_x": v["x"], "ind_theta_deg": v["theta_deg"],
        "transducer_x": 0.79,
        "P_median": float(ps.median()), "P_per_seed": list(ps.round(4)),
        "guards_all_pass": t.get("guards_all_pass"),
        "P_criterion": t.get("P_median"),
        "ind_frac_median": float((g["E_ind"] / g["E_inj"]).median()),
        "closure_max": float(g["closure_frac"].max()),
        "sat_max": float(g["sat_frac"].max()),
        "slope_max": float(g["slope_settled"].max()),
        "E_cap_end_median": float(g["E_cap_end"].median()),
        "is_winner": cell == WINNER})
grid = pd.DataFrame(grid_rows).sort_values("P_median", ascending=False)
grid.to_csv(f"{HD}/riemann_capacitor_grid.csv", index=False)
R["grid_ranking"] = grid.to_dict("records")

offp = persistence_of(off)
R["gc4_offcap_P_median"] = float(offp.median())
if WINNER and WINNER in crit["table"]:
    wp = np.array(crit["table"][WINNER]["P_per_seed"])
    pooled = float(np.std(np.concatenate([offp.to_numpy(), wp]), ddof=1))
    R["gc4_winner_P_median"] = float(np.median(wp))
    R["gc4_pooled_sd"] = pooled
    R["gc4_pass"] = bool(np.median(wp) - offp.median() > pooled)
else:
    R["gc4_pass"] = False

# ------------------------------------------------- GI4-analog (GC6) --
def dist_stats(df):
    if df.empty:
        return {"n": 0, "med_d_line": None, "med_d_off": None,
                "frac_line_015": None, "frac_off_015": None}
    return {"n": int(len(df)),
            "med_d_line": float(df["d_line"].median()),
            "med_d_off": float(df["d_off"].median()),
            "frac_line_015": float((df["d_line"] <= 0.15).mean()),
            "frac_off_015": float((df["d_off"] <= 0.15).mean())}


settled = cand[(cand["window"] == "SETTLED")
               & (cand["pickup"].isin(BASE12))]
dist_base, seed_meds = {}, {}
for (cell, arm), g in settled.groupby(["cell", "arm"]):
    dist_base[f"{cell}|{arm}"] = dist_stats(g)
    seed_meds[f"{cell}|{arm}"] = {
        str(sd): float(x["d_off"].median())
        for sd, x in g.groupby("seed")}
R["dist_base12_settled"] = dist_base


def pooled_sd(d):
    v = np.array(list(d.values()), dtype=float)
    return float(np.std(v, ddof=1)) if len(v) > 1 else 0.0


align = {}
for cell in ("OFFCAP", "PLACEBO", WINNER):
    if cell is None:
        continue
    for arm in ("LADDER", "DETUNE", "SHIFT"):
        k, ref = f"{cell}|{arm}", f"OFFCAP|{arm}"
        if k in dist_base and ref in dist_base \
                and dist_base[k]["med_d_off"] is not None:
            align[k] = {
                "delta_off": dist_base[k]["med_d_off"] - dist_base[ref]["med_d_off"],
                "delta_line": dist_base[k]["med_d_line"] - dist_base[ref]["med_d_line"],
                "med_d_off": dist_base[k]["med_d_off"],
                "med_d_line": dist_base[k]["med_d_line"],
                "n": dist_base[k]["n"]}
R["alignment_vs_offcap"] = align
if WINNER:
    ctrl_deltas = [abs(align[f"{WINNER}|{a}"]["delta_off"])
                   for a in ("DETUNE", "SHIFT") if f"{WINNER}|{a}" in align]
    pl_d = abs(align.get("PLACEBO|LADDER", {}).get("delta_off", 0.0))
    spread = max(ctrl_deltas + [pl_d]) if ctrl_deltas else pl_d
    wd = align.get(f"{WINNER}|LADDER", {}).get("delta_off")
    R["gc6_control_spread_off"] = float(spread)
    R["gc6_winner_delta_off"] = wd
    R["gc6_alignment_claim"] = bool(wd is not None and abs(wd) > spread)

# ------------------------------------------------------------- trap --
nsw = runs[runs["cell"].str.startswith("NSW_") & (runs["arm"] == "LADDER")]
trap_rows = []
for cell, g in nsw.groupby("cell"):
    d = float(cell.split("_d")[1])
    trap_rows.append({"d": d,
                      "ind_frac": float((g["E_ind"] / g["E_inj"]).median())})
trap = pd.DataFrame(trap_rows).sort_values("d")
if len(trap):
    ends = trap.iloc[[0, -1]]["ind_frac"].mean()
    imin = trap["ind_frac"].idxmin()
    R["trap_with_capacitor"] = {
        "d_star": float(trap.loc[imin, "d"]),
        "min_frac": float(trap.loc[imin, "ind_frac"]),
        "plateau": float(ends),
        "depth": float(trap.loc[imin, "ind_frac"] / ends),
        "curve": trap.to_dict("records")}
R["trap_iter5_G"] = {"d_star": 0.24, "depth": 0.164,
                     "note": "iter-5 G row, RIEMANN_NODE2D.md table"}

# ------------------------------------------------------- 0.3176 scan --
scan = {}
for cell in ("OFFCAP", "PLACEBO", WINNER):
    if cell is None:
        continue
    g = cand[(cand["cell"] == cell) & (cand["window"] == "SETTLED")
             & (cand["peak_freq"] >= 0.25) & (cand["peak_freq"] <= 0.40)]
    scan[cell] = {"n_peaks": int(len(g)),
                  "freqs": sorted(g["peak_freq"].round(4).unique().tolist()),
                  "max_amp": float(g["peak_amp"].max()) if len(g) else None,
                  "seeds": sorted(int(x) for x in g["seed"].unique())}
R["scan_03176"] = scan

# ------------------------------------------------------------ log CSV --
log_cols = ["run_id", "cell", "arm", "seed", "ind_x", "ind_theta_deg",
            "cap_x", "cap_theta_deg", "E_inj", "E_ind", "E_cap_end",
            "E_cap_flow_resid", "closure_frac", "sat_frac", "gate_frac",
            "slope_settled", "clamp_count", "pool_max",
            "demod_t0_A_in5_DRIVE_END_amp", "demod_t0_A_in5_SETTLED_amp"]
runs[log_cols].to_csv(f"{HD}/riemann_capacitor_log.csv", index=False)

json.dump(R, open(f"{HD}/riemann_capacitor_results.json", "w"), indent=1,
          default=str)
print("GC0 ind_frac", R["gc0_offcap_ind_frac_median"], "pass", R["gc0_pass"])
print("GC1 min recovery", R["gc1_min_recovery"])
print("GC2 max closure", R["gc2_max_closure_frac"], "pass", R["gc2_pass"])
print("GC3 slope/sat/clamps", R["gc3_max_slope_settled"],
      R["gc3_max_sat_frac_G"], R["gc3_clamp_count_total"], "pass", R["gc3_pass"])
print("GC4 winner", WINNER, "P off", R.get("gc4_offcap_P_median"),
      "P win", R.get("gc4_winner_P_median"), "pass", R.get("gc4_pass"))
print("GC5 unexplained cells", R["gc5_unexplained_by_cell"])
print("GC6 winner delta_off", R.get("gc6_winner_delta_off"),
      "spread", R.get("gc6_control_spread_off"),
      "claim", R.get("gc6_alignment_claim"))
print("TRAP", R.get("trap_with_capacitor", {}).get("d_star"),
      R.get("trap_with_capacitor", {}).get("depth"))
print("03176", {k: (v["n_peaks"], v["freqs"][:4]) for k, v in scan.items()})
print("DONE analyze")
