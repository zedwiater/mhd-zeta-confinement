"""
riemann_node2d_analyze.py — iteration-5 verdicts.
Reads riemann_node2d_{peaks_raw,series,runs}.csv + raw meta.
Pre-registration: ITERATION_LOG.md 2026-10-22 17:44 EDT (gates GN0-GN6,
node metric definitions, discriminator criterion — all stated there
BEFORE the grid ran; this script only evaluates them).

Writes:
- riemann_node2d_scatter.csv (every peak + audit5 + distances)
- riemann_node2d_log.csv     (per run: ledger/stability/demod log)
- riemann_node2d_grid.csv    (the full 2D (x, theta) metric grid)
- riemann_node2d_results.json (all verdict numbers)
Audit5: iter-4 audit4 (base 0.45/0.30/0.15 + flow 0.15) with the
inducer label generalized to the cell's tone list (+ pairwise tone
sums/differences for comb cells). Engine-expected != unique.
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

pk = pd.read_csv(f"{HD}/riemann_node2d_peaks_raw.csv")
runs = pd.read_csv(f"{HD}/riemann_node2d_runs.csv")
ser = pd.read_csv(f"{HD}/riemann_node2d_series.csv")
meta = json.load(open(f"{HD}/riemann_node2d_raw_meta.json"))
CELLS = meta["cells"]
TONES = {c: (v["tones"] if v else []) for c, v in CELLS.items()}


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


def audit5(fc, F, tones):
    """(u045, u030, u015_extended, label5)."""
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
    return u45, u30, u15, lab


# ------------------------------------------------------- scatter build --
def add_distances(df):
    f = df["peak_freq"].to_numpy()
    df["d_line"] = np.min(np.abs(f[:, None] - LAD[None, :]), axis=1)
    df["d_off"] = np.min(np.abs(f[:, None] - EIGA[None, :]), axis=1)
    df["d_offB"] = np.min(np.abs(f[:, None] - EIGB[None, :]), axis=1)
    return df


cand_mask = (~pk["is_drive"]) & (~pk["is_inducer"])
au = [audit5(r.peak_freq, S.ARMS[r.arm], TONES[r.cell])
      for r in pk[cand_mask].itertuples()]
pk.loc[cand_mask, ["u045", "u030", "u015x"]] = np.array(
    [(a, b, c) for a, b, c, _l in au], dtype=bool)
pk.loc[cand_mask, "audit5_label"] = np.array(
    [l for _a, _b, _c, l in au], dtype=object)
pk.loc[~cand_mask, "audit5_label"] = np.where(
    pk.loc[~cand_mask, "is_inducer"], "inducer_tone", "drive")
pk = add_distances(pk)
pk.to_csv(f"{HD}/riemann_node2d_scatter.csv", index=False)
cand = pk[cand_mask].copy()

R = {"tones": meta["tones"], "AIND": meta["AIND"],
     "profile_check": meta["profile_check"]}

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
        "n": len(cl), "eddy": int(sum(c["eddy"] for c in cl))}
    unexp_clusters[f"{cell}|{arm}|{win}"] = [
        c for c in cl if c["u015x"]]
R["clusters_reproducible_n"] = clusters
R["unexplained_clusters"] = {k: v for k, v in unexp_clusters.items() if v}
R["gn5_unexplained_cluster_cells"] = {
    k: len(v) for k, v in unexp_clusters.items() if v}
R["gn5_pass"] = not R["gn5_unexplained_cluster_cells"]

# ------------------------------------------------------ per-run helpers --
def settled_slope(run_id, col):
    g = ser[ser["run_id"] == run_id]
    w = g[(g["t"] >= S.WINDOWS["SETTLED"][0])
          & (g["t"] < S.WINDOWS["SETTLED"][1])]
    e = w[col].to_numpy()
    t = w["t"].to_numpy()
    ok = e > 1e-12
    if ok.sum() < 10:
        return float("nan")
    return float(np.polyfit(t[ok], np.log(e[ok]), 1)[0])


def global_rel_slope(run_id):
    g = ser[ser["run_id"] == run_id]
    w = g[(g["t"] >= S.WINDOWS["SETTLED"][0])
          & (g["t"] < S.WINDOWS["SETTLED"][1])]
    e = w["E_A"].to_numpy()
    t = w["t"].to_numpy()
    if len(e) < 10 or e.mean() <= 0:
        return float("nan")
    return float(np.polyfit(t, e, 1)[0] / e.mean())


def cell_power(cell, arm, pickup, win="SETTLED"):
    g = cand[(cand["cell"] == cell) & (cand["arm"] == arm)
             & (cand["pickup"] == pickup) & (cand["window"] == win)]
    return float((g["peak_amp"] ** 2).sum()) if len(g) else 0.0


def out_share_run(run_id):
    g = pk[(pk["run_id"] == run_id) & (pk["window"] == "SETTLED")]
    rin = g[g["pickup"].str.startswith("A_rin")]
    rout = g[g["pickup"].str.startswith("A_rout")]
    den = (rin["peak_amp"] ** 2).sum() + (rout["peak_amp"] ** 2).sum()
    if den <= 0:
        return float("nan")
    return float((rout["peak_amp"] ** 2).sum() / den)


off_pow_A_in5 = cell_power("OFF", "LADDER", "A_in5")
log_rows = []
for r in runs.itertuples():
    cell = CELLS[r.cell]
    row = {"run_id": r.run_id, "cell": r.cell, "arm": r.arm,
           "seed": r.seed,
           "tones": ";".join(f"{t:.4f}" for t in TONES[r.cell]),
           "x": (cell["x"] if cell else None),
           "theta_deg": (float(np.degrees(cell["theta"]))
                         if cell else None),
           "E_inj": r.E_inj, "E_ind": r.E_ind,
           "eind_0": r.eind_0, "eind_1": r.eind_1, "eind_2": r.eind_2,
           "ind_frac": r.E_ind / r.E_inj if r.E_inj > 0 else float("nan"),
           "closure_frac": r.closure_frac,
           "ext_frac": r.E_tr_ext / r.E_inj if r.E_inj > 0 else float("nan"),
           "sat_frac": r.sat_frac, "gate_frac": r.gate_frac,
           "pool_max": r.pool_max, "clamp_count": r.clamp_count,
           "settled_slope_A79": settled_slope(r.run_id, "eloc_A79"),
           "global_rel_slope_EA": global_rel_slope(r.run_id),
           "out_share_settled": out_share_run(r.run_id),
           "tau_A": getattr(r, "tau_A", float("nan")),
           "tau_B": getattr(r, "tau_B", float("nan"))}
    log_rows.append(row)
logdf = pd.DataFrame(log_rows)
logdf.to_csv(f"{HD}/riemann_node2d_log.csv", index=False)
R["log_rows"] = log_rows

# ---------------------------------------------------------- node sweep --
def sweep_table(prefix, ds):
    rows = []
    for d in ds:
        cell = f"{prefix}_d{d:.2f}"
        g = logdf[(logdf["cell"] == cell) & (logdf["arm"] == "LADDER")]
        rows.append({"d": d, "cell": cell,
                     "ind_frac_med": float(g["ind_frac"].median()),
                     "ind_frac_sd": float(g["ind_frac"].std(ddof=1)),
                     "per_seed": [float(v) for v in g["ind_frac"]]})
    return rows


def node_metrics(rows):
    d = np.array([r["d"] for r in rows])
    f = np.array([r["ind_frac_med"] for r in rows])
    i = int(np.argmin(f))
    plateau = float(np.mean([f[list(d).index(0.12)],
                             f[list(d).index(0.28)]]))
    thr = f[i] + 0.5 * (plateau - f[i])
    # contiguous span around argmin at/below threshold, interpolated
    lo = i
    while lo > 0 and f[lo - 1] <= thr:
        lo -= 1
    if lo < i:
        x0 = d[lo] + (d[lo + 1] - d[lo]) * (thr - f[lo]) / (f[lo + 1] - f[lo])
    else:
        x0 = d[lo]
    hi = i
    while hi < len(f) - 1 and f[hi + 1] <= thr:
        hi += 1
    if hi > i:
        x1 = d[hi] + (d[hi - 1] - d[hi]) * (thr - f[hi]) / (f[hi - 1] - f[hi])
    else:
        x1 = d[hi]
    return {"d_star": float(d[i]), "min_frac": float(f[i]),
            "plateau": plateau, "depth": float(f[i] / plateau),
            "width": float(x1 - x0), "span": [float(x0), float(x1)],
            "table": rows}


G_DS = [0.10, 0.12, 0.14, 0.16, 0.18, 0.19, 0.20, 0.21, 0.22, 0.24,
        0.26, 0.28]
ALT_DS = [0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28]
node = {}
node["G"] = node_metrics(sweep_table("G", G_DS))
for pfx in ("MID", "GM", "GP"):
    node[pfx] = node_metrics(sweep_table(pfx, ALT_DS))
R["node"] = node
g = node["G"]
disc = {"criterion": "SPATIAL iff all alt d* within +/-0.02 of G d* "
        "and widths within 2x; TRAP iff any alt d* shift >= 0.04 or "
        "node vanishes (depth>0.70) at an alt tone while G depth<=0.50"}
trap = False
spatial = True
for pfx in ("MID", "GM", "GP"):
    a = node[pfx]
    disc[pfx] = {"d_star": a["d_star"], "shift": a["d_star"] - g["d_star"],
                 "depth": a["depth"], "width": a["width"]}
    if abs(a["d_star"] - g["d_star"]) >= 0.04 \
            or (a["depth"] > 0.70 and g["depth"] <= 0.50):
        trap = True
    if abs(a["d_star"] - g["d_star"]) > 0.02 \
            or not (0.5 <= a["width"] / g["width"] <= 2.0):
        spatial = False
disc["verdict"] = "TRAP" if trap else ("SPATIAL" if spatial else "MIXED")
R["discriminator"] = disc

# ------------------------------------------------------------- 2D grid --
GRID_SRC = {}
for x in (0.59, 0.65, 0.71, 0.85, 0.91):
    for adeg in (0, 30, 60, 90, 120, 150, 180):
        if adeg == 0 and x == 0.59:
            GRID_SRC[(x, adeg)] = "G_d0.20"
        elif adeg == 0 and x == 0.65:
            GRID_SRC[(x, adeg)] = "G_d0.14"
        else:
            GRID_SRC[(x, adeg)] = f"G2D_x{x:.2f}_t{adeg:03d}"
grid_rows = []
for (x, adeg), cell in sorted(GRID_SRC.items()):
    lg = logdf[(logdf["cell"] == cell) & (logdf["arm"] == "LADDER")]
    rr = runs[(runs["cell"] == cell) & (runs["arm"] == "LADDER")]
    p = cell_power(cell, "LADDER", "A_in5")
    grid_rows.append({
        "x": x, "theta_deg": adeg, "cell": cell,
        "ind_frac_med": float(lg["ind_frac"].median()),
        "ind_frac_sd": float(lg["ind_frac"].std(ddof=1)),
        "demod_A_in5_SETTLED_med": float(
            rr["demod_t0_A_in5_SETTLED_amp"].median()),
        "demod_rin0_SETTLED_med": float(
            rr["demod_t0_rin0_SETTLED_amp"].median()),
        "settled_site_power_ratio_vsOFF": (
            float(p / off_pow_A_in5) if off_pow_A_in5 > 0 else None),
        "out_share_settled_med": float(lg["out_share_settled"].median()),
        "sat_frac_med": float(lg["sat_frac"].median()),
        "ext_frac_med": float(lg["ext_frac"].median())})
griddf = pd.DataFrame(grid_rows)
griddf.to_csv(f"{HD}/riemann_node2d_grid.csv", index=False)
R["grid"] = grid_rows

# node shape: x = 0.59 row across theta, theta = 0 column across x
thr_node = g["min_frac"] + 0.5 * (g["plateau"] - g["min_frac"])
row59 = griddf[griddf["x"] == 0.59].sort_values("theta_deg")
in_node_row = row59[row59["ind_frac_med"] <= thr_node]
shape = {"threshold": float(thr_node),
         "row_x0.59": [{"theta_deg": int(r.theta_deg),
                        "ind_frac_med": float(r.ind_frac_med)}
                       for r in row59.itertuples()],
         "n_angles_in_node_row59": int(len(in_node_row)),
         "col_theta0": [{"x": float(r.x),
                         "ind_frac_med": float(r.ind_frac_med)}
                        for r in griddf[griddf["theta_deg"] == 0]
                        .sort_values("x").itertuples()]}
if len(in_node_row) >= 5:
    shape["character"] = "LINE across theta at x=0.59 (ring in the " \
        "reduced pairing): node ind_frac at >=5 of 7 angles"
elif len(in_node_row) <= 2:
    shape["character"] = "POINT-like: node only near theta=0 at x=0.59"
else:
    shape["character"] = "PATCH/ARC: node at 3-4 of 7 angles at x=0.59"
R["node_shape"] = shape

# out-of-plane skew: per radius, share vs theta; mirror-pair asymmetry
skew = {}
for x, gx in griddf.groupby("x"):
    s = gx.sort_values("theta_deg")
    pairs = {}
    for a, b in ((30, 150), (60, 120)):
        va = s[s["theta_deg"] == a]["out_share_settled_med"]
        vb = s[s["theta_deg"] == b]["out_share_settled_med"]
        if len(va) and len(vb):
            pairs[f"{a}vs{b}"] = float(va.iloc[0] - vb.iloc[0])
    skew[f"x={x}"] = {
        "share_min": float(s["out_share_settled_med"].min()),
        "share_max": float(s["out_share_settled_med"].max()),
        "share_range": float(s["out_share_settled_med"].max()
                             - s["out_share_settled_med"].min()),
        "argmax_theta": int(s.loc[s["out_share_settled_med"].idxmax(),
                                  "theta_deg"]),
        "mirror_pair_diffs": pairs}
R["out_share_skew"] = skew

# ---------------------------------------------------------------- comb --
comb = {}
for cell in ("COMB_d0.20", "COMB_d0.10"):
    lg = logdf[(logdf["cell"] == cell) & (logdf["arm"] == "LADDER")]
    ent = {"ind_frac_total_med": float(lg["ind_frac"].median())}
    for k, nm in enumerate(("G", "MID", "H")):
        ent[f"frac_{nm}"] = float((lg[f"eind_{k}"] / lg["E_inj"]).median())
    comb[cell] = ent
ratios = {}
for nm in ("G", "MID", "H"):
    ratios[nm] = comb["COMB_d0.20"][f"frac_{nm}"] \
        / comb["COMB_d0.10"][f"frac_{nm}"]
spread = max(ratios.values()) / min(ratios.values())
comb["per_component_node_over_ref"] = ratios
comb["spread_max_over_min"] = float(spread)
comb["verdict"] = ("SELECTIVE" if spread >= 3.0 else
                   "BROADBAND" if spread <= 1.5 else "PARTIAL")
# single-tone cross-checks (same ratios from the P1 sweeps)
st = {}
for pfx, nm in (("G", "G"), ("MID", "MID")):
    tb = {r["d"]: r["ind_frac_med"] for r in node[pfx]["table"]}
    st[nm] = tb[0.20] / tb[0.10] if 0.10 in tb else None
comb["single_tone_node_over_ref"] = st
R["comb"] = comb

# ---------------------------------------------------------------- gates --
R["gn1_drive_recovery"] = meta["g1"]
R["gn1_pass"] = all(v["hits"] / v["total"] >= 0.95
                    for v in meta["g1"].values())
R["gn1_cells_below"] = {k: v["hits"] / v["total"]
                        for k, v in meta["g1"].items()
                        if v["hits"] / v["total"] < 0.95}
R["tone_rung_crowding"] = {
    nm: {"tone": float(t), "nearest_rung": float(
        LAD[np.argmin(np.abs(LAD - t))]),
        "distance": float(np.min(np.abs(LAD - t)))}
    for nm, t in meta["tones"].items()}
R["gn2_max_closure_frac"] = float(runs["closure_frac"].max())
R["gn2_pass"] = bool((runs["closure_frac"] <= 0.10).all())
R["gn3_max_sat_frac"] = float(runs["sat_frac"].max())
R["gn3_total_clamps"] = int(runs["clamp_count"].sum())
R["gn3_max_global_rel_slope"] = float(
    np.nanmax([row["global_rel_slope_EA"] for row in log_rows]))
R["gn3_pass"] = bool(R["gn3_max_sat_frac"] <= 0.30
                     and R["gn3_total_clamps"] == 0
                     and R["gn3_max_global_rel_slope"] < 0.005)
# GN0 run-level: check pair (0.59, 90deg) vs (0.71, 90deg)
def _indfrac_stats(cell):
    g_ = logdf[(logdf["cell"] == cell) & (logdf["arm"] == "LADDER")]
    return float(g_["ind_frac"].median()), float(g_["ind_frac"].std(ddof=1))
m1, s1 = _indfrac_stats("G2D_x0.59_t090")
m2, s2 = _indfrac_stats("G2D_x0.71_t090")
pooled = float(np.sqrt((s1 ** 2 + s2 ** 2) / 2.0))
R["gn0"] = {"profile_min_rel_L2": meta["profile_check"][
                "min_same_theta_rel_L2"],
            "profile_check_pair_rel_L2": meta["profile_check"][
                "check_pair_0.59_vs_0.71_th90_rel_L2"],
            "run_pair": {"med_0.59_t090": m1, "med_0.71_t090": m2,
                         "pooled_sd": pooled,
                         "diff_over_3sd": float(abs(m1 - m2) / (3 * pooled))
                         if pooled > 0 else None}}
R["gn0_pass"] = bool(meta["profile_check"]["min_same_theta_rel_L2"] > 0.05
                     and pooled > 0 and abs(m1 - m2) > 3 * pooled)
# OFF replication vs iteration 4 (B share of output, LADDER)
def b_share(df_runs):
    gg = df_runs[df_runs["arm"] == "LADDER"]
    out = (gg["E_diss_B"] + gg["E_rad_B"]).median()
    tot = (gg["E_diss_A"] + gg["E_diss_B"] + gg["E_rad_A"]
           + gg["E_rad_B"]).median()
    return float(out / tot)
R["b_share_off_iter5"] = b_share(runs[runs["cell"] == "OFF"])
R["b_share_off_iter4"] = 0.08583019239973681
# transducer stability summary
R["stability"] = {
    "sat_frac_max": R["gn3_max_sat_frac"],
    "sat_frac_median": float(runs["sat_frac"].median()),
    "ext_frac_median": float((runs["E_tr_ext"] / runs["E_inj"]).median()),
    "pool_max_max": float(runs["pool_max"].max()),
    "clamps_total": R["gn3_total_clamps"]}

json.dump(R, open(f"{HD}/riemann_node2d_results.json", "w"), indent=1)
print("GN0 pass:", R["gn0_pass"], R["gn0"]["run_pair"])
print("GN1 pass:", R["gn1_pass"], "below:", R["gn1_cells_below"])
print("GN2 max closure:", R["gn2_max_closure_frac"], "pass:", R["gn2_pass"])
print("GN3 pass:", R["gn3_pass"], "max sat:", R["gn3_max_sat_frac"],
      "clamps:", R["gn3_total_clamps"],
      "max rel slope:", R["gn3_max_global_rel_slope"])
print("GN5 pass:", R["gn5_pass"], "unexplained:", R["gn5_unexplained_cluster_cells"])
print("Node:", {k: (v["d_star"], round(v["depth"], 3), round(v["width"], 3))
                  for k, v in node.items()})
print("Discriminator:", disc["verdict"])
print("Node shape:", shape["character"], "in-node angles:",
      shape["n_angles_in_node_row59"])
print("Comb:", comb["verdict"], ratios)
print("B share OFF iter5:", R["b_share_off_iter5"], "vs iter4:",
      R["b_share_off_iter4"])
