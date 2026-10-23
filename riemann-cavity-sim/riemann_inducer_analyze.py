"""
riemann_inducer_analyze.py — iteration-4 verdicts.
Reads riemann_inducer_{peaks_raw,series,runs}.csv + raw meta, and the
iteration-3 OFF-config peaks (bit-identical system, same seeds) as the
inducer-OFF reference for the DETUNE/SHIFT arms (iteration 4's OFF
cell is LADDER-only; the restriction of all Delta contrasts to the 12
base channels is stated in the report).

Writes:
- riemann_inducer_scatter.csv  (every peak + extended audit + d_line /
  d_off / d_offB distances — the full-system scatter)
- riemann_inducer_log.csv      (per run: cell, arm, cap-binding
  fractions, stability metrics, local peak-table deltas)
- riemann_inducer_results.json (all verdict numbers)
Audit: iteration-3 audit3 verbatim (base at 0.45/0.30/0.15 + flow
extension at 0.15), extended with inducer targets at 0.15
(tone, 2tone, tone/2, tone+-f_i, 2tone+-f_i, tone+-k*s_i) — anything
so explained is engine-expected, never unique.
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

pk = pd.read_csv(f"{HD}/riemann_inducer_peaks_raw.csv")
runs = pd.read_csv(f"{HD}/riemann_inducer_runs.csv")
ser = pd.read_csv(f"{HD}/riemann_inducer_series.csv")
meta = json.load(open(f"{HD}/riemann_inducer_raw_meta.json"))
pk3 = pd.read_csv(f"{HD}/riemann_transducer_pickups.csv")
TONES = {c: (v["tone"] if v else 0.0) for c, v in meta["cells"].items()}


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


def ind_label(fc, F, tone):
    if tone <= 0:
        return None
    s = FP * np.arange(1, len(F) + 1)
    def hit(t):
        return t > 0.5 and abs(fc - t) <= 0.15
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
    return None


def audit4(fc, F, tone):
    """(u045, u030, u015_extended, label4)."""
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
        il = ind_label(fc, F, tone)
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
au = [audit4(r.peak_freq, S.ARMS[r.arm], TONES[r.cell])
      for r in pk[cand_mask].itertuples()]
pk.loc[cand_mask, ["u045", "u030", "u015x"]] = np.array(
    [(a, b, c) for a, b, c, _l in au], dtype=bool)
pk.loc[cand_mask, "audit4_label"] = np.array(
    [l for _a, _b, _c, l in au], dtype=object)
pk.loc[~cand_mask, "audit4_label"] = np.where(
    pk.loc[~cand_mask, "is_inducer"], "inducer_tone", "drive")
pk = add_distances(pk)
pk.to_csv(f"{HD}/riemann_inducer_scatter.csv", index=False)
cand = pk[cand_mask].copy()

# iter-3 reference = A79 config (transducer at A79, NO inducer — the
# exact iteration-4 base system; iter-3's OFF config had no transducer
# at all and is NOT the right reference. Caught by the replication
# check: iter-4 OFF B-share equals iter-3 A79's to all printed digits.)
ref3 = pk3[(pk3["config"] == "A79") & (pk3["modulation"] == "NONE")
           & (~pk3["is_drive"])].copy()
ref3 = add_distances(ref3)


def dist_stats(df):
    if df.empty:
        return {"n": 0, "med_d_line": None, "med_d_off": None,
                "frac_line_015": None, "frac_off_015": None}
    return {"n": int(len(df)),
            "med_d_line": float(df["d_line"].median()),
            "med_d_off": float(df["d_off"].median()),
            "frac_line_015": float((df["d_line"] <= 0.15).mean()),
            "frac_off_015": float((df["d_off"] <= 0.15).mean())}


R = {"tones": meta["tones"], "AIND": meta["AIND"]}
# scatter distributions per cell|arm|window (all 28 channels) + base-12
dist_all, dist_base = {}, {}
for (cell, arm, win), cell_df in cand.groupby(["cell", "arm", "window"]):
    dist_all[f"{cell}|{arm}|{win}"] = dist_stats(cell_df)
    dist_base[f"{cell}|{arm}|{win}"] = dist_stats(
        cell_df[cell_df["pickup"].isin(BASE12)])
R["dist_all_channels"] = dist_all
R["dist_base12"] = dist_base
# per-seed medians (SETTLED, base-12) for SD estimates
seed_meds = {}
for (cell, arm), cell_df in cand[cand["window"] == "SETTLED"].groupby(
        ["cell", "arm"]):
    b = cell_df[cell_df["pickup"].isin(BASE12)]
    seed_meds[f"{cell}|{arm}"] = {
        str(sd): float(g["d_off"].median())
        for sd, g in b.groupby("seed")}
    seed_meds[f"{cell}|{arm}_line"] = {
        str(sd): float(g["d_line"].median())
        for sd, g in b.groupby("seed")}
R["settled_per_seed_med_d_off"] = seed_meds
# iter-3 reference medians per arm (SETTLED, base-12)
ref_med = {}
for arm, g in ref3[ref3["window"] == "SETTLED"].groupby("arm"):
    ref_med[arm] = {"med_d_off": float(g["d_off"].median()),
                    "med_d_line": float(g["d_line"].median()),
                    "per_seed_off": {str(sd): float(x["d_off"].median())
                                     for sd, x in g.groupby("seed")},
                    "n": int(len(g))}
R["iter3_off_reference"] = ref_med

# alignment deltas (primary: d_off, SETTLED, base-12)
def pooled_sd(d):
    v = np.array(list(d.values()), dtype=float)
    return float(v.std(ddof=1)) if len(v) > 1 else 0.0

align = {}
for cell in meta["cells"]:
    for arm in ("LADDER", "DETUNE", "SHIFT"):
        k = f"{cell}|{arm}|SETTLED"
        if k not in dist_base or dist_base[k]["med_d_off"] is None:
            continue
        d_off = dist_base[k]["med_d_off"] - ref_med[arm]["med_d_off"]
        d_line = dist_base[k]["med_d_line"] - ref_med[arm]["med_d_line"]
        align[f"{cell}|{arm}"] = {"delta_off": d_off, "delta_line": d_line,
                                  "med_d_off": dist_base[k]["med_d_off"],
                                  "n": dist_base[k]["n"]}
R["alignment"] = align
# GI4 verdict for the sweep cells at d_ref + all LADDER cells vs controls
gi4 = {}
ref_sd = pooled_sd({**seed_meds.get("SWEEP_G_d0.10|LADDER", {}),
                    **ref_med["LADDER"]["per_seed_off"]})
for cell in meta["cells"]:
    kl, kd, ks = (f"{cell}|LADDER", f"{cell}|DETUNE", f"{cell}|SHIFT")
    if kl in align and kd in align and ks in align:
        gi4[cell] = {
            "aligns": bool(align[kl]["delta_off"] < align[kd]["delta_off"] - ref_sd
                           and align[kl]["delta_off"] < align[ks]["delta_off"] - ref_sd),
            "delta_off_LADDER": align[kl]["delta_off"],
            "delta_off_DETUNE": align[kd]["delta_off"],
            "delta_off_SHIFT": align[ks]["delta_off"]}
R["gi4_alignment"] = gi4
R["alignment_sd_ref"] = ref_sd

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
R["gi5_unexplained_cluster_cells"] = {
    k: len(v) for k, v in unexp_clusters.items() if v}

# ------------------------------------------------------ stability / log --
def settled_slope(run_id):
    g = ser[(ser["run_id"] == run_id)]
    w = g[(g["t"] >= S.WINDOWS["SETTLED"][0]) & (g["t"] < S.WINDOWS["SETTLED"][1])]
    e = w["eloc_A79"].to_numpy()
    t = w["t"].to_numpy()
    ok = e > 1e-12
    if ok.sum() < 10:
        return float("nan")
    return float(np.polyfit(t[ok], np.log(e[ok]), 1)[0])


# OFF (iter-4) SETTLED candidate power at pickups, for local deltas
def cell_power(cell, arm, pickup, win="SETTLED"):
    g = cand[(cand["cell"] == cell) & (cand["arm"] == arm)
             & (cand["pickup"] == pickup) & (cand["window"] == win)]
    return float((g["peak_amp"] ** 2).sum()) if len(g) else 0.0


off_pow = {(arm, pk_): cell_power("OFF", arm, pk_)
           for arm in ("LADDER",) for pk_ in ("A_in5", "A_rin0", "A_out")}
log_rows = []
for r in runs.itertuples():
    row = {"run_id": r.run_id, "cell": r.cell, "arm": r.arm,
           "seed": r.seed,
           "tone": TONES[r.cell],
           "E_inj": r.E_inj, "E_ind": r.E_ind,
           "ind_frac": r.E_ind / r.E_inj if r.E_inj > 0 else float("nan"),
           "closure_frac": r.closure_frac,
           "ext_frac": r.E_tr_ext / r.E_inj if r.E_inj > 0 else float("nan"),
           "sat_frac": r.sat_frac, "gate_frac": r.gate_frac,
           "pool_max": r.pool_max, "clamp_count": r.clamp_count,
           "settled_slope_A79": settled_slope(r.run_id),
           "tau_A": getattr(r, "tau_A", float("nan")),
           "tau_B": getattr(r, "tau_B", float("nan"))}
    for pk_ in ("A_in5", "A_rin0", "A_out"):
        p = cell_power(r.cell, r.arm, pk_) if r.arm == "LADDER" else float("nan")
        row[f"settled_pow_{pk_}"] = p
        row[f"settled_pow_ratio_{pk_}_vsOFF"] = (
            p / off_pow[("LADDER", pk_)]
            if r.arm == "LADDER" and off_pow[("LADDER", pk_)] > 0 else float("nan"))
    log_rows.append(row)
pd.DataFrame(log_rows).to_csv(f"{HD}/riemann_inducer_log.csv", index=False)
R["log_rows"] = log_rows

# ------------------------------------------------------------- purity --
pur = {}
for cell, v in meta["cells"].items():
    if v is None:
        continue
    g = runs[(runs["cell"] == cell) & (runs["arm"] == "LADDER")]
    ent = {"tone": v["tone"]}
    for pk_ in ("rin0", "A_in5", "A_out"):
        for win in ("EARLY", "DRIVE_END", "SETTLED"):
            amps = g[f"demod_{pk_}_{win}_amp"].to_numpy(dtype=float)
            phs = g[f"demod_{pk_}_{win}_phase"].to_numpy(dtype=float)
            ent[f"demod_{pk_}_{win}_amp_med"] = float(np.median(amps))
            ent[f"demod_{pk_}_{win}_phase_med"] = float(np.median(phs))
        tf = g["tpeak_rin0_DRIVE_END_freq"].replace("", np.nan).to_numpy(dtype=float)
        ent["tpeak_rin0_DRIVE_END_dev_med"] = (
            float(np.nanmedian(np.abs(tf - v["tone"])))
            if np.isfinite(tf).any() else None)
        tf = g["tpeak_rin0_SETTLED_freq"].replace("", np.nan).to_numpy(dtype=float)
        ent["tpeak_rin0_SETTLED_dev_med"] = (
            float(np.nanmedian(np.abs(tf - v["tone"])))
            if np.isfinite(tf).any() else None)
    # demod floor proxy: median peak amp at A_rin0 DRIVE_END per run
    fl = []
    for rid in g["run_id"]:
        gp = pk[(pk["run_id"] == rid) & (pk["pickup"] == "A_rin0")
                & (pk["window"] == "DRIVE_END")]
        fl.append(float(gp["peak_amp"].median()) if len(gp) else float("nan"))
    ent["rin0_driveend_peakamp_median"] = float(np.nanmedian(fl))
    pur[cell] = ent
R["tone_purity"] = pur

# GI3 at reference distance (SWEEP_G_d0.10, LADDER, tone G)
gref = runs[(runs["cell"] == "SWEEP_G_d0.10") & (runs["arm"] == "LADDER")]
gi3_seeds = []
for r in gref.itertuples():
    gp = pk[(pk["run_id"] == r.run_id) & (pk["pickup"] == "A_rin0")
            & (pk["window"] == "DRIVE_END")]
    floor = float(gp["peak_amp"].median()) if len(gp) else float("nan")
    amp = getattr(r, "demod_rin0_DRIVE_END_amp")
    ph_de = getattr(r, "demod_rin0_DRIVE_END_phase")
    ph_st = getattr(r, "demod_rin0_SETTLED_phase")
    drift = abs((ph_st - ph_de + 180.0) % 360.0 - 180.0)
    gi3_seeds.append({"seed": int(r.seed), "demod_amp": float(amp),
                      "floor_proxy": floor,
                      "ratio": float(amp / floor) if floor == floor and floor > 0 else None,
                      "phase_drift_deg": float(drift),
                      "pass": bool(amp >= 5 * floor and drift < 30.0)})
R["gi3_tone_hold"] = {"per_seed": gi3_seeds,
                      "pass": bool(sum(s["pass"] for s in gi3_seeds) >= 2)}

# --------------------------------------------------------------- gates --
R["gi1_drive_recovery"] = meta["g1"]
R["gi1_pass"] = all(v["hits"] / v["total"] >= 0.95
                    for v in meta["g1"].values())
R["gi2_max_closure_frac"] = float(runs["closure_frac"].max())
R["gi2_pass"] = bool((runs["closure_frac"] <= 0.10).all())
R["gi6_max_sat_frac"] = float(runs["sat_frac"].max())
R["gi6_total_clamps"] = int(runs["clamp_count"].sum())
slopes = [row["settled_slope_A79"] for row in log_rows]
R["gi6_settled_slopes_max"] = float(np.nanmax(slopes))
R["gi6_pass"] = bool(np.nanmax(slopes) < 0 and R["gi6_max_sat_frac"] <= 0.30)
# OFF replication vs iter-3 A79 exchange: B share of output, LADDER NONE
def b_share(df_runs):
    g = df_runs[(df_runs["arm"] == "LADDER")]
    out = (g["E_diss_B"] + g["E_rad_B"]).median()
    tot = (g["E_diss_A"] + g["E_diss_B"] + g["E_rad_A"] + g["E_rad_B"]).median()
    return float(out / tot)
R["b_share_off_iter4"] = b_share(runs[runs["cell"] == "OFF"])
runs3 = pd.read_csv(f"{HD}/riemann_transducer_runs.csv")
g3 = runs3[(runs3["config"] == "A79") & (runs3["arm"] == "LADDER")
           & (runs3["modulation"] == "NONE")]
R["b_share_a79_iter3"] = float(
    (g3["E_diss_B"] + g3["E_rad_B"]).median()
    / (g3["E_diss_A"] + g3["E_diss_B"] + g3["E_rad_A"] + g3["E_rad_B"]).median())

json.dump(R, open(f"{HD}/riemann_inducer_results.json", "w"), indent=1)
print("GI1 pass:", R["gi1_pass"], {k: f"{v['hits']}/{v['total']}" for k, v in meta["g1"].items()})
print("GI2 max closure:", R["gi2_max_closure_frac"], "pass:", R["gi2_pass"])
print("GI3 tone hold pass:", R["gi3_tone_hold"]["pass"], gi3_seeds)
print("GI5 unexplained cluster cells:", R["gi5_unexplained_cluster_cells"])
print("GI6 pass:", R["gi6_pass"], "max sat:", R["gi6_max_sat_frac"],
      "clamps:", R["gi6_total_clamps"], "max slope:", R["gi6_settled_slopes_max"])
print("GI4:", json.dumps(R["gi4_alignment"], indent=1))
print("B share OFF iter4:", R["b_share_off_iter4"], "vs iter3 A79:", R["b_share_a79_iter3"])
print("alignment deltas (SETTLED, base12):")
for k, v in sorted(R["alignment"].items()):
    print(f"  {k}: dOff={v['delta_off']:+.4f} dLine={v['delta_line']:+.4f} n={v['n']}")
