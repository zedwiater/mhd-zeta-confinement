"""
riemann_honing_analyze.py — ITERATION 7 per-round analysis.
Usage: python3 riemann_honing_analyze.py --round=N

Audits round N's peaks with audit6 (iter-6 construction; capacitor
targets from THAT round's kappa: {f_bare, f_ring} per
riemann_honing_rounds.json), appends them to
riemann_honing_scatter.csv, then evaluates the round against the
pre-registered anomaly criteria (a)-(d) and writes
riemann_honing_round_summary.json (one entry per round) +
prints a compact summary. Controls (OFF / PLACEBO) are resolved
at the focus inducer geometry from round N if present, else the
most recent earlier round with matching geometry/kappa/seeds.
"""

import itertools
import json
import sys
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

pk_all = pd.read_csv(f"{HD}/riemann_honing_peaks_raw.csv")
runs_all = pd.read_csv(f"{HD}/riemann_honing_runs.csv")
rounds_meta = json.load(open(f"{HD}/riemann_honing_rounds.json"))
F_CAPS_BY_MULT = {float(k): v for k, v in rounds_meta["f_caps"].items()}
iter6_runs = pd.read_csv(f"{HD}/riemann_capacitor_runs.csv")

SCATTER_F = f"{HD}/riemann_honing_scatter.csv"
SUMMARY_F = f"{HD}/riemann_honing_round_summary.json"


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


def cap_label(fc, F, tones, f_caps):
    def hit(t):
        return t > 0.5 and abs(fc - t) <= 0.15
    for fc_ in f_caps:
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


def audit6(fc, F, tones, f_caps):
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
        cl = cap_label(fc, F, tones, f_caps)
        if cl is not None:
            u15 = False
            lab = cl
    return u45, u30, u15, lab


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
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    for g in groups:
        sub = d.iloc[g]
        seeds = sorted(int(x) for x in sub["seed"].unique())
        lit = sorted(sub["pickup"].unique())
        out.append({
            "freq_median": float(sub["peak_freq"].median()),
            "amp_median": float(sub["peak_amp"].median()),
            "n_peaks": int(len(sub)), "seeds": seeds,
            "n_seeds": len(seeds), "lit": lit, "n_lit": len(lit),
            "eddy": bool(2 <= len(lit) <= 7),
            "u015x": bool(sub["u015x"].mean() > 0.5)})
    return out


def metrics_of(df):
    """Per-seed exchange metrics from a runs slice."""
    out = {}
    for seed, g in df.groupby("seed"):
        r = g.iloc[0]
        out[int(seed)] = {
            "ind_frac": r["E_ind"] / r["E_inj"],
            "B_share_out": (r["E_diss_B"] + r["E_rad_B"]) / max(
                r["E_diss_A"] + r["E_diss_B"] + r["E_rad_A"] + r["E_rad_B"],
                1e-12),
            "tr_ext_frac": r["E_tr_ext"] / r["E_inj"],
            "cap_end_share": r["E_cap_end"] / r["E_inj"],
            "sat_frac": r["sat_frac"],
            "closure_frac": r["closure_frac"],
            "slope_settled": r["slope_settled"]}
    return out


def medians(per_seed):
    if not per_seed:
        return {}
    keys = list(next(iter(per_seed.values())).keys())
    return {k: float(np.median([v[k] for v in per_seed.values()]))
            for k in keys}


# iter-6 baselines keyed by geometry (ind_x, ind_th, cap_x, cap_th)
def iter6_baselines():
    base = {}
    for cell, g in iter6_runs[iter6_runs["arm"] == "LADDER"].groupby("cell"):
        r0 = g.iloc[0]
        key = (round(float(r0["ind_x"]), 4), round(float(r0["ind_theta_deg"]), 2),
               None if pd.isna(r0["cap_x"]) else round(float(r0["cap_x"]), 4),
               None if pd.isna(r0["cap_theta_deg"])
               else round(float(r0["cap_theta_deg"]), 2))
        base[key] = {"cell": cell, "per_seed": metrics_of(g)}
    return base


BASE6 = iter6_baselines()


def geom_key(ind, cap):
    return (round(float(ind[0]), 4), round(float(ind[1]), 2),
            None if cap is None else round(float(cap[0]), 4),
            None if cap is None else round(float(cap[1]), 2))


def main():
    rnd = int([a.split("=", 1)[1] for a in sys.argv
               if a.startswith("--round=")][0])
    rmeta = {r["round"]: r for r in rounds_meta["rounds"]}[rnd]
    focus_name = rmeta["focus"]
    fcell = {c["name"]: c for c in rmeta["cells"]}[focus_name]
    focus_runs = runs_all[(runs_all["round"] == rnd)
                          & (runs_all["cell"] == focus_name)]

    # ------------------------------------------------ audit round peaks --
    pk = pk_all[pk_all["round"] == rnd].copy()
    cand_mask = (~pk["is_drive"]) & (~pk["is_inducer"])
    tones = [20.2147]
    au = []
    for r in pk[cand_mask].itertuples():
        f_caps = F_CAPS_BY_MULT[float(r.kappa_mult)]
        au.append(audit6(r.peak_freq, S.ARMS[r.arm], tones, f_caps))
    pk.loc[cand_mask, ["u045", "u030", "u015x"]] = np.array(
        [(a, b, c) for a, b, c, _l in au], dtype=bool)
    pk.loc[cand_mask, "audit6_label"] = np.array(
        [l for _a, _b, _c, l in au], dtype=object)
    pk.loc[~cand_mask, "audit6_label"] = np.where(
        pk.loc[~cand_mask, "is_inducer"], "inducer_tone", "drive")
    f = pk["peak_freq"].to_numpy()
    pk["d_line"] = np.min(np.abs(f[:, None] - LAD[None, :]), axis=1)
    pk["d_off"] = np.min(np.abs(f[:, None] - EIGA[None, :]), axis=1)
    pk["d_offB"] = np.min(np.abs(f[:, None] - EIGB[None, :]), axis=1)
    pk.to_csv(SCATTER_F, mode="a", header=not _exists(SCATTER_F),
              index=False)
    cand = pk[cand_mask].copy()

    # --------------------------------------- (a) unexplained clusters --
    unexp = {}
    for (cell, arm, win), g in cand.groupby(["cell", "arm", "window"]):
        cl = [c for c in chain_clusters(g)
              if c["n_seeds"] >= 2 and c["u015x"]]
        if cl:
            unexp[f"{cell}|{arm}|{win}"] = cl
    n_unexp_peaks = int(cand["u015x"].sum())

    # --------------------------------------- resolve controls ------------
    def find_runs(role, ind, cap, kappa_mult, seeds, upto):
        """Runs of a control role at geometry, latest round <= upto."""
        for rr in sorted(runs_all["round"].unique(), reverse=True):
            if rr > upto:
                continue
            g = runs_all[(runs_all["round"] == rr)
                         & (runs_all["role"] == role)]
            if g.empty:
                continue
            r0 = g.iloc[0]
            ck = (None if pd.isna(r0["cap_x"]) else float(r0["cap_x"]),
                  None if pd.isna(r0["cap_theta_deg"])
                  else float(r0["cap_theta_deg"]))
            want = (None, None) if cap is None else (float(cap[0]),
                                                     float(cap[1]))
            lad_seeds = sorted(int(x) for x in
                               g[g["arm"] == "LADDER"]["seed"].unique())
            if abs(float(r0["ind_x"]) - ind[0]) < 1e-6 \
                    and abs(float(r0["ind_theta_deg"]) - ind[1]) < 1e-6 \
                    and ck == want \
                    and abs(float(r0["kappa_mult"]) - kappa_mult) < 1e-9 \
                    and lad_seeds == sorted(seeds):
                return g, int(rr)
        return None, None

    f_ind = fcell["ind"]
    f_mult = float(fcell.get("kappa_mult", 1.0))
    f_seeds = sorted(int(x) for x in focus_runs["seed"].unique())
    off_runs, off_round = find_runs("off", f_ind, None, 1.0, f_seeds, rnd)
    pl_runs, pl_round = find_runs("placebo", f_ind, [0.13, 0], f_mult,
                                  f_seeds, rnd)

    scat_all = pd.read_csv(SCATTER_F)

    def settled_base12(cell_runs):
        if cell_runs is None:
            return None
        rids = set(cell_runs["run_id"])
        g = scat_all[scat_all["run_id"].isin(rids)]
        g = g[(~g["is_drive"]) & (~g["is_inducer"])
              & (g["window"] == "SETTLED") & (g["pickup"].isin(BASE12))]
        return g

    def dist_stats(g):
        if g is None or g.empty:
            return None
        return {"n": int(len(g)),
                "med_d_off": float(g["d_off"].median()),
                "med_d_line": float(g["d_line"].median())}

    def arm_stats(cell_runs, arm):
        if cell_runs is None:
            return None
        return dist_stats(settled_base12(cell_runs[cell_runs["arm"] == arm]))

    focus_L = arm_stats(focus_runs, "LADDER")
    off_L = arm_stats(off_runs, "LADDER")
    focus_D = arm_stats(focus_runs, "DETUNE")
    off_D = arm_stats(off_runs, "DETUNE")
    focus_S = arm_stats(focus_runs, "SHIFT")
    off_S = arm_stats(off_runs, "SHIFT")
    pl_L = arm_stats(pl_runs, "LADDER")

    def delta(a, b):
        return None if a is None or b is None else a["med_d_off"] - b["med_d_off"]

    d_focus = delta(focus_L, off_L)
    ctrls = {"focus_DETUNE": delta(focus_D, off_D),
             "focus_SHIFT": delta(focus_S, off_S),
             "placebo_LADDER": delta(pl_L, off_L)}
    spread = max([abs(v) for v in ctrls.values() if v is not None],
                 default=None)
    align_claim = bool(d_focus is not None and spread is not None
                       and abs(d_focus) > spread)
    d_line_focus = (None if focus_L is None or off_L is None
                    else focus_L["med_d_line"] - off_L["med_d_line"])

    # --------------------------------------- exchanges + (c) -------------
    focus_per_seed = metrics_of(focus_runs[focus_runs["arm"] == "LADDER"])
    focus_med = medians(focus_per_seed)
    key = geom_key(fcell["ind"], fcell.get("cap"))
    base = BASE6.get(key) if f_mult == 1.0 else None
    c_flags = []
    if base is not None:
        bp = base["per_seed"]
        for mk in ("ind_frac", "B_share_out", "tr_ext_frac",
                   "cap_end_share", "sat_frac"):
            bv = np.array([v[mk] for v in bp.values()], dtype=float)
            rv = np.array([v[mk] for v in focus_per_seed.values()],
                          dtype=float)
            pooled = float(np.sqrt((bv.std(ddof=1) ** 2
                                    + rv.std(ddof=1) ** 2) / 2.0)) \
                if len(bv) > 1 and len(rv) > 1 else 0.0
            dev = abs(float(np.median(rv)) - float(bv.mean()))
            if pooled > 0 and dev > 3.0 * pooled:
                c_flags.append({"metric": mk, "round_median":
                                float(np.median(rv)),
                                "baseline_mean": float(bv.mean()),
                                "pooled_sd": pooled, "dev_over_sd":
                                float(dev / pooled),
                                "baseline_cell": base["cell"]})
    # persistence continuity metric
    fr = focus_runs[focus_runs["arm"] == "LADDER"]
    ps = []
    for _seed, g in fr.groupby("seed"):
        r = g.iloc[0]
        de = r["demod_t0_A_in5_DRIVE_END_amp"]
        if de and de > 1e-9:
            ps.append(r["demod_t0_A_in5_SETTLED_amp"] / de)
    P_med = float(np.median(ps)) if ps else None

    # --------------------------------------- (d) stability ---------------
    rr = runs_all[runs_all["round"] == rnd]
    d_flags = []
    if rr["slope_settled"].max() > 0.005:
        d_flags.append(f"slope_settled max {rr['slope_settled'].max():.4f}")
    if int(rr["clamp_count"].sum()) > 0:
        d_flags.append(f"clamps {int(rr['clamp_count'].sum())}")
    if rr["closure_frac"].max() > 0.001:
        d_flags.append(f"closure max {rr['closure_frac'].max():.5f}")

    summary = {}
    if _exists(SUMMARY_F):
        summary = json.load(open(SUMMARY_F))
    entry = {
        "round": rnd, "focus": focus_name,
        "focus_geom": {"ind": fcell["ind"], "cap": fcell.get("cap"),
                       "kappa_mult": f_mult,
                       "f_bare": F_CAPS_BY_MULT[f_mult][0],
                       "f_ring": F_CAPS_BY_MULT[f_mult][1]},
        "seeds": f_seeds,
        "n_runs_round": int(len(rr)),
        "n_candidates": int(len(cand)),
        "n_unexplained_peaks": n_unexp_peaks,
        "exchanges_focus_median": focus_med,
        "cap_peak_share_median": float(fr["cap_peak_share"].median())
        if len(fr) else None,
        "P_median": P_med,
        "off_round_used": off_round, "placebo_round_used": pl_round,
        "alignment": {"delta_off_focus": d_focus,
                      "delta_line_focus": d_line_focus,
                      "controls": ctrls, "control_spread": spread,
                      "focus_stats": focus_L, "off_stats": off_L},
        "max_closure_frac": float(rr["closure_frac"].max()),
        "max_slope_settled": float(rr["slope_settled"].max()),
        "total_clamps": int(rr["clamp_count"].sum()),
        "anomaly": {"a_unexplained_clusters": unexp,
                    "b_alignment_claim": align_claim,
                    "c_exchange_flags": c_flags,
                    "d_stability_flags": d_flags},
        "candidate": bool(unexp or align_claim or c_flags or d_flags),
        "baseline_matched": base["cell"] if base else None,
    }
    summary[str(rnd)] = entry
    json.dump(summary, open(SUMMARY_F, "w"), indent=1, default=str)
    print(f"ROUND {rnd} focus={focus_name} runs={len(rr)} "
          f"cand={len(cand)} unexp_peaks={n_unexp_peaks}")
    print(f"  exchanges med: ind_frac={focus_med.get('ind_frac'):.5f} "
          f"B_share={focus_med.get('B_share_out'):.4f} "
          f"tr_ext={focus_med.get('tr_ext_frac'):.4f} "
          f"cap_peak={entry['cap_peak_share_median']:.2e} P={P_med}")
    print(f"  align: d_off={d_focus} d_line={d_line_focus} "
          f"spread={spread} claim={align_claim} "
          f"(off_r{off_round} pl_r{pl_round})")
    print(f"  ANOMALY a={list(unexp.keys())} b={align_claim} "
          f"c={[c['metric'] for c in c_flags]} d={d_flags} "
          f"-> candidate={entry['candidate']}")


def _exists(p):
    import os
    return os.path.exists(p)


if __name__ == "__main__":
    main()
