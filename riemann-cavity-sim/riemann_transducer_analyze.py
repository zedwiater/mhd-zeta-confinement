"""
riemann_transducer_analyze.py — iteration-3 verdicts.
Reads riemann_transducer_{pickups,transients,series,runs}.csv and the
iteration-2 products (pickups CSV + results JSON) as baselines.
Writes riemann_transducer_results.json and the final
riemann_transducer_log.csv (per run: site, arm, cap-binding, stability,
local peak-table deltas). Methods replicate iteration 2 exactly where a
comparison is made (leakage units fix, audit + flow extension at 0.15,
chain-clustering gap 0.45 stated here).
"""

import itertools
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import riemann_tp_cylinder_sim as S
import riemann_transducer_sim as TS

HD = S.HD
FP = 0.025
ARMS = S.ARMS
CONFIGS = ["OFF"] + TS.SITE_NAMES
SITE_OF = {nm: TS.SITES[nm] for nm in TS.SITE_NAMES}

pk = pd.read_csv(f"{HD}/riemann_transducer_pickups.csv")
trn = pd.read_csv(f"{HD}/riemann_transducer_transients.csv")
ser = pd.read_csv(f"{HD}/riemann_transducer_series.csv")
runs = pd.read_csv(f"{HD}/riemann_transducer_runs.csv")
meta = json.load(open(f"{HD}/riemann_transducer_raw_meta.json"))
pk2 = pd.read_csv(f"{HD}/riemann_tp_cylinder_pickups.csv")
res2 = json.load(open(f"{HD}/riemann_tp_cylinder_results.json"))

R = {"params": meta["TR"], "g1_by_config": meta["g1_by_config"]}

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

def audit3(fc, F, mod, rate):
    """(unexplained@0.45, unexplained@0.30, unexplained@0.15+flow)"""
    u45 = S.audit_label(fc, F, mod, rate, tol=0.45) == "UNEXPLAINED"
    u30 = S.audit_label(fc, F, mod, rate, tol=0.30) == "UNEXPLAINED"
    u15 = S.audit_label(fc, F, mod, rate, tol=0.15) == "UNEXPLAINED"
    if u15 and flow_label(fc, F) is not None:
        u15 = False
    return u45, u30, u15

cand = pk[~pk["is_drive"]].copy()
u = [audit3(r.peak_freq, ARMS[r.arm], r.modulation, r.mod_rate)
     for r in cand.itertuples()]
cand[["u045", "u030", "u015"]] = pd.DataFrame(u, index=cand.index)

# ------------------------------------------------------------- clusters --
def chain_clusters(df, tol=0.45):
    """Chain-cluster candidate peaks of one cell by frequency."""
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
            "u045": bool(sub["u045"].mean() > 0.5),
            "u030": bool(sub["u030"].mean() > 0.5),
            "u015": bool(sub["u015"].mean() > 0.5),
        })
    return out

clusters = {}
for (cfg, arm, mod, win), cell in cand.groupby(
        ["config", "arm", "modulation", "window"]):
    cl = [c for c in chain_clusters(cell) if c["n_seeds"] >= 2]
    clusters[f"{cfg}|{arm}|{mod}|{win}"] = cl
R["clusters_reproducible"] = {
    k: {"n": len(v), "eddy": int(sum(c["eddy"] for c in v)),
        "u045": int(sum(c["u045"] for c in v)),
        "u030": int(sum(c["u030"] for c in v)),
        "u015": int(sum(c["u015"] for c in v)),
        "list": v}
    for k, v in clusters.items()}

# neighbour eddy-lit fraction: over reproducible clusters of a cell, the
# fraction lighting each pickup; averaged over a config's neighbour set.
def neighbour_eddy(cfg, arm, mod):
    vals = {}
    for win in ("DRIVE_END", "SETTLED"):
        cl = clusters.get(f"{cfg}|{arm}|{mod}|{win}", [])
        nb = SITE_OF[cfg][3] if cfg in SITE_OF else []
        fr = []
        for p in nb:
            fr.append(np.mean([p in c["lit"] for c in cl]) if cl else np.nan)
        vals[win] = float(np.nanmean(fr)) if fr else None
    return vals

# ------------------------------------------------------- ledger / GT2 --
runs["closure_ok"] = runs["closure_frac"] <= 0.10
R["gt2_ledger"] = {
    "max_closure_frac": float(runs["closure_frac"].max()),
    "median_closure_frac": float(runs["closure_frac"].median()),
    "n_runs": int(len(runs)),
    "n_fail": int((~runs["closure_ok"]).sum()),
    "pool_balance_resid_max_abs": float(runs["pool_balance_resid"].abs().max()),
    "total_clamps": int(runs["clamp_count"].sum()),
}
R["gt1_sanity"] = {k: (v["hits"] / v["total"] if v["total"] else None)
                   for k, v in meta["g1_by_config"].items()}

# ------------------------------------------------------- stability GT3 --
def eloc_slope(g, site):
    st = g[(g["t"] >= 5.2) & (g["t"] <= 7.6)]
    y = np.log(np.maximum(st[f"eloc_{site}"].to_numpy(), 1e-300))
    if len(st) < 10 or np.all(st[f"eloc_{site}"] <= 0):
        return float("nan")
    return float(np.polyfit(st["t"].to_numpy(), y, 1)[0])

def late_ratio(g, site):
    a = g[(g["t"] >= 5.2) & (g["t"] < 6.0)][f"eloc_{site}"].mean()
    b = g[(g["t"] >= 7.2) & (g["t"] <= 8.0)][f"eloc_{site}"].mean()
    return float(b / a) if a > 0 else float("nan")

ser["key"] = (ser["config"] + "|" + ser["arm"] + "|" + ser["modulation"]
              + "|" + ser["seed"].astype(str))
run_slope, run_late = {}, {}
for key, g in ser.groupby("key"):
    cfg = g["config"].iloc[0]
    site = cfg if cfg in SITE_OF else None
    if site:
        run_slope[key] = eloc_slope(g, site)
        run_late[key] = late_ratio(g, site)
    else:  # OFF: slopes for every site (baseline)
        for nm in SITE_OF:
            run_slope[f"{key}@{nm}"] = eloc_slope(g, nm)
            run_late[f"{key}@{nm}"] = late_ratio(g, nm)

unstable = []
for r in runs.itertuples():
    if r.config == "OFF":
        continue
    key = f"{r.config}|{r.arm}|{r.modulation}|{r.seed}"
    off_slopes = [v for k, v in run_slope.items()
                  if k.startswith(f"OFF|{r.arm}|{r.modulation}|")
                  and k.endswith(f"@{r.config}")]
    sl = run_slope.get(key, float("nan"))
    lr = run_late.get(key, float("nan"))
    flags = []
    if off_slopes and sl == sl and sl > max(off_slopes) + 0.05:
        flags.append(f"slope {sl:.3f} > OFF max {max(off_slopes):.3f}+0.05")
    if lr == lr and lr > 1.2:
        flags.append(f"late growth ratio {lr:.2f} > 1.2")
    if not np.isfinite(r.E_tr_ext):
        flags.append("non-finite transducer energy")
    if flags:
        unstable.append({"run_id": int(r.run_id), "config": r.config,
                         "arm": r.arm, "modulation": r.modulation,
                         "seed": int(r.seed), "flags": flags,
                         "slope": sl, "late_ratio": lr})
R["gt3_unstable_runs"] = unstable

# ------------------------------------------- transducer metrics / site --
def cell_med(df, cfg, arm, mod, col):
    s = df[(df["config"] == cfg) & (df["arm"] == arm)
           & (df["modulation"] == mod)][col]
    return float(s.median()) if len(s) else None

def cell_range(df, cfg, arm, mod, col):
    s = df[(df["config"] == cfg) & (df["arm"] == arm)
           & (df["modulation"] == mod)][col]
    return [float(s.min()), float(s.max())] if len(s) else None

runs["ext_frac"] = runs["E_tr_ext"] / runs["E_inj"]
tr_metrics = {}
for cfg in CONFIGS:
    tr_metrics[cfg] = {}
    for arm in ARMS:
        for mod in ("NONE", "AM"):
            tr_metrics[cfg][f"{arm}|{mod}"] = {
                "ext_frac_med": cell_med(runs, cfg, arm, mod, "ext_frac"),
                "sat_frac_med": cell_med(runs, cfg, arm, mod, "sat_frac"),
                "gate_frac_med": cell_med(runs, cfg, arm, mod, "gate_frac"),
                "pool_max_med": cell_med(runs, cfg, arm, mod, "pool_max"),
                "reinj_net_med": cell_med(runs, cfg, arm, mod, "E_tr_reinj_net"),
                "closure_med": cell_med(runs, cfg, arm, mod, "closure_frac"),
            }
R["transducer_metrics"] = tr_metrics

# ------------------------------------------- local peak deltas at site --
def per_seed_stats(w, with_u015=False):
    ns, pws, us = [], [], []
    for _seed, g in w.groupby("seed"):
        ns.append(len(g))
        pws.append(float((g["peak_amp"] ** 2).sum()))
        if with_u015:
            us.append(int(g["u015"].sum()))
    if not ns:
        return {"n_med": 0.0, "pow_med": 0.0, "pow_range": None, "u015_med": 0.0}
    return {"n_med": float(np.median(ns)), "pow_med": float(np.median(pws)),
            "pow_range": [float(np.min(pws)), float(np.max(pws))],
            "u015_med": float(np.median(us)) if us else 0.0}

def site_peak_stats(cfg, arm, mod):
    pkname = SITE_OF[cfg][2] if cfg in SITE_OF else None
    if pkname is None:
        return None
    sub = cand[(cand["config"] == cfg) & (cand["arm"] == arm)
               & (cand["modulation"] == mod) & (cand["pickup"] == pkname)]
    return {win: per_seed_stats(sub[sub["window"] == win], with_u015=True)
            for win in ("DRIVE_END", "SETTLED")}

local_stats = {}
for cfg in CONFIGS:
    if cfg == "OFF":
        continue
    local_stats[cfg] = {f"{arm}|{mod}": site_peak_stats(cfg, arm, mod)
                        for arm in ARMS for mod in ("NONE", "AM")}
# OFF baseline at the same physical pickups
off_stats = {}
for cfg in SITE_OF:
    pkname = SITE_OF[cfg][2]
    sub = cand[(cand["config"] == "OFF") & (cand["pickup"] == pkname)]
    off_stats[cfg] = {}
    for arm in ARMS:
        for mod in ("NONE", "AM"):
            s2 = sub[(sub["arm"] == arm) & (sub["modulation"] == mod)]
            d = {}
            for win in ("DRIVE_END", "SETTLED"):
                d[win] = per_seed_stats(s2[s2["window"] == win])
            off_stats[cfg][f"{arm}|{mod}"] = d
R["local_peak_stats"] = local_stats
R["local_peak_stats_OFF_at_site_pickups"] = off_stats

# ------------------------------------------------------------- leakage --
def leakage_table(pkdf, arm, runsel):
    F = ARMS[arm]
    sub = pkdf[(pkdf["arm"] == arm) & (pkdf["modulation"] == "NONE")
               & (pkdf["window"] == "DRIVE_END")].copy()
    sub = sub[sub["run_id"].isin(runsel)]
    out_mask = sub["pickup_kind"] == "out"
    sub.loc[out_mask, "peak_amp"] = (
        sub.loc[out_mask, "peak_amp"]
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
            table[cav].append(round(float(np.median(ratios)), 4)
                              if ratios else None)
    return table

leak3 = {}
for cfg in CONFIGS:
    rids = set(runs[(runs["config"] == cfg) & (runs["modulation"] == "NONE")]["run_id"])
    leak3[cfg] = {arm: leakage_table(pk, arm, rids) for arm in ARMS}
R["leakage_by_config"] = leak3

def site_strong_components(cfg):
    cav, W = SITE_OF[cfg][0], SITE_OF[cfg][1]
    if cav == "A":
        xs = {"A79": 0.79, "A61": 0.61, "APORT": 0.97}[cfg]
        return cav, [n for n in range(1, 31)
                     if abs(np.sin(n * np.pi * xs)) >= 0.7]
    return cav, [i + 1 for i in range(30) if abs(W[i]) >= 0.7 * np.max(np.abs(W))]

leak_delta = {}
for cfg in SITE_OF:
    cav, comps = site_strong_components(cfg)
    base = res2["leakage_ratio_by_component"]["LADDER"][cav]
    cur = leak3[cfg]["LADDER"][cav]
    off = leak3["OFF"]["LADDER"][cav]
    ds = [cur[i - 1] - base[i - 1] for i in comps
          if cur[i - 1] is not None and base[i - 1] is not None]
    do = [off[i - 1] - base[i - 1] for i in comps
          if off[i - 1] is not None and base[i - 1] is not None]
    leak_delta[cfg] = {"cavity": cav, "components": comps,
                       "mean_delta_vs_iter2": float(np.mean(ds)) if ds else None,
                       "mean_delta_OFF_vs_iter2": float(np.mean(do)) if do else None,
                       "per_component_delta": {str(i): round(cur[i - 1] - base[i - 1], 4)
                                               for i in comps
                                               if cur[i - 1] is not None
                                               and base[i - 1] is not None}}
R["leakage_delta_site_strong"] = leak_delta

# ------------------------------------------------------------- pulling --
def pulling_proxy(pkdf, cfg):
    sub = pkdf[(pkdf["arm"] == "LADDER") & (pkdf["modulation"] == "NONE")
               & (pkdf["window"] == "SETTLED") & (pkdf["cavity"] == "A")
               & (pkdf["pickup_kind"] == "in")]
    if cfg is not None:
        sub = sub[sub["config"] == cfg]
    ratios = []
    for n in range(1, 8):
        fd = ARMS["LADDER"][n - 1]
        fe = fd + FP * n
        near = sub[np.abs(sub["peak_freq"] - fe) <= 1.5]
        if len(near):
            pkf = near.loc[near.groupby("run_id")["peak_amp"].idxmax()]
            off = (pkf["peak_freq"] - fd).median()
            ratios.append(float(off / (FP * n)))
    return {"per_component_ratio": ratios,
            "median_ratio": float(np.median(ratios)) if ratios else None}

R["pulling_proxy"] = {"iter2": pulling_proxy(pk2, None),
                      **{cfg: pulling_proxy(pk, cfg) for cfg in CONFIGS}}

# ------------------------------------------------------------ exchange --
def b_share(df):
    tot = df["E_diss_A"] + df["E_diss_B"] + df["E_rad_A"] + df["E_rad_B"]
    return (df["E_diss_B"] + df["E_rad_B"]) / tot

exch = {}
for cfg in CONFIGS:
    exch[cfg] = {}
    for arm in ARMS:
        s = runs[(runs["config"] == cfg) & (runs["arm"] == arm)
                 & (runs["modulation"] == "NONE")]
        exch[cfg][arm] = float(b_share(s).median()) if len(s) else None
R["exchange_B_share"] = exch
R["exchange_B_share_iter2"] = {"LADDER": 0.085, "DETUNE": 0.036, "SHIFT": 0.027}

# ------------------------------------------------- neighbour eddy frac --
R["neighbour_eddy_fraction"] = {
    cfg: {f"{arm}|{mod}": neighbour_eddy(cfg, arm, mod)
          for arm in ARMS for mod in ("NONE", "AM")}
    for cfg in CONFIGS if cfg != "OFF"}

# ------------------------------------------------------- final log CSV --
off_pow = {}
for cfg in SITE_OF:
    for arm in ARMS:
        for mod in ("NONE", "AM"):
            off_pow[(cfg, arm, mod)] = off_stats[cfg][f"{arm}|{mod}"]
log_rows = []
for r in runs.itertuples():
    row = {"run_id": int(r.run_id), "config": r.config,
           "site_cavity": SITE_OF[r.config][0] if r.config in SITE_OF else "",
           "site_pickup": SITE_OF[r.config][2] if r.config in SITE_OF else "",
           "arm": r.arm, "modulation": r.modulation, "seed": int(r.seed),
           "closure_frac": r.closure_frac,
           "E_ext": r.E_tr_ext, "ext_frac": r.E_tr_ext / r.E_inj,
           "E_reinj_net": r.E_tr_reinj_net, "E_tr_diss": r.E_tr_diss,
           "sat_frac": r.sat_frac, "gate_frac": r.gate_frac,
           "pool_max": r.pool_max, "pool_end": r.pool_end,
           "clamp_count": int(r.clamp_count),
           "pool_balance_resid": r.pool_balance_resid}
    if r.config in SITE_OF:
        key = f"{r.config}|{r.arm}|{r.modulation}|{r.seed}"
        st = local_stats[r.config][f"{r.arm}|{r.modulation}"]
        sub = cand[(cand["config"] == r.config) & (cand["arm"] == r.arm)
                   & (cand["modulation"] == r.modulation)
                   & (cand["seed"] == r.seed)
                   & (cand["pickup"] == SITE_OF[r.config][2])]
        row["site_cand_n_DE"] = int((sub["window"] == "DRIVE_END").sum())
        row["site_cand_pow_DE"] = float(
            (sub[sub["window"] == "DRIVE_END"]["peak_amp"] ** 2).sum())
        row["site_cand_pow_DE_delta_vs_OFFmed"] = (
            row["site_cand_pow_DE"]
            - off_pow[(r.config, r.arm, r.modulation)]["DRIVE_END"]["pow_med"])
        row["site_cand_n_SETTLED"] = int((sub["window"] == "SETTLED").sum())
        row["site_u015_DE_SETTLED"] = int(
            sub[sub["window"].isin(["DRIVE_END", "SETTLED"])]["u015"].sum())
        row["eloc_slope_settled"] = run_slope.get(key)
        row["eloc_late_ratio"] = run_late.get(key)
        row["unstable"] = any(u["run_id"] == int(r.run_id) for u in unstable)
    log_rows.append(row)
logdf = pd.DataFrame(log_rows)
logdf.to_csv(f"{HD}/riemann_transducer_log.csv", index=False)

# ------------------------------------- v2 additions (analysis fix) --
# GT4 placebo arm read at the SAME pickup as the site measurement
# (v1 compared the placebo config at its own B pickup — cross-cavity
# power scales differ ~1000x, which masked real effects; the gate text
# says "the PLACEBO seed range" for the effect, and the faithful reading
# is the placebo ARM's range at the same measurement point).
def pickup_stats_for_config(cfgname, pickup, arm, mod):
    sub = cand[(cand["config"] == cfgname) & (cand["arm"] == arm)
               & (cand["modulation"] == mod) & (cand["pickup"] == pickup)]
    return {win: per_seed_stats(sub[sub["window"] == win], with_u015=True)
            for win in ("DRIVE_END", "SETTLED")}

placebo_at_site = {}
for cfg in SITE_OF:
    placebo_at_site[cfg] = {
        f"{arm}|{mod}": pickup_stats_for_config(
            "PLACEBO", SITE_OF[cfg][2], arm, mod)
        for arm in ARMS for mod in ("NONE", "AM")}
R["placebo_at_site_pickups"] = placebo_at_site

# per-pickup SETTLED candidate-power ratio config/OFF (LADDER NONE):
# the spatial profile of any suppression — is it local to the site?
supp = {}
for cfg in CONFIGS:
    supp[cfg] = {}
    for pickup in S.CH_NAMES:
        def med_pow(df):
            vals = []
            for _seed, g in df.groupby("seed"):
                vals.append(float((g["peak_amp"] ** 2).sum()))
            return float(np.median(vals)) if vals else 0.0
        cur = med_pow(cand[(cand["config"] == cfg) & (cand["arm"] == "LADDER")
                           & (cand["modulation"] == "NONE")
                           & (cand["window"] == "SETTLED")
                           & (cand["pickup"] == pickup)])
        base = med_pow(cand[(cand["config"] == "OFF") & (cand["arm"] == "LADDER")
                            & (cand["modulation"] == "NONE")
                            & (cand["window"] == "SETTLED")
                            & (cand["pickup"] == pickup)])
        supp[cfg][pickup] = (cur / base) if base > 0 else None
R["settled_cand_power_ratio_vs_OFF_LADDER_NONE"] = supp

# cluster counts per cell (context for eddy fractions)
R["cluster_counts"] = {
    k: {"n": v["n"], "eddy": v["eddy"], "u045": v["u045"]}
    for k, v in R["clusters_reproducible"].items()}

# neighbour eddy fraction for the PLACEBO arm over each site's neighbour
# set (same-measurement-point placebo for the neighbour metric)
def neighbour_frac_for_arm(cfgname, neighbours, arm, mod, win):
    cl = clusters.get(f"{cfgname}|{arm}|{mod}|{win}", [])
    if not cl:
        return None
    return float(np.mean([np.mean([p in c["lit"] for c in cl])
                          for p in neighbours]))
R["neighbour_eddy_placebo_arm"] = {
    cfg: {f"{arm}|{mod}": {win: neighbour_frac_for_arm(
        "PLACEBO", SITE_OF[cfg][3], arm, mod, win)
        for win in ("DRIVE_END", "SETTLED")}
        for arm in ARMS for mod in ("NONE", "AM")}
    for cfg in SITE_OF}

# leakage attributable delta = config delta - OFF delta (seed-set control)
for cfg in SITE_OF:
    d = leak_delta[cfg]
    d["attributable_delta"] = (d["mean_delta_vs_iter2"]
                               - d["mean_delta_OFF_vs_iter2"])
R["leakage_delta_site_strong"] = leak_delta
R["pulling_per_component"] = {
    k: v["per_component_ratio"] for k, v in R["pulling_proxy"].items()}

# ------------------------------------------------------------- verdicts --
def gt4(val, off_range, pla_range):
    if val is None or off_range is None or pla_range is None:
        return False
    return val > max(off_range[1], pla_range[1]) or val < min(off_range[0], pla_range[0])

verdicts = {}
for cfg in SITE_OF:
    cell = "LADDER|NONE"
    eff = []
    # site pickup candidate power, DRIVE_END and SETTLED — placebo arm
    # read at the SAME pickup (v2 fix, see note above)
    for win in ("DRIVE_END", "SETTLED"):
        v = local_stats[cfg][cell][win]["pow_med"]
        o = off_stats[cfg][cell][win]["pow_range"]
        p = (placebo_at_site[cfg][cell][win]["pow_range"]
             if cfg != "PLACEBO" else None)
        if p is not None and gt4(v, o, p):
            eff.append(f"site_pickup_cand_power_{win}")
    # neighbour eddy fraction (SETTLED) vs OFF arm and PLACEBO arm,
    # both read over this site's neighbour set
    nv = R["neighbour_eddy_fraction"][cfg][cell]["SETTLED"]
    nb = SITE_OF[cfg][3]
    no = neighbour_frac_for_arm("OFF", nb, "LADDER", "NONE", "SETTLED")
    npla = R["neighbour_eddy_placebo_arm"][cfg][cell]["SETTLED"] \
        if cfg != "PLACEBO" else None
    if nv is not None and no is not None and abs(nv - no) > 0.15 \
            and (npla is None or abs(nv - npla) > 0.15):
        eff.append(f"neighbour_eddy_fraction_SETTLED ({nv:.2f} vs OFF "
                   f"{no:.2f}, placebo-arm {npla})")
    cl_cfg = R["clusters_reproducible"].get(f"{cfg}|LADDER|NONE|DRIVE_END",
                                            {"u015": None})
    verdicts[cfg] = {
        "unstable_runs": [u for u in unstable if u["config"] == cfg],
        "gt4_effects_LADDER_NONE": eff,
        "driveend_clusters_LADDER_NONE": cl_cfg,
        "ext_frac_LADDER_NONE": tr_metrics[cfg][cell]["ext_frac_med"],
        "leakage_delta_site_strong": leak_delta[cfg]["mean_delta_vs_iter2"],
        "pulling_median_ratio": R["pulling_proxy"][cfg]["median_ratio"],
    }
R["verdicts"] = verdicts
R["neighbour_eddy_OFF_baseline_note"] = (
    "OFF neighbour fractions are computed per config neighbour set in the "
    "verdict assembly; see RIEMANN_TRANSDUCER.md")

json.dump(R, open(f"{HD}/riemann_transducer_results.json", "w"), indent=1,
          default=str)
print("GT1:", R["gt1_sanity"])
print("GT2:", R["gt2_ledger"])
print("unstable runs:", len(unstable))
for cfg, v in verdicts.items():
    print(cfg, "| unstable:", len(v["unstable_runs"]),
          "| GT4 effects:", v["gt4_effects_LADDER_NONE"],
          "| ext_frac:", v["ext_frac_LADDER_NONE"],
          "| DE clusters u015:", v["driveend_clusters_LADDER_NONE"].get("u015"))
print("exchange B share:", R["exchange_B_share"])
print("pulling medians:", {k: v["median_ratio"]
                            for k, v in R["pulling_proxy"].items()})
print("leakage deltas:", {k: v["mean_delta_vs_iter2"]
                          for k, v in leak_delta.items()})
