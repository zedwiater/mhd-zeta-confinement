"""
riemann_reseed3d_analyze.py — iteration-8 verdicts.
Reads riemann_reseed3d_{peaks_raw,runs}.csv + ring json + rd npz.
Pre-registration: riemann_reseed3d_prereg.json (gates, A_max rule,
consistency definitions, resolution ladder fixed before the runs).
Writes riemann_reseed3d_{scatter,log}.csv + results.json.
"""

import glob
import itertools
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal.windows import blackmanharris

import riemann_tp_cylinder_sim as S

HD = S.HD
FP = 0.025
LAD = S.T_LADDER
EIGA = S.F_EIG_A
EIGB = S.F_EIG_B
FS = 1.0 / (S.DT * S.REC_EVERY)
TONE_G = 20.2147
BASE12Z = ([f"A_in{i+1}@z0.50" for i in range(5)] + ["A_out@z0.50"]
           + [f"B_in{i+1}" for i in range(5)] + ["B_out"])

ring = json.load(open(f"{HD}/riemann_reseed3d_ring.json"))
F_CAPS = sorted({10.10735, ring["f_ring_3d"]})
runs = pd.read_csv(f"{HD}/riemann_reseed3d_runs.csv")
pk = pd.read_csv(f"{HD}/riemann_reseed3d_peaks_raw.csv")
capres = json.load(open(f"{HD}/riemann_capacitor_results.json"))

MOD_OF_RUN = dict(zip(runs["run_id"], runs["mod"]))


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


def audit6(fc, F, tones, mod_name):
    mr = float(S.MEAN_GAP) if mod_name in ("AM", "FM") else 0.0
    u15 = S.audit_label(fc, F, mod_name, mr, tol=0.15) == "UNEXPLAINED"
    lab = "UNEXPLAINED"
    if not u15:
        lab = S.audit_label(fc, F, mod_name, mr, tol=0.15)
    if u15 and flow_label(fc, F) is not None:
        u15, lab = False, flow_label(fc, F)
    if u15:
        il = ind_label5(fc, F, tones)
        if il is not None:
            u15, lab = False, il
    if u15:
        cl = cap_label(fc, F, tones)
        if cl is not None:
            u15, lab = False, cl
    return u15, lab


cand_mask = (~pk["is_drive"]) & (~pk["is_inducer"])
au = [audit6(r.peak_freq, S.ARMS[r.arm], [TONE_G], MOD_OF_RUN[r.run_id])
      for r in pk[cand_mask].itertuples()]
pk.loc[cand_mask, "u015x"] = np.array([a for a, _l in au], dtype=bool)
pk.loc[cand_mask, "audit6_label"] = np.array([l for _a, l in au],
                                             dtype=object)
pk.loc[~cand_mask, "audit6_label"] = np.where(
    pk.loc[~cand_mask, "is_inducer"], "inducer_tone", "drive")
pk.loc[~cand_mask, "u015x"] = False
f = pk["peak_freq"].to_numpy()
pk["d_line"] = np.min(np.abs(f[:, None] - LAD[None, :]), axis=1)
pk["d_off"] = np.min(np.abs(f[:, None] - EIGA[None, :]), axis=1)
pk["d_offB"] = np.min(np.abs(f[:, None] - EIGB[None, :]), axis=1)
pk.to_csv(f"{HD}/riemann_reseed3d_scatter.csv", index=False)
cand = pk[cand_mask].copy()

# per-run settled distances (BASE12Z)
sett = cand[(cand["window"] == "SETTLED") & (cand["pickup"].isin(BASE12Z))]
per_run_dist = sett.groupby("run_id").agg(
    d_off_settled=("d_off", "median"),
    d_line_settled=("d_line", "median"),
    n_cand_settled=("d_off", "size")).reset_index()
runs = runs.merge(per_run_dist, on="run_id", how="left")

R = {"f_caps_audit": F_CAPS, "ring": ring,
     "iter6_refs": {"winner_ind_frac": 0.00685, "offcap_ind_frac": 0.016569,
                    "winner_P": capres.get("gc4_winner_P_median"),
                    "offcap_P": capres.get("gc4_offcap_P_median"),
                    "trap": capres.get("trap_with_capacitor")}}


def stats(df, cols):
    out = {}
    for c in cols:
        v = df[c].dropna().to_numpy(dtype=float) \
            if c in df.columns else np.array([])
        if len(v):
            out[c] = {"mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1))
                      if len(v) > 1 else 0.0,
                      "max_dev": float(np.max(np.abs(v - np.mean(v)))),
                      "n": int(len(v)),
                      "vals": [round(float(x), 6) for x in v]}
        else:
            out[c] = {"mean": float("nan"), "sd": float("nan"),
                      "max_dev": float("nan"), "n": 0, "vals": []}
    return out


# ------------------------------------------------------------------- G0 --
g0w = runs[runs["cell"] == "G0_WIN"]
g0o = runs[runs["cell"] == "G0_OFF"]
R["G0"] = {
    "winner_ind_frac_median": float(g0w["ind_frac"].median()),
    "winner_ind_frac_iter6": 0.00685,
    "winner_ratio": float(g0w["ind_frac"].median() / 0.00685),
    "offcap_ind_frac_median": float(g0o["ind_frac"].median()),
    "offcap_ind_frac_iter6": 0.016569,
    "offcap_ratio": float(g0o["ind_frac"].median() / 0.016569),
    "winner_P_median": float(g0w["P"].median()),
    "offcap_P_median": float(g0o["P"].median()),
    "max_closure": float(runs[runs["tag"] == "A"]["closure_frac"].max()),
}
R["G0"]["pass"] = bool(abs(R["G0"]["winner_ratio"] - 1) <= 0.02
                       and abs(R["G0"]["offcap_ratio"] - 1) <= 0.02
                       and abs(R["G0"]["winner_P_median"]
                               - R["iter6_refs"]["winner_P"]) <= 0.03)

# --------------------------------------------------------------- Phase A --
A = {}
uw = runs[runs["cell"] == "UNIF_WIN"]
A["zuniform_kz1"] = {
    "ind_frac_median": float(uw["ind_frac"].median()),
    "P_median": float(uw["P"].median()),
    "offset_vs_G0_ind_frac": float(uw["ind_frac"].median()
                                   / g0w["ind_frac"].median())}
capz = {}
for z in (0.1, 0.3, 0.5, 0.7, 0.9):
    g = runs[runs["cell"] == f"CAPZ_{z:.1f}"]
    capz[str(z)] = {"ind_frac_median": float(g["ind_frac"].median()),
                    "P_median": float(g["P"].median()),
                    "P_per_seed": [round(float(x), 4) for x in g["P"]],
                    "B_share_out_median": float(g["B_share_out"].median())}
A["cap_z_sweep"] = capz
A["cap_z_best"] = max(capz, key=lambda k: capz[k]["P_median"])
trap = {}
for dd in (0.10, 0.18, 0.24, 0.28, 0.40):
    g = runs[runs["cell"] == f"TRAP_d{dd:.2f}"]
    trap[str(dd)] = float(g["ind_frac"].median())
if trap:
    ends = np.mean([trap["0.1"], trap["0.4"]])
    dstar = min(trap, key=lambda k: trap[k])
    A["trap_z05"] = {"curve": trap, "d_star": float(dstar),
                     "depth": float(trap[dstar] / ends),
                     "plateau": float(ends)}
trapz = {}
for dd, zz in ((0.10, 0.25), (0.10, 0.75), (0.24, 0.25), (0.24, 0.75)):
    g = runs[runs["cell"] == f"TRAP_d{dd:.2f}_z{zz:.2f}"]
    trapz[f"d{dd}_z{zz}"] = float(g["ind_frac"].median())
A["trap_offplane"] = trapz
R["phaseA"] = A

# --------------------------------------------------------------- Phase B --
POST_METRICS = ["P", "demod_ST_amp", "B_share_post", "d_off_settled",
                "d_line_settled", "slope_settled"]
CTX_METRICS = ["ind_frac", "B_share_out"]
B = {"ensembles": {}}
win = runs[(runs["cell"] == "WIN3D") & (runs["arm"] == "LADDER")]
for a, g in win.groupby("a_reseed"):
    B["ensembles"][str(a)] = stats(g, POST_METRICS + CTX_METRICS)
base = B["ensembles"]["0.0"]
amax, amax_detail = 0.0, {}
for a in (0.1, 0.25, 0.5, 0.75, 1.0):
    e = B["ensembles"][str(a)]
    ok, det = True, {}
    for m in POST_METRICS:
        sd_ok = e[m]["sd"] <= 1.5 * max(base[m]["sd"], 1e-12)
        mean_ok = abs(e[m]["mean"] - base[m]["mean"]) \
            <= 1.5 * max(base[m]["sd"], 1e-12)
        det[m] = {"sd_ratio": float(e[m]["sd"] / max(base[m]["sd"], 1e-12)),
                  "mean_shift_in_base_sd":
                  float(abs(e[m]["mean"] - base[m]["mean"])
                        / max(base[m]["sd"], 1e-12)),
                  "pass": bool(sd_ok and mean_ok)}
        ok = ok and det[m]["pass"]
    amax_detail[str(a)] = det
    if ok:
        amax = a
B["A_max"] = amax
B["A_max_detail"] = amax_detail
ctrl = {}
for arm in ("DETUNE", "SHIFT"):
    for a in (0.0, 1.0):
        g = runs[(runs["cell"] == "WIN3D") & (runs["arm"] == arm)
                 & (runs["a_reseed"] == a)]
        ctrl[f"{arm}_A{a}"] = stats(g, POST_METRICS)
for a in (0.0, 1.0):
    g = runs[(runs["cell"] == "PLAC3D") & (runs["a_reseed"] == a)]
    ctrl[f"PLAC_A{a}"] = stats(g, POST_METRICS)
B["controls"] = ctrl
R["phaseB"] = B

# ------------------------------------------------------- dual-clock (B) --
def _cv(st):
    return float(st["sd"] / abs(st["mean"])) if st and st["mean"] != 0 else None


def paired_shifts(df, col):
    piv = df.pivot(index="seed", columns="a_reseed", values=col)
    out = {}
    if 0.0 not in piv.columns:
        return out
    base_sd = float(np.std(piv[0.0].to_numpy(dtype=float), ddof=1)) \
        if len(piv) > 1 else 0.0
    for a in piv.columns:
        if a == 0.0:
            continue
        d = (piv[a] - piv[0.0]).to_numpy(dtype=float)
        out[str(a)] = {
            "max_abs_paired_diff": float(np.max(np.abs(d))),
            "mean_shift": float(np.mean(d)),
            "mean_shift_in_base_sd": float(abs(np.mean(d))
                                           / max(base_sd, 1e-12))}
    return out


DC = {"clock1_metric": "demod_DE_amp (A_in5@z0.50, DRIVE_END window)",
      "clock2_metric": "demod_ST_amp (A_in5@z0.50, SETTLED window)",
      "per_A": {}}
for a, g in win.groupby("a_reseed"):
    st1 = stats(g, ["demod_DE_amp"])["demod_DE_amp"]
    st2 = stats(g, ["demod_ST_amp"])["demod_ST_amp"]
    DC["per_A"][str(a)] = {"clock1": st1, "clock2": st2,
                           "clock1_cv": _cv(st1), "clock2_cv": _cv(st2)}
DC["clock1_paired_shifts"] = paired_shifts(win, "demod_DE_amp")
DC["clock2_paired_shifts"] = paired_shifts(win, "demod_ST_amp")
_b0 = DC["per_A"]["0.0"]
DC["verdict_inputs"] = {
    "A0_cv_clock1": _b0["clock1_cv"], "A0_cv_clock2": _b0["clock2_cv"],
    "lower_cv_at_A0": "clock1" if (_b0["clock1_cv"] or 9) <= (_b0["clock2_cv"] or 9) else "clock2",
    "clock1_max_mean_shift_in_base_sd": max(
        (v["mean_shift_in_base_sd"] for v in DC["clock1_paired_shifts"].values()),
        default=0.0),
    "clock2_max_mean_shift_in_base_sd": max(
        (v["mean_shift_in_base_sd"] for v in DC["clock2_paired_shifts"].values()),
        default=0.0)}
B["dual_clock"] = DC
R["dual_clock"] = DC


def a1_absorption(df0, df1):
    """Phase-B rule evaluated at A=1 only: (pass, detail)."""
    if not len(df0) or not len(df1):
        return None, {}
    s0 = stats(df0, POST_METRICS)
    s1 = stats(df1, POST_METRICS)
    ok, det = True, {}
    for m in POST_METRICS:
        sd_ok = s1[m]["sd"] <= 1.5 * max(s0[m]["sd"], 1e-12)
        mean_ok = abs(s1[m]["mean"] - s0[m]["mean"]) \
            <= 1.5 * max(s0[m]["sd"], 1e-12)
        det[m] = bool(sd_ok and mean_ok)
        ok = ok and det[m]
    return bool(ok), det


def consistency_composite(df0, base_stats):
    ratios = {}
    for m in CMET_REF:
        if m in base_stats and m in stats(df0, [m]):
            ratios[m] = float(stats(df0, [m])[m]["sd"]
                              / max(base_stats[m]["sd"], 1e-12))
    return {"sd_ratios": ratios,
            "composite": float(np.median(list(ratios.values())))
            if ratios else None}


CMET_REF = ["ind_frac", "P", "B_share_out", "d_off_settled"]
WIN_A0_STATS = stats(win[win["a_reseed"] == 0.0], CMET_REF)

# --------------------------------------------------------------- Phase C --
CMET = ["ind_frac", "P", "B_share_out", "d_off_settled"]
C = {"cells": {"BASE": stats(win[win["a_reseed"] == 0.0], CMET)}}
for pname in ("STEP", "SLOWRAMP", "SHORT", "AM", "FM"):
    g = runs[(runs["cell"] == f"PROTO_{pname}") & (runs["tag"] == "C")]
    C["cells"][pname] = stats(g, CMET)
comp = {}
for pname, st in C["cells"].items():
    ratios = [st[m]["sd"] / max(C["cells"]["BASE"][m]["sd"], 1e-12)
              for m in CMET]
    comp[pname] = {"sd_ratios": {m: round(float(r), 4)
                                 for m, r in zip(CMET, ratios)},
                   "composite": float(np.nanmedian(ratios))}
C["composite"] = comp
C["winner"] = min(comp, key=lambda k: (np.inf if np.isnan(
        comp[k]["composite"]) else comp[k]["composite"]))
flags = {}
for pname, st in C["cells"].items():
    flags[pname] = {m: bool(st[m]["sd"] / abs(st[m]["mean"]) < 1e-3)
                    for m in CMET if st[m]["mean"] != 0}
C["seed_independence_flags"] = flags
R["phaseC"] = C
TD = {"dial1_need": "cross-seed consistency composite (Phase C/E tables)",
      "dial2_amplitude": {}}
TD["dial2_amplitude"]["BASE"] = {"A_max": amax, "source": "Phase B full sweep"}
wp_ = C["winner"]
if wp_ != "BASE":
    dw0 = runs[runs["cell"] == f"PROTO_{wp_}"]
    d0 = dw0[dw0["a_reseed"] == 0.0]
    d1 = dw0[dw0["a_reseed"] == 1.0]
    ok_, det_ = a1_absorption(d0, d1)
    TD["dial2_amplitude"][wp_] = {
        "A1_absorbed": ok_, "detail": det_,
        "A_max": 1.0 if ok_ else "<1 (unresolved; only A=1 tested)",
        "source": "Phase DW A=1 ensemble vs Phase C A=0 ensemble"}
for pname in ("STEP", "SLOWRAMP", "SHORT", "AM", "FM"):
    if pname not in TD["dial2_amplitude"]:
        TD["dial2_amplitude"][pname] = {
            "A_max": None,
            "source": "not a consistency finalist; amplitude not swept"}
R["two_dial_protocols"] = TD

# --------------------------------------------------------------- Phase D --
def win_peaks(x, fs=FS):
    x = x - np.mean(x)
    w = blackmanharris(len(x))
    spec = np.abs(np.fft.rfft(x * w)) * (2.0 / np.sum(w))
    freqs = np.fft.rfftfreq(len(x), d=1.0 / fs)
    med = np.median(spec) + 1e-15
    floor = max(8.0 * med, 0.002 * np.max(spec))
    idx, _ = find_peaks(spec, height=floor, prominence=4.0 * med,
                        distance=3)
    out = []
    for i in idx:
        if 0 < i < len(spec) - 1:
            a_, b_, c_ = spec[i - 1], spec[i], spec[i + 1]
            d = 0.5 * (a_ - c_) / (a_ - 2.0 * b_ + c_ + 1e-30)
        else:
            d = 0.0
        out.append((float(freqs[i] + d * (freqs[1] - freqs[0])),
                    float(spec[i])))
    return out


def load_rd():
    rd = {}
    for fn in glob.glob(f"{HD}/riemann_reseed3d_rd_*.npz"):
        z = np.load(fn)
        for k in z.files:
            rd[k] = z[k]
    return rd


def ens_peaks(rd, cell, proto, arm, a, ch="A_in5@z0.50"):
    """per-seed peak lists per resolution rung for one ensemble."""
    rungs = {0.9: [], 1.8: [], 3.6: []}
    phases_coh = []
    for seed in range(8):
        key = f"{cell}|{proto}|{arm}|A{a}|s{seed}|{ch}"
        if key not in rd:
            continue
        seg = rd[key].astype(float)          # t = 3.9 .. 7.9 s
        for L in rungs:
            n = int(round(L * FS))
            x = seg[-n:]
            rungs[L].append(win_peaks(x))
        # per-mode demod phase on the longest rung
        n = int(round(3.6 * FS))
        x = seg[-n:]
        tt = 4.3 + np.arange(n) / FS
        ph = []
        for fe in EIGA:
            zc = np.mean(x * np.exp(-2j * np.pi * fe * tt)) * 2.0
            ph.append(np.angle(zc))
        phases_coh.append(ph)
    return rungs, phases_coh


def line_test(rungs0, rungsA, df_tol_scale=1.0):
    """lines in >=6/8 of A-ensemble, <=1/8 of baseline, per rung."""
    found = {}
    for L, listsA in rungsA.items():
        lists0 = rungs0[L]
        if len(listsA) < 6 or not lists0:
            continue
        df = 1.0 / L * df_tol_scale
        fa = defaultdict(int)
        for lp in listsA:
            seen = set()
            for frq, _amp in lp:
                seen.add(round(frq / df))
            for b in seen:
                fa[b] += 1
        f0 = defaultdict(int)
        for lp in lists0:
            seen = set()
            for frq, _amp in lp:
                seen.add(round(frq / df))
            for b in seen:
                f0[b] += 1
        hits = []
        for b, cnt in fa.items():
            if cnt >= 6 and f0.get(b, 0) <= 1:
                fc = (b + 0.5) * df
                u15, lab = audit6(fc, S.T_LADDER, [TONE_G], "NONE")
                hits.append({"freq_bin_center": round(float(fc), 3),
                             "n_seeds_A": int(cnt),
                             "n_seeds_base": int(f0.get(b, 0)),
                             "audit6": lab, "u015x": bool(u15)})
        if hits:
            found[str(L)] = hits
    return found


def circ_sd(phases_list):
    if len(phases_list) < 2:
        return None
    P = np.array(phases_list)                # (seeds, modes)
    Rr = np.abs(np.mean(np.exp(1j * P), axis=0))
    sd = np.sqrt(np.maximum(-2.0 * np.log(np.maximum(Rr, 1e-12)), 0.0))
    return float(np.mean(sd))


rd = load_rd()
D = {"ensembles_found": sorted({k.rsplit("|", 1)[0] for k in rd})}
r0, ph0 = ens_peaks(rd, "WIN3D", "BASE", "LADDER", 0.0)
rA, phA = ens_peaks(rd, "WIN3D", "BASE", "LADDER", float(amax))
r1, ph1 = ens_peaks(rd, "WIN3D", "BASE", "LADDER", 1.0)
D["Amax_lines"] = line_test(r0, rA)
D["A1_lines"] = line_test(r0, r1)
D["phase_coherence_circsd"] = {"baseline": circ_sd(ph0),
                               "A_max": circ_sd(phA),
                               "A_1": circ_sd(ph1)}
wp = C["winner"]
if wp != "BASE":
    rw0, pwh0 = ens_peaks(rd, f"PROTO_{wp}", wp, "LADDER", 0.0)
    rw1, pwh1 = ens_peaks(rd, f"PROTO_{wp}", wp, "LADDER", 1.0)
    D["winner_proto_lines"] = line_test(rw0, rw1)
    D["winner_proto_circsd"] = {"baseline": circ_sd(pwh0),
                                "A_1": circ_sd(pwh1)}
else:
    D["winner_proto_lines"] = "winner is BASE; covered by A1 test"
# dual-clock on the Phase D ensembles (run-metric clocks; the
# resolution ladder itself is a settled-clock readout)
DD = {"note": "Clock1 = DRIVE_END demod (pre-reseed by construction: "
              "event fires at first step t>=T_OFF, DRIVE_END window "
              "ends at T_OFF). Clock2 = SETTLED demod. Ring-down "
              "resolution ladder = clock2 only.",
      "baseline_ensembles": DC["per_A"],
      "clock1_paired_shifts": DC["clock1_paired_shifts"],
      "clock2_paired_shifts": DC["clock2_paired_shifts"]}
if wp != "BASE":
    dw0 = runs[runs["cell"] == f"PROTO_{wp}"]
    DD["winner_protocol"] = {
        a: {"clock1_cv": _cv(stats(g, ["demod_DE_amp"])["demod_DE_amp"]),
            "clock2_cv": _cv(stats(g, ["demod_ST_amp"])["demod_ST_amp"])}
        for a, g in dw0.groupby("a_reseed")}
D["dual_clock"] = DD
R["phaseD"] = D

# --------------------------------------------------------------- Phase E --
E = {}
e1 = runs[(runs["cell"] == "E_INDTRAP") & (runs["arm"] == "LADDER")]
E1 = {"ensembles": {}}
for a, g in e1.groupby("a_reseed"):
    E1["ensembles"][str(a)] = stats(g, POST_METRICS + CTX_METRICS)
if "0.0" in E1["ensembles"]:
    b1 = E1["ensembles"]["0.0"]
    am1, det1 = 0.0, {}
    for a in (0.1, 0.25, 0.5, 0.75, 1.0):
        if str(a) not in E1["ensembles"]:
            continue
        e = E1["ensembles"][str(a)]
        ok = True
        for m in POST_METRICS:
            ok = ok and e[m]["sd"] <= 1.5 * max(b1[m]["sd"], 1e-12) \
                and abs(e[m]["mean"] - b1[m]["mean"]) \
                <= 1.5 * max(b1[m]["sd"], 1e-12)
        det1[str(a)] = bool(ok)
        if ok:
            am1 = a
    E1["A_max_adjusted_inducer"] = am1
    E1["A_max_detail_pass"] = det1
    E1["A_max_WIN3D"] = amax
E["E1_indtrap"] = E1

tc = runs[runs["cell"] == "E_TWOCAP"]
E["E2_twocap"] = {
    "A0": stats(tc[tc["a_reseed"] == 0.0], POST_METRICS + CTX_METRICS),
    "A1": stats(tc[tc["a_reseed"] == 1.0], POST_METRICS + CTX_METRICS),
    "closure_max": float(tc["closure_frac"].max()) if len(tc) else None,
    "E_cap1_end_median": float(tc["E_cap1_end"].median()) if len(tc) else None,
    "E_cap2_end_median": float(tc["E_cap2_end"].median()) if len(tc) else None,
    "E_cap_end_median": float(tc["E_cap_end"].median()) if len(tc) else None}

an = runs[runs["cell"] == "E_ANTI"]
E["E3_anti"] = {
    "LADDER_A0": stats(an[(an["arm"] == "LADDER") & (an["a_reseed"] == 0.0)],
                       POST_METRICS + CTX_METRICS),
    "LADDER_A1": stats(an[(an["arm"] == "LADDER") & (an["a_reseed"] == 1.0)],
                       POST_METRICS + CTX_METRICS),
    "arm_d_off_medians": {
        arm: float(an[(an["arm"] == arm) & (an["a_reseed"] == 0.0)]
                   ["d_off_settled"].median())
        for arm in ("LADDER", "DETUNE", "SHIFT") if len(an[an["arm"] == arm])},
    "std_d_off_median": float(win[win["a_reseed"] == 0.0]
                              ["d_off_settled"].median()),
    "std_arm_d_off_medians": {
        arm: float(runs[(runs["cell"] == "WIN3D") & (runs["arm"] == arm)
                        & (runs["a_reseed"] == 0.0)]["d_off_settled"].median())
        for arm in ("DETUNE", "SHIFT")},
    "closure_max": float(an["closure_frac"].max()) if len(an) else None}
at = runs[runs["cell"] == "E_ANTITRAP"]
E["E3_anti_trap"] = {
    f"d{dd}": float(at[np.isclose(at["ind_x"], 0.79 - dd)]["ind_frac"].median())
    for dd in (0.10, 0.24) if len(at)}
E["E3_anti_trap"]["std_trap_iter8A"] = A.get("trap_z05", {}).get("curve")

dw_ = runs[runs["cell"] == "E_DENSEWIN"]
E["E4_densewin"] = {
    "A0": stats(dw_[dw_["a_reseed"] == 0.0], POST_METRICS + CTX_METRICS),
    "A1": stats(dw_[dw_["a_reseed"] == 1.0], POST_METRICS + CTX_METRICS),
    "n_channels": 116, "base_channels": 40}
rdE0, _p = ens_peaks(rd, "E_DENSEWIN", "BASE", "LADDER", 0.0)
rdE1, _p = ens_peaks(rd, "E_DENSEWIN", "BASE", "LADDER", 1.0)
E["E4_dense_lines_A1"] = line_test(rdE0, rdE1)
if float(amax) not in (0.0, 1.0):
    rdEM, _p = ens_peaks(rd, "E_DENSEWIN", "BASE", "LADDER", float(amax))
    E["E4_dense_lines_Amax"] = line_test(rdE0, rdEM)
ring_c = cand[(cand["cell"] == "E_DENSEWIN") & (cand["window"] == "SETTLED")
              & (cand["pickup"].str.contains("rin|rout"))]
E["E4_ring_census"] = {
    "n_ring_candidate_peaks": int(len(ring_c)),
    "n_ring_unexplained": int(ring_c["u015x"].sum()) if len(ring_c) else 0,
    "by_A": {str(a): {"n": int(len(g)), "unexplained": int(g["u015x"].sum())}
             for a, g in ring_c.groupby(
                 ring_c["run_id"].map(dict(zip(runs["run_id"],
                                               runs["a_reseed"]))))}}
eruns = runs[runs["tag"] == "E"]
E["stability"] = {"max_closure_frac": float(eruns["closure_frac"].max())
                  if len(eruns) else None,
                  "max_slope_settled": float(eruns["slope_settled"].max())
                  if len(eruns) else None,
                  "total_clamps": int(eruns["clamp_count"].sum())
                  if len(eruns) else None,
                  "n_runs": int(len(eruns))}
# two-dial scoring for every Phase E adjustment
for key, cellname in (("E1_indtrap", "E_INDTRAP"), ("E2_twocap", "E_TWOCAP"),
                      ("E3_anti", "E_ANTI"), ("E4_densewin", "E_DENSEWIN")):
    g = runs[runs["cell"] == cellname]
    if cellname == "E_ANTI":
        g = g[g["arm"] == "LADDER"]
    if not len(g):
        continue
    g0 = g[g["a_reseed"] == 0.0]
    g1 = g[g["a_reseed"] == 1.0]
    ok_, det_ = a1_absorption(g0, g1)
    E[key]["two_dial"] = {
        "dial1_need_composite_vs_WIN3D": consistency_composite(
            g0, WIN_A0_STATS),
        "dial2_A1_absorbed": ok_, "dial2_detail": det_,
        "dial2_A_max": (E["E1_indtrap"].get("A_max_adjusted_inducer")
                        if key == "E1_indtrap" else
                        (1.0 if ok_ else "<1 (only A=1 tested)"))}
R["phaseE"] = E


# ------------------------------------------------------------- trail -----
def _fmt(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "none"
    return f"{v:g}" if isinstance(v, (int, float, np.floating)) else str(v)


def config_id(r):
    return ("kz{}|{}|ind({},{},{})|cap({},{},{})|cap2({},{},{})|"
            "tr(th{})|dense{}|ncaps{}").format(
        int(r["kz"]), r["protocol"],
        _fmt(r.get("ind_x")), _fmt(r.get("ind_theta_deg")),
        _fmt(r.get("ind_z")), _fmt(r.get("cap_x")),
        _fmt(r.get("cap_theta_deg")), _fmt(r.get("cap_z")),
        _fmt(r.get("cap2_x")), _fmt(r.get("cap2_theta_deg")),
        _fmt(r.get("cap2_z")), _fmt(r.get("tr_theta_deg")),
        int(bool(r.get("dense_pickups", False))),
        int(r.get("n_caps", 1)) if not pd.isna(r.get("n_caps", 1)) else 1)


tr = runs.copy()
tr["config_id"] = [config_id(r) for r in tr.to_dict("records")]
tr["ensemble_id"] = (tr["tag"] + ":" + tr["cell"] + ":" + tr["arm"]
                     + ":A" + tr["a_reseed"].astype(str))
tr["drive_rng_stream"] = 10000 + tr["seed"]
tr["reseed_rng_stream"] = np.where(tr["a_reseed"] > 0,
                                   20000 + tr["seed"], np.nan)
tr["reseed_event_time"] = np.where(tr["a_reseed"] > 0, tr["t_off"], np.nan)
tr["reseed_event_count"] = (tr["a_reseed"] > 0).astype(int)
trail_cols = ["run_id", "tag", "cell", "arm", "seed", "config_id",
              "ensemble_id", "protocol", "mod", "ramp", "t_off",
              "drive_rng_stream", "reseed_rng_stream",
              "reseed_event_time", "reseed_event_count", "a_reseed",
              "demod_DE_amp", "demod_ST_amp", "P", "B_share_post",
              "B_share_out", "ind_frac", "d_off_settled",
              "d_line_settled", "slope_settled", "closure_frac",
              "sat_frac", "clamp_count", "E_cap_end"]
tr = tr[[c for c in trail_cols if c in tr.columns]].rename(columns={
    "demod_DE_amp": "clock1_pre_event_demod_amp",
    "demod_ST_amp": "clock2_post_event_demod_amp",
    "P": "post_P", "B_share_post": "post_B_share_post",
    "d_off_settled": "post_d_off_settled",
    "d_line_settled": "post_d_line_settled",
    "slope_settled": "post_slope_settled"})
tr.to_csv(f"{HD}/riemann_reseed3d_trail.csv", index=False)
R["trail"] = {"path": "riemann_reseed3d_trail.csv",
              "n_rows": int(len(tr)),
              "ensembles": {eid: [int(x) for x in g["run_id"]]
                            for eid, g in tr.groupby("ensemble_id")}}


# ---------------------------------------------------- anomaly scan (all) --
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
        out.append({"freq_median": float(sub["peak_freq"].median()),
                    "n_seeds": len(seeds), "seeds": [int(x) for x in seeds],
                    "n_peaks": int(len(sub)),
                    "u015x": bool(sub["u015x"].mean() > 0.5)})
    return out


unexp = {}
for (cell, arm, win_), g in cand.groupby(["cell", "arm", "window"]):
    cl = [c for c in chain_clusters(g) if c["n_seeds"] >= 2 and c["u015x"]]
    if cl:
        unexp[f"{cell}|{arm}|{win_}"] = cl
R["unexplained_clusters"] = unexp
R["stability"] = {
    "max_closure_frac": float(runs["closure_frac"].max()),
    "max_slope_settled": float(runs["slope_settled"].max()),
    "max_sat_frac": float(runs["sat_frac"].max()),
    "total_clamps": int(runs["clamp_count"].sum()),
    "n_runs": int(len(runs)),
    "n_peaks": int(len(pk)),
    "n_candidates_audited": int(cand_mask.sum()),
}

log_cols = ["run_id", "tag", "cell", "arm", "seed", "kz", "protocol",
            "a_reseed", "ind_x", "ind_theta_deg", "ind_z", "cap_x",
            "cap_theta_deg", "cap_z", "cap2_x", "cap2_theta_deg", "cap2_z",
            "tr_theta_deg", "dense_pickups", "n_caps", "E_cap1_end",
            "E_cap2_end", "ind_frac", "P", "demod_DE_amp",
            "demod_ST_amp", "B_share_out", "B_share_post",
            "d_off_settled", "d_line_settled", "closure_frac",
            "sat_frac", "slope_settled", "clamp_count", "E_cap_end"]
runs[[c for c in log_cols if c in runs.columns]].to_csv(
    f"{HD}/riemann_reseed3d_log.csv", index=False)
json.dump(R, open(f"{HD}/riemann_reseed3d_results.json", "w"), indent=1,
          default=str)
print("G0 pass:", R["G0"]["pass"], R["G0"])
print("PhaseA cap_z_best:", A["cap_z_best"], "trap:", A.get("trap_z05"))
print("PhaseB A_max:", amax)
print("PhaseC winner:", C["winner"], comp)
print("PhaseD Amax lines:", D["Amax_lines"], "A1 lines:", D["A1_lines"])
print("PhaseD circsd:", D["phase_coherence_circsd"])
print("PhaseE:", json.dumps({k: E[k] for k in ("E2_twocap",) if k in E}, default=str)[:400])
print("E1 A_max adj ind:", E.get("E1_indtrap", {}).get("A_max_adjusted_inducer"), "vs WIN3D", amax)
print("E4 dense lines A1:", E.get("E4_dense_lines_A1"))
print("Unexplained clusters:", unexp)
print("Stability:", R["stability"])
print("DONE analyze")
