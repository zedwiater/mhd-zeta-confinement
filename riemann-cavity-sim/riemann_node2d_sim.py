"""
riemann_node2d_sim.py — RIEMANN ITERATION 5: NODE MAPPING + TRUE 2D
SPREAD (Zak's order 2026-10-22 17:39 EDT). Pre-registration:
ITERATION_LOG.md 2026-10-22 17:44 EDT — grids, tone confounds,
discriminator criterion and gates GN0-GN6 stated BEFORE any run.

BASE: iteration-2/3/4 system imported bit-identical (TP cavity A,
FP=0.025, cylinder cavity B, global ladder antenna, cubic mixing,
Q=30, RK4 DT=0.0004, record 8.0 s, windows EARLY/DRIVE_END/SETTLED).
TRANSDUCER at A79 with iteration-3 limits verbatim (G=250,
V_SAT=10.9, EGATE=5.4e3, eta=0.5) in EVERY run.

CHANGES vs iteration 4:
- Inducer generalized to a TONE LIST (comb-ready): force =
  sum_k amp_k * W * cos(2 pi f_k t), 0.5 s ramp, held the full 8.0 s.
  Per-component work E_ind_k integrated exactly; E_ind = sum_k.
- Cells are LITERAL (x, theta) pairs (the iter-4 PU/PD
  radius-independence degeneracy was a parameterization defect:
  angular cells derived x from a nominal distance label). Profile
  itself verbatim: bump sigma 0.06 on the 48-grid, W = (2/48) PHI^T b,
  W_n *= cos(n*theta). Theta in [0, pi] (cos even — stated).
- P1: densified node sweep (tone G, theta=0, steps <=0.02) + the
  same sweep at MID / G-5% / G+5% (frequency-dependence
  discriminator). P2: 2D grid x x theta, tone G. P3: comb
  {G, MID, H} at node + reference (amps 4500/sqrt(3) each).

TONES: G=20.2147; H=G*(t2/t1)=30.0646; MID=(G+H)/2=25.1397
(0.1288 from rung 3 — confound pre-stated); GM=0.95G=19.2040;
GP=1.05G=21.2254 (0.2034 from rung 2 — confound pre-stated).
A_IND_A=4500 every single-tone A cell; A_IND_B=900 (placebo).

Modes: --grid [--only=PREFIX ...] (substring match on cell names).
"""

import csv
import json
import sys
import time

import numpy as np

import riemann_tp_cylinder_sim as S
import riemann_transducer_sim as TS

HD = S.HD
M = S.M
DT = S.DT
T_OFF = S.T_OFF
T_TOTAL = S.T_TOTAL
WINDOWS = S.WINDOWS
ARMS = S.ARMS
SEEDS = [0, 1, 2]
SER_EVERY = 20

G_TR, VSAT, EGATE, ETA, TAU_ENV = 250.0, 10.9, 5400.0, 0.5, 0.25
TONE_G = 20.2147
TONE_H = TONE_G * float(S.T_LADDER[1] / S.T_LADDER[0])   # 30.0646
TONE_MID = 0.5 * (TONE_G + TONE_H)                       # 25.1397
TONE_GM = 0.95 * TONE_G                                  # 19.2040
TONE_GP = 1.05 * TONE_G                                  # 21.2254
AIND_A_DEF, AIND_B_DEF = 4500.0, 900.0

# ------------------------------------------------------------- rings --
RING_X = np.array([(0.79 + k / 8.0) % 1.0 for k in range(8)])
RING_IN_A = np.sin(np.pi * np.outer(RING_X, S.NIDX))                  # (8, M)
RING_OUT_A = RING_IN_A * np.cos(S.NIDX * 1.0)                         # (8, M)
CH28 = (S.CH_NAMES + [f"A_rin{k}" for k in range(8)]
        + [f"A_rout{k}" for k in range(8)])
_FS = 1.0 / (DT * S.REC_EVERY)
SITE_X = 0.79

# ---------------------------------------------------------- inducer --
_xgA = (np.arange(48) + 0.5) / 48.0

def a_ind_profile(x_ind, theta):
    b = np.exp(-((_xgA - x_ind) ** 2) / (2.0 * 0.06 ** 2))
    W = (2.0 / 48.0) * (S.PHI_G_A.T @ b)
    return W * np.cos(S.NIDX * theta)

B_IND_W = TS.b_profile(S.B_PICK[0])          # iteration-3 PLACEBO profile

def _acell(tones, x, theta, amp_scale=None):
    return {"tones": list(tones), "cav": "A", "x": x, "theta": theta,
            "amp_scale": amp_scale}

CELLS = {"OFF": None,
         "PLACEBO_B": {"tones": [TONE_G], "cav": "B", "x": None,
                       "theta": 0.0, "amp_scale": None}}
CELLS_ARMS = {c: ["LADDER"] for c in CELLS}

# P1 — densified node sweep, theta = 0, x = SITE_X - d
G_DS = [0.10, 0.12, 0.14, 0.16, 0.18, 0.19, 0.20, 0.21, 0.22, 0.24,
        0.26, 0.28]
ALT_DS = [0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28]
for d in G_DS:
    CELLS[f"G_d{d:.2f}"] = _acell([TONE_G], round(SITE_X - d, 4), 0.0)
for _pfx, _tone in (("MID", TONE_MID), ("GM", TONE_GM), ("GP", TONE_GP)):
    for d in ALT_DS:
        CELLS[f"{_pfx}_d{d:.2f}"] = _acell([_tone], round(SITE_X - d, 4), 0.0)
for c in CELLS:
    CELLS_ARMS.setdefault(c, ["LADDER"])
CELLS_ARMS["G_d0.10"] = ["LADDER", "DETUNE", "SHIFT"]
CELLS_ARMS["G_d0.20"] = ["LADDER", "DETUNE", "SHIFT"]

# P2 — true 2D grid, tone G; (0.59,0) = G_d0.20 and (0.65,0) = G_d0.14
# are NOT duplicated (analysis reads the sweep runs for those points).
GRID_RADII = [0.59, 0.65, 0.71, 0.85, 0.91]
GRID_ANGLES_DEG = [0, 30, 60, 90, 120, 150, 180]
for x in GRID_RADII:
    for adeg in GRID_ANGLES_DEG:
        if adeg == 0 and x in (0.59, 0.65):
            continue
        nm = f"G2D_x{x:.2f}_t{adeg:03d}"
        CELLS[nm] = _acell([TONE_G], x, np.radians(adeg))
        CELLS_ARMS[nm] = ["LADDER"]

# P3 — comb {G, MID, H} at node + reference
_SQ3 = 1.0 / np.sqrt(3.0)
for _nm, _x in (("COMB_d0.20", 0.59), ("COMB_d0.10", 0.69)):
    CELLS[_nm] = _acell([TONE_G, TONE_MID, TONE_H], _x, 0.0,
                        amp_scale=[_SQ3] * 3)
    CELLS_ARMS[_nm] = ["LADDER"]

AIND = {"A": AIND_A_DEF, "B": AIND_B_DEF}

def cell_amps(cell):
    base = AIND[cell["cav"]]
    if cell.get("amp_scale"):
        return [base * s for s in cell["amp_scale"]]
    return [base] * len(cell["tones"])

def ind_profile(cell):
    if cell is None:
        return None, None
    if cell["cav"] == "A":
        return a_ind_profile(cell["x"], cell["theta"]), "A"
    return B_IND_W, "B"

# GN0 (vector part): min same-theta profile separation over grid radii
def profile_check():
    out = {}
    worst = 1e9
    for adeg in GRID_ANGLES_DEG:
        th = np.radians(adeg)
        Ws = [a_ind_profile(x, th) for x in GRID_RADII]
        for i in range(len(Ws)):
            for j in range(i + 1, len(Ws)):
                rel = (np.linalg.norm(Ws[i] - Ws[j])
                       / max(np.linalg.norm(Ws[i]), np.linalg.norm(Ws[j])))
                worst = min(worst, float(rel))
        out[f"theta_{adeg}"] = True
    chk = (np.linalg.norm(a_ind_profile(0.59, np.pi / 2)
                          - a_ind_profile(0.71, np.pi / 2))
           / np.linalg.norm(a_ind_profile(0.71, np.pi / 2)))
    return {"min_same_theta_rel_L2": worst,
            "check_pair_0.59_vs_0.71_th90_rel_L2": float(chk)}

# ---------------------------------------------------------------- sim --
def run_node(drive_freqs, seed, cell):
    rng = np.random.default_rng(10_000 + seed)
    phases = rng.uniform(0.0, 2.0 * np.pi, M)
    F_A, F_B = S.F_EIG_A, S.F_EIG_B
    omA, omB = 2.0 * np.pi * F_A, 2.0 * np.pi * F_B
    gA = F_A / S.Q * (1.0 + S.RAD_FRAC)
    gB = F_B / S.Q * (1.0 + S.RAD_FRAC)
    gA_rad = F_A / S.Q * S.RAD_FRAC
    gB_rad = F_B / S.Q * S.RAD_FRAC
    stiffA, stiffB = omA ** 2, omB ** 2
    kappa = S.KAPPA_FRAC * stiffA
    omd = 2.0 * np.pi * drive_freqs
    force = 2.0 * (drive_freqs / S.Q) * omd * S.A_TARGET
    cavT, W_T = TS.SITES["A79"][0], TS.SITES["A79"][1]
    W_I, cavI = ind_profile(cell)
    tones = cell["tones"] if cell else []
    amps = cell_amps(cell) if cell else []
    oms = [2.0 * np.pi * f for f in tones]
    K = len(tones)

    n_steps = int(round(T_TOTAL / DT))
    n_rec = int(round(T_TOTAL / (DT * S.REC_EVERY)))
    n_es = int(round(T_TOTAL / (DT * S.ESAMP_EVERY)))
    rec = np.empty((28, n_rec))
    esamp = np.empty((n_es, 5))
    emode = np.empty((n_es, 2 * M))
    series = []
    inj_m = np.zeros(M); dissA = np.zeros(M); dissB = np.zeros(M)
    radA_m = np.zeros(M); radB_m = np.zeros(M)
    winE = {w: (np.zeros(M), np.zeros(M)) for w in WINDOWS}
    ri = ei = 0
    E0 = None
    E_indk = np.zeros(K)
    E_ext = E_reinj = E_tr_diss = 0.0
    sat_time = gate_time = 0.0
    clamp_count = 0
    clamp_E = 0.0
    pool_max = 0.0

    def drive(t):
        if t >= T_OFF:
            return np.zeros(M)
        return min(1.0, t / 0.5) * force * np.cos(omd * t + phases)

    def ind_forces(t):
        if cell is None:
            return None
        ramp = min(1.0, t / 0.5)
        return [amps[k] * W_I * np.cos(oms[k] * t) * ramp
                for k in range(K)]

    def tr_forces(va, vb, env, pool):
        v = va if cavT == "A" else vb
        s = float(W_T @ v)
        f1 = -G_TR * VSAT * np.tanh(s / VSAT) * W_T
        p_ext = -float(f1 @ v)
        sig = pool / (pool + EGATE) if pool > 0.0 else 0.0
        f2 = G_TR * sig * env * W_T
        p_reinj = float(f2 @ v)
        return f1 + f2, p_ext, p_reinj, s

    def rhs(y, t):
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        ua = S.PHI_G_A @ qa
        ub = S.S_B @ qb
        nla = S.BETA * (S.PHI_G_A.T @ (ua ** 3)) * (2.0 / 48.0)
        nlb = S.BETA * (S.S_B.T @ (ub ** 3)) * (2.0 / 64.0)
        aa = -2.0 * gA * va - stiffA * qa + kappa * (qb - qa) + nla + drive(t)
        ab = -2.0 * gB * vb - stiffB * qb + kappa * (qa - qb) + nlb
        F, p_ext, p_reinj, s = tr_forces(va, vb, env, pool)
        if cavT == "A":
            aa = aa + F
        else:
            ab = ab + F
        Fk = ind_forces(t)
        if Fk is not None:
            Fi = np.sum(Fk, axis=0)
            if cavI == "A":
                aa = aa + Fi
            else:
                ab = ab + Fi
        denv = (abs(s) - env) / TAU_ENV
        dpool = ETA * p_ext - p_reinj
        return np.concatenate([va, aa, vb, ab, [denv, dpool]])

    y = np.zeros(4 * M + 2)
    for i in range(n_steps):
        t = i * DT
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        dv = drive(t)
        inj_m += dv * va * DT
        dissA += 2.0 * (F_A / S.Q) * va ** 2 * DT
        dissB += 2.0 * (F_B / S.Q) * vb ** 2 * DT
        radA_m += 2.0 * gA_rad * va ** 2 * DT
        radB_m += 2.0 * gB_rad * vb ** 2 * DT
        Fk = ind_forces(t)
        if Fk is not None:
            vv = va if cavI == "A" else vb
            for k in range(K):
                E_indk[k] += float(Fk[k] @ vv) * DT
        F, p_ext, p_reinj, s_act = tr_forces(va, vb, env, pool)
        E_ext += p_ext * DT
        E_reinj += p_reinj * DT
        E_tr_diss += (1.0 - ETA) * p_ext * DT
        if abs(s_act) >= 2.0 * VSAT:
            sat_time += DT
        sig_now = pool / (pool + EGATE) if pool > 0.0 else 0.0
        if sig_now < 0.5:
            gate_time += DT
        eA_m = 0.5 * (va ** 2 + stiffA * qa ** 2)
        eB_m = 0.5 * (vb ** 2 + stiffB * qb ** 2)
        if E0 is None:
            E0 = (eA_m.sum(), eB_m.sum())
        for wname, (w0, w1) in WINDOWS.items():
            if w0 <= t < w1:
                winE[wname][0][:] += (va ** 2 + stiffA * qa ** 2) * DT
                winE[wname][1][:] += (vb ** 2 + stiffB * qb ** 2) * DT
        if i % S.ESAMP_EVERY == 0 and ei < n_es:
            pab = float(np.sum(kappa * (qa - qb) * vb))
            pba = float(np.sum(kappa * (qb - qa) * va))
            esamp[ei] = (t, eA_m.sum(), eB_m.sum(), pab, pba)
            emode[ei] = np.concatenate([eA_m, eB_m])
            ei += 1
        if i % SER_EVERY == 0:
            row = {"t": round(t, 4), "E_A": float(eA_m.sum()),
                   "E_B": float(eB_m.sum()), "pool": float(pool),
                   "env": float(env), "s_act": float(s_act),
                   "E_ind_cum": float(E_indk.sum())}
            for nm, (cv, W, _pk, _nb, cw) in TS.SITES.items():
                vv = va if cv == "A" else vb
                ee = eA_m if cv == "A" else eB_m
                row[f"s_{nm}"] = float(W @ vv)
                row[f"eloc_{nm}"] = float(cw @ ee)
            series.append(row)
        if i % S.REC_EVERY == 0 and ri < n_rec:
            rec[:5, ri] = S.PHI_A_IN @ qa
            rec[5, ri] = np.sum(np.sqrt(2.0 * gA_rad) * va * S.PHI_A_AP)
            rec[6:11, ri] = S.PHI_B_IN @ qb
            rec[11, ri] = np.sum(np.sqrt(2.0 * gB_rad) * vb * S.PHI_B_AP)
            rec[12:20, ri] = RING_IN_A @ qa
            rec[20:28, ri] = RING_OUT_A @ qa
            ri += 1
        k1 = rhs(y, t)
        k2 = rhs(y + 0.5 * DT * k1, t + 0.5 * DT)
        k3 = rhs(y + 0.5 * DT * k2, t + 0.5 * DT)
        k4 = rhs(y + DT * k3, t + DT)
        y = y + (DT / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if y[4 * M + 1] < 0.0:
            clamp_count += 1
            clamp_E += -y[4 * M + 1]
            y[4 * M + 1] = 0.0
        if y[4 * M] < 0.0:
            y[4 * M] = 0.0
        pool_max = max(pool_max, float(y[4 * M + 1]))
    qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
    pool_end = float(y[4 * M + 1])
    Eend = (0.5 * (va ** 2 + stiffA * qa ** 2).sum(),
            0.5 * (vb ** 2 + stiffB * qb ** 2).sum())
    rms = np.sqrt(np.mean(rec ** 2, axis=1, keepdims=True)) + 1e-12
    rec = rec + rng.standard_normal(rec.shape) * S.NOISE_FRAC * rms

    E_ind = float(E_indk.sum())
    exch = {
        "E_inj": float(inj_m.sum()), "E_ind": E_ind,
        "E_diss_A": float(dissA.sum()), "E_diss_B": float(dissB.sum()),
        "E_rad_A": float(radA_m.sum()), "E_rad_B": float(radB_m.sum()),
        "E_stored0_A": float(E0[0]), "E_stored0_B": float(E0[1]),
        "E_stored_end_A": float(Eend[0]), "E_stored_end_B": float(Eend[1]),
        "E_tr_ext": float(E_ext), "E_tr_reinj_net": float(E_reinj),
        "E_tr_diss": float(E_tr_diss), "pool_end": pool_end,
        "pool_max": float(pool_max), "clamp_count": int(clamp_count),
        "clamp_E": float(clamp_E),
        "sat_frac": float(sat_time / T_TOTAL),
        "gate_frac": float(gate_time / T_TOTAL),
    }
    exch["closure_resid"] = (exch["E_inj"] + exch["E_ind"]) - (
        (exch["E_stored_end_A"] - exch["E_stored0_A"])
        + (exch["E_stored_end_B"] - exch["E_stored0_B"])
        + exch["E_diss_A"] + exch["E_diss_B"]
        + exch["E_rad_A"] + exch["E_rad_B"]
        + exch["E_tr_diss"] + pool_end)
    exch["closure_frac"] = (abs(exch["closure_resid"])
                            / (exch["E_inj"] + abs(exch["E_ind"]))
                            if exch["E_inj"] > 0 else float("nan"))
    exch["pool_balance_resid"] = pool_end - (ETA * E_ext - E_reinj)
    return {"rec": rec, "esamp": esamp[:ei], "emode": emode[:ei],
            "inj_m": inj_m, "radA_m": radA_m, "radB_m": radB_m,
            "winE": winE, "exch": exch, "series": series,
            "E_indk": E_indk}


# ------------------------------------------------------------- purity --
def demod(series, f_tone, wname):
    w0, w1 = WINDOWS[wname]
    i0, i1 = int(w0 * _FS), int(w1 * _FS)
    x = series[i0:i1]
    tt = np.arange(i0, i1) / _FS
    z = np.mean(x * np.exp(-2j * np.pi * f_tone * tt)) * 2.0
    return float(abs(z)), float(np.degrees(np.angle(z)))


# --------------------------------------------------------------- driver --
def _flags():
    fl = {"only": []}
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            fl["only"].append(a.split("=", 1)[1])
        elif a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            fl[k] = v
        elif a.startswith("--"):
            fl[a[2:]] = True
    return fl


def main():
    fl = _flags()
    AIND["A"] = float(fl.get("aind-a", AIND_A_DEF))
    AIND["B"] = float(fl.get("aind-b", AIND_B_DEF))
    only = fl["only"]
    t0 = time.time()

    peak_rows, ser_rows, run_rows = [], [], []
    meta = {"runs": [], "AIND": dict(AIND),
            "tones": {"G": TONE_G, "H": TONE_H, "MID": TONE_MID,
                      "GM": TONE_GM, "GP": TONE_GP},
            "profile_check": profile_check(),
            "cells": {c: (None if v is None else
                          {"tones": v["tones"], "cav": v["cav"],
                           "x": v["x"], "theta": v["theta"],
                           "amps": cell_amps(v)})
                      for c, v in CELLS.items()}}
    run_id = 0
    g1 = {}
    for cell_name, cell in CELLS.items():
        if only and not any(o in cell_name for o in only):
            continue
        for arm in CELLS_ARMS[cell_name]:
            freqs = ARMS[arm]
            key = f"{cell_name}|{arm}"
            g1[key] = [0, 0]
            for seed in SEEDS:
                run_id += 1
                res = run_node(freqs, seed, cell)
                tones = cell["tones"] if cell else []
                for ci, ch in enumerate(CH28):
                    if ci < 6:
                        cav, kind = "A", ("out" if ch.endswith("out") else "in")
                    elif ci < 12:
                        cav, kind = "B", ("out" if ch.endswith("out") else "in")
                    elif ci < 20:
                        cav, kind = "A", "ring_in"
                    else:
                        cav, kind = "A", "ring_out"
                    for wname in WINDOWS:
                        for pf, pa in S.window_peaks(res["rec"][ci], wname):
                            is_ind = bool(tones) and bool(
                                np.min(np.abs(np.array(tones) - pf)) <= 0.15)
                            is_drv = (not is_ind) and bool(
                                np.min(np.abs(freqs - pf)) <= S.PEAK_TOL)
                            peak_rows.append({
                                "run_id": run_id, "cell": cell_name,
                                "arm": arm, "seed": seed, "cavity": cav,
                                "pickup": ch, "pickup_kind": kind,
                                "window": wname,
                                "peak_freq": round(pf, 4),
                                "peak_amp": round(pa, 6),
                                "is_drive": is_drv, "is_inducer": is_ind,
                                "audit_label": "inducer" if is_ind else
                                ("drive" if is_drv else
                                 S.audit_label(pf, freqs, "NONE", 0.0)),
                            })
                seen = {r["peak_freq"] for r in peak_rows
                        if r["run_id"] == run_id and r["cavity"] == "A"
                        and r["pickup_kind"] == "in"
                        and r["window"] == "DRIVE_END"}
                for fd in freqs:
                    g1[key][1] += 1
                    if any(abs(sv - fd) <= S.PEAK_TOL for sv in seen):
                        g1[key][0] += 1
                ts, _delays = S.transient_scalars(res)
                prow = {"run_id": run_id, "cell": cell_name, "arm": arm,
                        "seed": seed, **res["exch"], **ts}
                for k in range(3):
                    prow[f"eind_{k}"] = (round(float(res["E_indk"][k]), 2)
                                         if k < len(tones) else 0.0)
                if cell:
                    for k, f_tone in enumerate(tones):
                        for pkname, chidx in (("rin0", 12), ("A_in5", 4),
                                              ("A_out", 5)):
                            for wname in WINDOWS:
                                da, dp = demod(res["rec"][chidx], f_tone,
                                               wname)
                                prow[f"demod_t{k}_{pkname}_{wname}_amp"] = \
                                    round(da, 6)
                                prow[f"demod_t{k}_{pkname}_{wname}_phase"] = \
                                    round(dp, 2)
                run_rows.append(prow)
                for r in res["series"]:
                    ser_rows.append({"run_id": run_id, "cell": cell_name,
                                     "arm": arm, "seed": seed, **r})
                meta["runs"].append({"run_id": run_id, "cell": cell_name,
                                     "arm": arm, "seed": seed})
                print(f"run {run_id} {cell_name} {arm} seed{seed} "
                      f"closure={res['exch']['closure_frac']:.5f} "
                      f"Eind={res['exch']['E_ind']:.1f} "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
    for name, rows in (("riemann_node2d_peaks_raw.csv", peak_rows),
                       ("riemann_node2d_series.csv", ser_rows),
                       ("riemann_node2d_runs.csv", run_rows)):
        flds = []
        for r in rows:
            for k in r.keys():
                if k not in flds:
                    flds.append(k)
        with open(f"{HD}/{name}", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=flds)
            w.writeheader(); w.writerows(rows)
    meta["g1"] = {k: {"hits": v[0], "total": v[1]} for k, v in g1.items()}
    json.dump(meta, open(f"{HD}/riemann_node2d_raw_meta.json", "w"),
              indent=1)
    print(f"DONE runs={run_id} peaks={len(peak_rows)} "
          f"elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
