"""
riemann_capacitor_sim.py — RIEMANN ITERATION 6: PROXIMAL CAPACITORS
(Zak's order 2026-10-22 19:29 EDT). Pre-registration:
ITERATION_LOG.md 2026-10-22 19:3x EDT — capacitor model, calibration
plan, grids, best-standing criterion, gates GC0-GC6 stated BEFORE
any run.

BASE: iteration-5 engine (TP cavity A FP=0.025 + cylinder B, global
ladder antenna, cubic mixing, Q=30, RK4 DT=0.0004, record 8.0 s,
windows EARLY/DRIVE_END/SETTLED). TRANSDUCER at A79 with
iteration-3 limits verbatim (G=250, V_SAT=10.9, EGATE=5.4e3,
eta=0.5) in EVERY run. Inducer: tone list verbatim from iter-5;
reference cell = tone G at (x=0.69, theta=0) (iter-5 d=0.10 ref).

NEW - CAPACITOR (cavity A only): state q; dq/dt = KAPPA*(W_C . v_A);
modal force on A = -KAPPA*q*W_C. W_C = iter-5 bump profile
(sigma 0.06, cos(n*theta)) NORMALIZED to unit modal L2, so the bare
resonance f_bare = KAPPA/(2 pi) is placement-independent.
E_cap = q^2/2 is a ledger STORAGE line (field power -KAPPA*q*(W_C.v)
is identically dE_cap/dt). KAPPA = 2*pi*G/2 (f_bare = tone
subharmonic 10.10735 Hz). Calibration (--calib): quiet ring test
(q(0)=300, no drive/inducer) -> FFT of q(t) = in-situ f_ring;
audit targets use {f_bare, f_ring}. Driven sensing run logs the
E_cap share (split-gain contingency pre-registered in the log).

PHASES (main): 1) OFFCAP / PLACEBO (arms L/D/S) + CAP grid
(LADDER); criterion pick; 2) REFINE 3x3 around winner + winner
under DETUNE/SHIFT + NODESWEEP (winner cap, inducer d sweep).
Modes: --calib | (default full) [--only=PREFIX ...]
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
AIND = {"A": AIND_A_DEF, "B": AIND_B_DEF}

# ---------------------------------------------------------- capacitor --
F_BARE = TONE_G / 2.0                    # 10.10735 Hz, tone subharmonic
KAPPA_C = 2.0 * np.pi * F_BARE           # 63.5035
DRIVE_SCALE = 1.0                        # 0.0 only in calibration ring

# ------------------------------------------------------------- rings --
RING_X = np.array([(0.79 + k / 8.0) % 1.0 for k in range(8)])
RING_IN_A = np.sin(np.pi * np.outer(RING_X, S.NIDX))                  # (8, M)
RING_OUT_A = RING_IN_A * np.cos(S.NIDX * 1.0)                         # (8, M)
CH28 = (S.CH_NAMES + [f"A_rin{k}" for k in range(8)]
        + [f"A_rout{k}" for k in range(8)])
_FS = 1.0 / (DT * S.REC_EVERY)
SITE_X = 0.79

# ---------------------------------------------------------- profiles --
_xgA = (np.arange(48) + 0.5) / 48.0

def a_ind_profile(x_ind, theta):
    b = np.exp(-((_xgA - x_ind) ** 2) / (2.0 * 0.06 ** 2))
    W = (2.0 / 48.0) * (S.PHI_G_A.T @ b)
    return W * np.cos(S.NIDX * theta)

def cap_profile(x_cap, theta):
    W = a_ind_profile(x_cap, theta)
    n = np.linalg.norm(W)
    return W / n if n > 0 else W

B_IND_W = TS.b_profile(S.B_PICK[0])          # iteration-3 PLACEBO profile

def _cap(x, theta_deg):
    return {"x": float(x), "theta": float(np.radians(theta_deg))}

def _cell(x, theta, cap=None, tones=None):
    return {"tones": list(tones) if tones else [TONE_G], "cav": "A",
            "x": float(x), "theta": float(theta), "amp_scale": None,
            "cap": cap}

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

# ------------------------------------------------------------- cells --
REF_X, REF_TH = 0.69, 0.0
CAP_RADII = [0.59, 0.65, 0.71, 0.79, 0.85]
CAP_ANGLES = [0, 30, 60, 90, 120, 150, 180]
NODE_DS = [0.10, 0.14, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.34, 0.40]

def stage1_cells():
    cells, arms = {}, {}
    cells["OFFCAP"] = _cell(REF_X, REF_TH, None)
    arms["OFFCAP"] = ["LADDER", "DETUNE", "SHIFT"]
    cells["PLACEBO"] = _cell(REF_X, REF_TH, _cap(0.13, 0))
    arms["PLACEBO"] = ["LADDER", "DETUNE", "SHIFT"]
    for x in CAP_RADII:
        for adeg in CAP_ANGLES:
            nm = f"CAP_x{x:.2f}_t{adeg:03d}"
            cells[nm] = _cell(REF_X, REF_TH, _cap(x, adeg))
            arms[nm] = ["LADDER"]
    return cells, arms

def refine_cells(winner_cell, existing):
    wx, wt = winner_cell["cap"]["x"], np.degrees(winner_cell["cap"]["theta"])
    cells, arms = {}, {}
    for dx in (-0.03, 0.0, 0.03):
        for dt_ in (-15.0, 0.0, 15.0):
            if dx == 0.0 and dt_ == 0.0:
                continue
            x = round(min(0.94, max(0.53, wx + dx)), 4)
            th = wt + dt_
            th = abs(th) if th < 0 else (360.0 - th if th > 180.0 else th)
            nm = f"RFC_x{x:.2f}_t{int(round(th)):03d}"
            if nm in existing or nm in cells:
                continue
            cells[nm] = _cell(REF_X, REF_TH, _cap(x, th))
            arms[nm] = ["LADDER"]
    return cells, arms

def nodesweep_cells(winner_cell):
    cells, arms = {}, {}
    wc = winner_cell["cap"]
    for d in NODE_DS:
        nm = f"NSW_d{d:.2f}"
        cells[nm] = _cell(round(SITE_X - d, 4), 0.0,
                          {"x": wc["x"], "theta": wc["theta"]})
        arms[nm] = ["LADDER"]
    return cells, arms

# ---------------------------------------------------------------- sim --
def run_cap(drive_freqs, seed, cell):
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
    force = 2.0 * (drive_freqs / S.Q) * omd * S.A_TARGET * DRIVE_SCALE
    cavT, W_T = TS.SITES["A79"][0], TS.SITES["A79"][1]
    W_I, cavI = ind_profile(cell)
    tones = cell["tones"] if cell else []
    amps = cell_amps(cell) if cell else []
    oms = [2.0 * np.pi * f for f in tones]
    K = len(tones)
    cap = cell.get("cap") if cell else None
    W_C = cap_profile(cap["x"], cap["theta"]) if cap else None
    has_cap = cap is not None
    NEXT = 3 if has_cap else 2

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
    E_cap_flow = 0.0
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
        q = y[4 * M + 2] if has_cap else 0.0
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
        if has_cap:
            dq = KAPPA_C * float(W_C @ va)
            aa = aa - KAPPA_C * q * W_C
            return np.concatenate([va, aa, vb, ab, [denv, dpool, dq]])
        return np.concatenate([va, aa, vb, ab, [denv, dpool]])

    y = np.zeros(4 * M + NEXT)
    if has_cap:
        y[4 * M + 2] = float(cell.get("q0", 0.0))
    for i in range(n_steps):
        t = i * DT
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        q = y[4 * M + 2] if has_cap else 0.0
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
        if has_cap:
            E_cap_flow += KAPPA_C * q * float(W_C @ va) * DT
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
                   "E_ind_cum": float(E_indk.sum()),
                   "q_cap": float(q), "E_cap": float(0.5 * q * q)}
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
    q_end = float(y[4 * M + 2]) if has_cap else 0.0
    E_cap_end = 0.5 * q_end * q_end
    Eend = (0.5 * (va ** 2 + stiffA * qa ** 2).sum(),
            0.5 * (vb ** 2 + stiffB * qb ** 2).sum())
    es = esamp[:ei]
    st = (es[:, 0] >= WINDOWS["SETTLED"][0]) & (es[:, 0] < WINDOWS["SETTLED"][1])
    Es = es[st, 1] + es[st, 2]
    slope_settled = float(np.polyfit(es[st, 0], np.log(np.maximum(Es, 1e-12)), 1)[0]) \
        if st.sum() > 10 else float("nan")
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
        "E_cap_end": float(E_cap_end), "E_cap_flow": float(E_cap_flow),
        "E_cap_flow_resid": float(E_cap_end - E_cap_flow),
        "slope_settled": slope_settled,
    }
    exch["closure_resid"] = (exch["E_inj"] + exch["E_ind"]) - (
        (exch["E_stored_end_A"] - exch["E_stored0_A"])
        + (exch["E_stored_end_B"] - exch["E_stored0_B"])
        + exch["E_diss_A"] + exch["E_diss_B"]
        + exch["E_rad_A"] + exch["E_rad_B"]
        + exch["E_tr_diss"] + pool_end + E_cap_end)
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


# ------------------------------------------------------------- coords --
def coordinates():
    return {
        "transducer": {"cavity": "A", "x": SITE_X,
                       "profile": "bump sigma 0.06 on 48-grid (iter-3 A79)"},
        "inducer_reference": {"cavity": "A", "x": REF_X, "theta_deg": 0.0},
        "capacitor": {"element": "q-state, KAPPA=2pi*G/2",
                      "f_bare": F_BARE, "kappa": float(KAPPA_C)},
        "pickups_A_in": [{"name": f"A_in{i + 1}", "cavity": "A",
                          "x": float(x)}
                         for i, x in enumerate(S.A_PICK_X)],
        "pickup_A_out": {"name": "A_out", "cavity": "A",
                         "x": float(S.A_APER_X), "kind": "aperture"},
        "pickups_B_in": [{"name": f"B_in{i + 1}", "cavity": "B",
                          "r": float(p[0]), "theta": float(p[1]),
                          "z": float(p[2])}
                         for i, p in enumerate(S.B_PICK)],
        "pickup_B_out": {"name": "B_out", "cavity": "B",
                         "r": float(S.B_APER[0]),
                         "theta": float(S.B_APER[1]),
                         "z": float(S.B_APER[2]), "kind": "aperture"},
        "rings_A": [{"name": f"A_rin{k}", "cavity": "A",
                     "x": float(RING_X[k]), "kind": "ring_in"}
                    for k in range(8)]
                   + [{"name": f"A_rout{k}", "cavity": "A",
                       "x": float(RING_X[k]), "kind": "ring_out"}
                      for k in range(8)],
    }


def cell_record(cell):
    if cell is None:
        return None
    cap = cell.get("cap")
    return {"tones": cell["tones"], "cav": cell["cav"], "x": cell["x"],
            "theta": cell["theta"],
            "theta_deg": float(np.degrees(cell["theta"])),
            "amps": cell_amps(cell),
            "cap": (None if cap is None else
                    {"x": cap["x"], "theta": cap["theta"],
                     "theta_deg": float(np.degrees(cap["theta"])),
                     "kappa": float(KAPPA_C), "f_bare": F_BARE})}


# --------------------------------------------------------------- driver --
class Runner:
    def __init__(self):
        self.peak_rows, self.ser_rows, self.run_rows = [], [], []
        self.meta_runs = []
        self.g1 = {}
        self.run_id = 0
        self.t0 = time.time()
        self.cells = {}

    def run_set(self, cells, arms_map):
        for cell_name, cell in cells.items():
            self.cells[cell_name] = cell
            for arm in arms_map[cell_name]:
                freqs = ARMS[arm]
                key = f"{cell_name}|{arm}"
                self.g1[key] = [0, 0]
                for seed in SEEDS:
                    self.run_id += 1
                    rid = self.run_id
                    res = run_cap(freqs, seed, cell)
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
                                self.peak_rows.append({
                                    "run_id": rid, "cell": cell_name,
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
                    seen = {r["peak_freq"] for r in self.peak_rows
                            if r["run_id"] == rid and r["cavity"] == "A"
                            and r["pickup_kind"] == "in"
                            and r["window"] == "DRIVE_END"}
                    for fd in freqs:
                        self.g1[key][1] += 1
                        if any(abs(sv - fd) <= S.PEAK_TOL for sv in seen):
                            self.g1[key][0] += 1
                    ts, _delays = S.transient_scalars(res)
                    prow = {"run_id": rid, "cell": cell_name, "arm": arm,
                            "seed": seed, **res["exch"], **ts}
                    cap = cell.get("cap") if cell else None
                    prow["ind_x"] = cell["x"] if cell else None
                    prow["ind_theta_deg"] = (float(np.degrees(cell["theta"]))
                                             if cell else None)
                    prow["cap_x"] = cap["x"] if cap else None
                    prow["cap_theta_deg"] = (float(np.degrees(cap["theta"]))
                                             if cap else None)
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
                    self.run_rows.append(prow)
                    for r in res["series"]:
                        self.ser_rows.append({"run_id": rid, "cell": cell_name,
                                              "arm": arm, "seed": seed, **r})
                    self.meta_runs.append({"run_id": rid, "cell": cell_name,
                                           "arm": arm, "seed": seed})
                    print(f"run {rid} {cell_name} {arm} seed{seed} "
                          f"closure={res['exch']['closure_frac']:.5f} "
                          f"Eind={res['exch']['E_ind']:.1f} "
                          f"Ecap={res['exch']['E_cap_end']:.1f} "
                          f"elapsed={time.time() - self.t0:.0f}s", flush=True)

    def dump(self, cells_all):
        for name, rows in (("riemann_capacitor_peaks_raw.csv", self.peak_rows),
                           ("riemann_capacitor_series.csv", self.ser_rows),
                           ("riemann_capacitor_runs.csv", self.run_rows)):
            flds = []
            for r in rows:
                for k in r.keys():
                    if k not in flds:
                        flds.append(k)
            with open(f"{HD}/{name}", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=flds)
                w.writeheader(); w.writerows(rows)
        meta = {"runs": self.meta_runs, "AIND": dict(AIND),
                "tones": {"G": TONE_G, "H": TONE_H, "MID": TONE_MID,
                          "GM": TONE_GM, "GP": TONE_GP},
                "capacitor": {"kappa": float(KAPPA_C), "f_bare": F_BARE},
                "coordinates": coordinates(),
                "cells": {c: cell_record(v) for c, v in cells_all.items()},
                "g1": {k: {"hits": v[0], "total": v[1]}
                       for k, v in self.g1.items()}}
        json.dump(meta, open(f"{HD}/riemann_capacitor_raw_meta.json", "w"),
                  indent=1)


def persistence(run_rows, cell_name):
    """Per-seed P = SETTLED/DRIVE_END demod amp of tone 0 at A_in5."""
    out = {}
    for r in run_rows:
        if r["cell"] == cell_name and r["arm"] == "LADDER":
            de = r.get("demod_t0_A_in5_DRIVE_END_amp", 0.0)
            st = r.get("demod_t0_A_in5_SETTLED_amp", 0.0)
            out[r["seed"]] = (st / de if de and de > 1e-9 else float("nan"),
                              de, r)
    return out


def pick_winner(runner, cand_names):
    off = persistence(runner.run_rows, "OFFCAP")
    off_de = np.median([v[1] for v in off.values()]) if off else 0.0
    table = {}
    for nm in cand_names:
        per = persistence(runner.run_rows, nm)
        if len(per) < 3:
            continue
        ps, ok = [], True
        for seed, (p, de, r) in sorted(per.items()):
            guard = (de >= 0.5 * off_de and r["sat_frac"] <= 0.30
                     and r["closure_frac"] <= 0.001
                     and r["slope_settled"] <= 0.005)
            ps.append(p)
            ok = ok and bool(guard)
        table[nm] = {"P_median": float(np.median(ps)),
                     "P_per_seed": [float(x) for x in ps],
                     "guards_all_pass": bool(ok)}
    passing = {k: v for k, v in table.items() if v["guards_all_pass"]}
    pool = passing if passing else table
    winner = max(pool, key=lambda k: pool[k]["P_median"]) if pool else None
    return winner, table, {"OFFCAP_P_median":
                           float(np.median([v[0] for v in off.values()]))
                           if off else None,
                           "OFFCAP_DE_median": float(off_de)}


# ---------------------------------------------------------- calibration --
def calib():
    global DRIVE_SCALE
    out = {"f_bare": F_BARE, "kappa": float(KAPPA_C)}
    # ring test: no drive, no inducer tones, q(0) = 300
    DRIVE_SCALE = 0.0
    cell = {"tones": [], "cav": "A", "x": REF_X, "theta": 0.0,
            "amp_scale": None, "cap": _cap(0.69, 0), "q0": 300.0}
    res = run_cap(S.T_LADDER, 0, cell)
    DRIVE_SCALE = 1.0
    ser = res["series"]
    tt = np.array([r["t"] for r in ser])
    qq = np.array([r["q_cap"] for r in ser])
    m = (tt >= 1.0) & (tt <= 7.5)
    x = qq[m] - qq[m].mean()
    fs = 1.0 / (tt[1] - tt[0])
    spec = np.abs(np.fft.rfft(x * np.hanning(len(x))))
    fr = np.fft.rfftfreq(len(x), d=1.0 / fs)
    band = (fr >= 2.0) & (fr <= 40.0)
    f_ring = float(fr[band][np.argmax(spec[band])])
    out["ring"] = {"f_ring": f_ring, "q_max": float(np.max(np.abs(qq))),
                   "q_end": float(qq[-1]),
                   "E_cap_end": res["exch"]["E_cap_end"],
                   "closure_frac": res["exch"]["closure_frac"]}
    # driven sensing run at the reference geometry
    cell2 = _cell(REF_X, REF_TH, _cap(0.69, 0))
    res2 = run_cap(S.T_LADDER, 0, cell2)
    ex = res2["exch"]
    out["sensing"] = {"cell": "REF + cap(0.69,0)", "E_inj": ex["E_inj"],
                      "E_ind": ex["E_ind"], "E_cap_end": ex["E_cap_end"],
                      "E_cap_flow": ex["E_cap_flow"],
                      "E_cap_flow_resid": ex["E_cap_flow_resid"],
                      "E_cap_share_of_inj": ex["E_cap_end"] / ex["E_inj"],
                      "closure_frac": ex["closure_frac"],
                      "sat_frac": ex["sat_frac"]}
    json.dump(out, open(f"{HD}/riemann_capacitor_calib.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


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
    if fl.get("calib"):
        calib()
        return
    AIND["A"] = float(fl.get("aind-a", AIND_A_DEF))
    AIND["B"] = float(fl.get("aind-b", AIND_B_DEF))
    only = fl["only"]
    runner = Runner()
    cells1, arms1 = stage1_cells()
    if only:
        cells1 = {k: v for k, v in cells1.items()
                  if any(o in k for o in only)}
        arms1 = {k: arms1[k] for k in cells1}
    runner.run_set(cells1, arms1)
    all_cells = dict(cells1)
    runner.dump(all_cells)
    cand = [c for c in cells1 if c.startswith("CAP_")]
    winner, table, offinfo = pick_winner(runner, cand)
    print("CRITERION:", json.dumps({"winner": winner, **offinfo}), flush=True)
    json.dump({"winner": winner, "table": table, **offinfo},
              open(f"{HD}/riemann_capacitor_criterion.json", "w"), indent=1)
    if winner is not None:
        wcell = cells1[winner]
        cells2, arms2 = refine_cells(wcell, set(all_cells))
        arms2 = {k: arms2[k] for k in cells2}
        runner.run_set(cells2, arms2)
        all_cells.update(cells2)
        # re-pick over the enlarged candidate pool (pre-registered:
        # winner = best P over grid+refine under the same guards)
        cand2 = cand + list(cells2)
        winner2, table2, _ = pick_winner(runner, cand2)
        if winner2 is not None:
            winner, table = winner2, table2
            wcell = all_cells[winner]
            json.dump({"winner": winner, "table": table, **offinfo},
                      open(f"{HD}/riemann_capacitor_criterion.json", "w"),
                      indent=1)
            print("CRITERION FINAL:", winner, flush=True)
        cells3 = {winner: wcell}
        arms3 = {winner: ["DETUNE", "SHIFT"]}
        runner.run_set(cells3, arms3)
        cells4, arms4 = nodesweep_cells(wcell)
        runner.run_set(cells4, arms4)
        all_cells.update(cells4)
    runner.dump(all_cells)
    print(f"DONE runs={runner.run_id} peaks={len(runner.peak_rows)} "
          f"elapsed={time.time() - runner.t0:.0f}s")


if __name__ == "__main__":
    main()
