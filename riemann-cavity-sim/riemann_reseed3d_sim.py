"""
riemann_reseed3d_sim.py — RIEMANN ITERATION 8: THIRD DIMENSION +
RE-SEED CONSISTENCY (Zak's order 2026-10-22 ~23:38 EDT).
PRE-REGISTRATION: riemann_reseed3d_prereg.json (written before any
run; definitions, gates, A_max rule, pivot, cut order fixed there).

BASE: iteration-6/7 engine (riemann_capacitor_sim.run_cap) — TP
cavity A (FP=0.025) + cylinder B, transducer A79 in EVERY run
(iter-3 limits), inducer tone G, capacitor kappa = 2 pi G/2.
Base geometry = iter-6 winner: inducer (0.69, 0deg), cap (0.79, 15deg).

PHASE A (3D): cavity A gains axial coordinate z in [0,1]. Mode n
gets axial wavenumber m_n = n mod 3, envelope chi_n(z) =
sqrt(2) cos(pi m_n z) (1 for m_n=0), shapes on the (48 x 16) grid
Phi[(x,z),n] = sin(pi n x) chi_n(z); cubic mixing projected on
that grid. Eigenfrequencies unchanged (the ladder IS the
spectrum; axial structure is spatial only — stated convention).
KZ=0 collapses chi to 1 and uses iteration-6 builders verbatim:
the z-uniform limit, gated by G0 against iter-6 headline numbers.
Element bumps are (x,theta,z) Gaussians (sigma_x 0.06, sigma_z
0.12); pickups on 3 z-planes {0.25, 0.5, 0.75}; 40 channels.

PHASE B (re-seed): SEED = rng(10000+seed) drive-phase draw +
readout noise (engine convention). RE-SEED EVENT at drive-off:
per mode of both cavities, a_n = omega q + i v; phase rotated a
fraction A_reseed toward a fresh uniform draw (rng 20000+seed);
per-mode energies exactly invariant; A=0 bitwise baseline.
Ensembles of 8 seeds; A_max per the pre-registered rule.

PHASE C: ladder-drive protocol grid (duration / ramp / AM / FM
at iteration-1 depth) x 8 seeds; maximal-consistency input =
argmin median SD ratio vs BASE.

PHASE D data: ring-down segments (t >= 3.9 s) of channels
[A_in5@z.5, A_in5@z.25, A_out@z.5, B_in1] + full A_in5@z.5 saved
to riemann_reseed3d_ringdown.npz for the resolution ladder.

Modes: --ring3d | --phase=A|B|C|DW (append-only CSV outputs)
"""

import csv
import json
import os
import sys
import time

import numpy as np

import riemann_tp_cylinder_sim as S
import riemann_transducer_sim as TS
import riemann_capacitor_sim as C

HD = S.HD
M = S.M
DT = S.DT
T_TOTAL = S.T_TOTAL
WINDOWS = S.WINDOWS
NIDX = S.NIDX
SEEDS8 = list(range(8))

G_TR, VSAT, EGATE, ETA, TAU_ENV = C.G_TR, C.VSAT, C.EGATE, C.ETA, C.TAU_ENV
TONE_G = C.TONE_G
F_BARE = C.F_BARE
KAPPA_C = float(C.KAPPA_C)
AIND_A = C.AIND_A_DEF
SER_EVERY = 20
_FS = 1.0 / (DT * S.REC_EVERY)

# ------------------------------------------------------------- 3D grid --
NZ = 16
_ZG = (np.arange(NZ) + 0.5) / NZ
_XG = (np.arange(48) + 0.5) / 48.0
M_AX = (NIDX % 3).astype(float)                    # axial wavenumbers
SIGMA_X, SIGMA_Z = 0.06, 0.12
Z_PLANES = [0.25, 0.5, 0.75]
Z_PLANES_DENSE = [0.125, 0.25, 0.5, 0.75, 0.875]


def chi(z, kz):
    """(len(z), M) axial envelopes; kz=0 -> all ones (collapsed)."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    if not kz:
        return np.ones((len(z), M))
    out = np.ones((len(z), M))
    mpos = M_AX > 0
    out[:, mpos] = np.sqrt(2.0) * np.cos(
        np.pi * np.outer(z, M_AX[mpos]))
    return out


_CHI_G = chi(_ZG, 1)                               # (NZ, M) kz=1 grid
PHI3 = (np.sin(np.pi * np.outer(_XG, NIDX))[:, None, :]
        * _CHI_G[None, :, :]).reshape(48 * NZ, M)  # (768, M)
_GS3 = 2.0 / (48 * NZ)


def profile3(x_c, theta, z_c, z_uniform, kz, normalize=False):
    """Element modal profile on the 3D grid (iter-6 conventions)."""
    if not kz:
        W = C.a_ind_profile(x_c, theta)
    else:
        bx = np.exp(-((_XG - x_c) ** 2) / (2.0 * SIGMA_X ** 2))
        if z_uniform:
            bz = np.ones(NZ)
        else:
            bz = np.exp(-((_ZG - z_c) ** 2) / (2.0 * SIGMA_Z ** 2))
        b = (bx[:, None] * bz[None, :]).reshape(-1)
        W = _GS3 * (PHI3.T @ b) * np.cos(NIDX * theta)
    if normalize:
        n = np.linalg.norm(W)
        return W / n if n > 0 else W
    return W


def transducer_W(kz, z_uniform, theta_deg=0.0, z_c=0.5):
    if not kz:
        return TS.SITES["A79"][1]
    bx = np.exp(-((_XG - 0.79) ** 2) / (2.0 * SIGMA_X ** 2))
    bz = (np.ones(NZ) if z_uniform else
          np.exp(-((_ZG - z_c) ** 2) / (2.0 * SIGMA_Z ** 2)))
    b = (bx[:, None] * bz[None, :]).reshape(-1)
    W = _GS3 * (PHI3.T @ b)
    if theta_deg:
        W = W * np.cos(NIDX * np.radians(theta_deg))
    return W


# ------------------------------------------------------------- pickups --
def pickup_mats(kz, dense=False):
    """Channel layout. Base: 40 ch (A_in 5x x 3z, A_out 3z, B 6,
    rings 16 at z=0.5). Dense (Phase E): 84 ch (A_in 5x x 5z,
    A_out 5z, B 6, rings 16 x 3z = 48)."""
    names, rowsA, rowsB = [], [], []
    zplanes = Z_PLANES_DENSE if dense else Z_PLANES
    sinA = np.sin(np.pi * np.outer(S.A_PICK_X, NIDX))
    for z in zplanes:
        cz = chi([z], kz)[0]
        for i in range(5):
            names.append(f"A_in{i+1}@z{z:.2f}")
            rowsA.append(sinA[i] * cz)
            rowsB.append(None)
    sinAP = np.sin(np.pi * S.A_APER_X * NIDX)
    for z in zplanes:
        cz = chi([z], kz)[0]
        names.append(f"A_out@z{z:.2f}")
        rowsA.append(sinAP * cz)
        rowsB.append(None)
    for i in range(5):
        names.append(f"B_in{i+1}")
        rowsA.append(None)
        rowsB.append(S.PHI_B_IN[i])
    names.append("B_out")
    rowsA.append(None)
    rowsB.append(S.PHI_B_AP)
    ring_z = zplanes if dense else [0.5]
    for z in ring_z:
        cz = chi([z], kz)[0]
        for k in range(8):
            names.append(f"A_rin{k}@z{z:.2f}" if dense else f"A_rin{k}")
            rowsA.append(C.RING_IN_A[k] * cz)
            rowsB.append(None)
        for k in range(8):
            names.append(f"A_rout{k}@z{z:.2f}" if dense else f"A_rout{k}")
            rowsA.append(C.RING_OUT_A[k] * cz)
            rowsB.append(None)
    return names, rowsA, rowsB


CH40, _RA, _RB = pickup_mats(1)
CH40_0, _RA0, _RB0 = pickup_mats(0)
CH84, _RAD, _RBD = pickup_mats(1, dense=True)
CH84_0, _RAD0, _RBD0 = pickup_mats(0, dense=True)
RD_DENSE_CH = ["A_rin0@z0.50", "A_rout0@z0.50", "A_in1@z0.50"]
BASE12Z = ([f"A_in{i+1}@z0.50" for i in range(5)] + ["A_out@z0.50"]
           + [f"B_in{i+1}" for i in range(5)] + ["B_out"])
RD_CH = ["A_in5@z0.50", "A_in5@z0.25", "A_out@z0.50", "B_in1"]

PROTOCOLS = {
    "BASE":     {"mod": "NONE", "ramp": 0.5, "t_off": 4.0},
    "STEP":     {"mod": "NONE", "ramp": 0.0, "t_off": 4.0},
    "SLOWRAMP": {"mod": "NONE", "ramp": 2.0, "t_off": 4.0},
    "SHORT":    {"mod": "NONE", "ramp": 0.5, "t_off": 2.0},
    "AM":       {"mod": "AM",   "ramp": 0.5, "t_off": 4.0},
    "FM":       {"mod": "FM",   "ramp": 0.5, "t_off": 4.0},
}
MOD_RATE = float(S.MEAN_GAP)


# ------------------------------------------------------------------ sim --
def run8(drive_freqs, seed, cfg):
    """cfg keys: kz, ind{x,theta_deg,z,z_uniform}|None,
    cap{x,theta_deg,z,z_uniform}|None, protocol{name,mod,ramp,t_off},
    a_reseed. Returns rec/exch/peaks-ready dict (iter-6 conventions
    + post-off ledger + optional ring-down capture)."""
    kz = int(cfg.get("kz", 1))
    proto = cfg["protocol"]
    mod_name, ramp_tau, t_off = proto["mod"], proto["ramp"], proto["t_off"]
    a_reseed = float(cfg.get("a_reseed", 0.0))
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
    _tr = cfg.get("tr") or {}
    W_T = transducer_W(kz, cfg.get("tr_z_uniform", False), _tr.get("theta_deg", 0.0), _tr.get("z", 0.5))
    ind = cfg.get("ind")
    cap = cfg.get("cap")
    W_I = (profile3(ind["x"], np.radians(ind["theta_deg"]), ind["z"],
                    ind.get("z_uniform", False), kz) if ind else None)
    cap_specs = list(cfg["caps"]) if cfg.get("caps") else \
        ([cap] if cap is not None else [])
    W_Cs = [profile3(cs["x"], np.radians(cs["theta_deg"]), cs["z"],
                     cs.get("z_uniform", False), kz, normalize=True)
            for cs in cap_specs]
    W_C = W_Cs[0] if W_Cs else None
    has_cap = len(W_Cs) > 0
    NEXT = 2 + len(W_Cs)
    PHI = PHI3 if kz else S.PHI_G_A
    GS = _GS3 if kz else (2.0 / 48.0)
    dense_pk = bool(cfg.get("dense_pickups", False))
    if dense_pk:
        names, rowsA, rowsB = (CH84, _RAD, _RBD) if kz else (CH84_0, _RAD0, _RBD0)
    else:
        names, rowsA, rowsB = (CH40, _RA, _RB) if kz else (CH40_0, _RA0, _RB0)
    i_off = int(round(t_off / DT))
    rng_r = np.random.default_rng(20_000 + seed) if a_reseed > 0 else None

    n_steps = int(round(T_TOTAL / DT))
    n_rec = int(round(T_TOTAL / (DT * S.REC_EVERY)))
    n_es = int(round(T_TOTAL / (DT * S.ESAMP_EVERY)))
    rec = np.empty((len(pickup_mats(kz, bool(cfg.get("dense_pickups", False)))[0]), n_rec))
    esamp = np.empty((n_es, 5))
    emode = np.empty((n_es, 2 * M))
    series = []
    inj_m = np.zeros(M); dissA = np.zeros(M); dissB = np.zeros(M)
    radA_m = np.zeros(M); radB_m = np.zeros(M)
    dAp = dBp = rAp = rBp = 0.0
    winE = {w: (np.zeros(M), np.zeros(M)) for w in WINDOWS}
    ri = ei = 0
    E0 = None
    E_ind = 0.0
    E_ext = E_reinj = E_tr_diss = 0.0
    E_cap_flow = 0.0
    sat_time = gate_time = 0.0
    clamp_count = 0
    clamp_E = 0.0
    pool_max = 0.0
    reseeded = False

    def rampf(t):
        return 1.0 if ramp_tau <= 0 else min(1.0, t / ramp_tau)

    def drive(t):
        if t >= t_off:
            return np.zeros(M)
        base = omd * t + phases
        r = rampf(t)
        if mod_name == "AM":
            return r * force * np.cos(base) * (
                1.0 + S.AM_DEPTH * np.cos(2.0 * np.pi * MOD_RATE * t))
        if mod_name == "FM":
            return r * force * np.cos(
                base + (S.FM_DEV / MOD_RATE)
                * np.sin(2.0 * np.pi * MOD_RATE * t))
        return r * force * np.cos(base)

    def ind_force(t):
        if ind is None:
            return None
        return AIND_A * W_I * np.cos(2.0 * np.pi * TONE_G * t) \
            * min(1.0, t / 0.5)

    def tr_forces(va, env, pool):
        s = float(W_T @ va)
        f1 = -G_TR * VSAT * np.tanh(s / VSAT) * W_T
        p_ext = -float(f1 @ va)
        sig = pool / (pool + EGATE) if pool > 0.0 else 0.0
        f2 = G_TR * sig * env * W_T
        p_reinj = float(f2 @ va)
        return f1 + f2, p_ext, p_reinj, s

    def rhs(y, t):
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        q = y[4 * M + 2] if has_cap else 0.0
        ua = PHI @ qa
        ub = S.S_B @ qb
        nla = S.BETA * (PHI.T @ (ua ** 3)) * GS
        nlb = S.BETA * (S.S_B.T @ (ub ** 3)) * (2.0 / 64.0)
        aa = -2.0 * gA * va - stiffA * qa + kappa * (qb - qa) + nla \
            + drive(t)
        ab = -2.0 * gB * vb - stiffB * qb + kappa * (qa - qb) + nlb
        F, p_ext, p_reinj, s = tr_forces(va, env, pool)
        aa = aa + F
        Fi = ind_force(t)
        if Fi is not None:
            aa = aa + Fi
        denv = (abs(s) - env) / TAU_ENV
        dpool = ETA * p_ext - p_reinj
        if has_cap:
            qs = y[4 * M + 2: 4 * M + 2 + len(W_Cs)]
            dqs = [KAPPA_C * float(Wc @ va) for Wc in W_Cs]
            fcap = np.zeros(M)
            for Wc, qv in zip(W_Cs, qs):
                fcap = fcap + KAPPA_C * qv * Wc
            aa = aa - fcap
            return np.concatenate([va, aa, vb, ab, [denv, dpool], dqs])
        return np.concatenate([va, aa, vb, ab, [denv, dpool]])

    def do_reseed(y):
        psiA = rng_r.uniform(0.0, 2.0 * np.pi, M)
        psiB = rng_r.uniform(0.0, 2.0 * np.pi, M)
        for sl, om, psi in ((slice(0, M), omA, psiA),
                            (slice(2 * M, 3 * M), omB, psiB)):
            q_ = y[sl]
            v_ = y[sl.start + M: sl.stop + M]
            a_re = om * q_
            a_im = v_
            mag = np.hypot(a_re, a_im)
            ang = np.arctan2(a_im, a_re)
            dpsi = (psi - ang + np.pi) % (2.0 * np.pi) - np.pi
            na = ang + a_reseed * dpsi
            y[sl] = mag * np.cos(na) / om
            y[sl.start + M: sl.stop + M] = mag * np.sin(na)
        return y

    y = np.zeros(4 * M + NEXT)
    if has_cap and cfg.get("q0"):
        y[4 * M + 2] = float(cfg["q0"])
    for i in range(n_steps):
        t = i * DT
        if a_reseed > 0 and not reseeded and i >= i_off:
            y = do_reseed(y)
            reseeded = True
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        qs = y[4 * M + 2: 4 * M + 2 + len(W_Cs)] if has_cap else np.zeros(0)
        q = float(qs[0]) if has_cap else 0.0
        dv = drive(t)
        inj_m += dv * va * DT
        dissA += 2.0 * (F_A / S.Q) * va ** 2 * DT
        dissB += 2.0 * (F_B / S.Q) * vb ** 2 * DT
        radA_m += 2.0 * gA_rad * va ** 2 * DT
        radB_m += 2.0 * gB_rad * vb ** 2 * DT
        if t >= t_off:
            dAp += float(np.sum(2.0 * (F_A / S.Q) * va ** 2)) * DT
            dBp += float(np.sum(2.0 * (F_B / S.Q) * vb ** 2)) * DT
            rAp += float(np.sum(2.0 * gA_rad * va ** 2)) * DT
            rBp += float(np.sum(2.0 * gB_rad * vb ** 2)) * DT
        Fi = ind_force(t)
        if Fi is not None:
            E_ind += float(Fi @ va) * DT
        F, p_ext, p_reinj, s_act = tr_forces(va, env, pool)
        E_ext += p_ext * DT
        E_reinj += p_reinj * DT
        E_tr_diss += (1.0 - ETA) * p_ext * DT
        if has_cap:
            E_cap_flow += float(sum(KAPPA_C * qv * float(Wc @ va)
                                    for Wc, qv in zip(W_Cs, qs))) * DT
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
            series.append({"t": round(t, 4), "E_A": float(eA_m.sum()),
                           "E_B": float(eB_m.sum()), "pool": float(pool),
                           "env": float(env), "s_act": float(s_act),
                           "E_ind_cum": float(E_ind),
                           "q_cap": float(q), "E_cap": float(0.5 * float(qs @ qs)) if has_cap else 0.0})
        if i % S.REC_EVERY == 0 and ri < n_rec:
            portA = np.sqrt(2.0 * gA_rad)
            portB = np.sqrt(2.0 * gB_rad)
            for ci in range(len(names)):
                if rowsA[ci] is not None:
                    if names[ci].startswith("A_out"):
                        rec[ci, ri] = np.sum(portA * va * rowsA[ci])
                    else:
                        rec[ci, ri] = rowsA[ci] @ qa
                else:
                    if names[ci] == "B_out":
                        rec[ci, ri] = np.sum(portB * vb * rowsB[ci])
                    else:
                        rec[ci, ri] = rowsB[ci] @ qb
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
    qs_end = y[4 * M + 2: 4 * M + 2 + len(W_Cs)] if has_cap else np.zeros(0)
    q_end = float(qs_end[0]) if has_cap else 0.0
    E_cap_end = 0.5 * float(qs_end @ qs_end) if has_cap else 0.0
    Eend = (0.5 * (va ** 2 + stiffA * qa ** 2).sum(),
            0.5 * (vb ** 2 + stiffB * qb ** 2).sum())
    es = esamp[:ei]
    st = (es[:, 0] >= WINDOWS["SETTLED"][0]) & (es[:, 0] < WINDOWS["SETTLED"][1])
    Es = es[st, 1] + es[st, 2]
    slope_settled = float(np.polyfit(es[st, 0],
                                     np.log(np.maximum(Es, 1e-12)), 1)[0]) \
        if st.sum() > 10 else float("nan")
    rms = np.sqrt(np.mean(rec ** 2, axis=1, keepdims=True)) + 1e-12
    rec = rec + rng.standard_normal(rec.shape) * S.NOISE_FRAC * rms

    exch = {
        "E_inj": float(inj_m.sum()), "E_ind": float(E_ind),
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
        "E_cap1_end": float(0.5 * qs_end[0] ** 2) if len(qs_end) > 0 else 0.0,
        "E_cap2_end": float(0.5 * qs_end[1] ** 2) if len(qs_end) > 1 else 0.0,
        "n_caps": int(len(W_Cs)),
        "slope_settled": slope_settled,
        "dissA_post": dAp, "dissB_post": dBp,
        "radA_post": rAp, "radB_post": rBp,
        "reseeded": bool(reseeded),
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
    tot_out = exch["E_diss_A"] + exch["E_diss_B"] \
        + exch["E_rad_A"] + exch["E_rad_B"]
    exch["B_share_out"] = (exch["E_diss_B"] + exch["E_rad_B"]) \
        / max(tot_out, 1e-12)
    tot_post = dAp + dBp + rAp + rBp
    exch["B_share_post"] = (dBp + rBp) / max(tot_post, 1e-12)
    exch["ind_frac"] = exch["E_ind"] / exch["E_inj"] \
        if exch["E_inj"] > 0 else float("nan")
    return {"rec": rec, "names": names, "esamp": esamp[:ei],
            "emode": emode[:ei], "winE": winE, "exch": exch,
            "series": series}


# ------------------------------------------------------------- readout --
def demod(series, f_tone, wname):
    w0, w1 = WINDOWS[wname]
    i0, i1 = int(w0 * _FS), int(w1 * _FS)
    x = series[i0:i1]
    tt = np.arange(i0, i1) / _FS
    z = np.mean(x * np.exp(-2j * np.pi * f_tone * tt)) * 2.0
    return float(abs(z)), float(np.degrees(np.angle(z)))


WIN_IND = {"x": 0.69, "theta_deg": 0.0, "z": 0.5}
WIN_CAP = {"x": 0.79, "theta_deg": 15.0, "z": 0.5}
PLAC_CAP = {"x": 0.13, "theta_deg": 0.0, "z": 0.5}

PEAKS_F = f"{HD}/riemann_reseed3d_peaks_raw.csv"
RUNS_F = f"{HD}/riemann_reseed3d_runs.csv"
STATE_F = f"{HD}/riemann_reseed3d_state.json"

PEAK_FIELDS = ["run_id", "tag", "cell", "arm", "seed", "cavity",
               "pickup", "pickup_kind", "window", "peak_freq",
               "peak_amp", "is_drive", "is_inducer", "audit_label"]


def _append(path, fields, rows):
    fields = list(fields)
    if os.path.exists(path):
        with open(path, newline="") as fh:
            rd = list(csv.DictReader(fh))
            old_fields = list(rd[0].keys()) if rd else fields
        if set(old_fields) != set(fields):
            fields = old_fields + [f for f in fields if f not in old_fields]
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields,
                                   extrasaction="ignore")
                w.writeheader()
                for r in rd:
                    w.writerow({k: r.get(k) for k in fields})
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


class Driver:
    def __init__(self):
        self.state = (json.load(open(STATE_F)) if os.path.exists(STATE_F)
                      else {"next_run_id": 1, "run_fields": None})
        self.t0 = time.time()
        self.rd = {}

    def save_state(self):
        json.dump(self.state, open(STATE_F, "w"), indent=1)

    def run(self, tag, cell_name, cfg, arm, seeds):
        freqs = S.ARMS[arm]
        proto = cfg["protocol"]
        names = ((CH84 if cfg.get("kz", 1) else CH84_0)
                 if cfg.get("dense_pickups") else
                 (CH40 if cfg.get("kz", 1) else CH40_0))
        i_dem = names.index("A_in5@z0.50")
        run_rows, peak_rows = [], []
        for seed in seeds:
            rid = self.state["next_run_id"]
            self.state["next_run_id"] += 1
            res = run8(freqs, seed, cfg)
            ex = res["exch"]
            ts, _dl = S.transient_scalars(res)
            de, _p1 = demod(res["rec"][i_dem], TONE_G, "DRIVE_END")
            st, ph_s = demod(res["rec"][i_dem], TONE_G, "SETTLED")
            ind = cfg.get("ind") or {}
            cap = cfg.get("cap") or {}
            cap2 = cfg.get("cap2") or {}
            tr = cfg.get("tr") or {}
            row = {"run_id": rid, "tag": tag, "cell": cell_name,
                   "arm": arm, "seed": seed, "kz": cfg.get("kz", 1),
                   "protocol": proto["name"], "mod": proto["mod"],
                   "ramp": proto["ramp"], "t_off": proto["t_off"],
                   "a_reseed": cfg.get("a_reseed", 0.0),
                   "ind_x": ind.get("x"), "ind_theta_deg": ind.get("theta_deg"),
                   "ind_z": ind.get("z"),
                   "cap_x": cap.get("x"), "cap_theta_deg": cap.get("theta_deg"),
                   "cap_z": cap.get("z"),
                   "cap2_x": cap2.get("x"),
                   "cap2_theta_deg": cap2.get("theta_deg"),
                   "cap2_z": cap2.get("z"),
                   "tr_theta_deg": tr.get("theta_deg", 0.0),
                   "dense_pickups": bool(cfg.get("dense_pickups", False)),
                   **ex, **ts,
                   "demod_DE_amp": de, "demod_ST_amp": st,
                   "demod_ST_phase": ph_s,
                   "P": (st / de if de > 1e-9 else float("nan"))}
            run_rows.append(row)
            for ci, ch in enumerate(names):
                cav = "B" if ch.startswith("B_") else "A"
                kind = ("out" if "out" in ch else
                        ("ring" if "rin" in ch or "rout" in ch else "in"))
                for wname in WINDOWS:
                    for pf, pa in S.window_peaks(res["rec"][ci], wname):
                        is_ind = bool(abs(pf - TONE_G) <= 0.15)
                        is_drv = (not is_ind) and bool(
                            np.min(np.abs(freqs - pf)) <= S.PEAK_TOL)
                        peak_rows.append({
                            "run_id": rid, "tag": tag, "cell": cell_name,
                            "arm": arm, "seed": seed, "cavity": cav,
                            "pickup": ch, "pickup_kind": kind,
                            "window": wname, "peak_freq": round(pf, 4),
                            "peak_amp": round(pa, 6), "is_drive": is_drv,
                            "is_inducer": is_ind,
                            "audit_label": "inducer" if is_ind else
                            ("drive" if is_drv else
                             S.audit_label(pf, freqs, proto["mod"],
                                           MOD_RATE if proto["mod"] != "NONE"
                                           else 0.0))})
            # ring-down capture (t >= 3.9 s) for Phase D
            i39 = int(3.9 * _FS)
            key = f"{cell_name}|{proto['name']}|{arm}|A{cfg.get('a_reseed', 0.0)}|s{seed}"
            for chn in RD_CH:
                self.rd[f"{key}|{chn}"] = \
                    res["rec"][names.index(chn), i39:].astype(np.float32)
            if cfg.get("dense_pickups"):
                for chn in RD_DENSE_CH:
                    if chn in names:
                        self.rd[f"{key}|{chn}"] = res["rec"][
                            names.index(chn), i39:].astype(np.float32)
            print(f"run {rid} {tag} {cell_name} {arm} s{seed} "
                  f"A={cfg.get('a_reseed', 0.0)} "
                  f"closure={ex['closure_frac']:.5f} "
                  f"ind_frac={ex['ind_frac']:.5f} P={row['P']:.3f} "
                  f"elapsed={time.time() - self.t0:.0f}s", flush=True)
        if self.state["run_fields"] is None:
            self.state["run_fields"] = []
        for r in run_rows:
            for k in r.keys():
                if k not in self.state["run_fields"]:
                    self.state["run_fields"].append(k)
        _append(RUNS_F, self.state["run_fields"], run_rows)
        _append(PEAKS_F, PEAK_FIELDS, peak_rows)
        self.save_state()
        return run_rows

    def dump_rd(self, phase):
        if self.rd:
            np.savez_compressed(
                f"{HD}/riemann_reseed3d_rd_{phase}.npz", **self.rd)
            self.rd = {}


def cfg_of(kz, ind, cap, proto_name, a_reseed=0.0, tr_zunif=False,
           caps=None, tr=None, dense_pickups=False):
    d = {"kz": kz, "ind": ind, "cap": cap,
         "protocol": {"name": proto_name, **PROTOCOLS[proto_name]},
         "a_reseed": a_reseed, "tr_z_uniform": tr_zunif,
         "caps": caps, "tr": tr, "dense_pickups": dense_pickups}
    if caps:
        d["cap"] = caps[0]
        d["cap2"] = caps[1] if len(caps) > 1 else None
    return d


def zunif(d):
    return {**d, "z_uniform": True} if d else None


# --------------------------------------------------------------- phases --
def phase_A(drv):
    b = "BASE"
    drv.run("A", "G0_WIN", cfg_of(0, zunif(WIN_IND), zunif(WIN_CAP), b),
            "LADDER", [0, 1, 2])
    drv.run("A", "G0_OFF", cfg_of(0, zunif(WIN_IND), None, b),
            "LADDER", [0, 1, 2])
    drv.run("A", "UNIF_WIN", cfg_of(1, zunif(WIN_IND), zunif(WIN_CAP), b,
                                    tr_zunif=True), "LADDER", [0, 1, 2])
    for z in (0.1, 0.3, 0.5, 0.7, 0.9):
        drv.run("A", f"CAPZ_{z:.1f}",
                cfg_of(1, WIN_IND, {**WIN_CAP, "z": z}, b),
                "LADDER", [0, 1, 2])
    for dd in (0.10, 0.18, 0.24, 0.28, 0.40):
        drv.run("A", f"TRAP_d{dd:.2f}",
                cfg_of(1, {**WIN_IND, "x": round(0.79 - dd, 4)}, WIN_CAP, b),
                "LADDER", [0, 1, 2])
    for dd, zz in ((0.10, 0.25), (0.10, 0.75), (0.24, 0.25), (0.24, 0.75)):
        drv.run("A", f"TRAP_d{dd:.2f}_z{zz:.2f}",
                cfg_of(1, {**WIN_IND, "x": round(0.79 - dd, 4), "z": zz},
                       WIN_CAP, b), "LADDER", [0, 1, 2])


def phase_B(drv):
    for a in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
        drv.run("B", "WIN3D", cfg_of(1, WIN_IND, WIN_CAP, "BASE", a),
                "LADDER", SEEDS8)
    for arm in ("DETUNE", "SHIFT"):
        for a in (0.0, 1.0):
            drv.run("B", "WIN3D", cfg_of(1, WIN_IND, WIN_CAP, "BASE", a),
                    arm, [0, 1, 2, 3])
    for a in (0.0, 1.0):
        drv.run("B", "PLAC3D", cfg_of(1, WIN_IND, PLAC_CAP, "BASE", a),
                "LADDER", [0, 1, 2, 3])
    drv.dump_rd("B")


def phase_C(drv):
    for pname in ("STEP", "SLOWRAMP", "SHORT", "AM", "FM"):
        drv.run("C", f"PROTO_{pname}",
                cfg_of(1, WIN_IND, WIN_CAP, pname, 0.0), "LADDER", SEEDS8)
    drv.dump_rd("C")



# ------------------------------------------------------- Phase E (scope
# addition, Zak 2026-10-22 ~23:54 EDT; pre-stated in prereg addendum
# riemann_reseed3d_prereg.json["phaseE_addendum"] before any E run) --
A_VALS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
ANTI_TR = {"theta_deg": 180.0, "z": 0.5}
ANTI_IND = {"x": 0.69, "theta_deg": 180.0, "z": 0.5}
ANTI_CAP = {"x": 0.79, "theta_deg": 195.0, "z": 0.5}
CAP2_ZOFF = {"x": 0.79, "theta_deg": 15.0, "z": 0.30}
IND_TRAP_Z = {"x": 0.55, "theta_deg": 0.0, "z": 0.25}


def phase_E(drv):
    b = "BASE"
    # E1 adjusted inducer: trap standoff d=0.24 AND z-offset 0.25,
    # full re-seed sweep x 8 seeds (A_max picture vs WIN3D)
    for a in A_VALS:
        drv.run("E", "E_INDTRAP", cfg_of(1, IND_TRAP_Z, WIN_CAP, b, a,
                                          dense_pickups=True),
                "LADDER", SEEDS8)
    # E2 two capacitors: winner cap + same (x,theta) partner at
    # z=0.30 (axial offset -0.20); A in {0,1} x 8 seeds
    for a in (0.0, 1.0):
        drv.run("E", "E_TWOCAP",
                cfg_of(1, WIN_IND, None, b, a,
                       caps=[WIN_CAP, CAP2_ZOFF], dense_pickups=True),
                "LADDER", SEEDS8)
    # E3 anti-mirror layout: (x,theta)->(x,theta+180), z->1-z applied
    # to transducer, inducer, capacitor (all z=0.5 -> 0.5)
    for a in (0.0, 1.0):
        drv.run("E", "E_ANTI",
                cfg_of(1, ANTI_IND, ANTI_CAP, b, a, tr=ANTI_TR,
                       dense_pickups=True), "LADDER", SEEDS8)
    for arm in ("DETUNE", "SHIFT"):
        drv.run("E", "E_ANTI",
                cfg_of(1, ANTI_IND, ANTI_CAP, b, 0.0, tr=ANTI_TR,
                       dense_pickups=True), arm, [0, 1, 2, 3])
    for dd in (0.10, 0.24):
        drv.run("E", "E_ANTITRAP",
                cfg_of(1, {**ANTI_IND, "x": round(0.79 - dd, 4)},
                       ANTI_CAP, b, 0.0, tr=ANTI_TR, dense_pickups=True),
                "LADDER", [0, 1, 2])
    # E4 dense-array standard winner (Phase D dense follow)
    for a in (0.0, 1.0):
        drv.run("E", "E_DENSEWIN",
                cfg_of(1, WIN_IND, WIN_CAP, b, a, dense_pickups=True),
                "LADDER", SEEDS8)
    drv.dump_rd("E")


def phase_EDMAX(drv, a):
    drv.run("E", "E_DENSEWIN",
            cfg_of(1, WIN_IND, WIN_CAP, "BASE", a, dense_pickups=True),
            "LADDER", SEEDS8)
    drv.dump_rd("EDMAX")


def phase_DW(drv, proto_name):
    drv.run("DW", f"PROTO_{proto_name}",
            cfg_of(1, WIN_IND, WIN_CAP, proto_name, 1.0),
            "LADDER", SEEDS8)
    drv.dump_rd("DW")


def ring3d():
    """3D capacitor ring test: q(0)=300, drive silenced (t_off=0),
    no inducer, cap at (0.69, 0deg, z=0.5) -> in-situ f_ring."""
    cfg = {"kz": 1, "ind": None,
           "cap": {"x": 0.69, "theta_deg": 0.0, "z": 0.5},
           "protocol": {"name": "RING", "mod": "NONE", "ramp": 0.5,
                        "t_off": 0.0},
           "a_reseed": 0.0, "q0": 300.0}
    res = run8(S.T_LADDER, 0, cfg)
    ser = res["series"]
    tt = np.array([r["t"] for r in ser])
    qq = np.array([r["q_cap"] for r in ser])
    m = (tt >= 1.0) & (tt <= 7.5)
    x = qq[m] - qq[m].mean()
    fs = 1.0 / (tt[1] - tt[0])
    spec = np.abs(np.fft.rfft(x * np.hanning(len(x))))
    fr = np.fft.rfftfreq(len(x), d=1.0 / fs)
    band = (fr >= 2.0) & (fr <= 80.0)
    f_ring = float(fr[band][np.argmax(spec[band])])
    out = {"f_bare": F_BARE, "f_ring_3d": f_ring,
           "q_max": float(np.max(np.abs(qq))), "q_end": float(qq[-1]),
           "closure_frac": res["exch"]["closure_frac"]}
    json.dump(out, open(f"{HD}/riemann_reseed3d_ring.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def main():
    args = sys.argv[1:]
    if "--ring3d" in args:
        ring3d()
        return
    drv = Driver()
    for a in args:
        if a.startswith("--phase="):
            ph = a.split("=", 1)[1]
            if ph == "A":
                phase_A(drv)
            elif ph == "B":
                phase_B(drv)
            elif ph == "C":
                phase_C(drv)
            elif ph == "E":
                phase_E(drv)
            elif ph.startswith("EDMAX:"):
                phase_EDMAX(drv, float(ph.split(":", 1)[1]))
            elif ph.startswith("DW:"):
                phase_DW(drv, ph.split(":", 1)[1])
    print(f"DRIVER DONE elapsed={time.time() - drv.t0:.0f}s")


if __name__ == "__main__":
    main()
