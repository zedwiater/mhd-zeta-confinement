"""
riemann_tp_cylinder_sim.py — RIEMANN TP + CYLINDER-FLOW DUAL CAVITY (iter 2)
(Zak's order 2026-10-22 14:35 EDT. Pre-registration: ITERATION_LOG.md
2026-10-22 14:38 EDT — geometry, grid and gates stated BEFORE this run.)

WHAT CHANGED vs run 1 (riemann_dual_cavity_sim.py): the two generic
cavities are replaced by (A) a TP cavity — toroidal modes with a poloidal
flow — and (B) a cylinder-flow cavity — cylindrical Bessel-root modes with
an axial flow. Pickups now sit INSIDE and OUTSIDE each cavity (the outside
channel is a radiation-port readout). The run records transients (drive
attack, drive-end, ring-down) and a full energy ledger, not just a steady
window. Source is the Riemann ladder, same controls as run 1.

GEOMETRY / APPROXIMATIONS (stated, per pre-registration):
- A (TP): toroidal modes n=1..30, shapes sin(n*pi*x) on the unwrapped
  toroidal coordinate x in (0,1) (run-1 basis). Plasma-frame eigenfreqs =
  the 30-zero ladder. POLOIDAL FLOW = rigid poloidal rotation, first
  order: lab-frame resonance shifts by df_n = FP * n (FP = 0.025). No
  flow shear, no flow-induced harmonic coupling (neglected, stated).
- B (cylinder): modes (m, nr, p), m=0..4, nr=1..4, p=0..3,
  f_raw = sqrt(j_{m,nr}^2 + (p*pi)^2) — Bessel-root modal approximation,
  a = L = 1 (labelled as such). Lowest 30 taken, rescaled so the set
  median equals the ladder median. Shapes J_m(j_{m,nr} r) cos(m*theta)
  cos(p*pi*z), RMS-normalized over the NL sample points to 1/sqrt(2)
  (the sine-basis convention). AXIAL FLOW = co-flow Doppler shift
  df/f = FA * axial_share (FA = 0.02, axial_share = (p*pi)^2 / f_raw^2);
  counter-flow branch neglected (stated).
- Arms change the DRIVE (lab antenna) set only: LADDER / DETUNE x1.01 /
  SHIFT +mean_gap/2. The cavity pair is identical in every arm.
- Coupling A<->B per mode index, stiffness 2% (run-1 convention — a modal
  approximation of aperture coupling). Cubic mixing per cavity in its own
  basis, force = BETA * (S^T u^3) * (2/N) (run-1 formula; for B, S is the
  RMS-normalized Bessel shape matrix on 64 volume sample points).
- Radiation port: gamma_rad = RAD_FRAC * gamma_m extra damping on every
  mode; modal radiated energy exact (int 2*gamma_rad v^2 dt). Exterior
  pickup = port series sum_m sqrt(2*gamma_rad,m) v_m phi_m(aperture).
- Ledger: E_inj (drive work, A), E_diss (viscous), E_rad (port),
  E_stored (modal, cubic potential neglected — small at these amplitudes;
  closure residual reported, gate G2).

PROTOCOL: drive on 0 -> T_OFF=4.0 (0.5 ramp), off, record to 8.0.
Windows: EARLY [0.6,2.0], DRIVE_END [2.2,4.0], SETTLED [5.2,7.6].
Grid: 3 arms x 3 modulations (NONE, AM@mean_gap d=0.5, FM@mean_gap
dev=0.5) x 5 seeds = 45 runs. Q=30, BETA=10, A_TARGET=1.2, DT=0.0004,
RK4 — run-1 values. Fixed per-mode detune rng(7) (run-1 draw for A,
extended for B). SCALE-FREE: ladder t values are dimensionless freqs.
Sim-internal claims only.
"""

import csv
import itertools
import json
import sys
import time

import numpy as np
from scipy.signal import find_peaks
from scipy.signal.windows import blackmanharris
from scipy.special import jn_zeros, jv

HD = "/home/hatch/workspace/harness-dev"
LADDER_JSON = "/home/hatch/workspace/riemann-tp/edge_frequencies.json"

# ------------------------------------------------------------ parameters --
M = 30
Q = 30.0
KAPPA_FRAC = 0.02
BETA = 10.0
A_TARGET = 1.2
RAD_FRAC = 0.25              # radiation port share of modal damping
FP = 0.025                   # TP poloidal flow: df_n = FP * n (pre-reg)
for _a in sys.argv:          # flow-sweep override (log 2026-10-22 14:47 EDT)
    if _a.startswith("--fp="):
        FP = float(_a.split("=")[1])
FA_FLOW = 0.02               # cylinder axial flow: df/f = FA * axial_share
DT = 0.0004
T_OFF = 4.0                  # drive off
T_TOTAL = 8.0
REC_EVERY = 2                # channel record -> fs = 1250
ESAMP_EVERY = 10             # energy sample every 10 steps
NOISE_FRAC = 0.002
PEAK_TOL = 0.45
SEEDS = [0, 1, 2, 3, 4]
WINDOWS = {"EARLY": (0.6, 2.0), "DRIVE_END": (2.2, 4.0), "SETTLED": (5.2, 7.6)}
for _a in sys.argv:          # fine-window pass override (log 14:47 EDT)
    if _a.startswith("--windows="):
        WINDOWS = {}
        for _part in _a.split("=", 1)[1].split(";"):
            _nm, _rng = _part.split(":")
            _a0, _b0 = _rng.split(",")
            WINDOWS[_nm] = (float(_a0), float(_b0))
AM_DEPTH = 0.5
FM_DEV = 0.5

_lad = json.load(open(LADDER_JSON))
T_LADDER = np.array([z["t"] for z in _lad["zeros"]], dtype=float)
MEAN_GAP = float(_lad["mean_gap"])

ARMS = {
    "LADDER": T_LADDER.copy(),
    "DETUNE": T_LADDER * 1.01,
    "SHIFT": T_LADDER + MEAN_GAP / 2.0,
}
MODS = [("NONE", 0.0), ("AM", MEAN_GAP), ("FM", MEAN_GAP)]

# fixed per-mode detune (run-1 convention): first 30 draws = cavity A
_det = np.random.default_rng(7).uniform(-0.0015, 0.0015, 2 * M)
DET_A, DET_B = _det[:M], _det[M:]

# ------------------------------------------------------- cavity A (TP) --
NIDX = np.arange(1, M + 1)
F_EIG_A = (T_LADDER + FP * NIDX) * (1.0 + DET_A)     # lab-frame resonances
A_PICK_X = [0.13, 0.29, 0.47, 0.61, 0.79]
A_APER_X = 0.97
_xgA = (np.arange(48) + 0.5) / 48.0
PHI_G_A = np.sin(np.pi * np.outer(_xgA, NIDX))                    # (48, M)
PHI_A_IN = np.sin(np.pi * np.outer(A_PICK_X, NIDX))               # (5, M)
PHI_A_AP = np.sin(np.pi * A_APER_X * NIDX)                        # (M,)

# ------------------------------------------------- cavity B (cylinder) --
def build_cylinder():
    combos = []
    for m_ in range(0, 5):
        roots = jn_zeros(m_, 4)
        for nr in range(1, 5):
            j = roots[nr - 1]
            for p in range(0, 4):
                f_raw = np.sqrt(j ** 2 + (p * np.pi) ** 2)
                combos.append((f_raw, m_, nr, p, j))
    combos.sort(key=lambda c: c[0])
    sel = combos[:M]
    raw = np.array([c[0] for c in sel])
    scale = np.median(T_LADDER) / np.median(raw)
    freqs = raw * scale
    axial_share = np.array([(c[3] * np.pi) ** 2 / c[0] ** 2 for c in sel])
    freqs = freqs * (1.0 + FA_FLOW * axial_share) * (1.0 + DET_B)
    return sel, freqs, axial_share, scale

B_MODES, F_EIG_B, B_AXIAL_SHARE, B_SCALE = build_cylinder()

def bessel_shape(m_, j, p, r, th, z):
    return jv(m_, j * r) * np.cos(m_ * th) * (1.0 if p == 0 else np.cos(p * np.pi * z))

# NL sample points: 64 fixed volume samples (rng(11))
_rngs = np.random.default_rng(11)
_s_r = 0.98 * np.sqrt(_rngs.uniform(0, 1, 64))
_s_th = _rngs.uniform(0, 2 * np.pi, 64)
_s_z = _rngs.uniform(0, 1, 64)
S_B = np.column_stack([bessel_shape(c[1], c[4], c[3], _s_r, _s_th, _s_z)
                       for c in B_MODES])                        # (64, M)
_rms = np.sqrt(np.mean(S_B ** 2, axis=0))
S_B = S_B / _rms * (1.0 / np.sqrt(2.0))                           # RMS -> 1/sqrt(2)
B_SCALE_NORM = _rms * np.sqrt(2.0)        # physical shape = normalized / this

B_PICK = [(0.35, 0.4, 0.3), (0.6, 1.5, 0.6), (0.85, 2.8, 0.45),
          (0.5, 4.2, 0.7), (0.75, 5.5, 0.25)]
B_APER = (0.9, 0.0, 0.5)
PHI_B_IN = np.array([[bessel_shape(c[1], c[4], c[3], *pt) for c in B_MODES]
                     for pt in B_PICK]) / B_SCALE_NORM           # (5, M)
PHI_B_AP = np.array([bessel_shape(c[1], c[4], c[3], *B_APER)
                     for c in B_MODES]) / B_SCALE_NORM           # (M,)

CH_NAMES = ([f"A_in{i+1}" for i in range(5)] + ["A_out"] +
            [f"B_in{i+1}" for i in range(5)] + ["B_out"])

# ------------------------------------------------------------------ sim --
def run_once(drive_freqs, mod_name, mod_rate, seed):
    rng = np.random.default_rng(10_000 + seed)
    phases = rng.uniform(0.0, 2.0 * np.pi, M)
    omA, omB = 2.0 * np.pi * F_EIG_A, 2.0 * np.pi * F_EIG_B
    gA = F_EIG_A / Q * (1.0 + RAD_FRAC)          # total damping per mode
    gB = F_EIG_B / Q * (1.0 + RAD_FRAC)
    gA_rad = F_EIG_A / Q * RAD_FRAC
    gB_rad = F_EIG_B / Q * RAD_FRAC
    stiffA, stiffB = omA ** 2, omB ** 2
    kappa = KAPPA_FRAC * stiffA
    omd = 2.0 * np.pi * drive_freqs
    force = 2.0 * (drive_freqs / Q) * omd * A_TARGET

    n_steps = int(round(T_TOTAL / DT))
    n_rec = int(round(T_TOTAL / (DT * REC_EVERY)))
    n_es = int(round(T_TOTAL / (DT * ESAMP_EVERY)))
    rec = np.empty((12, n_rec))
    esamp = np.empty((n_es, 5))                  # t, E_A, E_B, P_AB, P_BA
    emode = np.empty((n_es, 2 * M))              # per-mode energies
    inj_m = np.zeros(M); dissA = np.zeros(M); dissB = np.zeros(M)
    radA_m = np.zeros(M); radB_m = np.zeros(M)
    winE = {w: (np.zeros(M), np.zeros(M)) for w in WINDOWS}   # sum (v^2+w^2 q^2) dt
    ri = ei = 0
    E0 = None

    def drive(t):
        if t >= T_OFF:
            return np.zeros(M)
        base = omd * t + phases
        ramp = min(1.0, t / 0.5)
        if mod_name == "AM":
            return ramp * force * np.cos(base) * (1.0 + AM_DEPTH * np.cos(2.0 * np.pi * mod_rate * t))
        if mod_name == "FM":
            return ramp * force * np.cos(base + (FM_DEV / mod_rate) * np.sin(2.0 * np.pi * mod_rate * t))
        return ramp * force * np.cos(base)

    def rhs(y, t):
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:]
        ua = PHI_G_A @ qa
        ub = S_B @ qb
        nla = BETA * (PHI_G_A.T @ (ua ** 3)) * (2.0 / 48.0)
        nlb = BETA * (S_B.T @ (ub ** 3)) * (2.0 / 64.0)
        aa = -2.0 * gA * va - stiffA * qa + kappa * (qb - qa) + nla + drive(t)
        ab = -2.0 * gB * vb - stiffB * qb + kappa * (qa - qb) + nlb
        return np.concatenate([va, aa, vb, ab])

    y = np.zeros(4 * M)
    for i in range(n_steps):
        t = i * DT
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:]
        # bookkeeping at step start (rectangle rule; closure gate absorbs it)
        dv = drive(t)
        inj_m += dv * va * DT
        dissA += 2.0 * (F_EIG_A / Q) * va ** 2 * DT
        dissB += 2.0 * (F_EIG_B / Q) * vb ** 2 * DT
        radA_m += 2.0 * gA_rad * va ** 2 * DT
        radB_m += 2.0 * gB_rad * vb ** 2 * DT
        eA_m = 0.5 * (va ** 2 + stiffA * qa ** 2)
        eB_m = 0.5 * (vb ** 2 + stiffB * qb ** 2)
        if E0 is None:
            E0 = (eA_m.sum(), eB_m.sum())
        for wname, (w0, w1) in WINDOWS.items():
            if w0 <= t < w1:
                winE[wname][0][:] += (va ** 2 + stiffA * qa ** 2) * DT
                winE[wname][1][:] += (vb ** 2 + stiffB * qb ** 2) * DT
        if i % ESAMP_EVERY == 0 and ei < n_es:
            pab = float(np.sum(kappa * (qa - qb) * vb))
            pba = float(np.sum(kappa * (qb - qa) * va))
            esamp[ei] = (t, eA_m.sum(), eB_m.sum(), pab, pba)
            emode[ei] = np.concatenate([eA_m, eB_m])
            ei += 1
        if i % REC_EVERY == 0 and ri < n_rec:
            rec[:5, ri] = PHI_A_IN @ qa
            rec[5, ri] = np.sum(np.sqrt(2.0 * gA_rad) * va * PHI_A_AP)
            rec[6:11, ri] = PHI_B_IN @ qb
            rec[11, ri] = np.sum(np.sqrt(2.0 * gB_rad) * vb * PHI_B_AP)
            ri += 1
        k1 = rhs(y, t)
        k2 = rhs(y + 0.5 * DT * k1, t + 0.5 * DT)
        k3 = rhs(y + 0.5 * DT * k2, t + 0.5 * DT)
        k4 = rhs(y + DT * k3, t + DT)
        y = y + (DT / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:]
    Eend = (0.5 * (va ** 2 + stiffA * qa ** 2).sum(),
            0.5 * (vb ** 2 + stiffB * qb ** 2).sum())
    rms = np.sqrt(np.mean(rec ** 2, axis=1, keepdims=True)) + 1e-12
    rec = rec + rng.standard_normal(rec.shape) * NOISE_FRAC * rms

    exch = {
        "E_inj": float(inj_m.sum()), "E_diss_A": float(dissA.sum()),
        "E_diss_B": float(dissB.sum()), "E_rad_A": float(radA_m.sum()),
        "E_rad_B": float(radB_m.sum()),
        "E_stored0_A": float(E0[0]), "E_stored0_B": float(E0[1]),
        "E_stored_end_A": float(Eend[0]), "E_stored_end_B": float(Eend[1]),
    }
    exch["closure_resid"] = exch["E_inj"] - (
        (exch["E_stored_end_A"] - exch["E_stored0_A"])
        + (exch["E_stored_end_B"] - exch["E_stored0_B"])
        + exch["E_diss_A"] + exch["E_diss_B"]
        + exch["E_rad_A"] + exch["E_rad_B"])
    exch["closure_frac"] = (abs(exch["closure_resid"]) / exch["E_inj"]
                            if exch["E_inj"] > 0 else float("nan"))
    return {"rec": rec, "esamp": esamp[:ei], "emode": emode[:ei],
            "inj_m": inj_m, "radA_m": radA_m, "radB_m": radB_m,
            "winE": winE, "exch": exch}


# --------------------------------------------------------------- readout --
_FS = 1.0 / (DT * REC_EVERY)

def window_peaks(series, wname):
    w0, w1 = WINDOWS[wname]
    i0, i1 = int(w0 * _FS), int(w1 * _FS)
    x = series[i0:i1] - np.mean(series[i0:i1])
    w = blackmanharris(len(x))
    spec = np.abs(np.fft.rfft(x * w)) * (2.0 / np.sum(w))
    freqs = np.fft.rfftfreq(len(x), d=1.0 / _FS)
    med = np.median(spec) + 1e-15
    floor = max(8.0 * med, 0.002 * np.max(spec))
    idx, _ = find_peaks(spec, height=floor, prominence=4.0 * med, distance=3)
    out = []
    for i in idx:
        if 0 < i < len(spec) - 1:
            a, b, c = spec[i - 1], spec[i], spec[i + 1]
            d = 0.5 * (a - c) / (a - 2.0 * b + c + 1e-30)
        else:
            d = 0.0
        out.append((float(freqs[i] + d * (freqs[1] - freqs[0])), float(spec[i])))
    return out

def audit_label(fc, drive_freqs, mod_name, mod_rate, tol=PEAK_TOL):
    F = drive_freqs
    def hit(target):
        return target > 0.5 and abs(fc - target) <= tol
    for fi in F:
        if hit(2.0 * fi):
            return "2fi"
        if hit(3.0 * fi):
            return "3fi"
        if hit(fi / 2.0):
            return "fi/2"
    for fi, fj in itertools.combinations(F, 2):
        if hit(fi + fj):
            return "fi+fj"
        if hit(abs(fi - fj)):
            return "fi-fj"
        if hit(2.0 * fi + fj) or hit(2.0 * fj + fi):
            return "2fi+fj"
        if hit(abs(2.0 * fi - fj)) or hit(abs(2.0 * fj - fi)):
            return "2fi-fj"
    if mod_name in ("AM", "FM") and mod_rate > 0:
        for fi in F:
            for k in (1, 2):
                if hit(fi + k * mod_rate) or hit(fi - k * mod_rate):
                    return f"sideband_fi+-{k}r"
    return "UNEXPLAINED"


# ------------------------------------------------------------- transients --
def transient_scalars(res):
    """Per-cavity t50 + settled ring-down tau from the energy samples."""
    es = res["esamp"]
    t = es[:, 0]
    out = {}
    de = (t >= WINDOWS["DRIVE_END"][0]) & (t < WINDOWS["DRIVE_END"][1])
    st = (t >= WINDOWS["SETTLED"][0]) & (t < WINDOWS["SETTLED"][1])
    for ci, cav in enumerate(("A", "B")):
        E = es[:, 1 + ci]
        ref = E[de].mean()
        above = np.where(E >= 0.5 * ref)[0]
        out[f"t50_{cav}"] = float(t[above[0]]) if len(above) else float("nan")
        Es = E[st]; ts = t[st]
        ok = Es > 1e-9
        if ok.sum() > 10:
            slope = np.polyfit(ts[ok], np.log(Es[ok]), 1)[0]
            out[f"tau_{cav}"] = float(-1.0 / slope) if slope < 0 else float("nan")
        else:
            out[f"tau_{cav}"] = float("nan")
    # per-mode delay = t50 difference (pilot adjustment, log 14:44 EDT:
    # energy cross-correlation degenerates on monotone attack envelopes).
    # t50_m = first t (<= T_OFF) where modal energy crosses 50% of its
    # DRIVE_END-window mean; delay = t50_B - t50_A per mode index.
    em = res["emode"]
    tt = es[:, 0]
    de_m = (tt >= WINDOWS["DRIVE_END"][0]) & (tt < WINDOWS["DRIVE_END"][1])
    atk_i = np.where(tt <= T_OFF)[0]
    delays = []
    for m_ in range(M):
        t50s = []
        for sl in (slice(0, M), slice(M, 2 * M)):
            E = em[:, sl][:, m_]
            ref = E[de_m].mean()
            hit = atk_i[E[atk_i] >= 0.5 * ref] if ref > 0 else []
            t50s.append(float(tt[hit[0]]) if len(hit) else float("nan"))
        delays.append(t50s[1] - t50s[0]
                      if t50s[0] == t50s[0] and t50s[1] == t50s[1]
                      else float("nan"))
    return out, delays


# --------------------------------------------------------------- driver --
def main():
    pilot = "--pilot" in sys.argv
    outprefix = "riemann_tp_cylinder"
    cells_want = None
    for _a in sys.argv:
        if _a.startswith("--outprefix="):
            outprefix = _a.split("=", 1)[1]
        if _a.startswith("--cells="):
            cells_want = {tuple(p.split(":")) for p in _a.split("=", 1)[1].split(",")}
    t0 = time.time()
    if pilot:
        res = run_once(ARMS["LADDER"], "NONE", 0.0, 0)
        print("B modes (m,nr,p) f_eig_B:")
        for i, c in enumerate(B_MODES):
            print(f"  B{i+1}: m={c[1]} nr={c[2]} p={c[3]} f={F_EIG_B[i]:.2f} "
                  f"axial_share={B_AXIAL_SHARE[i]:.2f}")
        print(f"B scale={B_SCALE:.3f}  A eig range {F_EIG_A[0]:.2f}..{F_EIG_A[-1]:.2f}")
        print("exchange:", json.dumps(res["exch"], indent=1))
        ts, delays = transient_scalars(res)
        print("transients:", ts)
        print("delays (s), first 10 modes:", [round(d, 3) for d in delays[:10]])
        wA, wB = res["winE"]["DRIVE_END"]
        print("DRIVE_END modal energy A (per mode):",
              [round(float(x), 1) for x in wA])
        print("DRIVE_END modal energy B (per mode):",
              [round(float(x), 1) for x in wB])
        pk = window_peaks(res["rec"][0], "DRIVE_END")
        print(f"A_in1 DRIVE_END peaks: {len(pk)}; top:",
              [(round(f, 2), round(a, 3)) for f, a in
               sorted(pk, key=lambda z: -z[1])[:8]])
        print(f"elapsed={time.time()-t0:.1f}s")
        return

    peak_rows, tran_rows, es_rows, meta = [], [], [], {"runs": []}
    run_id = 0
    g1_hits = g1_tot = 0
    for arm, freqs in ARMS.items():
        for mod_name, mod_rate in MODS:
            if cells_want is not None and (arm, mod_name) not in cells_want:
                continue
            for seed in SEEDS:
                run_id += 1
                res = run_once(freqs, mod_name, mod_rate, seed)
                rt = time.time() - t0
                for ci, ch in enumerate(CH_NAMES):
                    cav = "A" if ci < 6 else "B"
                    kind = "out" if ch.endswith("out") else "in"
                    for wname in WINDOWS:
                        for pf, pa in window_peaks(res["rec"][ci], wname):
                            is_drv = bool(np.min(np.abs(freqs - pf)) <= PEAK_TOL)
                            peak_rows.append({
                                "run_id": run_id, "arm": arm,
                                "modulation": mod_name, "mod_rate": mod_rate,
                                "seed": seed, "cavity": cav, "pickup": ch,
                                "pickup_kind": kind, "window": wname,
                                "peak_freq": round(pf, 4),
                                "peak_amp": round(pa, 6), "is_drive": is_drv,
                                "audit_label": "drive" if is_drv else
                                audit_label(pf, freqs, mod_name, mod_rate),
                            })
                if mod_name == "NONE":
                    seen = {r["peak_freq"] for r in peak_rows
                            if r["run_id"] == run_id and r["cavity"] == "A"
                            and r["pickup_kind"] == "in"
                            and r["window"] == "DRIVE_END"}
                    for fd in freqs:
                        g1_tot += 1
                        if any(abs(s - fd) <= PEAK_TOL for s in seen):
                            g1_hits += 1
                ts, delays = transient_scalars(res)
                wA_e, wB_e = res["winE"].get("EARLY", (np.zeros(M), np.zeros(M)))
                wA_d, wB_d = res["winE"]["DRIVE_END"]
                wA_s, wB_s = res["winE"]["SETTLED"]
                for m_ in range(M):
                    tran_rows.append({
                        "run_id": run_id, "arm": arm, "modulation": mod_name,
                        "mod_rate": mod_rate, "seed": seed, "mode_idx": m_ + 1,
                        "f_drive": round(float(freqs[m_]), 4),
                        "f_eig_A": round(float(F_EIG_A[m_]), 4),
                        "f_eig_B": round(float(F_EIG_B[m_]), 4),
                        "E_A_early": round(float(wA_e[m_]), 3),
                        "E_B_early": round(float(wB_e[m_]), 3),
                        "E_A_driveend": round(float(wA_d[m_]), 3),
                        "E_B_driveend": round(float(wB_d[m_]), 3),
                        "E_A_settled": round(float(wA_s[m_]), 3),
                        "E_B_settled": round(float(wB_s[m_]), 3),
                        "transfer_ratio_driveend":
                            round(float(wB_d[m_] / wA_d[m_]), 5)
                            if wA_d[m_] > 0 else float("nan"),
                        "delay_s": round(delays[m_], 4)
                        if delays[m_] == delays[m_] else "",
                        "E_inj_mode": round(float(res["inj_m"][m_]), 3),
                        "E_rad_A_mode": round(float(res["radA_m"][m_]), 3),
                        "E_rad_B_mode": round(float(res["radB_m"][m_]), 3),
                        "esc_frac_A_mode":
                            round(float(res["radA_m"][m_] / res["inj_m"][m_]), 5)
                            if res["inj_m"][m_] > 0 else "",
                    })
                for _er in res["esamp"]:
                    es_rows.append({"run_id": run_id, "arm": arm,
                                    "modulation": mod_name, "seed": seed,
                                    "t": round(float(_er[0]), 3),
                                    "E_A": round(float(_er[1]), 3),
                                    "E_B": round(float(_er[2]), 3),
                                    "P_AB": round(float(_er[3]), 3),
                                    "P_BA": round(float(_er[4]), 3)})
                meta["runs"].append({"run_id": run_id, "arm": arm,
                                     "modulation": mod_name,
                                     "mod_rate": mod_rate, "seed": seed,
                                     **res["exch"], **ts})
                print(f"run {run_id}/45 {arm} {mod_name} seed{seed} "
                      f"closure={res['exch']['closure_frac']:.4f} "
                      f"elapsed={rt:.0f}s", flush=True)
    for name, rows in ((f"{outprefix}_pickups.csv", peak_rows),
                       (f"{outprefix}_transients.csv", tran_rows),
                       (f"{outprefix}_energy.csv", es_rows)):
        with open(f"{HD}/{name}", "w", newline="") as fh:
            wcsv = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wcsv.writeheader()
            wcsv.writerows(rows)
    meta["g1_drive_hits"] = g1_hits
    meta["g1_drive_total"] = g1_tot
    meta["FP"] = FP
    meta["F_EIG_A"] = [round(float(x), 4) for x in F_EIG_A]
    meta["F_EIG_B"] = [round(float(x), 4) for x in F_EIG_B]
    meta["B_MODES"] = [{"m": int(c[1]), "nr": int(c[2]), "p": int(c[3])}
                       for c in B_MODES]
    json.dump(meta, open(f"{HD}/{outprefix}_raw_meta.json", "w"), indent=1)
    print(f"DONE runs={run_id} peaks={len(peak_rows)} "
          f"G1={g1_hits}/{g1_tot} elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
