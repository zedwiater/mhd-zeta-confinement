"""
riemann_dual_cavity_sim.py — RIEMANN-LADDER DUAL-CAVITY SIM
(Zak's order 2026-10-22 13:21 EDT; follow-the-data amendment 13:24 EDT.
Pre-registration: ITERATION_LOG.md entry 2026-10-22 13:26 EDT — gates and
sweep grid stated there BEFORE this run.)

WHAT THIS IS: two coupled modal cavities (A driven, B coupled/undriven —
the transducer-series second-cavity picture), driven with standing waves
whose frequencies are the first-30 Riemann-zero ladder (dimensionless t
values from ~/workspace/riemann-tp/edge_frequencies.json). The sim is
SCALE-FREE: t values are used directly as dimensionless frequencies. No
Hz claim is made; any physical scaling is a free bench choice (see the
test plan), not a sim output.

PHYSICS REUSED (lab stock, not reinvented):
- audio_cavity_sim.py modal-oscillator conventions: per-mode damping
  gamma_m = f_m / Q, fixed per-mode detune (rng(7) uniform, as there),
  RK4 stepping.
- Cubic mixing as the intermodulation source (the lab's superlinear
  activation, here as a physical-space u^3 term projected back onto the
  modes — pseudospectral; this is what generates fi+-fj, 2fi+-fj, 3fi).
- transducer_sim.py dual-cavity bookkeeping picture: cavity B receives
  energy only through the coupling (residual second cavity).
Mode shapes phi_m(x) = sin(m*pi*x), m = 1..30 — a loaded cavity, whose
eigenfrequencies need not be harmonic; the ladder supplies the freqs.

SWEEP GRID (pre-stated, mirrors the log):
  arms       LADDER (t), DETUNE (1.01 t), SHIFT (t + mean_gap/2)
  modulation NONE | AM d=0.5 r in {mean_gap/2, mean_gap}
                   | FM dev=0.5 r in {mean_gap/2, mean_gap}
  seeds      5 per cell (drive phases + pickup noise)  -> 75 runs
  pickups    6 per cavity, x = .13 .29 .47 .61 .79 .93  (12 total)

READOUT: per-pickup Hann-windowed FFT of the steady window; peaks logged
to riemann_dual_cavity_pickups.csv (plot-ready vs the ladder). Candidate
= peak > 0.45 from every driven frequency. INTERMOD AUDIT (tol 0.45)
against fi+fj, |fi-fj|, 2fi, 3fi, fi/2, 2fi+fj, 2fi-fj, fi+-r, fi+-2r.
EDDY = candidate cluster in >=3/5 seeds at consistent freq, localized to
2..7 of 12 pickups. VERDICT gate G3: ladder unexplained+eddy rate must
beat BOTH controls by > 2 pooled SD, else NULL. Sim-internal claims only.
"""

import csv
import itertools
import json
import sys
import time

import numpy as np
from scipy.signal import find_peaks
from scipy.signal.windows import blackmanharris

HD = "/home/hatch/workspace/harness-dev"
LADDER_JSON = "/home/hatch/workspace/riemann-tp/edge_frequencies.json"

# ------------------------------------------------------------ parameters --
M = 30                       # modes per cavity (ladder length)
Q = 30.0                     # modal quality factor, gamma_m = f_m / Q
KAPPA_FRAC = 0.02            # inter-cavity stiffness coupling (fraction)
# Cubic strength: pilot 1 (BETA=0.30, A_TARGET=0.35) produced ZERO intermod
# peaks above threshold in all 12 pickups — NL generalized force ~1e-3 of
# the linear restoring force at those amplitudes. Intermod response scales
# as BETA * A^3, so the pilot-2 values below raise it ~100x into the
# 1-10%-of-drive band where products are visible but drives still dominate.
# Duffing self-shift at these values: |df| ~ 0.002 << PEAK_TOL — safe.
BETA = 10.0                  # cubic mixing coefficient (pilot-tuned; see log 13:31 EDT)
A_TARGET = 1.2               # per-mode resonant drive amplitude target
for _a in sys.argv:
    if _a.startswith("--beta="):
        BETA = float(_a.split("=")[1])
    if _a.startswith("--atarget="):
        A_TARGET = float(_a.split("=")[1])
DT = 0.0004                  # RK4 step (dimensionless time)
T_TOTAL = 6.0                # drive on throughout; record from T_REC
T_REC = 2.5                  # steady-state record window start
REC_EVERY = 2                # record every 2nd step -> fs = 1250
NOISE_FRAC = 0.002           # pickup measurement noise (frac of RMS)
N_GRID = 48                  # spatial grid for the u^3 projection
PICKUP_X = [0.13, 0.29, 0.47, 0.61, 0.79, 0.93]
PEAK_TOL = 0.45              # drive/audit match tolerance (freq units)
SEEDS = [0, 1, 2, 3, 4]

_lad = json.load(open(LADDER_JSON))
T_LADDER = np.array([z["t"] for z in _lad["zeros"]], dtype=float)
MEAN_GAP = float(_lad["mean_gap"])          # 3.0063

ARMS = {
    "LADDER": T_LADDER.copy(),
    "DETUNE": T_LADDER * 1.01,
    "SHIFT": T_LADDER + MEAN_GAP / 2.0,
}
MODS = [  # (name, rate) — rate 0.0 for NONE
    ("NONE", 0.0),
    ("AM", MEAN_GAP / 2.0),
    ("AM", MEAN_GAP),
    ("FM", MEAN_GAP / 2.0),
    ("FM", MEAN_GAP),
]
AM_DEPTH = 0.5
FM_DEV = 0.5

# mode shapes on the NL grid and at the pickups (loaded cavity)
_xg = (np.arange(N_GRID) + 0.5) / N_GRID
_midx = np.arange(1, M + 1)
PHI_G = np.sin(np.pi * np.outer(_xg, _midx))            # (N_GRID, M)
PHI_P = np.sin(np.pi * np.outer(PICKUP_X, _midx))       # (n_pick, M)
_GRID_W = 1.0 / N_GRID                                  # trapezoid-ish weights

# fixed per-mode detune, lab convention (audio_cavity_sim rng(7)); here a
# small static frequency jitter per mode, drawn once, same for all runs
DETUNE_FRAC = np.random.default_rng(7).uniform(-0.0015, 0.0015, M)


# ------------------------------------------------------------------ sim --
def run_once(freqs, mod_name, mod_rate, seed):
    """One dual-cavity run. Returns pickup time series (12, n_rec)."""
    rng = np.random.default_rng(10_000 + seed)
    phases = rng.uniform(0.0, 2.0 * np.pi, M)
    f = freqs * (1.0 + DETUNE_FRAC)                      # actual mode freqs
    om = 2.0 * np.pi * f
    gam = f / Q
    stiff = om ** 2
    kappa = KAPPA_FRAC * stiff
    # per-mode drive force: resonant steady amplitude ~= A_TARGET
    force = 2.0 * gam * om * A_TARGET

    n_steps = int(round(T_TOTAL / DT))
    rec_steps = int(round((T_TOTAL - T_REC) / (DT * REC_EVERY)))
    rec = np.empty((2 * len(PICKUP_X), rec_steps))
    ri = 0

    def drive(t):
        base = 2.0 * np.pi * freqs * t + phases
        if mod_name == "AM":
            return force * np.cos(base) * (1.0 + AM_DEPTH * np.cos(2.0 * np.pi * mod_rate * t))
        if mod_name == "FM":
            return force * np.cos(base + (FM_DEV / mod_rate) * np.sin(2.0 * np.pi * mod_rate * t))
        return force * np.cos(base)

    def rhs(y, t):
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:]
        ua = PHI_G @ qa
        ub = PHI_G @ qb
        nla = BETA * (PHI_G.T @ (ua ** 3)) * _GRID_W * 2.0
        nlb = BETA * (PHI_G.T @ (ub ** 3)) * _GRID_W * 2.0
        ramp = min(1.0, t / 0.5)
        aa = -2.0 * gam * va - stiff * qa + kappa * (qb - qa) + nla + ramp * drive(t)
        ab = -2.0 * gam * vb - stiff * qb + kappa * (qa - qb) + nlb
        return np.concatenate([va, aa, vb, ab])

    y = np.zeros(4 * M)
    for i in range(n_steps):
        t = i * DT
        if t >= T_REC and i % REC_EVERY == 0:
            rec[:len(PICKUP_X), ri] = PHI_P @ y[:M]
            rec[len(PICKUP_X):, ri] = PHI_P @ y[2 * M:3 * M]
            ri += 1
        k1 = rhs(y, t)
        k2 = rhs(y + 0.5 * DT * k1, t + 0.5 * DT)
        k3 = rhs(y + 0.5 * DT * k2, t + 0.5 * DT)
        k4 = rhs(y + DT * k3, t + DT)
        y = y + (DT / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    # seeded pickup measurement noise
    rms = np.sqrt(np.mean(rec ** 2, axis=1, keepdims=True)) + 1e-12
    rec = rec + rng.standard_normal(rec.shape) * NOISE_FRAC * rms
    return rec


# --------------------------------------------------------------- readout --
_FS = 1.0 / (DT * REC_EVERY)

def pickup_peaks(series):
    """Hann-windowed FFT peak pick for one pickup series -> (freq, amp)."""
    x = series - np.mean(series)
    w = blackmanharris(len(x))  # -92 dB sidelobes: kills the Hann leakage skirt
    spec = np.abs(np.fft.rfft(x * w)) * (2.0 / np.sum(w))
    freqs = np.fft.rfftfreq(len(x), d=1.0 / _FS)
    med = np.median(spec) + 1e-15
    floor = max(8.0 * med, 0.002 * np.max(spec))
    idx, props = find_peaks(spec, height=floor, prominence=4.0 * med, distance=3)
    out = []
    for i in idx:
        # parabolic sub-bin interpolation
        if 0 < i < len(spec) - 1:
            a, b, c = spec[i - 1], spec[i], spec[i + 1]
            d = 0.5 * (a - c) / (a - 2.0 * b + c + 1e-30)
        else:
            d = 0.0
        out.append((float(freqs[i] + d * (freqs[1] - freqs[0])), float(spec[i])))
    return out


def audit_label(fc, drive_freqs, mod_name, mod_rate):
    """Classify a candidate frequency against the intermod audit set.
    Returns the explaining rule, or 'UNEXPLAINED'."""
    F = drive_freqs
    def hit(target):
        return target > 0.5 and abs(fc - target) <= PEAK_TOL
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


# --------------------------------------------------------------- driver --
def main():
    pilot = "--pilot" in sys.argv
    csv_name = "riemann_dual_cavity_pickups.csv"
    for _a in sys.argv:  # pass-2 overrides (log 2026-10-22 13:47 EDT)
        if _a.startswith("--ttotal="):
            globals()["T_TOTAL"] = float(_a.split("=")[1])
        if _a.startswith("--trec="):
            globals()["T_REC"] = float(_a.split("=")[1])
        if _a.startswith("--csv="):
            csv_name = _a.split("=")[1]
    t0 = time.time()
    rows = []          # CSV rows (every detected peak)
    drive_hits = 0
    drive_total = 0
    run_id = 0
    cells = [(a, mn, mr) for a in ARMS for (mn, mr) in MODS]
    if "--cells=none" in sys.argv:
        cells = [(a, "NONE", 0.0) for a in ARMS]
    if pilot:
        cells = [("LADDER", "NONE", 0.0)]
    for arm, mod_name, mod_rate in cells:
        freqs = ARMS[arm]
        for seed in SEEDS:
            run_id += 1
            series = run_once(freqs, mod_name, mod_rate, seed)
            for cav_i, cav in enumerate(("A", "B")):
                for p_i, px in enumerate(PICKUP_X):
                    s = series[cav_i * len(PICKUP_X) + p_i]
                    for pf, pa in pickup_peaks(s):
                        is_drv = bool(np.min(np.abs(freqs - pf)) <= PEAK_TOL)
                        lab = "drive" if is_drv else audit_label(pf, freqs, mod_name, mod_rate)
                        rows.append({
                            "run_id": run_id, "arm": arm, "modulation": mod_name,
                            "mod_rate": mod_rate, "seed": seed, "cavity": cav,
                            "pickup_idx": p_i + 1, "pickup_x": px,
                            "peak_freq": round(pf, 4), "peak_amp": round(pa, 6),
                            "is_drive": is_drv, "audit_label": lab,
                        })
            if mod_name == "NONE":
                # G1 bookkeeping: driven freqs seen at >= 1 pickup this run
                seen = {r["peak_freq"] for r in rows if r["run_id"] == run_id}
                for fd in freqs:
                    drive_total += 1
                    if any(abs(s - fd) <= PEAK_TOL for s in seen):
                        drive_hits += 1
            if pilot:
                pk = [r for r in rows if r["run_id"] == run_id]
                drv = [r for r in pk if r["is_drive"]]
                cand = [r for r in pk if not r["is_drive"]]
                unx = [r for r in cand if r["audit_label"] == "UNEXPLAINED"]
                print(f"PILOT peaks={len(pk)} drive={len(drv)} cand={len(cand)} "
                      f"unexplained={len(unx)}")
                top = sorted(pk, key=lambda r: -r["peak_amp"])[:12]
                for r in top:
                    print(f"  cav{r['cavity']} p{r['pickup_idx']} f={r['peak_freq']:.2f} "
                          f"amp={r['peak_amp']:.4f} {r['audit_label']}")
                return
    with open(f"{HD}/{csv_name}", "w", newline="") as fh:
        wcsv = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(rows)
    print(f"runs={run_id} peaks_logged={len(rows)} "
          f"G1_drive_recovery={drive_hits}/{drive_total} "
          f"({drive_hits / max(drive_total,1):.3f}) elapsed={time.time()-t0:.1f}s")
    json.dump({"rows_file": csv_name,
               "g1_drive_hits": drive_hits, "g1_drive_total": drive_total},
              open(f"{HD}/riemann_dual_cavity_raw_meta.json", "w"), indent=1)


if __name__ == "__main__":
    main()
