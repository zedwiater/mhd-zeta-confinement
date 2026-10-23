"""
riemann_tp_cylinder_phase.py — leakage PHASE matrix, focused pass
(logged 2026-10-22 14:56 EDT; second delta item 1). NONE cells only,
3 arms x 5 seeds, FP = 0.025 (main-grid value). For each drive component
and cavity: complex DFT (Blackman-Harris window, DRIVE_END window) at the
exact drive frequency, interior channels vs exterior (port) channel:
  H = X_out / complex_mean(X_in)   -> |H|, arg(H)
  amp_ratio = |X_out| / mean_in |X_in|   (matches the amplitude leakage)
  interior phase spread = circular std of arg(X_in,k) across pickups.
Writes riemann_tp_cylinder_leakage_phase.csv and merges a summary into
riemann_tp_cylinder_results.json under "leakage_phase".
"""

import csv
import json

import numpy as np
from scipy.signal.windows import blackmanharris

import riemann_tp_cylinder_sim as S

HD = S.HD
FS = 1250.0
W0, W1 = S.WINDOWS["DRIVE_END"]

rows = []
for arm, freqs in S.ARMS.items():
    for seed in S.SEEDS:
        res = S.run_once(freqs, "NONE", 0.0, seed)
        rec = res["rec"]
        i0, i1 = int(W0 * FS), int(W1 * FS)
        seg = rec[:, i0:i1]
        t = np.arange(seg.shape[1]) / FS
        w = blackmanharris(seg.shape[1])
        for cav, chs_in, ch_out in (("A", range(0, 5), 5),
                                    ("B", range(6, 11), 11)):
            Xin = []
            for ch in chs_in:
                x = (seg[ch] - seg[ch].mean()) * w
                Xin.append(np.array([np.sum(x * np.exp(-2j * np.pi * f * t))
                                     for f in freqs]))
            Xin = np.array(Xin)                       # (5, M) complex
            xo = (seg[ch_out] - seg[ch_out].mean()) * w
            Xout = np.array([np.sum(xo * np.exp(-2j * np.pi * f * t))
                             for f in freqs])
            Xin_mean = Xin.mean(axis=0)
            # units fix (same as amplitude leakage): port series is
            # velocity-weighted; divide by 2*pi*f for displacement
            # equivalent. Positive real factor — arg(H) unaffected.
            H = (Xout / (2.0 * np.pi * freqs)) / Xin_mean
            ph_in = np.angle(Xin)
            spread = np.sqrt(-2.0 * np.log(np.abs(np.exp(1j * ph_in).mean(axis=0))))
            for m_ in range(S.M):
                rows.append({
                    "arm": arm, "seed": seed, "cavity": cav,
                    "mode_idx": m_ + 1, "f_drive": round(float(freqs[m_]), 4),
                    "amp_ratio": round(float(np.abs(Xout[m_])
                                             / np.mean(np.abs(Xin[:, m_]))), 4),
                    "H_amp": round(float(np.abs(H[m_])), 4),
                    "H_phase_deg": round(float(np.degrees(np.angle(H[m_]))), 2),
                    "phase_spread_in_deg": round(float(np.degrees(spread[m_])), 2),
                })
        print(f"{arm} seed{seed} done", flush=True)

with open(f"{HD}/riemann_tp_cylinder_leakage_phase.csv", "w", newline="") as fh:
    wcsv = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    wcsv.writeheader()
    wcsv.writerows(rows)

import pandas as pd
df = pd.DataFrame(rows)
summ = {}
for arm in ("LADDER", "DETUNE", "SHIFT"):
    summ[arm] = {}
    for cav in ("A", "B"):
        s = df[(df["arm"] == arm) & (df["cavity"] == cav)]
        med = s.groupby("mode_idx").median(numeric_only=True)
        summ[arm][cav] = {
            "H_amp_by_mode": [round(float(x), 3) for x in med["H_amp"]],
            "H_phase_deg_by_mode": [round(float(x), 1) for x in med["H_phase_deg"]],
            "phase_spread_deg_by_mode": [round(float(x), 1)
                                         for x in med["phase_spread_in_deg"]],
            "H_phase_band_median_deg": {
                "n1_10": float(np.median(med["H_phase_deg"].iloc[:10])),
                "n11_20": float(np.median(med["H_phase_deg"].iloc[10:20])),
                "n21_30": float(np.median(med["H_phase_deg"].iloc[20:]))},
            "H_amp_median": float(med["H_amp"].median()),
        }
try:
    R = json.load(open(f"{HD}/riemann_tp_cylinder_results.json"))
except FileNotFoundError:
    R = {}
R["leakage_phase"] = summ
json.dump(R, open(f"{HD}/riemann_tp_cylinder_results.json", "w"), indent=1)
print(json.dumps({a: {c: {"H_amp_median": summ[a][c]["H_amp_median"],
                          "phase_bands": summ[a][c]["H_phase_band_median_deg"]}
                      for c in ("A", "B")} for a in summ}, indent=1))
