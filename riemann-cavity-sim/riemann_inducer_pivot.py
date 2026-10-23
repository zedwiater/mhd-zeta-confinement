#!/usr/bin/env python3
"""
Iteration 4 PIVOT (pre-stated in ITERATION_LOG 2026-10-22 17:26 EDT):
re-measure gate GI3 with the instrument as pre-registered — coherent
demod amplitude at tone G at A_rin0 vs the MEDIAN OF THE WINDOW
SPECTRUM (2..120 Hz, Blackman, amplitude units) at the same pickup —
for SWEEP_G_d0.10 / LADDER, seeds 0-2. Deterministic rerun of the
grid configuration; window spectra are the only new output.
GATE (unchanged from pre-registration): ratio >= 5 AND demod phase
drift DRIVE_END->SETTLED < 30 deg in >= 2/3 seeds.
"""
import json
import numpy as np
import riemann_inducer_sim as I
import riemann_tp_cylinder_sim as S

I.AIND["A"] = 4500.0
I.AIND["B"] = 900.0
cell = I.CELLS["SWEEP_G_d0.10"]
tone = I.TONE_G
FS = I._FS
CH = {"A_rin0": 12, "A_out": 5}


def spec_median(x):
    x = x - x.mean()
    n = len(x)
    w = np.blackman(n)
    X = np.fft.rfft(x * w) * (2.0 / np.sum(w))
    fr = np.fft.rfftfreq(n, d=1.0 / FS)
    band = (fr >= 2.0) & (fr <= 120.0)
    return float(np.median(np.abs(X[band])))


out = []
for seed in (0, 1, 2):
    run = I.run_ind(S.ARMS["LADDER"], seed, cell)
    res = {"seed": seed}
    for win in ("DRIVE_END", "SETTLED"):
        w0, w1 = S.WINDOWS[win]
        i0, i1 = int(w0 * FS), int(w1 * FS)
        for pk, ci in CH.items():
            series = run["rec"][ci]
            med = spec_median(series[i0:i1])
            amp, ph = I.demod(series, tone, win)
            res[f"{pk}_{win}_spec_median"] = med
            res[f"{pk}_{win}_demod_amp"] = float(amp)
            res[f"{pk}_{win}_demod_phase"] = float(ph)
            if pk == "A_rin0":
                res[f"{win}_ratio"] = float(amp / med) if med > 0 else None
    dph = abs((res["A_rin0_SETTLED_demod_phase"]
               - res["A_rin0_DRIVE_END_demod_phase"] + 180.0) % 360.0 - 180.0)
    res["phase_drift_deg"] = float(dph)
    res["pass"] = bool(res["DRIVE_END_ratio"] >= 5.0 and dph < 30.0)
    out.append(res)
    print(f"seed {seed}: DE ratio={res['DRIVE_END_ratio']:.2f} "
          f"SET ratio={res['SETTLED_ratio']:.2f} drift={dph:.2f} "
          f"pass={res['pass']}")
npass = sum(r["pass"] for r in out)
print(f"GI3 (as registered): {'PASS' if npass >= 2 else 'FAIL'} ({npass}/3 seeds)")
json.dump({"cell": "SWEEP_G_d0.10", "arm": "LADDER", "tone": tone,
           "gi3_pass": bool(npass >= 2), "rows": out},
          open(f"{I.HD}/riemann_inducer_pivot.json", "w"), indent=1)
