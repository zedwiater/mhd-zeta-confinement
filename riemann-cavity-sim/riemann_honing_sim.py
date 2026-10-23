"""
riemann_honing_sim.py — RIEMANN ITERATION 7: "HONING ANY ANOMALY"
(Zak's order 2026-10-22 21:55 EDT). Pre-registration:
ITERATION_LOG.md 2026-10-22 21:5x EDT — round definition,
standing measurement set, anomaly criteria (a)-(d), confirmation
rule, 15-round hard cap — stated BEFORE round 1.

Engine: iteration 6 verbatim (riemann_capacitor_sim.Runner /
run_cap). This driver only adds: per-cell KAPPA multiplier
(capacitance as a tuning knob; f_bare = mult * G/2), a ring
test re-measuring f_ring whenever a new kappa multiplier
appears, fresh-seed support for confirmation rounds, and
append-only honing outputs. Iteration-6 files are NEVER written.

Usage: python3 riemann_honing_sim.py --spec=round_spec.json
Spec: {"round": N, "why": "...", "focus": "<cell name>",
 "cells": [{"name": ..., "role": focus|off|placebo|sweep,
   "ind": [x, theta_deg], "cap": [x, theta_deg] | null,
   "kappa_mult": 1.0, "arms": ["LADDER", ...],
   "seeds": [0, 1, 2]}]}
Outputs (append): riemann_honing_peaks_raw.csv,
riemann_honing_runs.csv, riemann_honing_rounds.json.
"""

import csv
import json
import os
import sys

import numpy as np

import riemann_capacitor_sim as C
import riemann_tp_cylinder_sim as S

HD = S.HD
KAPPA_BASE = float(C.KAPPA_C)          # 63.5035, f_bare = G/2
F_BARE_BASE = float(C.F_BARE)          # 10.10735
PEAKS_F = f"{HD}/riemann_honing_peaks_raw.csv"
RUNS_F = f"{HD}/riemann_honing_runs.csv"
ROUNDS_F = f"{HD}/riemann_honing_rounds.json"

PEAK_FIELDS = ["run_id", "round", "cell", "role", "arm", "seed",
               "kappa_mult", "cavity", "pickup", "pickup_kind",
               "window", "peak_freq", "peak_amp", "is_drive",
               "is_inducer", "audit_label"]


def ring_f_ring(kappa_mult):
    """Quiet ring test at the given kappa (iter-6 calib method)."""
    C.KAPPA_C = KAPPA_BASE * kappa_mult
    C.DRIVE_SCALE = 0.0
    cell = {"tones": [], "cav": "A", "x": 0.69, "theta": 0.0,
            "amp_scale": None, "cap": C._cap(0.69, 0), "q0": 300.0}
    res = C.run_cap(S.T_LADDER, 0, cell)
    C.DRIVE_SCALE = 1.0
    ser = res["series"]
    tt = np.array([r["t"] for r in ser])
    qq = np.array([r["q_cap"] for r in ser])
    m = (tt >= 1.0) & (tt <= 7.5)
    x = qq[m] - qq[m].mean()
    fs = 1.0 / (tt[1] - tt[0])
    spec = np.abs(np.fft.rfft(x * np.hanning(len(x))))
    fr = np.fft.rfftfreq(len(x), d=1.0 / fs)
    band = (fr >= 2.0) & (fr <= 80.0)
    return float(fr[band][np.argmax(spec[band])])


def append_rows(path, fields, rows):
    new_file = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new_file:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def existing_run_fields():
    if os.path.exists(RUNS_F):
        with open(RUNS_F) as fh:
            return next(csv.reader(fh))
    return None


def main():
    spec_path = None
    for a in sys.argv[1:]:
        if a.startswith("--spec="):
            spec_path = a.split("=", 1)[1]
    spec = json.load(open(spec_path))
    rnd = spec["round"]

    rounds_meta = (json.load(open(ROUNDS_F)) if os.path.exists(ROUNDS_F)
                   else {"rounds": [], "f_caps": {}})
    f_caps = rounds_meta.setdefault("f_caps", {})
    for cs in spec["cells"]:
        mult = float(cs.get("kappa_mult", 1.0))
        key = str(mult)
        if key not in f_caps:
            if mult == 1.0:
                f_caps[key] = [F_BARE_BASE, 33.6716]  # iter-6 calib
            else:
                fr = ring_f_ring(mult)
                f_caps[key] = [F_BARE_BASE * mult, fr]
                print(f"RING kappa x{mult}: f_bare={F_BARE_BASE * mult:.4f} "
                      f"f_ring={fr:.4f}", flush=True)

    offset = 0
    if os.path.exists(RUNS_F):
        import pandas as pd
        offset = int(pd.read_csv(RUNS_F, usecols=["run_id"])["run_id"].max())

    runner = C.Runner()
    runner.run_id = offset
    run_fields = existing_run_fields()
    all_run_rows = []
    for cs in spec["cells"]:
        mult = float(cs.get("kappa_mult", 1.0))
        C.KAPPA_C = KAPPA_BASE * mult
        C.SEEDS = list(cs.get("seeds", [0, 1, 2]))
        cap = C._cap(cs["cap"][0], cs["cap"][1]) if cs.get("cap") else None
        cell = C._cell(cs["ind"][0], np.radians(cs["ind"][1]), cap)
        name = cs["name"]
        n_pk, n_rn, n_ser = (len(runner.peak_rows), len(runner.run_rows),
                             len(runner.ser_rows))
        runner.run_set({name: cell}, {name: list(cs["arms"])})
        # tag new rows
        for r in runner.peak_rows[n_pk:]:
            r["round"] = rnd
            r["role"] = cs["role"]
            r["kappa_mult"] = mult
        # cap peak share from series (max E_cap / E_inj of the run)
        peak_share = {}
        for sr in runner.ser_rows[n_ser:]:
            rid = sr["run_id"]
            peak_share[rid] = max(peak_share.get(rid, 0.0), sr["E_cap"])
        for r in runner.run_rows[n_rn:]:
            r["round"] = rnd
            r["role"] = cs["role"]
            r["kappa_mult"] = mult
            r["f_bare"] = F_BARE_BASE * mult
            ps = peak_share.get(r["run_id"], 0.0)
            r["cap_peak_share"] = (ps / r["E_inj"]) if r["E_inj"] else 0.0
            all_run_rows.append(r)
        append_rows(PEAKS_F, PEAK_FIELDS, runner.peak_rows[n_pk:])
        runner.peak_rows = []
        runner.ser_rows = []
    if run_fields is None:
        run_fields = []
        for r in all_run_rows:
            for k in r.keys():
                if k not in run_fields:
                    run_fields.append(k)
    append_rows(RUNS_F, run_fields, all_run_rows)
    rounds_meta["rounds"].append(
        {"round": rnd, "why": spec.get("why", ""), "focus": spec["focus"],
         "cells": spec["cells"]})
    json.dump(rounds_meta, open(ROUNDS_F, "w"), indent=1)
    print(f"ROUND {rnd} DONE runs={len(all_run_rows)} "
          f"run_id<={runner.run_id}", flush=True)


if __name__ == "__main__":
    main()
