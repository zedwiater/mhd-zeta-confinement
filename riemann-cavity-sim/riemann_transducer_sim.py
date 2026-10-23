"""
riemann_transducer_sim.py — RIEMANN ITERATION 3: LOCAL TRANSDUCERS AT
EDDY SITES (Zak's order 2026-10-22 16:11 EDT, correcting a Gemini
proposal). Pre-registration: ITERATION_LOG.md 2026-10-22 16:17 EDT —
model, sites, grid and gates stated BEFORE any run.

NOT ORDERED (rejected): no drive masking to eddy sites, no eddy
injection, no removal of TP flow. The iteration-2 system is untouched:
TP cavity A + cylinder cavity B, FP = 0.025, GLOBAL Riemann-ladder
antenna drive, LADDER/DETUNE/SHIFT arms. Geometry, constants, windows,
readout and audit are IMPORTED from riemann_tp_cylinder_sim so the
cavities are bit-identical to iteration 2.

TRANSDUCER (one per run, at one site; full statement in the log):
- Profile: bump b on the cavity sample grid (A: 48-pt x-grid, Gaussian
  sigma 0.06; B: 64 volume samples, Gaussian in Cartesian distance,
  sigma 0.35). Modal profile W = (2/N) PHI^T b (engine convention).
  Sense s(t) = W . v_cav.
- Extraction: F1 = -G * V_SAT * tanh(s/V_SAT) * W (saturating damper);
  p_ext = G*V_SAT*tanh(s/V_SAT)*s >= 0 by construction. Force cap
  G*V_SAT*|W|. eta = 0.5 of p_ext feeds the pool; (1-eta) is transducer
  dissipation (ledger line).
- Re-injection: envelope state env (tau 0.25) of |s|; F2 = +G * sigma(e)
  * env * W with pool state e, de/dt = eta*p_ext - p_reinj,
  sigma(e) = e/(e+E_GATE), p_reinj = F2 . v exact signed. Pool clamped
  >= 0 at step boundaries (clamps counted).
- Ledger: E_inj = dStored + E_diss + E_rad + E_tr_diss + dPool + resid.

SITES: OFF (none) / A79 (A x=0.79) / A61 (A x=0.61) / APORT (A x=0.97) /
BCROSS (B at B_in4 point) / PLACEBO (B at B_in1 point). Provenance in
the pre-registration. Sensing is passive at ALL sites in every run.

Modes: --calib0 (G=0 sense runs, LADDER NONE seeds 0-2, series only),
--calib1 --config=X (one actuated run), --grid (full pre-registered
grid). Params via --g= --vsat= --egate= (G2 = G, ETA = 0.5 fixed).
"""

import csv
import json
import sys
import time

import numpy as np

sys.argv = [a for a in sys.argv
            if not a.startswith(("--fp=", "--windows="))]
import riemann_tp_cylinder_sim as S

HD = S.HD
M = S.M
DT = S.DT
T_OFF = S.T_OFF
T_TOTAL = S.T_TOTAL
WINDOWS = S.WINDOWS
ARMS = S.ARMS
SEEDS = [0, 1, 2]
MODS3 = [("NONE", 0.0), ("AM", S.MEAN_GAP)]
SER_EVERY = 20                     # series point every 20 steps (8 ms)

ETA = 0.5
TAU_ENV = 0.25
TR = {"G": 0.0, "V_SAT": 1.0, "EGATE": 1.0}

# ---------------------------------------------------------------- sites --
_xgA = (np.arange(48) + 0.5) / 48.0

def a_profile(x_s, sigma=0.06):
    b = np.exp(-((_xgA - x_s) ** 2) / (2.0 * sigma ** 2))
    return (2.0 / 48.0) * (S.PHI_G_A.T @ b)

_sx = S._s_r * np.cos(S._s_th)
_sy = S._s_r * np.sin(S._s_th)
_sz = S._s_z

def b_profile(pt, sigma=0.35):
    r, th, z = pt
    dx = _sx - r * np.cos(th)
    dy = _sy - r * np.sin(th)
    dz = _sz - z
    b = np.exp(-((dx ** 2 + dy ** 2 + dz ** 2) / (2.0 * sigma ** 2)))
    return (2.0 / 64.0) * (S.S_B.T @ b)

# site -> (cavity, W, nearest pickup, neighbour pickups)
SITES = {
    "A79":     ("A", a_profile(0.79), "A_in5",
                ["A_in1", "A_in2", "A_in3", "A_in4", "A_out"]),
    "A61":     ("A", a_profile(0.61), "A_in4",
                ["A_in1", "A_in2", "A_in3", "A_in5", "A_out"]),
    "APORT":   ("A", a_profile(0.97), "A_out",
                ["A_in1", "A_in2", "A_in3", "A_in4", "A_in5"]),
    "BCROSS":  ("B", b_profile(S.B_PICK[3]), "B_in4",
                ["B_in1", "B_in2", "B_in3", "B_in5", "B_out"]),
    "PLACEBO": ("B", b_profile(S.B_PICK[0]), "B_in1",
                ["B_in2", "B_in3", "B_in4", "B_in5", "B_out"]),
}
SITE_NAMES = list(SITES.keys())
for _nm, (_cav, _W, _pk, _nb) in SITES.items():
    _c = _W ** 2
    SITES[_nm] = (_cav, _W, _pk, _nb, _c / _c.sum())   # + eloc weights

# ------------------------------------------------------------------ sim --
def run_tr(drive_freqs, mod_name, mod_rate, seed, config):
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
    G, VS, EG = TR["G"], TR["V_SAT"], TR["EGATE"]
    act = None
    if config != "OFF":
        act = SITES[config]              # (cav, W, pk, nb, elocw)

    n_steps = int(round(T_TOTAL / DT))
    n_rec = int(round(T_TOTAL / (DT * S.REC_EVERY)))
    n_es = int(round(T_TOTAL / (DT * S.ESAMP_EVERY)))
    rec = np.empty((12, n_rec))
    esamp = np.empty((n_es, 5))
    emode = np.empty((n_es, 2 * M))
    series = []
    inj_m = np.zeros(M); dissA = np.zeros(M); dissB = np.zeros(M)
    radA_m = np.zeros(M); radB_m = np.zeros(M)
    winE = {w: (np.zeros(M), np.zeros(M)) for w in WINDOWS}
    ri = ei = 0
    E0 = None
    # transducer bookkeeping
    E_ext = E_reinj = E_tr_diss = 0.0
    sat_time = gate_time = 0.0
    clamp_count = 0
    clamp_E = 0.0
    pool_max = 0.0

    def drive(t):
        if t >= T_OFF:
            return np.zeros(M)
        base = omd * t + phases
        ramp = min(1.0, t / 0.5)
        if mod_name == "AM":
            return ramp * force * np.cos(base) * (
                1.0 + S.AM_DEPTH * np.cos(2.0 * np.pi * mod_rate * t))
        if mod_name == "FM":
            return ramp * force * np.cos(
                base + (S.FM_DEV / mod_rate) * np.sin(2.0 * np.pi * mod_rate * t))
        return ramp * force * np.cos(base)

    def tr_forces(va, vb, env, pool):
        """Returns (F1+F2 modal force on actuated cavity, cav, p_ext, p_reinj, s)."""
        if act is None or G == 0.0:
            return None, None, 0.0, 0.0, 0.0
        cav, W = act[0], act[1]
        v = va if cav == "A" else vb
        s = float(W @ v)
        f1 = -G * VS * np.tanh(s / VS) * W
        p_ext = -float(f1 @ v)
        sig = pool / (pool + EG) if pool > 0.0 else 0.0
        f2 = G * sig * env * W
        p_reinj = float(f2 @ v)
        return f1 + f2, cav, p_ext, p_reinj, s

    def rhs(y, t):
        qa, va, qb, vb = y[:M], y[M:2 * M], y[2 * M:3 * M], y[3 * M:4 * M]
        env, pool = y[4 * M], y[4 * M + 1]
        ua = S.PHI_G_A @ qa
        ub = S.S_B @ qb
        nla = S.BETA * (S.PHI_G_A.T @ (ua ** 3)) * (2.0 / 48.0)
        nlb = S.BETA * (S.S_B.T @ (ub ** 3)) * (2.0 / 64.0)
        aa = -2.0 * gA * va - stiffA * qa + kappa * (qb - qa) + nla + drive(t)
        ab = -2.0 * gB * vb - stiffB * qb + kappa * (qa - qb) + nlb
        denv = dpool = 0.0
        F, cav, p_ext, p_reinj, s = tr_forces(va, vb, env, pool)
        if F is not None:
            if cav == "A":
                aa = aa + F
            else:
                ab = ab + F
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
        F, cav, p_ext, p_reinj, s_act = tr_forces(va, vb, env, pool)
        if F is not None:
            E_ext += p_ext * DT
            E_reinj += p_reinj * DT
            E_tr_diss += (1.0 - ETA) * p_ext * DT
            if abs(s_act) >= 2.0 * VS:
                sat_time += DT
            sig_now = pool / (pool + EG) if pool > 0.0 else 0.0
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
                   "env": float(env), "s_act": float(s_act)}
            for nm, (cv, W, _pk, _nb, cw) in SITES.items():
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

    exch = {
        "E_inj": float(inj_m.sum()), "E_diss_A": float(dissA.sum()),
        "E_diss_B": float(dissB.sum()), "E_rad_A": float(radA_m.sum()),
        "E_rad_B": float(radB_m.sum()),
        "E_stored0_A": float(E0[0]), "E_stored0_B": float(E0[1]),
        "E_stored_end_A": float(Eend[0]), "E_stored_end_B": float(Eend[1]),
        "E_tr_ext": float(E_ext), "E_tr_reinj_net": float(E_reinj),
        "E_tr_diss": float(E_tr_diss), "pool_end": pool_end,
        "pool_max": float(pool_max), "clamp_count": int(clamp_count),
        "clamp_E": float(clamp_E),
        "sat_frac": float(sat_time / T_TOTAL), "gate_frac": float(gate_time / T_TOTAL),
    }
    exch["closure_resid"] = exch["E_inj"] - (
        (exch["E_stored_end_A"] - exch["E_stored0_A"])
        + (exch["E_stored_end_B"] - exch["E_stored0_B"])
        + exch["E_diss_A"] + exch["E_diss_B"]
        + exch["E_rad_A"] + exch["E_rad_B"]
        + exch["E_tr_diss"] + pool_end)
    exch["closure_frac"] = (abs(exch["closure_resid"]) / exch["E_inj"]
                            if exch["E_inj"] > 0 else float("nan"))
    exch["pool_balance_resid"] = pool_end - (ETA * E_ext - E_reinj)
    return {"rec": rec, "esamp": esamp[:ei], "emode": emode[:ei],
            "inj_m": inj_m, "radA_m": radA_m, "radB_m": radB_m,
            "winE": winE, "exch": exch, "series": series}


# --------------------------------------------------------------- driver --
def _flags():
    fl = {}
    for a in sys.argv[1:]:
        if a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            fl[k] = v
        elif a.startswith("--"):
            fl[a[2:]] = True
    return fl


def main():
    fl = _flags()
    TR["G"] = float(fl.get("g", 0.0))
    TR["V_SAT"] = float(fl.get("vsat", 1.0))
    TR["EGATE"] = float(fl.get("egate", 1.0))
    t0 = time.time()

    if "calib0" in fl:
        rows = []
        for seed in SEEDS:
            res = run_tr(ARMS["LADDER"], "NONE", 0.0, seed, "A79")
            for r in res["series"]:
                rows.append({"seed": seed, **r})
            print(f"calib0 seed{seed} closure={res['exch']['closure_frac']:.5f}",
                  flush=True)
        with open(f"{HD}/riemann_transducer_calib0_series.csv", "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"calib0 done {time.time()-t0:.0f}s")
        return

    if "calib1" in fl:
        cfg = fl.get("config", "A79")
        res = run_tr(ARMS["LADDER"], "NONE", 0.0, 0, cfg)
        print(f"calib1 config={cfg} G={TR['G']} V_SAT={TR['V_SAT']} "
              f"EGATE={TR['EGATE']}")
        print(json.dumps(res["exch"], indent=1))
        print(f"E_ext/E_inj = {res['exch']['E_tr_ext']/res['exch']['E_inj']:.4f}")
        print(f"elapsed={time.time()-t0:.0f}s")
        return

    # full grid
    configs = ["OFF"] + SITE_NAMES
    peak_rows, tran_rows, ser_rows, run_rows = [], [], [], []
    meta = {"runs": [], "TR": dict(TR), "ETA": ETA, "TAU_ENV": TAU_ENV,
            "sites": {nm: {"cavity": SITES[nm][0],
                           "nearest_pickup": SITES[nm][2],
                           "neighbours": SITES[nm][3]}
                      for nm in SITE_NAMES}}
    run_id = 0
    g1 = {}                                # config -> [hits, tot]
    for config in configs:
        g1[config] = [0, 0]
        for arm, freqs in ARMS.items():
            for mod_name, mod_rate in MODS3:
                for seed in SEEDS:
                    run_id += 1
                    res = run_tr(freqs, mod_name, mod_rate, seed, config)
                    for ci, ch in enumerate(S.CH_NAMES):
                        cav = "A" if ci < 6 else "B"
                        kind = "out" if ch.endswith("out") else "in"
                        for wname in WINDOWS:
                            for pf, pa in S.window_peaks(res["rec"][ci], wname):
                                is_drv = bool(np.min(np.abs(freqs - pf)) <= S.PEAK_TOL)
                                peak_rows.append({
                                    "run_id": run_id, "config": config,
                                    "arm": arm, "modulation": mod_name,
                                    "mod_rate": mod_rate, "seed": seed,
                                    "cavity": cav, "pickup": ch,
                                    "pickup_kind": kind, "window": wname,
                                    "peak_freq": round(pf, 4),
                                    "peak_amp": round(pa, 6), "is_drive": is_drv,
                                    "audit_label": "drive" if is_drv else
                                    S.audit_label(pf, freqs, mod_name, mod_rate),
                                })
                    if mod_name == "NONE":
                        seen = {r["peak_freq"] for r in peak_rows
                                if r["run_id"] == run_id and r["cavity"] == "A"
                                and r["pickup_kind"] == "in"
                                and r["window"] == "DRIVE_END"}
                        for fd in freqs:
                            g1[config][1] += 1
                            if any(abs(sv - fd) <= S.PEAK_TOL for sv in seen):
                                g1[config][0] += 1
                    ts, delays = S.transient_scalars(res)
                    wA_e, wB_e = res["winE"].get("EARLY", (np.zeros(M), np.zeros(M)))
                    wA_d, wB_d = res["winE"]["DRIVE_END"]
                    wA_s, wB_s = res["winE"]["SETTLED"]
                    for m_ in range(M):
                        tran_rows.append({
                            "run_id": run_id, "config": config, "arm": arm,
                            "modulation": mod_name, "mod_rate": mod_rate,
                            "seed": seed, "mode_idx": m_ + 1,
                            "f_drive": round(float(freqs[m_]), 4),
                            "f_eig_A": round(float(S.F_EIG_A[m_]), 4),
                            "f_eig_B": round(float(S.F_EIG_B[m_]), 4),
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
                        })
                    for r in res["series"]:
                        ser_rows.append({"run_id": run_id, "config": config,
                                         "arm": arm, "modulation": mod_name,
                                         "seed": seed, **r})
                    run_rows.append({"run_id": run_id, "config": config,
                                     "arm": arm, "modulation": mod_name,
                                     "mod_rate": mod_rate, "seed": seed,
                                     **res["exch"], **ts})
                    meta["runs"].append({"run_id": run_id, "config": config,
                                         "arm": arm, "modulation": mod_name,
                                         "seed": seed})
                    print(f"run {run_id}/108 {config} {arm} {mod_name} "
                          f"seed{seed} closure={res['exch']['closure_frac']:.5f} "
                          f"Eext={res['exch']['E_tr_ext']:.3f} "
                          f"elapsed={time.time()-t0:.0f}s", flush=True)
    for name, rows in (("riemann_transducer_pickups.csv", peak_rows),
                       ("riemann_transducer_transients.csv", tran_rows),
                       ("riemann_transducer_series.csv", ser_rows),
                       ("riemann_transducer_runs.csv", run_rows)):
        with open(f"{HD}/{name}", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    meta["g1_by_config"] = {k: {"hits": v[0], "total": v[1]}
                            for k, v in g1.items()}
    json.dump(meta, open(f"{HD}/riemann_transducer_raw_meta.json", "w"), indent=1)
    print(f"DONE runs={run_id} peaks={len(peak_rows)} "
          f"G1={meta['g1_by_config']} elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
