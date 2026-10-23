# RIEMANN_HONING — Iteration 7: "honing any anomaly"

Zak's order (2026-10-22 21:55 EDT): reiterate tuning of element
positioning in relation, record and track exchanges and
interactions, no more than 15 rounds, sweeping and adjusting as
needed; run against the Riemann line with controls maintained.
Pre-registration (anomaly criteria (a)–(d), round definition,
confirmation rule) is in ITERATION_LOG.md, stated before round 1.

Base: iteration-6 engine verbatim (TP cavity A + cylinder B,
transducer A79, inducer tone G, capacitor q-state element,
RK4, 3 seeds/round). Start point: the iteration-6 winner —
inducer (0.69, 0°), capacitor (0.79, 15°), κ×1 (f_bare = G/2).

**Verdict: NO ANOMALY in 15/15 rounds.** 231 runs, 398,482
peaks, 60,450 audited candidates. The one candidate (round 14)
failed its confirmation round and is characterized in full
below — its point estimate replicated; its anomaly status did
not.

## The 15 rounds

| # | focus (inducer / capacitor / κ) | knob moved | ind_frac | P | Δ_off vs spread | flags |
|---|---|---|---|---|---|---|
| 1 | (0.69,0°) / (0.79,15°) / ×1 | baseline = iter-6 winner | 0.00685 | 0.8536 | +0.0038 / 0.0449 | — |
| 2 | (0.69,0°) / (0.79,7.5°) / ×1 | cap θ −7.5° (+sweep 22.5°) | 0.00745 | 0.8512 | +0.0022 / 0.0276 | — |
| 3 | (0.69,0°) / (0.76,15°) / ×1 | cap x −0.03 (+sweep 0.82) | 0.00652 | 0.8523 | +0.0055 / 0.0496 | — |
| 4 | (0.66,0°) / (0.79,15°) / ×1 | inducer x −0.03 (+sweep 0.72) | 0.00498 | 0.8383 | +0.0029 / 0.0175 | — |
| 5 | (0.55,0°) / (0.79,15°) / ×1 | inducer onto trap d=0.24 | 0.00127 | 0.6836 | +0.0016 / 0.0152 | — |
| 6 | (0.69,45°) / (0.79,15°) / ×1 | inducer θ 45° (+sweep 90°) | 0.00031 | 0.8765* | +0.0020 / 0.0157 | — |
| 7 | (0.69,0°) / (0.79,15°) / ×2 | κ×2 (f_bare = G) | 0.00165 | 0.9782* | +0.0224 / 0.2468 | — |
| 8 | (0.69,0°) / (0.79,15°) / ×0.5 | κ×0.5 (f_bare = G/4) | 0.01281 | 0.7753 | −0.0013 / 0.0070 | — |
| 9 | (0.69,0°) / (0.55,0°) / ×1 | cap to trap radius, on-axis | 0.01678 | 0.7615 | +0.0083 / 0.0173 | — |
| 10 | (0.79,15°) / (0.79,15°) / ×1 | inducer+cap co-located at eddy | 0.00680 | 0.8600 | +0.0040 / 0.0453 | — |
| 11 | (0.69,0°) / (0.79,180°) / ×1 | cap anti-phase | 0.00658 | 0.8066 | +0.0033 / 0.0422 | — |
| 12 | (0.55,0°) / (0.79,15°) / ×2 | trap + κ×2 stacked | 0.00127 | 0.9371* | +0.0207 / 0.2048 | — |
| 13 | (0.69,0°) / (0.85,30°) / ×1 | cap wide toward aperture | 0.01050 | 0.7777 | +0.0224 / 0.1639 | — |
| 14 | (0.69,0°) / (0.79,15°) / ×4 | κ×4 (f_bare = 2G) | 0.00035 | 0.8455* | **+0.1744 / 0.1652** | **(b) candidate** |
| 15 | = R14 at fresh seeds 3/4/5 | confirmation | 0.00035 | 1.1488* | +0.1690 / 0.2191 | not confirmed |

\* P in a suppression/extinction regime (DRIVE_END demod below
the iter-6 guard floor 0.240): a ratio of residuals, not
standing. Honest standing points are the guard-passing ones —
R1 (0.8536 @ demod 0.258) and R10 (0.8600 @ demod 0.279).

Every round carried the full control set at the focus geometry
(OFF under all three arms, placebo capacitor, focus DETUNE +
SHIFT), the full audit6 interaction audit (with κ-matched
capacitor targets; f_ring re-measured by ring test at each new
κ: 33.36 Hz @×0.5, 33.67 @×1, 16.30 @×2, 17.22 @×4), and
per-run ledger/stability accounting. Coordinates, knob history
with one-line reasons, exchanges, and decisions per round:
`riemann_honing_log.csv`.

## The round-14 candidate, characterized

At κ×4 (f_bare = 40.43 Hz), the winner geometry's SETTLED
base-12 scatter sits at med d_off 0.6219 vs OFF's 0.4475:
Δ_off = +0.1744, exceeding that round's control spread 0.1652
(DETUNE +0.082, SHIFT +0.165, placebo +0.108); Δ_line +0.1540
as well. The confirmation round (fresh seeds 3/4/5, same
everything) reproduced the displacement almost exactly
(Δ_off +0.1690) — **but its own control spread grew to
0.2191**, so the pre-registered exceedance test fails. The
displacement is real and seed-robust; it is also a property
of the entire κ×4 regime: with the driven tone extinguished
(ind_frac 0.00035, DRIVE_END demod 0.0188), every scatter in
the regime — focus, arms, placebo — is displaced and wild,
and the focus is not special within it. Not an anomaly; the
control-spread criterion did the job it was registered to do.

## Bounds (the null, quantified)

- **(a) Spectral uniqueness:** zero unexplained peaks — not
  merely zero reproducible clusters — in 60,450 audited
  candidates across 15 rounds and 5 κ/geometry families.
  Sixth configuration family running; the cumulative
  uniqueness null strengthens.
- **(b) Riemann alignment:** at κ×1, focus |Δ_off| ≤ 0.0083
  in 8 of 10 rounds (max 0.0224 at R13, 7× inside its
  spread); Δ_line tracks Δ_off and never leads it. No
  geometry tuned in this search aligns scatter with the
  ladder or the eigen-offsets beyond control noise.
- **(c) Exchanges:** at the two iter-6-matched geometries
  (winner, trap) honing reproduced baselines digit-exactly
  (ind_frac 0.00685, 0.00127). No metric exceeded 3 pooled SD
  anywhere it could be tested. **B-share of output is
  geometry-blind: 0.0857–0.0881 in all 15 rounds** — relative
  positioning of inducer/capacitor does not move cross-cavity
  exchange at all in this engine.
- **(d) Stability/ledger:** max closure 0.00056 over 231
  runs; every settled slope negative (max −0.0518); zero
  pool clamps, including the κ×4 extinction regime.

## What the honing mapped (engine-expected, none anomalous)

- κ is a pure loading dial at the winner geometry: ind_frac
  0.0128 / 0.0069 / 0.0017 / 0.00035 at ×0.5 / ×1 / ×2 / ×4,
  monotone, with storage share rising to 9.5e-04 at ×4.
- The trap is a hard floor: κ×2 stacked on it (R12) changes
  ind_frac not at all (0.00127 = trap-only R5).
- Inducer angle is the strongest positional knob: θ=45°
  nearly decouples the profile (OFF acceptance 3,661 vs
  143,763 at θ=0 — a cos(nθ) modal-overlap zero, present
  cap-free), while θ=90° couples at ~half the θ=0 strength.
- The capacitor acts only when proximal (R9: at the trap
  radius on-axis it is inert, ind_frac 0.01678 ≈ OFFCAP).
- A DETUNE-arm interaction at the wide cap geometry (R13):
  the (0.85, 30°) capacitor pulls detuned-drive scatter
  toward the eigen-offsets (control Δ −0.1639) while leaving
  the LADDER condition unaligned — knob interaction in a
  control arm, recorded, not claimed.

## Limitations

- Modal toy, sim-internal claims only; no turbulence, no
  bench implication. "Anomaly" means anomalous *within this
  engine's audit/alignment/control construction*.
- Criterion (c) could only be tested at geometries with an
  iteration-5/6 baseline (winner, trap, OFF, placebo); κ≠1
  and novel-angle geometries are unmatched by construction
  and were covered by criteria (a), (b), (d) instead.
- 3 seeds per cell: pooled-SD and spread estimates are
  small-sample; the R14/R15 pair shows the spread itself
  carries ±30% seed noise at extreme κ.
- Alignment uses SETTLED base-12 medians (the standing
  GI4/GC6 construction); per-seed alignment distributions
  are in `riemann_honing_scatter.csv` for re-analysis.

## Files

- `riemann_honing_sim.py` — driver (iter-6 engine + κ knob,
  ring re-measurement, fresh-seed support, append-only I/O)
- `riemann_honing_analyze.py` — per-round audit6 + criteria
- `riemann_honing_log.csv` — the tracking log: one row per
  round (coordinates, knobs + why, exchanges, audit,
  alignment vs controls, decision)
- `riemann_honing_results.json` — all verdict numbers
- `riemann_honing_round_summary.json` — per-round detail
- `riemann_honing_scatter.csv` — every peak, audit6-labelled,
  with d_line / d_off / d_offB
- `riemann_honing_peaks_raw.csv`, `riemann_honing_runs.csv`,
  `riemann_honing_rounds.json`, `honing_spec_r*.json` — raw
  data + the round specs as run
