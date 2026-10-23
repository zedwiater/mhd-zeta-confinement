# RIEMANN_CAPACITOR — Iteration 6: proximal capacitors in the eddy field

Sim-internal claims only; modal toy (M=30 modes/cavity, RK4, DT=0.0004,
T=8 s, drive off at 4 s). Base system is iteration 5 unchanged (TP flow,
iter-3 transducers, inducer tone G at the (0.69, 0°) reference).
Pre-registration, calibration, and the one logged pivot are in
ITERATION_LOG.md (2026-10-22). 174 runs, 3 seeds, 302,599 peaks
(refinement skipped 4 cells already covered by the coarse grid).

## The capacitor element

A lossless reactive storage element (cavity A only), placed at (x, θ):

- charge state q; **dq/dt = KAPPA · (W_C · v_A)** — charges from the local
  modal field through profile W_C;
- back-force on the modes **−KAPPA · q · W_C** — spring-like, lossless;
- W_C = the iter-5 bump profile family (σ = 0.06, cos(nθ)) sampled on the
  28-channel ring geometry, **normalized to unit modal L2 norm**, so the
  bare resonance f_bare = KAPPA/(2π) is placement-independent;
- **E_cap = q²/2 is a storage line in the ledger**: field power
  −KAPPA·q·(W_C·v) ≡ dE_cap/dt identically, so closure carries it as a
  stored-energy term, not a source or sink.

One capacitance, per Zak's "inducer is subharmonic of the transducers'
mid-range released frequencies" framing: **f_bare = G/2 = 10.10735 Hz,
KAPPA = 2π·f_bare = 63.5064**.

## Calibration (in the open, before the grid)

- **Ring test** (no drive, no inducer, q(0) = 300, cap at (0.69, 0°)):
  q_max = 300.0 — the element never exceeds its initial charge; passive
  in situ. In-situ ring frequency of q(t): **f_ring = 33.6716 Hz** — set
  by the W-weighted modal content, not by f_bare (KAPPA² ≪ modal
  stiffness). Audit6 therefore tests both f_bare and f_ring.
- **Sensing run** (reference geometry, cap co-located with the inducer):
  closure 0.000445; the capacitor **halves accepted inducer energy**
  (E_ind 143,763 → 74,072, −48.5%) and cuts sat_frac 0.146 → 0.091.
  Its stored energy is small and brief (peak share of E_inj 9.3e-5 at
  t = 1.44 s; end share 1.5e-6).
- **Contingency disposition (the one pivot, logged):** the storage share
  fell below the pre-registered 0.01% floor, formally triggering the
  split-gain remedy — but that remedy is a change of variables (only
  a·b enters the closed-loop dynamics) and would change nothing.
  KAPPA stayed as designed; the element engages through the **loading
  channel** (impedance at the injection point), which is the physical
  action under test. Gates and criterion unchanged.
- First calibration attempt produced non-conservative artifacts with a
  one-sided-coupling signature and was discarded; two independent
  re-runs reproduce the clean trajectory digit-for-digit (see log).

## Placement search — best standing position

Criterion (pre-registered): maximize median-over-seeds
P = demod_amp_G(SETTLED, A_in5) / demod_amp_G(DRIVE_END, A_in5),
subject to per-seed guards (DRIVE_END demod ≥ 0.5× OFFCAP median;
sat ≤ 0.30; closure ≤ 0.001; settled slope ≤ +0.005/s).
Reference points: **OFFCAP P = 0.7448; PLACEBO (0.13, 0°) P = 0.7804.**

Coarse grid: radii {0.59, 0.65, 0.71, 0.79, 0.85} × θ {0…180° step 30°},
35 cells, then a 3×3 refinement around the stage-1 winner.
Full ranking with coordinates: `riemann_capacitor_grid.csv`.

Top of the criterion pool:

| cell | cap (x, θ°) | P median | guards |
|---|---|---|---|
| CAP_x0.71_t090 | (0.71, 90°) | 0.8732 | **FAIL — excluded** |
| CAP_x0.79_t090 | (0.79, 90°) | 0.8540 | **FAIL — excluded** |
| **RFC_x0.79_t015** | **(0.79, 15°)** | **0.8536** | **PASS — winner** |
| RFC_x0.82_t015 | (0.82, 15°) | 0.8492 | PASS |
| CAP_x0.79_t000 | (0.79, 0°) | 0.8434 | PASS |
| RFC_x0.76_t000 | (0.76, 0°) | 0.8399 | PASS |
| CAP_x0.85_t000 | (0.85, 0°) | 0.8267 | PASS |

**The guards earned their keep.** The θ=90° ridge posts the highest raw
ratios by *suppressing* the driven tone (DRIVE_END demod 0.164 vs the
0.240 floor) — persistence-by-suppression, a ratio artifact. Excluded
by the pre-registered floor, exactly as designed.

**Best standing position: the eddy site itself — x = 0.79 (the A79
transducer/eddy center), θ = 15°.** P per-seed 0.854 / 0.843 / 0.916;
margin over OFFCAP (+0.109) exceeds the pooled SD (GC4 PASS). The
winner trades acceptance for persistence (ind_frac 0.00685 vs OFFCAP
0.01657): a capacitor parked on the eddy loads the injection, and what
does get in, stands longer.

## Subspectra — full audit (audit6)

Audit6 = base intermod audit (0.45 / 0.30 / 0.15) + flow-shift targets
(0.15) + inducer targets (0.15) + **capacitor targets at 0.15 built from
both f_bare and f_ring**: f_cap, 2f_cap, f_cap/2, tone ± k·f_cap,
f_i ± k·f_cap (k = 1, 2) → labels cap_res / cap_harm / cap_mix / cap_sb.

- 56,238 candidate peaks audited; label counts are dominated by base
  terms (2f_i−f_j 32,258; f_i−f_j 15,711; f_i/2 4,378), flow_adv 1,558.
- **Capacitor spectral footprint: NONE.** Zero candidates needed a cap_*
  label — no line at the capacitor's resonance, harmonics, or sidebands
  that base/flow/inducer terms didn't already explain. A lossless
  linear element creates no frequencies; it loads, detunes, stores
  briefly. Its action is impedance, not emission.
- **Uniqueness NULL, 5th configuration (GC5 PASS):** 2 single-seed
  UNEXPLAINED peaks total; **zero reproducible unexplained clusters**
  in any cell, including OFFCAP, PLACEBO, and the winner.

## Ledger & stability (GC0–GC3)

- GC0 PASS: OFFCAP ind_frac median **0.0165693** vs iter-5 target
  0.016569 — the cap-free path replicates iteration 5 exactly.
- GC2 PASS: max closure **0.000553** across all 174 runs with E_cap as
  a storage line; flow-identity resid max 2.80 absolute.
- GC3 PASS: max settled slope −0.0505 /s (everything decays), max
  sat_frac 0.2356, zero pool clamps.

## Riemann comparison (GC6) — NULL

Winner SETTLED base-12 scatter vs the ladder and the TP natural offsets,
GI4 method: med d_off 0.4513 (winner) vs 0.4475 (OFFCAP), Δ = +0.0038 —
wrong direction and 12× inside the control spread 0.0449 (winner
DETUNE/SHIFT arms + placebo). med d_line 0.7211 vs 0.7099. **No
alignment claim.**

## Trap re-measurement — the capacitor does not move the node

Node sweep at the winner capacitor, tone G, θ = 0 (accepted fraction):
d = 0.10 → 0.00685; 0.18 → 0.00252; 0.22 → 0.00146; **0.24 → 0.00127
(min)**; 0.28 → 0.00154; 0.40 → 0.00602. **d\* stays 0.24** (iteration 5:
0.24); depth 0.197 vs 0.164 (marginally shallower); the whole curve
sits ~2.4× lower — global loading, node geometry intact.

## The 0.3176 Hz exception

**Vanishes.** Zero SETTLED peaks in [0.25, 0.40] Hz in the OFFCAP,
PLACEBO, and winner scatters alike — absent even with no capacitor
under this iteration's G-tone reference condition (iteration 5's
sighting was under tone GP at d = 0.14). No capacitor effect claimed.

## Coordinates (every configuration; full table in the log CSV / JSON)

- Transducer: cavity A, x = 0.79 (bump σ = 0.06). Inducer (reference):
  A, x = 0.69, θ = 0°. Capacitor: per cell in
  `riemann_capacitor_grid.csv` (winner: x = 0.79, θ = 15°; placebo:
  x = 0.13, θ = 0°).
- Pickups: A_in1–6 at the sim's A probe xs; A_out aperture x = 0.97;
  B_in1–6 + B_out (cylinder probes); 16 ring channels A_rin0–7 /
  A_rout0–7 at x = {0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.92, 0.97}.
  The complete static table is in `riemann_capacitor_raw_meta.json`
  (`coordinates`), and per-run element coordinates in
  `riemann_capacitor_log.csv`.

## Limitations

- f_ring was measured at one placement (0.69, 0°); it is
  placement-dependent in principle. Audit6 also carried f_bare, and no
  cap label was needed anywhere, so this does not drive the null.
- P is a ratio at one pickup (A_in5, per the standing demod convention);
  persistence at other pickups was recorded but not criterion-ranked.
- One capacitance value only (second value cut, as pre-authorized);
  the d = 0.24 second inducer condition ran only inside NODESWEEP.
- Modal toy: no turbulence, no bench claim. Positives would be bench
  candidates; there are no positives here beyond the placement effect.

## Files

- `riemann_capacitor_sim.py` — engine + phased driver + `--calib`
- `riemann_capacitor_analyze.py` — audit6, clusters, GC verdicts
- `riemann_capacitor_results.json` — all verdict numbers, trap curve,
  alignment table, 0.3176 scan
- `riemann_capacitor_grid.csv` — placement search: coordinates,
  criterion values, ranking (51 cap-bearing cells incl. node sweep)
- `riemann_capacitor_scatter.csv` — every peak, audit6-labelled
- `riemann_capacitor_log.csv` — per-run ledger/stability/coordinates
- `riemann_capacitor_calib.json`, `riemann_capacitor_criterion.json`,
  `riemann_capacitor_raw_meta.json`, `*_peaks_raw.csv`, `*_series.csv`,
  `*_runs.csv` — supporting raw data
