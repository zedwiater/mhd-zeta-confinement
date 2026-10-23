# RIEMANN_TRANSDUCER — local transducers at eddy sites (iteration 3)

**Order:** Zak, 2026-10-22 16:11 EDT — "just add transducer in eddy
regions and set proper limits on local transduction and induction from
eddy subspectra, try different locations each iteration and monitor
transduction stability and its effects on local." The order CORRECTED a
Gemini proposal: **no drive masking to eddy sites, no eddy injection,
TP flow stays on.** The iteration-2 system is bit-identical (geometry
imported from `riemann_tp_cylinder_sim.py`): TP cavity A + cylinder
cavity B, FP = 0.025, global Riemann-ladder antenna drive,
LADDER/DETUNE/SHIFT arms. Pre-registration: ITERATION_LOG.md 16:17 EDT;
calibration logged 16:24 EDT — all before the grid.

**HONESTY BOUNDARY: sim-internal only. No claim about the physical
world. The positives below are properties of this engine's cavities
with a local nonlinear element added — bench candidates at most.**

---

## 1. The transducer (model + hard limits, as pre-registered)

One local element per run, at one site. Profile: bump on the cavity
sample grid (A: 48-pt x-grid, Gaussian σ = 0.06; B: 64 volume samples,
Gaussian in Cartesian distance, σ = 0.35); modal profile
W = (2/N)·Φᵀb (the engine's projection convention). Sense
s(t) = W·v (site-profile velocity).

- **Extraction (saturating damper):** F1 = −G·V_SAT·tanh(s/V_SAT)·W.
  p_ext = G·V_SAT·tanh(s/V_SAT)·s ≥ 0 by construction. **Limits:**
  force cap G·V_SAT·|W|; conversion efficiency η = 0.5 (pool gains
  η·p_ext; (1−η)·p_ext is transducer dissipation, its own ledger line).
- **Re-injection (rectified-envelope push):** envelope state of |s|
  (τ = 0.25 s); F2 = +G·σ(e)·env·W, gated by the pool state e,
  σ(e) = e/(e+E_GATE), de/dt = η·p_ext − p_reinj, p_reinj exact signed.
  Pool clamped ≥ 0 at step boundaries — **0 clamps in 108 runs.**
- **Ledger:** E_inj = Δstored + E_diss + E_rad + E_tr_diss + Δpool +
  residual. Every transduced joule is drawn from the local field
  through F1; there is no other source.
- **Calibration (instrument, logged 16:24):** V_SAT = 10.9 (pooled
  driven-phase median |s| from G = 0 sensing runs), G = 250, G2 = G,
  E_GATE = 5.4e3. Measured operating points (LADDER NONE medians):
  extraction 1.42% (A79) / 1.47% (A61) / **2.83% (APORT)** / 0.44%
  (BCROSS) / 0.50% (PLACEBO) of E_inj.

**Sites** (provenance: run-1 + iter-2 audits): **A79** — A, x = 0.79
(run-1 A5, co-lit in the only unexplained+eddy cluster of runs 1–2;
iter-2 A_in5 eddy-lit 22, tied top). **A61** — A, x = 0.61 (run-1 A4,
top eddy-lit pickup of run 1; iter-2 A_in4 eddy-lit 22). **APORT** —
A, x = 0.97, the radiation aperture (iter-2 components 14/28 trapped
there, leakage 0.0). **BCROSS** — B at the B_in4 point (top B eddy-lit
pickup of iter-2; B is the sink of crossed modes 15/17/18).
**PLACEBO** — B at the B_in1 point (non-eddy control: quietest B
interior site, iter-2 SETTLED candidates 36 = lowest anywhere).
**OFF** — no transducer (iter-2 replication). Sensing is passive at
all sites in every run; actuation at one.

Grid: 6 configs × 3 arms × {NONE, AM@mean_gap} × 3 seeds = **108 runs**
(FM dropped, seeds 3 not 5 — stated in the pre-registration).
78,029 peaks logged.

## 2. Gates

| gate | statement | result |
|---|---|---|
| GT1 sanity | drive recovery ≥ 95% per config | **PASS — 0.985 in every config** (266/270; identical to iter-2's 0.984 — the 4 misses are a readout property, transducer-independent) |
| GT2 ledger | closure WITH transducer line ≤ 10% every run | **PASS — median 0.045%, max 0.055%, 108/108**; pool-balance residual ≤ 10.4 on pools up to 1.2e5; 0 pool clamps |
| GT3 stability | no runaway: growth slope / late growth / blowup | **PASS in substance** — all 108 runs decay (SETTLED E_loc slopes −1.5…−3.4, late ratios ≤ 0.05, pools bounded). Criterion (a) as literally written trips in 17/90 actuated runs — on **decay-slowing**, dissected in §3 |
| GT4 locality | effect counts only beyond OFF **and** PLACEBO seed ranges, same cell, same pickup | per site, §4: **A79 PASS, A61 PASS, APORT PASS, BCROSS marginal PASS, PLACEBO none (by construction)** |
| GT5 audit | unexplained at 0.45/0.30/0.15+flow, reproducible ≥ 2/3 seeds | **ZERO unexplained clusters in every cell of every config** — the uniqueness null survives iteration 3 |

OFF-arm replication check: exchange shares reproduce iter-2 to the
printed precision (8.50/3.64/2.69% vs 8.5/3.6/2.7) and OFF pulling
proxies are identical to iter-2's — the added machinery changed
nothing it wasn't aimed at.

## 3. Transduction stability (per site)

Operating points (LADDER NONE medians; full table in results JSON):

| site | ext % | sat_frac | gate_frac | pool_max | reinj net | decay-slope Δ vs OFF |
|---|---|---|---|---|---|---|
| A79 | 1.42 | 0.089 | 0.105 | 6.1e4 | −5.4 | none (in OFF range) |
| A61 | 1.47 | 0.097 | 0.104 | 6.3e4 | −8.1 | none (in OFF range) |
| APORT | 2.83 | **0.223** | 0.071 | 1.2e5 | +16.4 | **+0.23 /s** (τ 0.445→0.497 s) |
| BCROSS | 0.44 | 0.001 | 0.199 | 1.9e4 | −2.0 | +0.05 /s (3 of 6 LADDER runs) |
| PLACEBO | 0.50 | 0.000 | 0.171 | 2.2e4 | +0.3 | +0.08 /s (all LADDER runs) |

- **No runaway, no oscillation, no blowup in any run.** Every SETTLED
  local-energy slope is strongly negative; late-window (7.2–8.0 s)
  local energy is 0–5% of its 5.2–6.0 s level everywhere.
- **GT3 criterion (a) dissection (reported, not hidden):** the
  pre-stated tripwire compared the actuated SETTLED log-slope against
  max(OFF slopes)+0.05. It fired 17 times — always because the local
  decay got *slower*, never faster: APORT −2.01 vs OFF −2.245, PLACEBO
  −1.47 vs −1.56, BCROSS −1.82 vs −1.87. The mechanism is the pool
  return: the envelope push trickles stored energy back into a
  decaying tail. It fires at the placebo too, so the slowing is a
  **generic transducer effect, not an eddy-site effect**; it is
  largest where the transducer meets the exit channel (APORT) and
  where the local field is weakest relative to the return (B sites).
  The gate's intent was growth/runaway; what it caught is the opposite
  sign of danger. Recorded as a criterion artifact with the numbers
  kept.
- **Limits binding:** the force cap binds materially only in cavity A
  (sat_frac 0.09–0.22; APORT most — the port sense signal is the
  strongest in the system). At B sites the cap essentially never binds
  (≤ 0.001) — B's field is too weak to saturate a cap calibrated on
  the pooled median. The pool gate binds only early (gate_frac
  0.07–0.20 ≈ the charge-up phase); after that the pool is never the
  constraint.
- **The pool is a one-way accumulator at these settings:**
  re-injection net work over a full run is ≈ 0 (−8…+16 against
  E_ext ≈ 1.2e5 at A sites) — the envelope push does positive and
  negative work in near-equal measure while the pool keeps η·E_ext.
  As built, the element is a **saturating local sink with a
  time-redistributing push-back**, not a recirculating converter.
  Stated plainly because "transduction" in Zak's order implies
  conversion; the conversion channel exists and is ledger-exact, but
  its net throughput is ~zero in this parameter regime.

## 4. Local effects (the spatial story)

**SETTLED candidate-power ratio vs OFF** (LADDER NONE, median per
pickup; 1.00 = no change):

| config | A_in1 | A_in2 | A_in3 | A_in4 | A_in5 | A_out | B_in1–5 | B_out |
|---|---|---|---|---|---|---|---|---|
| A79 | 0.28 | 0.20 | 0.17 | 0.14 | **0.07** | 0.40 | 0.74–1.58 | 0.94 |
| A61 | 0.24 | 0.22 | 0.19 | **0.09** | 0.14 | 0.40 | 0.49–0.85 | 0.87 |
| APORT | 0.02 | 0.01 | 0.01 | 0.01 | 0.01 | **0.03** | 0.13–0.40 | 0.37 |
| BCROSS | 0.98 | 0.99 | 0.99 | 0.99 | 1.00 | 0.97 | 0.92–1.08 | 0.92 |
| PLACEBO | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.83–1.07 | 0.85 |

- **A-cavity transducers suppress ring-down candidate power
  cavity-wide, graded by distance from the site** (strongest at the
  site's own pickup: A79 → 14× down at A_in5; A61 → 11× at A_in4).
  The force is modal, so "local" was always going to be graded, not
  binary — the gradient is the locality evidence. DRIVE_END effects
  are much weaker (forced response; the damper works on the free
  oscillation, same drive/ring-down split as iter-2's pulling).
- **APORT is the extreme:** 30–100× suppression across all of A and
  3–8× into B. A transducer on the port damps the channel the
  ring-down lives in (A_out carried the largest SETTLED candidate
  population of any channel in iter-2). Its own pickup: candidate
  power 0.19 vs OFF [5.45, 5.64] — GT4 PASS by 29×.
- **GT4 (same-pickup placebo contrast, LADDER NONE):** A79 SETTLED
  site power 1.08e-5 vs OFF range [1.48e-4, 1.53e-4] and placebo-arm
  range [1.48e-4, 1.53e-4] — the placebo arm at the same pickup is
  indistinguishable from OFF, so the suppression is placement, not
  "a transducer anywhere". A61 identical in shape (1.18e-5 vs
  ~1.38e-4). BCROSS passes at its own pickup by ~3% over both range
  edges (3.69e-7 vs 3.57e-7 / 3.55e-7) — **marginal, negligible
  absolute size**, and its spatial profile is flat (~1.0): reported
  as a borderline pass, not an effect to build on.
- **Eddy concentration:** surviving SETTLED clusters at A sites light
  *fewer* neighbouring pickups (neighbour eddy-lit fraction: A79 0.60,
  A61 0.65, APORT 0.42 vs OFF 0.80/0.82/0.78 and placebo-arm
  0.80/0.82/0.78 — GT4 PASS for all three A sites) while cluster
  counts stay similar (11/11/8 vs OFF 11) and the eddy *share* of
  clusters rises (APORT 8/8). Reading: the transducers don't create
  or destroy the ring-down cluster population; they shrink its power
  and tighten its spatial footprint. (Caveat: cluster sets are
  re-derived per config; fractions compare different populations of
  similar size.)
- **B sites do none of this.** BCROSS and PLACEBO spatial profiles
  are ~1.0 everywhere; their eddy fractions (0.49 / 0.56) sit between
  OFF and placebo-arm baselines with no GT4 pass on the neighbour
  metric. At 0.44–0.50% extraction, a B transducer is a spectator.

## 5. Global effects

- **Exchange:** A-site transducers move the B share by ≤ +0.09 pp
  (noise). B-site transducers reduce it by ~0.3 pp (LADDER: BCROSS
  8.19%, PLACEBO 8.20% vs OFF 8.50%) — identical at eddy and placebo
  sites, i.e. placement-independent: the B dissipation channel
  (outside the share's denominator) skims B-bound energy before it
  dissipates/radiates. Bookkeeping-visible, small.
- **Leakage: NULL, with the seed-set trap named.** Raw iter-3 vs
  iter-2 deltas at site-strong components look nonzero (A79 −0.073,
  A61 +0.092) — and the OFF arm shows the *same* deltas (−0.071,
  +0.092): they are an artifact of iter-3 using seeds 0–2 against
  iter-2's 0–4. Attributable deltas (config − OFF): A79 −0.003,
  A61 −0.000, APORT −0.011, B sites 0.000. **No transducer changes
  the leakage ratio anywhere.**
- **Pulling (proxy — single-FP offset ratio, like-for-like code on
  all arms):** B-site transducers leave A's ring-down offsets
  *identical* to OFF (per-component ratios equal to 2 decimals).
  A-interior transducers suppress the tracked ring-down peaks so
  components ≥ 3 (A79) / ≥ 4 (A61) leave the search window — the
  proxy's median over the survivors drops; that is the suppression
  of §4 re-measured, not a slope change. APORT's surviving low-n
  peaks sit *further* past the drive frequency (ratios 6.26/5.73/
  4.92/4.55 vs OFF 6.16/5.46/4.65/3.77 for n = 1–4, then peaks
  vanish) — a proxy-level shift, mechanism unresolved, reported as
  an observation, not a measurement of the pulling slope.
- **Subspectra:** zero reproducible unexplained clusters at any
  tolerance in any config, arm, or window; site-pickup unexplained
  medians 0 everywhere. Adding a saturating nonlinear element at the
  eddy sites produced **no new spectral structure** — it removed
  power from the existing (fully audited) structure. The iter-1/2
  uniqueness null is intact.

## 6. Analysis bug found and fixed (transparency)

The first analysis pass implemented GT4's placebo contrast by
reading the PLACEBO config at *its own* B pickup — cross-cavity
candidate powers differ ~1000×, which silently masked the A-site
suppressions. Fixed to the faithful reading of the gate (placebo
*arm* at the same measurement point); all §4 numbers are from the
corrected pass. The bug was in analysis code only; no sim output
changed. (The v1 GT4 had also flagged BCROSS SETTLED — it survives
the fix, marginally, §4.)

## 7. Follow-the-data pivot: not used

Every pre-stated question received a discriminating answer from the
grid + the two controls; nothing was left ambiguous that a focused
pass could decide. The allowance is reported unused, per iter-1
precedent of logging the decision either way.

## 8. Limitations

- One transducer family: saturating damper + envelope push-back,
  η = 0.5, one gain for all sites. B sites are under-driven by that
  fairness rule (their sense signals run ~3× below A's against a cap
  calibrated on the pooled median) — a per-site gain normalization is
  a *different* experiment, not a fix to this one.
- The pool's near-zero net return is a property of this
  re-injection law (a sign-definite envelope push against an
  oscillating field); a phase-aware re-injection could recirculate
  for real. Not tested here.
- Modal forces are cavity-global by construction; spatial locality
  is measured (pickup gradient), not imposed.
- Seeds 3 (not 5): cluster reproducibility is ≥ 2/3; seed ranges are
  narrower than iter-2's, which is why GT4 contrasts lean on the
  OFF/placebo ranges of the *same* seed set throughout.
- Audit stops at 3rd-order combinations + flow terms, as before.
- Scale-free, sim-internal: even the §4 suppression needs the bench
  plan's frozen scale + ABBA protocol before it means anything
  physical.

## 9. Files

- `riemann_transducer_sim.py` — sim (model + gates in header);
  `riemann_transducer_analyze.py` — verdict code (incl. v2 GT4 fix).
- `riemann_transducer_pickups.csv` — 78,029 peaks (run, config, arm,
  modulation, seed, cavity, pickup, kind, window, freq, amp, drive
  flag, audit label). `riemann_transducer_transients.csv` — per-mode
  window energies / transfer ratios / delays (iter-2 schema + config).
  `riemann_transducer_series.csv` — per-run series (t, E_A, E_B,
  pool, env, s and E_loc at all sites). `riemann_transducer_runs.csv`
  — per-run ledger + transducer scalars. `riemann_transducer_log.csv`
  — final per-run log (site, arm, cap-binding fractions, stability
  metrics, local peak-table deltas vs OFF).
- `riemann_transducer_results.json` — all verdict numbers.
- `riemann_transducer_calib0_series.csv`, `riemann_transducer_grid.log`.
- `ITERATION_LOG.md` — pre-registration 16:17, calibration 16:24,
  verdict entry 16:47 EDT.
