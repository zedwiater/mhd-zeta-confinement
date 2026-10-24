# Riemann Cavity Simulation Program — Iterations 1–13

A simulation study, begun 2026-10-22, of two coupled modal cavities driven by a
"ladder" of tones whose frequencies are the imaginary parts of the first 30
nontrivial zeros of the Riemann zeta function, used as dimensionless numbers.
One cavity carries a toroidal–poloidal (TP) flow field; the other is a plain
cylinder. Every run is compared against two control arms — the same ladder
detuned by ×1.01, and the ladder shifted by half a gap — so that anything the
ladder appears to do has to beat its own controls before it counts.

**Scope, stated plainly:** this is a modal toy model (30 modes per cavity,
RK4 integration, cubic mixing as the only nonlinearity). Every result below
is internal to the simulation. Nothing here is evidence about any physical
system, and nothing here proves, disproves, or bears on the Riemann
Hypothesis itself. The related conjecture in the root of this repository is
stated formally in Lean 4; formal statement is not proof, and no result in
this folder should be read as one.

## Method discipline

- Gates and pass/fail criteria were written down before each run
  (see `ITERATION_LOG.md`), and gates that failed are reported as failed.
- An extended intermodulation audit (sum/difference tones, harmonics,
  subharmonics, modulation and flow sidebands) is applied before any peak
  may be called "unexplained."
- Energy ledgers must close for a run to count.
- Null results are kept, not rounded up.

## The first eight iterations

1. **Dual cavity** (`riemann_dual_cavity_*`, `RIEMANN_DUAL_CAVITY.md`) —
   Ladder-driven cavity pair, three modulation modes. **Verdict: NULL.**
   Of 14,047 candidate peaks, zero survived the audit unexplained in any
   arm; the ladder arm was in fact the cleanest of the three.

2. **TP flow + cylinder** (`riemann_tp_cylinder_*`, `riemann_fine_*`,
   `riemann_fp0_*`, `riemann_fp015_*`, `RIEMANN_TP_CYLINDER.md`) — Added
   toroidal–poloidal flow to cavity A. Subspectra NULL again, but the
   exchange between cavities is resonance-selected: the cylinder takes
   8.5% of output power under the ladder drive vs 3.6% / 2.7% in the
   controls, and the crossed modes are exactly the ones nearest a cylinder
   eigenfrequency. Settled peaks track the flow-shifted eigenfrequencies
   (frequency pulling, slopes ≈ +1). One pre-stated correlation gate
   (Spearman ρ ≥ 0.5) was missed as written (ρ = 0.389) and is recorded
   as missed.

3. **Transducers at eddy sites** (`riemann_transducer_*`,
   `RIEMANN_TRANSDUCER.md`) — Local saturating elements placed where the
   flow sheds eddies. Uniqueness NULL survives a third configuration.
   The element is a saturating sink with time-redistributing push-back,
   not a recirculating converter. Site transducers suppress ring-down
   candidate power graded by distance (11–100×); a placebo element at a
   quiet site changes nothing.

4. **Inducers in the eddy** (`riemann_inducer_*`, `RIEMANN_INDUCER.md`) —
   A driven element holding one standing tone near a transducer.
   Coupling vs distance is modal, not monotone: at standoff d = 0.20 the
   inducer nearly decouples (0.30% accepted energy) while the farthest
   placement couples ~5× better — a decoupling node, confirmed by an
   independent pickup ring. Scatter-to-line alignment: NULL, blind to
   which arm was running.

5. **Node mapping + true 2D spread** (`riemann_node2d_*`,
   `RIEMANN_NODE2D.md`) — The node was densified and mapped in two
   dimensions. **Verdict: the node is a frequency-dependent trap, not a
   fixed place** — its position migrates with tone frequency (0.20 at
   −5%, ~0.24 at the tone itself, 0.28 at +5%). A three-tone comb shows
   the node is a notch tuned to one component that reshuffles energy
   into the others rather than removing it. Uniqueness took its one
   dent: a single reproducible noise-floor cluster at 0.3176 Hz, with a
   named envelope-artifact suspect 0.0007 Hz away.

6. **Proximal capacitors** (`riemann_capacitor_*`,
   `RIEMANN_CAPACITOR.md`) — Lossless reactive storage elements swept
   over a 2D grid around the inducer and adjusted to their best standing
   position. Winner: the eddy site itself (x = 0.79, θ = 15°) — it
   trades acceptance for persistence (what gets in, stands longer). The
   capacitor produces no spectral footprint of its own, does not move
   the trap (node stays at 0.24), and the Riemann-line comparison is
   NULL. The 0.3176 Hz exception from iteration 5 vanishes under these
   conditions.

7. **Honing** (`riemann_honing_*`, `RIEMANN_HONING.md`) — Fifteen
   rounds of tuning the elements' positions in relation (inducer and
   capacitor micro-steps, the trap distance, capacitance \u00d70.5/\u00d72/\u00d74),
   hunting any anomaly under criteria pre-registered before round 1,
   with the Riemann-line comparison and full controls at every round.
   **Verdict: NO ANOMALY in 15/15 rounds** — 231 runs, 398,482 peaks,
   and zero unexplained *peaks* (not merely zero clusters) across
   60,450 audited candidates. The one candidate (round 14, capacitance
   \u00d74: alignment shift +0.1744 against a control spread of 0.1652)
   failed its fresh-seed confirmation round — the shift replicated
   (+0.1690) but the spread grew past it (0.2191), so it is a property
   of the \u00d74 extinction regime, not of the tuned geometry. Mapped en
   route, all engine-expected: capacitance is a pure loading dial; the
   trap is a hard floor; an inducer at \u03b8 = 45\u00b0 nearly decouples without
   any capacitor; cavity B's share of output is geometry-blind
   (0.0857\u20130.0881 in all 15 rounds).

8. **Third dimension + re-seed consistency** (`riemann_reseed3d_*`,
   `RIEMANN_RESEED3D.md`) — The cavities gained an axial dimension z,
   and the study turned to seed dependence: how consistent outcomes
   are across re-seeds, and how large a re-seed perturbation the
   system absorbs. The 3D engine reproduces the iteration-6 winner
   exactly in the collapsed (z-uniform) limit (accepted-fraction
   ratios 0.9995 / 1.00002). In true 3D the winner geometry
   transforms: the mid-plane is an axial node (acceptance
   0.00685 → 0.00019), the best-standing capacitor splits off-plane
   to z ≈ 0.3 / 0.7, and the trap survives but deepens (depth
   0.227). **Re-seed absorption: A_max = 0.0** — the baseline
   ensemble is so consistent (post-off B-share CV 0.26%) that even
   a 10% re-seed fails the pre-registered absorption rule on that
   metric; absorption is metric-split — the driven tone never sees
   the re-seed at any amplitude, while post-off exchange and
   ring-down slope see it fully, identically in all control arms.
   The best input protocol (STEP) is only marginally more
   consistent than baseline (composite 0.952); no protocol reaches
   seed-independence. The re-seed leaves no amplitude lines at any
   ring-down resolution — its signature is phase dispersion
   (circular SD 1.393 → 1.651 at full re-seed). The dual-clock
   ordering inverts for this readout: the settled clock is the
   tighter channel (CV 0.40% vs 10.9%). Anti-mirror layout: NULL.
   No anomaly: zero reproducible unexplained clusters in 170,893
   audited candidates (1,037,542 peaks); the energy ledger closes
   in all 281 runs. (Iteration 8's `*_results.json`, `*_trail.csv`,
   `*_runs.csv` and the updated `ITERATION_LOG.md` exceed the size
   limits of the upload routes available here and remain in the
   lab record; the report and run log in this folder state their
   numbers.)

## Pattern across the first eight

Every element added — transducer, inducer, capacitor — acts on
*acceptance and persistence*: how much energy gets in, and how long it
stands. None of them composes new spectral identity. The system
redistributes; it does not compose. The ladder drive has never produced
a peak the audit could not explain.

## Progress note — 2026-10-23 (iterations 9–13, and a trial now running)

Iterations 9–13 ran on 2026-10-23. Their full reports are written and held in
the lab record; they join this folder on the next bulk push (several files
exceed the per-file limits of the upload routes available here). Headlines —
sim-internal, under the same scope note as everything above:

9. **The plate** — A driven plate locks to the room mode at ratio 1 across a
   >2.4 Hz range; a subharmonic drive does not lock. A pickup placed at a
   phase node removes measurable phase diffusion for the protected mode,
   down to the no-noise floor. A corrective feedback loop is cheap and
   stable but does not correct — recorded as a null.

10. **The R-local key** — The plate's response signature verifies *locally*:
    all nine displacement families (geometry, room, protocol, drive family)
    are rejected with perfect separation. But a plain replay of the published
    drive verifies too — the signature authenticates a configuration, not an
    instance; there was no secret in the design.

11. **Key honing** — Minting the drive phases into the key (fixed φ\*) kills
    the replay (distance 16.45 → 34.62 against a threshold of 31.07). An
    adaptive optimizer with query access still crosses; widening enrollment
    opens a pocket (more enrollment made the key weaker); and with zero
    environmental wander both key systems collapse — the wander excitation
    is part of what the signature runs on.

12. **Wave-compute honing** — A nested subharmonic eddy (G/2) *anti-phases*
    the primary mode's phase residuals (correlation −0.9922, confirmed at
    fresh seeds) at roughly 1% of the primary's energy cost: a near-free
    phase mirror via cubic mixing. Drift-based separation through ratio
    statistics is a null (ratios cancel common-mode wander). Maximum wander
    band measured: the lock holds up to a peak wander between 0.625 and
    1.250 Hz.

13. **The wrapper** — A verifier that receives only a coarse, quantized
    location approximation of the generator: replay now fails at every
    layer, and the adaptive optimizer cannot cross through the coarse
    channel (quantization flattens the gradient it climbs; its best found
    point sits at true distance ~129). The coarse channel alone cannot
    verify everything — four axes are invisible to it, one of them
    inverted — and bonded observers pool insurance, not information: the
    pool ceiling is the coarse ceiling. A configuration's identity is now
    specified as a triple — (macro harmonic noise band, informational
    resolution, maximum environmental wander band): the *channel bands of
    MHD individuality*.

**Now running:** a five-day hands-off evolution trial of the stable 3D
cavity (started 2026-10-23), checkpointed against restarts, with daily
check-ins scored against evolution criteria registered before it started.
Final report due 2026-10-28.

**Implicated / pending:** iterations 9–13 reports, scripts and results join
this folder on the next bulk push (the large files need a route without the
per-file cap). The Lean 4 scaffold in the root of this repository is
unchanged by any of this — these are simulation results about a toy model,
and the scope note at the top applies to every line of this note.

## Files

- `*_sim.py` / `*_analyze.py` / `*_pivot*.py` — engines, audits, analyses.
- `RIEMANN_*.md` — per-iteration reports, written before any summary.
- `ITERATION_LOG.md` — the full pre-registration and verdict log.
- `*_results.json`, `*_grid.csv`, `*_log.csv`, `*_runs.csv`,
  `*_raw_meta.json`, calibration/criterion JSONs — run data.
- **Not included:** the multi-megabyte raw CSVs (`*_peaks_raw.csv`,
  `*_series.csv`, `*_scatter.csv`, full pickup/energy dumps, 1.3–79 MB
  each) exceed the 1 MB per-file limit of the upload route used here.
  They exist in the lab record; the reports and results files above are
  derived from them and state their numbers.
