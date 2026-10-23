# Riemann Cavity Simulation Program — Iterations 1–7

A simulation study, run 2026-10-22, of two coupled modal cavities driven by a
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

## The seven iterations

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

## Pattern across all seven

Every element added — transducer, inducer, capacitor — acts on
*acceptance and persistence*: how much energy gets in, and how long it
stands. None of them composes new spectral identity. The system
redistributes; it does not compose. The ladder drive has never produced
a peak the audit could not explain.

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
