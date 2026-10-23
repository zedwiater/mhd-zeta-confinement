# RIEMANN_RESEED3D — Iteration 8: third dimension + re-seed consistency

Zak's orders (2026-10-22 ~23:38 EDT; scope addition ~23:54;
refinement 2026-10-23 ~00:02): add a third dimension; define
re-seed transitions and consistencies, up to a maximal-consistency
input that makes re-seeding unnecessary or a bound on the maximum
re-seed amplitude; follow seed vs re-seed ring-down at increasing
resolution through the modulation playbook. Phase E (adjustable
inducers, two capacitors, denser pickups, anti-mirror layout) and
the tracker-trail / dual-clock / two-dial refinements integrated
as labelled sections. Pre-registration:
`riemann_reseed3d_prereg.json` (incl. Phase E addendum), fixed
before any run of each phase. Sim-internal claims only; modal toy.

281 runs, 1,037,542 peaks, 170,893 audited candidates.
**Zero reproducible unexplained clusters anywhere** (audit6:
base + flow + inducer + capacitor targets at f_bare = 10.10735
and the 3D ring f_ring = 22.448 Hz). Ledger: max closure
6.14e-04 over all 281 runs (gate 1e-3), zero pool clamps, every
settled slope negative (max −0.060), max sat_frac 0.147.

Definitions (pre-registered). SEED = rng(10000+seed) drive-phase
draw + readout noise, ensembles of 8. RE-SEED EVENT at drive-off
(first step t ≥ T_OFF): per mode of both cavities, the phase of
a_n = ω_n q_n + i v_n is rotated a fraction A_reseed toward a
fresh uniform draw; **per-mode energies exactly invariant**;
A = 0 is bitwise the baseline. CONSISTENCY = SD / max-deviation
across the 8-seed ensemble. A_max = largest A with every
post-reseed metric's SD ≤ 1.5× baseline SD and |mean shift|
≤ 1.5× baseline SD.

## Phase A — third dimension

Construction: cavity A modes gain axial envelopes χ_n(z) =
√2 cos(π m_n z), m_n = n mod 3, on a 48×16 (x,z) grid;
eigenfrequencies unchanged (the ladder is the spectrum; axial
structure is spatial only). KZ=0 collapses χ ≡ 1 and reproduces
the 2D engine verbatim.

**G0 — PASS, essentially exact.** Collapsed winner ind_frac
median 0.006846 vs iteration-6 0.00685 (ratio 0.9995); OFFCAP
0.016569 vs 0.016569 (ratio 1.00002); P median 0.8536 vs 0.8536.

What 3D does to the winner geometry (all numbers medians,
3 seeds):

- Acceptance collapses: inducer ind_frac 0.00685 (2D) →
  1.92e-04 (3D, z = 0.5 continuation). With z-uniform element
  profiles it falls further to 4.88e-05 (offset 0.0071×): a
  z-uniform profile couples only to the m_n = 0 modes (⟨χ_m⟩ = 0
  for m ≥ 1) — the exact mapping, stated as pre-registered.
- The mid-plane is an axial node for anything carrying a
  θ/z-symmetric bump: a z = 0.5 bump is symmetric about 0.5,
  χ_{m=1} is antisymmetric → overlap ≈ 0 with all ten m = 1
  modes. The A_in5@z0.5 pickup is likewise blind to m = 1 modes
  (DRIVE_END demod 0.110 vs the iteration-6 guard floor 0.240):
  **all kz1 winner-geometry P values are ratios in a suppressed
  regime**, flagged as in iteration 7, not standing persistence
  in the iteration-6 sense.
- Capacitor z-sweep (winner x, θ): acceptance peaks off-plane —
  z = 0.3 / 0.7: ind_frac 7.0e-04, P ≈ 1.187; z = 0.5 (node):
  1.9e-04; z = 0.1 / 0.9: 5.0e-04. Symmetric under z ↔ 1−z as
  the χ structure requires. **The best-standing capacitor
  position does not survive at mid-plane; it splits into the
  symmetric off-plane pair z ≈ 0.3 / 0.7.**
- Trap: survives, deepens, changes shape. Mid-plane curve:
  d = 0.10 → 1.90e-04, 0.18 → 1.05e-04, 0.24 → 5.12e-05,
  0.28 → 3.17e-05 (min), 0.40 → 8.95e-05; depth 0.227
  (2D: min at 0.24, depth 0.164). But at d = 0.10, moving the
  inducer off-plane reaches the floor immediately
  (z = 0.25: 2.65e-05; z = 0.75: 2.16e-05 — below the mid-plane
  minimum). The 2D trap point is a 3D basin: floor ≈ 2–3e-05,
  reachable along z at small standoff.

## Phase B — re-seed transitions and the amplitude bound

Winner geometry, BASE protocol, A ∈ {0, 0.1, 0.25, 0.5, 0.75, 1}
× 8 seeds. Baseline (A = 0) ensemble: B_share_post
0.10525 ± 0.00028 (CV 0.26%), slope_settled −1.0624 ± 0.0013,
demod_ST 0.11223 ± 0.00045, P 1.0348 ± 0.108, ind_frac
1.918e-04 ± 2.7e-06.

**A_max = 0.0.** The pre-registered rule fails at the first
nonzero amplitude, on one metric: **B_share_post**. At A = 0.1
its SD inflates 9.3× and its mean shifts 4.5 baseline SD; every
other metric passes at A = 0.1. The failure then spreads:
slope_settled joins at A = 0.25 (SD 2.9×, shift 3.3 SD); by
A = 1, B_share_post is 0.1264 ± 0.0132 (mean shift +77 SD),
slope SD 17×, and demod_ST fails marginally (SD 1.74×, shift
1.58 SD vs the 1.5 line). P passes at every A (SD ratio
0.95–1.0, shift ≤ 0.06 SD) — the persistence ratio does not see
the re-seed at all.

Mechanism (engine-internal, stated as such): the re-seed holds
every modal energy fixed and changes only phases. Post-off
cross-cavity exchange runs through the coupling κ(q_a − q_b)
and the cubic mixing, both phase-sensitive, so the ring-down
redistributes between the cavities along a different trajectory
for every re-phased seed — while the baseline ensemble is
nearly deterministic in B_share_post because the drive-phase
draw barely perturbs the settled modal-energy distribution.
The driven tone is immune because the inducer never turns off.

Controls replicate the signature exactly: B_share_post A = 0 →
A = 1 means — DETUNE 0.0407 → 0.0836, SHIFT 0.0363 → 0.0774,
placebo geometry 0.1062 → 0.1333, each with SD inflation of the
same order. The signature is arm- and geometry-independent.

**Answer to the brief's question**: re-seed necessity is NOT
prevented at the winner geometry, and the maximum absorbed
amplitude is **zero** — not because outcomes explode, but
because the baseline is so consistent that any phase
re-randomization of the residual field is immediately visible
in the post-off exchange ledger. The metrics split cleanly:
exchange/slope/scatter metrics carry the signature; driven-tone
metrics (P, and demod_ST until A = 1) absorb it.

## Phase C — maximal input consistency

Consistency composite (median SD ratio vs BASE over ind_frac,
P, B_share_out, d_off_settled; 8 seeds, A = 0):

| protocol | composite | notes |
|---|---|---|
| BASE | 1.000 | — |
| STEP (ramp 0) | **0.952** | winner, marginal: wins on B_share_out (0.85) + d_off (0.90); ind_frac SD 1.44× worse |
| FM @ mean-gap | 0.966 | uniform small gains |
| AM @ mean-gap | 1.008 | inert |
| SLOWRAMP (2 s) | 1.359 | ind_frac SD 9.4× worse |
| SHORT (T_OFF 2) | 2.627 | ind_frac SD 32× worse; settled scatter empty at the picker floor (d_off undefined) |

**Maximal-consistency input = STEP, nominally** — but the
honest statement is that no protocol approaches
seed-independence on the headline metrics: P's spread is
protocol-blind (SD ratio 0.95–1.01 everywhere), and the only
metric that is seed-independent at all (SD/|mean| < 1e-3) is
B_share_out — **already at BASE, in every protocol**. STEP at
A = 1 (Phase DW): full re-seed NOT absorbed (only P passes) —
the consistency finalist does not move the amplitude dial
either (see Two-dial section).

## Phase D — ring-down at increasing resolution

Windows ending 7.9 s, L = 0.9 / 1.8 / 3.6 s (Δf = 1.11 / 0.556 /
0.278 Hz), Blackman-Harris peak-picking, A_in5@z0.50.
Amplitude-line test (lines in ≥6/8 reseeded seeds, ≤1/8
baseline, audit6 before any claim):

- BASE A = 1 vs A = 0: **no lines at any rung.**
- STEP A = 1 vs A = 0: **no lines at any rung.**
- Dense array A = 1 vs A = 0 (Phase E4): **no lines**; ring
  census 10,101 settled candidate peaks on the dense ring
  channels, 0 unexplained.

Phase-domain (heterodyne) readout: per-mode cross-seed
circular SD at F_EIG_A on the 3.6 s rung rises from 1.393
(baseline) to 1.651 (A = 1) under BASE, and 1.392 → 1.651
under STEP — the re-seed IS visible in phase dispersion
(+19%), exactly the channel it acts on, while the amplitude
spectrum stays clean at ×4 resolution. The pre-stated pivot
(A_max = 1) did not trigger; the phase readout is reported as
the secondary channel regardless.

## Dual-clock capture (refinement §2)

Clock 1 = DRIVE_END intensity readout (demod_DE at A_in5@z0.5);
Clock 2 = SETTLED readout (demod_ST). Phase B ensembles,
paired by seed:

| | Clock 1 | Clock 2 |
|---|---|---|
| CV across seeds, A = 0 | 0.109 | **0.0040** |
| max paired |Δ| vs A = 0, A = 0.1 → 1 | 8.1e-09 → **1.7e-07** | 1.3e-04 → 2.1e-03 |
| mean shift at A = 1 (base SD) | 4.2e-07 | 1.58 |
| CV at A = 1 | 0.109 (unchanged) | 0.0071 |

Explicit answers:

1. **Re-seed visibility differs absolutely between clocks.**
   Clock 1 carries no signature at any amplitude (paired
   differences ≤ 1.7e-07 = numerical zero) — by construction:
   the event fires at the first step t ≥ T_OFF and the
   DRIVE_END window ends at T_OFF, so Clock 1 is a pre-event
   control. It validates the pairing; it cannot detect.
2. **All re-seed signal is on Clock 2**, graded monotonically
   in A (paired demod diffs 1.3e-04 → 2.1e-03), plus the
   exchange/slope signature of Phase B.
3. **Consistency verdict also differs, against the prior lab
   ordering.** For this readout (tone-G demod amplitude), the
   settled clock is ~27× tighter across seeds than the
   drive-end clock (CV 0.40% vs 10.9%): P's cross-seed spread
   is inherited from its Clock-1 denominator. The lab's earlier
   dualclock result (drive-end = higher fidelity, lower
   flapping) was established for broadband intensity on the
   harness telemetry — it does **not** transfer to the
   tone-demod readout. Same pattern under STEP (C1 CV 0.1099
   at both A; C2 0.0040 → 0.0071).

## Two-dial tuning outcomes (refinement §3)

Every adjustment in this iteration scored on exactly two
dials — (i) re-seed NEED: cross-seed consistency / approach to
seed-independence; (ii) re-seed AMPLITUDE: A_max absorbed.

| adjustment | dial (i) need | dial (ii) amplitude | verdict |
|---|---|---|---|
| STEP ramp | composite 0.952 (best) | A = 1 not absorbed | marginal on (i), inert on (ii) |
| FM | 0.966 | not finalist-swept | slight on (i) |
| AM | 1.008 | — | inert |
| SLOWRAMP | 1.359 | — | harmful on (i) |
| SHORT | 2.627 | — | harmful on (i) |
| E1 adjusted inducer (trap + z-offset) | composite 0.874 but P-SD 2.83× worse | A_max 0.0, unchanged | moves neither dial favourably |
| E2 second capacitor | 0.976 | A = 1 not absorbed | inert |
| E3 anti-mirror | 1.057 | A = 1 not absorbed | inert |
| E4 dense pickups | 0.999 | unchanged by construction | inert (readout only) |

**Bottom line: nothing moved the amplitude dial.** A_max = 0
at the winner geometry and at the adjusted-inducer geometry;
full re-seed fails absorption under STEP, two capacitors, and
the anti-mirror layout alike — always through the same
exchange/slope metrics. The need dial moves at most ~5%
(STEP). If a future iteration wants A_max > 0, the target is
named by the data: the post-off exchange channel
(B_share_post), not the tone channel.

## Phase E — scope addition

**E1 Adjusting inducers.** The adjustment set is Phase A's
TRAP/TRAPZ cells (standoffs × z-offsets, 3 seeds) plus an
8-seed re-seed sweep at the trap standoff d = 0.24 with
z-offset 0.25 (inducer (0.55, 0°, 0.25)). A_max there = 0.0,
same binding metric (B_share_post fails first). Moving the
inducer in z does not change the re-seed picture.

**E2 Two capacitors.** Winner cap (0.79, 15°, z 0.5) + partner
at the same (x, θ), z = 0.30. Ledger carries both storage
lines (E_cap ends, medians: q1 0.638, q2 0.186 — the off-plane
partner ends with ~29% of the primary's stored energy);
closure max 5.19e-04, gate held. Headline metrics ≈ single-cap
within seed noise; A = 1 not absorbed; dial-1 composite 0.976
(inert).

**E3 Anti-mirror layout** (named condition): (x, θ) → (x,
θ + 180°), z → 1 − z applied to transducer (0.79, 180°, 0.5),
inducer (0.69, 180°, 0.5), capacitor (0.79, 195°, 0.5).
**Null result.** Settled d_off medians: anti LADDER 0.4620 vs
standard 0.4575 (Δ +0.0045, ~100× inside the layouts' arm
spread: DETUNE 1.062 vs 1.037, SHIFT 0.439 vs 0.438); trap
points match the standard layout (d 0.10: 1.878e-04 vs
1.901e-04; d 0.24: 4.98e-05 vs 5.12e-05); exchanges,
persistence and B-shares within seed noise of standard in
every arm. The coherent (−1)^n sign flip of every element
profile leaves the system invariant at this resolution.
Consistency ensemble (8 seeds): composite 1.057 vs standard —
inert on both dials.

**E4 Denser pickups.** Dense array = 116 channels (A_in 25 =
5 x-positions × 5 z-planes {0.125…0.875}, A_out 5, B 6, ring
channels 80 = 8 positions × in/out × 5 z) vs base 40 (15/3/6/16,
3 z-planes). All Phase E focus configs ran dense. Dense
standard-winner ensembles (A ∈ {0,1} × 8) replicate the base
array's run metrics exactly where they must (dynamics are
pickup-independent; ind_frac/B_share ratios 1.000) and the
Phase D line test on dense A_in5@z0.50 finds no lines —
**the denser array changes no Phase D conclusion.**

## Tracker trail (refinement §1)

`riemann_reseed3d_trail.csv` — 281 rows (one per run), 50
ensembles: run seed, drive/reseed rng streams, protocol,
config_id (full 3D coordinates of inducer, both capacitors,
transducer angle, dense flag, κ count), ensemble_id, the
re-seed event (sim-time = T_OFF, A_reseed) and pre/post
headline metrics (Clock-1 demod as pre-event; Clock-2 demod,
P, B_share_post, d_off/d_line, slope as post). Every Phase B/C/E
number in this report is traceable to trail rows; ensemble →
run_id mapping is in `riemann_reseed3d_results.json["trail"]`.

## Limitations

- Modal toy, sim-internal claims only. "Re-seed" is a phase
  re-randomization at fixed modal energies — one specific
  transition class, pre-registered as such.
- A_max = 0 is a statement about this rule (1.5× tolerance)
  against an unusually tight baseline: B_share_post's baseline
  CV is 0.26%, so the bound is as much a property of baseline
  consistency as of re-seed strength. Absolute effect sizes
  are reported throughout for that reason.
- P at the 3D winner geometry is a suppressed-regime ratio
  (see Phase A); the A_max rule uses it as one metric among
  six and it is never the binding one.
- The runtime restarted twice during Phase E; lost and
  re-run: the first E_INDTRAP A = 0 ring-down segments (run
  rows intact) and the first E_DENSEWIN attempt (re-run in
  full). No run-metric rows were lost; rd_E contains the
  E_DENSEWIN ensembles used by the Phase D dense follow.
- SHORT's settled scatter is empty on the base-12 channels at
  the standing picker floor; its d_off/d_line are undefined,
  not zero.

## Files

- `riemann_reseed3d_sim.py` — 3D engine, re-seed event,
  protocol grid, Phase E drivers (+ `riemann_reseed3d_phaseE_resume.py`,
  `riemann_reseed3d_densewin_resume.py` — restart resumes)
- `riemann_reseed3d_analyze.py` — audit6, consistency/A_max,
  dual-clock, two-dial, Phase D ladder, trail builder
- `riemann_reseed3d_results.json` — all verdict numbers
- `riemann_reseed3d_trail.csv` — seed/re-seed lineage (281 rows)
- `riemann_reseed3d_log.csv` — per-run ledger + 3D coordinates
- `riemann_reseed3d_scatter.csv` — every peak, audit6-labelled
- `riemann_reseed3d_{runs,peaks_raw}.csv`, `riemann_reseed3d_rd_{B,C,DW,E}.npz`,
  `riemann_reseed3d_ring.json`, `riemann_reseed3d_prereg.json` — raw data + pre-registration
