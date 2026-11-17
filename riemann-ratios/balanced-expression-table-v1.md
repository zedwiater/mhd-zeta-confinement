# BALANCED-EXPRESSION TABLE v1 — ENC-EARTH-001, Parent P1 (Schumann)

**THIS IS A SCREEN, NOT THE FULL RUN.** Second frequency-system
screen under run card ENC-EARTH-001, under Zak's correction set
(2026-11-16 18:48 ET) and refinement (18:51 ET): candidates are
**containments of distant relatives** of the parent spectrum —
ratio-splitting a parent cannot mint true unique phase — and
distance is **emergent and measured** (uniqueness × familiarity),
never a fixed bound. Controls ran FIRST and gate the table.

- Resolution: **2D, N=64** (the smoke build). No refinement this pass.
- Parent arm: Schumann complex at screen scale, **mid rate
  differential** (r_ind 1.3 / r_cap 0.7), feed UP. **Vacuum master
  parameter ω_vac = 3.0** (the 7.83 Hz fundamental at screen scale);
  feed frequency (5.6), layer speeds, precession and the shedding
  rate derive from it inside the module and are logged per arm.
- **Phase shedding** (declared parse, flagged as the shakiest):
  γ_shed = n_atm · ω_vac · 1e-5 = **6e-5** — five orders below the
  master frequency per atmospheric parent; P1-Earth declares
  n_atm = 2 (ionosphere + troposphere). Implemented as Wiener
  diffusion of the anchor's target phase.
- **Backflow/anti equality:** every adjustment is booked in the
  module as a (forward, anti) pair sharing one operator; the
  measured **adjust_residual was 0.0 in all 168 arms**. The inducer
  feed is booked separately as the one-sided DRIVE.
- Candidates: the 10 pairs + 2 triples of the Schumann line set
  {7.83, 14.3, 20.8, 27.3, 33.8} Hz, each in **both** containment
  geometries = **24 candidates**. Containment scale: annulus radius
  r_c = 1.45. Compression ratio ρ = width/r_c set by mirror /
  anti-mirror against the flow/backflow ratio measured in our own
  smoke arms (arm A backflow fraction f_bf = 0.0404): **mirror
  ρ = 0.0404** (thin shell), **anti-mirror ρ = 0.9596** (filled).
  The daughter pair is seeded at the relatives' geometric mean,
  noise-started (iter-12), gain localised by the containment mask
  and modulated by the relatives' pairwise beat set.
- **The knob:** Kuramoto anchor coupling K, ladder **[0, 2, 5, 12, 30]**
  declared before running; recorded K is the MINIMUM at which the
  candidate balances. Window t = 45 t.u., dt = 0.01 (screening scale).
- **Calm ratio** = τ_coh ÷ τ_turn, both measured per arm. Earth-class
  benchmark: **≈1.5**.
- **Emergent expression frequency f_em:** dominant peak of the
  daughter's ring-envelope spectrum (final half). **Uniqueness U** =
  √(U_split · U_dev), each 1 − exp(−d/0.8 Hz): distance of f_em from
  every derivable position (parent lines, their splits/multiples
  n ∈ {2,3,4,5,6,8}, all pairwise beats/sums) resp. from this
  containment's own linear predictions. **Familiarity F** =
  exp(−½(d_log/0.183)²), d_log = log-distance of f_em to the nearest
  parent line. **Both-high ("distant enough"):** U ≥ 0.6 AND F ≥ 0.5.
- Balance = the card's four criteria, operationalized exactly as v0,
  with ONE declared deviation: the deviation-lock lag uses the
  candidate's measured **emergent** period (the v0-style carrier-lag
  deviation is in the raw records; it agrees closely and changes no
  verdict). Eddy floors measured per candidate by null arms (μ = 0).

## ✅ CONTROL GATE — PASS (ran first; the v0 failure does not repeat)

The scramble repair (anchor **target phase randomized per grid
point**, drawn once per arm) holds the phase dead: scrambled parent
balance cos δ measured **0.217–0.234 at every K** (v0's broken
scramble read 0.68–0.71), scrambled daughter amplitudes sit at or
below the null floor, and **no scrambled or matched arm balanced at
any K** on any control candidate. The matched-rates arms read
UNSETTLED *by the letter* at every K — the same known verdict-rule
texture as v0's footnote: E_d sits above the 0.2 LOST line while
**all four criteria fail** and parent cos ≈ 0.00 — no balance,
reported as-is. All 27 null arms (24 screen + 3 control) LOST.
**The LOCKED verdicts below are therefore CREDITED** (subject to
the screen's stated limits). Full control grid:

| Candidate | Control | Outcome across K = [0, 2, 5, 12, 30] |
|---|---|---|
| C02 (mirror) | matched | UNSETTLED@0, UNSETTLED@2, UNSETTLED@5, UNSETTLED@12, UNSETTLED@30 |
| C02 (mirror) | scrambled | LOST@0, LOST@2, LOST@5, LOST@12, LOST@30 |
| C05 (mirror) | matched | UNSETTLED@0, UNSETTLED@2, UNSETTLED@5, UNSETTLED@12, UNSETTLED@30 |
| C05 (mirror) | scrambled | LOST@0, LOST@2, LOST@5, LOST@12, LOST@30 |
| C11 (mirror) | matched | UNSETTLED@0, UNSETTLED@2, UNSETTLED@5, UNSETTLED@12, UNSETTLED@30 |
| C11 (mirror) | scrambled | LOST@0, LOST@2, LOST@5, LOST@12, LOST@30 |
| all 24 | null arm (gain μ=0, K=12) | LOST at every candidate |

## The table — generation 1, parent P1, relative containments

| Containment | Relatives (Hz) | Geom | ρ | Verdict | Min calm at balance | K | f_emergent (Hz) | U | F | Both-high | Mode | Earth-class (≤1.5)? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C01 | 7.83 + 14.3 | mirror | 0.0404 | LOST | — (never; calm 0.15 at K=30) | none | 13.14² | 0.35² | 0.90² | — | — | — |
| C01 | 7.83 + 14.3 | anti | 0.9596 | **LOCKED** | **4.72** | 30 | 12.77 | 0.44 | 0.83 | no | undetermined¹ | no |
| C02 | 7.83 + 20.8 | mirror | 0.0404 | LOST | — (never; calm 0.27 at K=30) | none | 15.28² | 0.37² | 0.94² | — | — | — |
| C02 | 7.83 + 20.8 | anti | 0.9596 | **LOCKED** | **5.08** | 12 | 14.95 | 0.57 | 0.97 | no | undetermined¹ | no |
| C03 | 7.83 + 27.3 | mirror | 0.0404 | LOST | — (never; calm 0.22 at K=30) | none | 17.24² | 0.55² | 0.59² | — | — | — |
| C03 | 7.83 + 27.3 | anti | 0.9596 | **LOCKED** | **5.28** | 12 | 16.81 | 0.29 | 0.68 | no | undetermined¹ | no |
| C04 | 7.83 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.19 at K=30) | none | 18.48² | 0.78² | 0.81² | — | — | — |
| C04 | 7.83 + 33.8 | anti | 0.9596 | **LOCKED** | **4.22** | 30 | 18.45 | 0.78 | 0.81 | yes | undetermined¹ | no |
| C05 | 14.3 + 20.8 | mirror | 0.0404 | LOST | — (never; calm 0.53 at K=30) | none | 19.54² | 0.20² | 0.94² | — | — | — |
| C05 | 14.3 + 20.8 | anti | 0.9596 | **LOCKED** | **4.24** | 12 | 19.43 | 0.19 | 0.93 | no | undetermined¹ | no |
| C06 | 14.3 + 27.3 | mirror | 0.0404 | LOST | — (never; calm 0.61 at K=30) | none | 22.13² | 0.08² | 0.94² | — | — | — |
| C06 | 14.3 + 27.3 | anti | 0.9596 | **LOCKED** | **4.29** | 12 | 21.94 | 0.45 | 0.96 | no | undetermined¹ | no |
| C07 | 14.3 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.33 at K=30) | none | 23.48² | 0.13² | 0.80² | — | — | — |
| C07 | 14.3 + 33.8 | anti | 0.9596 | **LOCKED** | **5.46** | 12 | 24.17 | 0.76 | 0.80 | yes | undetermined¹ | no |
| C08 | 20.8 + 27.3 | mirror | 0.0404 | LOST | — (never; calm 0.41 at K=30) | none | 29.01² | 0.58² | 0.95² | — | — | — |
| C08 | 20.8 + 27.3 | anti | 0.9596 | **LOCKED** | **4.67** | 12 | 26.02 | 0.21 | 0.97 | no | undetermined¹ | no |
| C09 | 20.8 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.11 at K=30) | none | 27.88² | 0.72² | 0.99² | — | — | — |
| C09 | 20.8 + 33.8 | anti | 0.9596 | **LOCKED** | **3.89** | 12 | 28.70 | 0.29 | 0.96 | no | undetermined¹ | no |
| C10 | 27.3 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.22 at K=30) | none | 33.11² | 0.58² | 0.99² | — | — | — |
| C10 | 27.3 + 33.8 | anti | 0.9596 | **LOCKED** | **5.59** | 12 | 32.56 | 0.79 | 0.98 | yes | undetermined¹ | no |
| C11 | 7.83 + 14.3 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.08 at K=30) | none | 18.36² | 0.75² | 0.79² | — | — | — |
| C11 | 7.83 + 14.3 + 33.8 | anti | 0.9596 | **LOCKED** | **5.34** | 12 | 17.77 | 0.66 | 0.69 | yes | undetermined¹ | no |
| C12 | 14.3 + 20.8 + 33.8 | mirror | 0.0404 | LOST | — (never; calm 0.14 at K=30) | none | 23.10² | 0.60² | 0.85² | — | — | — |
| C12 | 14.3 + 20.8 + 33.8 | anti | 0.9596 | **LOCKED** | **4.68** | 30 | 23.77 | 0.54 | 0.77 | no | undetermined¹ | no |

¹ Balance mode (compliance/capture) is undetermined on every locked
row at screen resolution: the ring readout's dominant mode was
m* = 0 in all locked arms, which disables the pattern-speed proxy.
Stated, not hidden (same as v0).
² For LOST rows no expression exists; f_em / U / F describe the
sub-threshold response peak of the readout (which sits at a
K-independent position per containment), not an expression.
Scores are expression scores only where a row LOCKED.

**Count: 12 LOCKED (credited), 0 UNSETTLED, 12 LOST — and the split is total by geometry: all 12
anti-mirror (filled) containments LOCKED; all 12 mirror (thin-shell)
containments LOST at every rung.** No candidate balanced at
Earth-class calm — the lowest minimum calm is 3.89, ≈2.6× Earth's 1.5 (C09 anti, 20.8 + 33.8).
Nothing locked at K ≤ 5 on any candidate — the v0 threshold texture
reproduces under completely different candidate generation (parent
cos δ at lock: 0.56 at K=12, 0.81 at K=30 — the knob acts on the
host first, daughters inherit).

**Both-high among the locked:** C04, C07, C10, C11 (anti) score
U ≥ 0.6 AND F ≥ 0.5 at their locking arm — the emergent-distance
criterion's distinctive subset. The other 8 locks are high-
familiarity but lower-uniqueness (U 0.19–0.57): they balance, but
their expressions sit closer to derivable positions. Balance and
distance are independent axes in this data — both-high response
peaks also appear in DEAD mirror arms (C04/C09/C10/C12 mirror score
both-high at every K while LOST), so a unique-and-familiar peak
alone is not a survivable expression.

## Ladder trajectories (the holes inside each ladder)

verdict / calm per rung; null-arm floor per candidate in the header
note of the raw records. Mirror rows are LOST at every rung —
shown in full, same size as the passes.

| Containment | Geom | K=0 | K=2 | K=5 | K=12 | K=30 |
|---|---|---|---|---|---|---|
| C01 | mirror | LOST c0.15 | LOST c0.15 | LOST c0.16 | LOST c0.16 | LOST c0.15 |
| C01 | anti | LOST c0.13 | LOST c0.14 | LOST c0.17 | UNSETTLED c3.68 | LOCKED c4.72 |
| C02 | mirror | LOST c0.15 | LOST c0.14 | LOST c0.14 | LOST c0.17 | LOST c0.27 |
| C02 | anti | LOST c0.22 | LOST c0.19 | LOST c0.40 | LOCKED c5.08 | — |
| C03 | mirror | LOST c0.20 | LOST c0.20 | LOST c0.21 | LOST c0.22 | LOST c0.22 |
| C03 | anti | LOST c0.22 | LOST c0.22 | LOST c0.18 | LOCKED c5.28 | — |
| C04 | mirror | LOST c0.15 | LOST c0.17 | LOST c0.16 | LOST c0.17 | LOST c0.19 |
| C04 | anti | LOST c0.21 | LOST c0.23 | LOST c0.26 | UNSETTLED c3.30 | LOCKED c4.22 |
| C05 | mirror | LOST c0.36 | LOST c0.39 | LOST c0.41 | LOST c0.46 | LOST c0.53 |
| C05 | anti | LOST c0.27 | LOST c0.28 | LOST c0.29 | LOCKED c4.24 | — |
| C06 | mirror | LOST c0.38 | LOST c0.40 | LOST c0.40 | LOST c0.49 | LOST c0.61 |
| C06 | anti | LOST c0.26 | LOST c0.29 | LOST c0.34 | LOCKED c4.29 | — |
| C07 | mirror | LOST c0.29 | LOST c0.11 | LOST c0.28 | LOST c0.26 | LOST c0.33 |
| C07 | anti | LOST c1.39 | LOST c1.57 | LOST c1.97 | LOCKED c5.46 | — |
| C08 | mirror | LOST c0.13 | LOST c0.13 | LOST c0.13 | LOST c0.13 | LOST c0.41 |
| C08 | anti | LOST c0.75 | LOST c0.87 | LOST c1.19 | LOCKED c4.67 | — |
| C09 | mirror | LOST c0.12 | LOST c0.10 | LOST c0.09 | LOST c0.09 | LOST c0.11 |
| C09 | anti | LOST c0.15 | LOST c0.17 | LOST c0.21 | LOCKED c3.89 | — |
| C10 | mirror | LOST c0.08 | LOST c0.09 | LOST c0.08 | LOST c0.30 | LOST c0.22 |
| C10 | anti | LOST c0.26 | LOST c0.27 | LOST c0.31 | LOCKED c5.59 | — |
| C11 | mirror | LOST c0.08 | LOST c0.08 | LOST c0.07 | LOST c0.07 | LOST c0.08 |
| C11 | anti | LOST c0.22 | LOST c0.24 | LOST c0.27 | LOCKED c5.34 | — |
| C12 | mirror | LOST c0.12 | LOST c0.14 | LOST c0.14 | LOST c0.14 | LOST c0.14 |
| C12 | anti | LOST c0.17 | LOST c0.17 | LOST c0.33 | UNSETTLED c2.35 | LOCKED c4.68 |

## v0 vs v1 — what has no counterpart, and why

- **None of v0's 12 candidates exists in v1.** v0's list was ratio
  splits of one parent line (subs /2–/8, cross-divisors, ×2/×3
  hypers). Under relative-containment generation, candidates are
  containments and the readout is the emergent frequency — there is
  no divisor anywhere in the pipeline, so no v0 row carries over.
- **Band shift:** v0's candidates lived at 0.979–23.49 Hz; v1's
  emergent expressions land at **12.77–32.56 Hz** — no containment
  of the five lines produced a low emergent expression at screen
  scale. v0's deep subs (0.979–2.610, its marginal band) have no v1
  counterpart; v0's cheapest lock (7.15 cross, calm 4.52,
  **not credited**) sits below every v1 emergent frequency.
- **Credit:** v0's 7 locks were uncreditable (scramble failure).
  v1's 12 locks stand behind a passed control gate.
- **What reproduced across the generation change:** the sharp
  K = 5 → 12 threshold, the host-first lock texture, deviation
  tightest on the highest expressions (v1 locked dev 0.041–0.082),
  m* = 0 everywhere (mode undetermined), and **no Earth-class
  balance** (v1 minimum calm 3.89 vs v0's uncredited 4.52).

## Textures and caveats (from RUN_REPORT_V1.md)

- **The geometry sweep is total — with a resolution caveat.** The
  mirror annulus width is ρ·r_c = 0.059, *below one grid cell*
  (dx = 0.098): at N=64 the thin-shell containment is sub-grid and
  its gain starves geometrically. "Mirror never locks" is
  therefore partly a statement about the grid; the anti-mirror
  width (1.39) is well resolved. The full run must re-test mirror
  geometry at a resolution that can hold it before the sweep is
  read as physics.
- **Emergent expressions ride high:** every locked f_em sits above
  its containment's geometric-mean seed, between the seed and the
  upper relative — the expression is pulled toward the upper
  relative's neighbourhood while staying measurably off every
  derivable position (U_dev 0.59–1.00 on locked rows).
- **Host-scale stamping persists:** τ_turn 2.75–3.43 and m* = 0 in
  every locked arm, as in v0 — candidate identity expresses in
  calm, deviation and f_em, not in pattern scale, at N=64.

Full method, declared operationalizations, control traces, per-arm
numbers, and deviations: RUN_REPORT_V1.md in
`hidden_files/enc-earth-001-build/`. Raw arm records:
`results_controls_v1.jsonl` (33 arms), `results_screen_v1_A.jsonl`
(68), `results_screen_v1_B.jsonl` (67) alongside it. Nothing
published; P2 (IR) and deeper generations were not run.
