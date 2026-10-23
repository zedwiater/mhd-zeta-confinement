# Riemann Iteration 5 — The Node, Densified: a Frequency-Dependent Trap, Mapped in 2D

Zak's brief (2026-10-22, incorporating fine-tuned Gemini suggestions):
densify the d = 0.20 decoupling node from iteration 4; test whether it
is a spatial wave node or a resonance trap (frequency dependence);
break iteration 4's reduced-pairing degeneracy for a true 2D
(x, theta) spread; comb test if budget allows. All three priorities
ran; nothing was cut. Gates GN0–GN6 and the node-vs-trap
discriminator were pre-registered in ITERATION_LOG.md (17:44 EDT)
before the grid; the one follow-the-data pivot (bracketing the two
edge minima) was pre-stated at 18:26 EDT before its runs.

System: iteration 4's engine (TP cavity A with flow + cylinder B,
global 30-rung Riemann-ladder drive, transducer at A79 with
iteration-3 limits, in-plane + out-of-plane rings, 3 seeds/cell),
with the inducer generalized: per-tone amplitude bookkeeping is
exact, and placements are literal (x, theta) pairs. 240 main runs +
21 pivot runs; 421,639 peaks. OFF replication is exact for the fifth
straight configuration (B-share 0.08583019239973681, digit-for-digit
vs iteration 4). Sim-internal claims only; modal toy.

## Priority 1 — the node is a trap, not a place

Accepted inducer energy fraction (E_ind/E_inj, medians, theta = 0,
LADDER arm) vs standoff distance d from the transducer site:

| d | G (20.2147) | GM = G−5% (19.2040) | GP = G+5% (21.2254) | MID (25.1397) |
|---|---|---|---|---|
| 0.02 | — | — | — | 3.867% |
| 0.05 | — | — | — | 1.864% |
| 0.08 | — | — | — | 0.544% |
| 0.10 | 1.657% | — | — | 0.135% |
| 0.12 | 1.331% | 0.369% | 6.400% | **0.123%** |
| 0.14 | 1.013% | 0.255% | 4.928% | 0.487% |
| 0.16 | 0.725% | 0.168% | 3.635% | 1.246% |
| 0.18 | 0.484% | 0.117% | 2.547% | 2.392% |
| 0.20 | 0.301% | **0.099%** | 1.691% | 3.894% |
| 0.22 | 0.181% | 0.115% | 1.034% | 5.618% |
| 0.24 | **0.127%** | 0.161% | 0.559% | 7.263% |
| 0.26 | 0.138% | 0.232% | 0.270% | 8.493% |
| 0.28 | 0.214% | 0.324% | **0.136%** | 9.055% |
| 0.30 | — | — | 0.159% | — |
| 0.34 | — | — | 0.704% | — |
| 0.40 | — | — | 2.911% | — |

Node metrics (pre-registered definitions; plateau = mean of the
sweep-end fractions; width = half-depth span, bracketing
interpolation — see disclosure note):

| tone | d* | min | plateau | depth | width | span |
|---|---|---|---|---|---|---|
| GM 19.204 | 0.20 | 0.099% | 0.346% | 0.287 | 0.110 | [0.147, 0.257] |
| G 20.215 | 0.24 | 0.127% | 0.772% | 0.164 | ≥ 0.097 | [0.183, >0.28] |
| GP 21.225 | 0.28 | 0.136% | 3.268% | 0.042 | 0.167 | [0.200, 0.367] |
| MID 25.140 | 0.12 | 0.123% | 4.589% | 0.027 | 0.137 | [0.043, 0.179] |

VERDICT — RESONANCE TRAP (the pre-registered discriminator reads
TRAP; the pivot then bracketed both edge minima, hardening it):
- The minimum MIGRATES with tone frequency: GM bottoms at d = 0.20,
  G in a broad basin centred ~0.24, GP at 0.28–0.30 — monotone
  outward across the ±5% pair. A fixed spatial node would sit at
  one d for all tones (criterion: ±0.02). It does not.
- MID's minimum sits at d ≈ 0.10–0.12 — off the GM→G→GP sequence
  entirely. MID lives in a different modal neighbourhood (0.129
  above rung 3, pre-disclosed confound) and couples 10–70× harder
  overall (plateau 4.59%, up to 9.06% accepted, vs G's 0.77%).
- Iteration 4's "node at d = 0.20" was real but tone-specific: at
  G's own frequency the basin is broad and its floor is at ~0.24;
  0.20 is where GM — not G — bottoms out. The d = 0.20 reading was
  the coarse 5-point sweep aliasing a wide, frequency-sliding basin.
- Depths are deep everywhere (0.03–0.29 of plateau): each tone has
  SOME standoff where it nearly decouples. What moves is where.
  Mechanism available in-model: the inducer radiates into a
  standing/flow-shifted modal field; the standoff at which its
  radiation destructively cancels at the transducer is set by the
  tone's wavelength against the local modal structure — i.e., an
  interference minimum that must move with frequency. The sim
  cannot distinguish "trap" from "moving interference null" beyond
  this; the discriminator's operational question (does d* move
  with f?) is answered unambiguously: yes, by up to 0.18 in d
  (≈ 1.5 basin widths) across a 31% frequency span.

Controls: placebo (in B, matched energy) remains inert (site demod
0.0083 vs 0.358 active at the same geometry — 43×). DETUNE/SHIFT
arms ran at d = 0.10 and d = 0.20 for tone G (their ind_frac values
are in the log CSV); the node analysis above is the LADDER arm, as
pre-registered.

## Priority 2 — true 2D: the degeneracy is dead; the node is a band

GN0 (degeneracy-fix check) PASSES at both levels:
- Profile level: two same-theta/different-radius placements now
  produce different drive profiles — minimum pairwise relative L2
  across all same-theta radius pairs = 0.610 (gate: > 0.05).
- Run level: (x=0.59, theta=90°) vs (x=0.71, theta=90°) accepted
  fractions 0.683% vs 2.263% — separation 79× the 3×pooled-SD line.
Iteration 4's radius-independent PU/PD pairing does not survive
into this engine (the placements here are literal (x, theta)
pairs, not projections of a shared pairing table).

Node shape at tone G (accepted fraction, %):

theta row at the node radius x = 0.59:
theta:    0°     30°    60°    90°    120°   150°   180°
frac:   0.301  0.107  0.165  0.683  0.202  0.401  1.287

radius column at theta = 0:
x:      0.59   0.65   0.71   0.85   0.91
frac:   0.301  1.013  1.964  1.883  0.869

Characterization: at x = 0.59 the accepted fraction stays at or
below G's half-depth threshold (0.449%) at 5 of 7 angles — in the
reduced pairing this is a LINE (an x-localized band extended in
theta), not a point and not a spiral arm (no monotone theta drift
of the minimum). But the band is strongly theta-modulated (12×
across the row): deepest at theta = 30° (0.107%), breached at
theta = 90° (0.683%) and gone at theta = 180° (1.287% — above even
the theta=0 plateau). So: an annular decoupling band whose depth
is set by angle to the flow, deepest ~30° off-axis, absent
downstream (180°). Site demod at A_in5 tracks the accepted
fraction cell-by-cell (0.069 at 30° → 0.323 at 180°); settled site
power ratio vs OFF stays 1.01–1.19 — the inducer remains a weak
perturbation at the site even where it couples best.

Out-of-plane skew (SETTLED ring-power share rout/(rin+rout), A):
share spans 0.157–0.326 across the 35-cell grid — leakage into the
out-of-plane ring is real and placement-dependent. Per radius, the
theta-range of the share is 0.137 (x=0.59), 0.044 (0.65), 0.057
(0.71), 0.144 (0.85), 0.147 (0.91). The maximum sits at theta = 60°
for x ≤ 0.65 and at 120° for x ≥ 0.71. Mirror-pair asymmetries
(theta vs 180°−theta): 30-vs-150 differs by only −0.011…−0.021 at
every radius (consistent small sign); 60-vs-120 flips sign with
radius (+0.036/+0.044 inside, −0.039/−0.071/−0.063 outside). So
the skew exists but is not a single flow-aligned bias: the
out-of-plane leakage pattern reorganizes across the eddy radius,
with the mirror asymmetry reversing through it. (Caveat: theta is
sampled on [0°, 180°]; the pairing is even in theta by
construction, so this is a statement about the reduced geometry.)

## Priority 3 — comb: the node is a G-notch, and MID inverts

Comb {G, MID, H} simultaneously, per-component amplitude matched
(4500/sqrt(3); per-component E_inj equal). Per-component accepted
fractions (net work per component — cross-component velocity terms
included, see caveat):

| component | at d = 0.10 (ref) | at d = 0.20 (node) | ratio |
|---|---|---|---|
| G | 0.554% | 0.0998% | 0.180 |
| MID | 0.0116% | 1.208% | 103.9 |
| H | 0.529% | 0.544% | 1.028 |
| total | 1.091% | 1.795% | 1.64 |

VERDICT — SELECTIVE (pre-registered criterion: spread 577× >> 3×):
- G is stripped at the node by exactly its single-tone factor
  (0.180 in comb vs 0.181 single-tone node/reference) — the notch
  is a property of the G component, undisturbed by the comb.
- H passes untouched (1.03×).
- MID inverts: at the reference distance its net accepted work is
  ~zero (0.0116%); at the node it becomes the dominant accepted
  component (1.208%). The comb's total accepted energy is HIGHER
  at the "node" than at the reference — the node is not a general
  decoupling point under multi-tone drive; it is a single-tone
  notch with a component reshuffle.
- Caveat (stated in pre-registration terms): per-component work
  E_ind_k = ∫F_k·v dt shares the total velocity v, so components
  are not independent acceptances; the claim is about net work
  per component, which is the bench-measurable quantity.

## Gates

- GN0 degeneracy fix: PASS (above).
- GN1 drive recovery ≥ 0.95: FAIL as counted, one cell of 225 —
  G_d0.10|DETUNE = 0.9444 (85/90). Same peak-crowding readout
  class as iterations 3–4 (a marginal 5-miss cell); tabulated in
  the results JSON, not waived.
- GN2 ledger closure: PASS — max 0.055%, median 0.045%, 240/240
  runs, now including per-tone inducer lines.
- GN3 transducer hygiene: FAIL as counted. Max saturation duty
  0.590 (MID_d0.28) vs the 0.30 gate. Mechanism, not mystery: the
  iteration-3 limits were calibrated at G-scale extraction; MID
  and GP couple 10–70× harder and drive the transducer into
  saturation duty in proportion to their accepted energy (top-sat
  cells are exactly the top-coupling cells: MID 0.42–0.59 at
  d ≥ 0.22, GP_d0.12 0.51). Every G-tone cell passes (max 0.175).
  Zero clamp events in 240 runs; pool bounded (max 2.87e5 in the
  hardest MID run); global relative energy slope positive in
  exactly 1 of 240 runs (+0.0022/s, COMB_d0.20 seed 0, bounded,
  closure 0.04%) — no runaway anywhere. The transducer is being
  overdriven by the stronger tones, not destabilized.
- GN4 placebo contrast: PASS (43×, above).
- GN5 uniqueness audit: FAIL as counted — by a single cluster.
  Of 421,639 peaks (53,336 audit candidates under the extended
  audit5: ladder, intermods, flow products, inducer tones incl.
  comb pairwise sums/differences), exactly ONE reproducible
  unexplained cluster: GP_d0.14, SETTLED, 0.3176 Hz, median amp
  0.0031 (noise-floor class), lit at 12 pickups (not an eddy),
  present in 3/3 seeds. It is NOT explained here. Two post-hoc
  observations, labelled as such and not counted: (i) the
  transducer envelope timescale gives 1/(4π·TAU_ENV) = 0.31831,
  within 0.0007 of the cluster frequency; (ii) the audit's hit
  test structurally excludes sub-0.5 targets, so any sub-0.5
  product is unauditable in this pipeline by construction. The
  uniqueness null has now survived five configurations with this
  single noise-floor exception, in the one cell whose tone sits
  closest above a rung boundary region (GP is 0.203 from rung 2).
- GN6 pivot discipline: the one pivot was the pre-stated
  bracketing extension (21 runs); it confirmed TRAP.

## Disclosures (analysis defects, per lab precedent)

1. First-pass node widths used the wrong interpolation segment at
   the half-depth crossings (extrapolated the inner segment,
   inflating widths and placing G's right wall at a fictitious
   0.343). Corrected with bracketing interpolation; both versions
   are in riemann_node2d_results.json ("node_widths_corrected"
   supersedes). d*, minima, and depths were unaffected.
2. G's basin right wall remains unbracketed (still below
   half-depth at the sweep edge d = 0.28); its width is reported
   as a lower bound. GP's bracketing runs (to d = 0.40) belong to
   GP only.

## Limitations

- Modal toy: positions are in the reduced (x, theta) pairing with
  cos(n·theta) weights; theta ∈ [0, pi] by evenness — no claim
  about a physical coil's full 3D field.
- MID and GP carry pre-disclosed rung-proximity confounds (0.129
  and 0.203 from rungs 3 and 2); their absolute coupling scales
  are inflated by modal neighbourhood, though the node POSITION
  measurement is a within-tone comparison and stands.
- 3 seeds; medians are tight (per-seed ind_frac spreads < ~30% of
  the median in sweep cells) but no formal intervals are claimed.
- The comb's per-component split is net work, not independent
  acceptance (cross-terms shared through v).
- The transducer limits are iteration-3 calibrations; a
  MID-scale recalibration would change sat_frac numbers (GN3) but
  not accepted-energy fractions, which are computed from the
  inducer force, upstream of the transducer.

## Files

riemann_node2d_sim.py (engine + grid driver; engine diff vs
iteration 4 confined to the inducer/pickup/analysis surface, core
integrate() byte-identical), riemann_node2d_analyze.py (audit5 +
all metrics), riemann_node2d_pivot.py (bracketing extension),
riemann_node2d_results.json, riemann_node2d_grid.csv (35-cell 2D
metrics), riemann_node2d_log.csv (per-run log, 240 rows),
riemann_node2d_pivot_runs.csv (21 rows), raw CSVs
(peaks/series/runs).
