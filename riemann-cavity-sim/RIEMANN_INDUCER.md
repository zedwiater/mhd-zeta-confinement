# Riemann Iteration 4 — Inducers in the Eddy

Zak's order 2026-10-22 16:53 EDT. Pre-registration: `ITERATION_LOG.md`
2026-10-22 16:58 EDT (tone provenance, placements, grid, gates GI1–GI7
stated before any run). Calibration logged 17:04 EDT. One
follow-the-data pivot (GI3 re-measurement with the registered
instrument), pre-stated 17:26 EDT. Two harness bugs found and
disclosed (driver CSV writer; analysis reference config), neither in
the physics. Sim-internal claims only.

**Headline:** the uniqueness null survives a fourth configuration —
with a tone actively held inside the eddy, the extended audit still
explains every reproducible cluster in all 17 cells (GI5). The scatter
does **not** align to the Riemann line or the TP offset set beyond the
controls (GI4 null). The iteration's positive result is spatial:
proximity response is **modal, not monotone** — at distance d = 0.20
the inducer sits in a coupling node and effectively decouples, while
the farthest point tested (d = 0.40) couples ~5× better. The held
tone itself is real and phase-stable at the transducer site by the
pre-registered instrument (GI3 pass, pivot-verified).

## System (unchanged base + one addition)

Iteration-2/3 system imported bit-identical: TP cavity A (flow
FP = 0.025; lab eigenfrequencies = ladder + FP·n) exchange-coupled
(2% of A's stiffness) to cylinder cavity B; global Riemann-ladder
antenna (NONE modulation); DETUNE (×1.01) and SHIFT (= A eigens)
control arms; Q = 30, RK4, DT = 0.0004, record 8.0 s, antenna off at
4.0 s. **Transducer at A79 in every run**, iteration-3 limits verbatim
(saturating damper + envelope push-back, G = 250, V_SAT = 10.9,
EGATE = 5.4e3, η = 0.5). Replication anchor: iteration-4 OFF's B-share
of radiated output = 0.08583019239973681, equal to iteration-3 A79's
to all printed digits, and its SETTLED scatter medians reproduce
iteration-3 A79's exactly (Δ = +0.0000/+0.0000, §4) — the base system
is bit-identical.

**Inducer:** local driven element F = A_IND·W_ind·cos(2π f_tone t),
0.5 s ramp, **held the full 8.0 s** (SETTLED = antenna-off,
inducer-on). A placements: σ = 0.06 bump on the 48-point grid,
W = (2/48)Φᵀb, poloidal projection W_n *= cos(n·θ_ind) (reduced modal
pairing — position enters through the bump center x and one poloidal
phase; stated approximation). B placement: iteration-3 PLACEBO volume
profile at B_in1. Inducer work E_ind integrated exactly, its own
ledger line: E_inj + E_ind = Δstored + diss + rad + tr_diss + Δpool +
resid.

**Calibration (pre-stated band: E_ind/E_inj ∈ [1%, 3%]):**
A_IND_A = 4500 → 1.67% (E_ind ≈ 1.44e5); A_IND_B = 900 → 1.66%
(E_ind ≈ 1.43e5). Matched accepted energies = fair placebo. Accepted
fraction varies by placement thereafter (that variation is a result,
§3). Calibration closure 0.044–0.045%.

## 1. Tone provenance (Zak's rule: the tiniest reproducible eddy subspectrum)

**Ground unit G = 20.2147.** Rule (pre-registered): smallest
freq_median among reproducible (≥2/3 seeds) eddy clusters (lit at
2–7 pickups) in iteration 3's logged cluster table, LADDER / NONE /
SETTLED, pooled over iteration 3's six transducer configurations.
Winning instances: the A79 configuration's 20.2147 cluster (lit at
B_in1/B_in2/B_in4, 3/3 seeds) and the A61 configuration's 20.2144 —
the same component at two transducer sites. Provenance stated whole:
the component sits 0.097 from cylinder B's fundamental eigenfrequency
(20.1177, a ring-down-settled B mode) **and** within 0.0044 of a
2fᵢ−fⱼ intermod target — fully audit-explained in iteration 3. The
claim is only "the smallest released eddy component in the table" —
the base unit of the local series, not a new object.

**Alternatives logged, not run:** 30.9244 (APORT configuration,
smallest there), 33.4008 (A79), 43.9287 (A79, smallest A-interior-lit).
**Rejected family:** 36 near-DC clusters (smallest 0.5185), all in
EARLY/DRIVE_END — driven-phase envelope artifacts (iteration-3
precedent), none survive into SETTLED.

**Second tone H = G × (t₂/t₁) = 20.2147 × 1.487262 = 30.0646** ("first
ladder-ratio harmonic"), run as separate cells. Stated confound:
H sits 0.3603 from ladder rung 4 (30.424876), inside the drive-peak
tolerance — purity for H rests on coherent demodulation, and H's
behavior cannot be fully separated from rung-4 neighborhood effects.

## 2. Grid

17 cells × 3 seeds = 51 runs, 89,857 peaks, 16,311 non-drive/
non-inducer candidates. Pickups: the 12 base channels + in-plane ring
A_rin0..7 at x_k = frac(0.79 + k/8) + out-of-plane ring A_rout0..7
(same x_k, poloidal readout phase 1.0 rad — rational multiples of π
were avoided because they blind whole mode families). 28 channels.

- **Proximity sweep** (tone G, θ = 0, x = 0.79 − d): d ∈ {0.02, 0.05,
  0.10, 0.20, 0.40} → x = 0.77 / 0.74 / 0.69 / 0.59 / 0.39. Reference
  distance d_ref = 0.10.
- **Position iterations** at d ∈ {0.02, 0.10}: TP (Δx = +d, θ = 0),
  PU (Δx = 0, θ = +π/2), PD (Δx = 0, θ = π). The −d/θ = 0 position is
  the sweep point itself (reused, not re-run). θ = −π/2 was replaced
  by π pre-run: cos(nθ) is even in θ, so −π/2 would duplicate +π/2
  mode-for-mode (substitution logged in the pre-registration).
- **Controls:** OFF (transducer only); SWEEP_G_d0.10 in DETUNE and
  SHIFT arms; PLACEBO_B (tone G in cavity B at B_in1, matched E_ind);
  H_d0.02, H_d0.10.

## 3. Proximity and position verdicts (GI7)

LADDER medians over 3 seeds. "Site power ratio" = SETTLED candidate
power at A_in5 relative to OFF, same seed.

| cell | E_ind/E_inj | ext. share | sat_frac | demod @ site (DE/SET) | site power ratio |
|---|---|---|---|---|---|
| OFF | — | 0.0142 | 0.089 | — | 1.000 |
| d = 0.02 | 2.50% | 0.0324 | 0.180 | 0.644 / 0.519 | 1.596 |
| d = 0.05 | 2.33% | 0.0308 | 0.172 | 0.606 / 0.480 | 1.466 |
| d = 0.10 | 1.66% | 0.0255 | 0.146 | 0.480 / 0.358 | 1.210 |
| **d = 0.20** | **0.30%** | **0.0153** | **0.095** | **0.152 / 0.093** | **1.019** |
| d = 0.40 | 1.54% | 0.0258 | 0.145 | 0.461 / 0.363 | 1.194 |
| PLACEBO_B | 1.66% | 0.0142 | 0.089 | 0.037 / 0.008 | 0.998 |
| H d = 0.02 | 0.34% | 0.0155 | 0.094 | 0.130 / 0.050 | 1.030 |
| H d = 0.10 | 1.66% | 0.0227 | 0.118 | 0.280 / 0.266 | 1.372 |

**State change: d = 0.20 (x = 0.59) is a decoupling node.** Accepted
inducer energy collapses to 0.30% (≈ 1/8 of the adjacent placements),
tone amplitude at the transducer site drops ~4× vs d = 0.02, site
power ratio returns to 1.019 (≈ OFF), saturation duty returns to
0.095 (≈ OFF's 0.089). One step farther, d = 0.40 couples again
(1.54%, ratio 1.194). Proximity in this cavity is therefore **not a
distance law** — it is the modal overlap of the placement with the
eddy's field pattern, with nodes. Cross-confirmation from the ring:
in SETTLED at d = 0.10, tone-flagged peak amplitudes around the
in-plane ring run 0.007 (A_rin6, x = 0.54 — nearest pickup to the
node) to 0.444 (A_rin4) — the standing pattern itself has its null in
the same neighborhood.

Otherwise the sweep grades smoothly: closer → more accepted energy,
more exterior radiation (share roughly doubles at d = 0.02), a busier
transducer (sat_frac 0.180 vs 0.089), more settled site power.

**Positions matter less than the node.** At fixed distance, TP/PU/PD
move the site demod by ≤ ~±25% (d = 0.10: TP 0.426, PU 0.543,
PD 0.546, sweep-point 0.480) and the scatter deltas by ~0.001–0.003 —
no angular position changes state. Disclosed degeneracy: in the
reduced pairing the PU/PD placements depend only on θ, not on the
nominal radius, so the d = 0.02 and d = 0.10 PU (resp. PD) cells are
the *same configuration* run twice — their outputs are bit-identical
(n = 198/198, deltas equal), which doubles as a determinism check on
the whole pipeline.

**H tone behaves differently from G:** weak at the site (demod
0.13–0.28) but the largest exterior demod of any cell (A_out
14.0–30.7 vs G's 5.7–11.3) — H couples outward to the global
structure rather than into the eddy neighborhood, consistent with the
stated rung-4 proximity confound. Its site demod is flat across
ring-down at d = 0.10 (0.280 DRIVE_END → 0.266 SETTLED): held.

**Placebo contrast (decisive):** tone G driven in cavity B at matched
E_ind leaves cavity A untouched — A-site demod 0.037/0.008 (vs
0.480/0.358 for the eddy placement at d = 0.10), site power ratio
0.998, exterior share identical to OFF, B's transducer extraction
unchanged (~0.5% class, iteration 3). The tone's effects are
placement-specific, not a generic consequence of injecting 1.44e5
energy units somewhere. (Caveat from the pre-registration: PLACEBO_B
is a non-eddy site in the *other* cavity; "distance to the A79
transducer" is undefined across cavities — the contrast is
eddy-vs-elsewhere, not a distance point.) The held tone is not
contained, though: tone-flagged peaks appear at both exterior ports
(A_out amp 4.27, B_out 0.91, SETTLED medians at d = 0.10) — a local
inducer in this system radiates globally.

## 4. Scatter vs the Riemann line and the natural TP offsets (GI4)

Per Zak's ask, after every configuration the full-system scatter
(89,857 peaks) was compared against (a) the **Riemann line**:
distance from each candidate peak to the nearest ladder rung, and
(b) the **natural offsets**: distance to the nearest flow-shifted A
eigenfrequency (iteration 2's offset set). Primary window SETTLED
(the only window whose peaks are cavity-settled structure; DRIVE_END
medians sit at ~4.25–4.64 from both references in every arm — the
driven-phase scatter is rung-interleaved by construction and
discriminates nothing).

Medians, SETTLED, base-12 channels, candidates (pooled seeds):

| cell / arm | n | med d_line | med d_off | Δ_off vs iter-3 A79 | Δ_line |
|---|---|---|---|---|---|
| OFF / LADDER | 202 | 0.7232 | 0.4498 | +0.0000 | +0.0000 |
| d0.02 / LADDER | 194 | 0.7099 | 0.4453 | −0.0045 | −0.0134 |
| d0.10 / LADDER | 194 | 0.7099 | 0.4475 | −0.0023 | −0.0133 |
| d0.40 / LADDER | 197 | 0.7101 | 0.4486 | −0.0012 | −0.0132 |
| d0.10 / DETUNE | 31 | — | — | −0.0056 | −0.1590 |
| d0.10 / SHIFT | 168 | — | — | −0.0038 | −0.0054 |
| PLACEBO_B / LADDER | 137 | 0.7164 | 0.4472 | −0.0026 | −0.0068 |
| H_d0.10 / LADDER | 219 | 0.7205 | 0.4499 | +0.0001 | −0.0028 |

Reference (iteration-3 A79, SETTLED): med d_off LADDER 0.4498 /
DETUNE 0.9192 / SHIFT 0.4484; per-seed SD of the reference median
(pooled) = 0.00157. Alignment rule (pre-registered): a cell aligns
only if Δ_off(LADDER) < Δ_off(both controls) − SD.

**Verdict: NO ALIGNMENT (GI4 null).** At the reference distance,
Δ_off = LADDER −0.0023, DETUNE −0.0056, SHIFT −0.0038 — the ladder
arm is the *least* shifted of the three, the wrong direction for the
claim, and all shifts are the same small size as the placebo's
(−0.0026). What the inducer does to the settled scatter is real but
tiny and arm-blind: every configuration, controls and placebo
included, nudges the median candidate ~0.001–0.006 closer to the
eigen-offset set (and the G cells ~0.013 closer to the line). The
Δ_line shift is **not** audit bookkeeping: zero ind_-labelled
candidates appear in the SETTLED pool at d = 0.10 — the shift lives in
the ordinary peak population. DETUNE's Δ_line (−0.159) rests on
n = 31 SETTLED candidates and is not a stable median; flagged, not
weighed. Reproducible-cluster counts per cell (SETTLED, LADDER)
11–14 vs OFF 13, eddy-lit 1–5 vs OFF 4 — the inducer neither creates
nor destroys cluster families.

**GI5: zero reproducible unexplained clusters in all 17 cells** under
the extended audit (iteration-3 audit + flow targets + inducer
targets: tone, 2tone, tone/2, tone±fᵢ, 2tone±fᵢ, tone±k·sᵢ). Fourth
consecutive uniqueness null — now against an actively driven eddy.

## 5. Tone purity (GI3)

As pre-registered, purity = coherent demodulation at the tone (peak
finding is not the instrument). The demod was evaluated at A_rin0 —
which sits at x = 0.79, the A79 transducer site itself (ring pickup 0
coincides with A_in5's point and readout phase; the sweep grades
distance *to that site*, so this is the graded reference point).

- **First pass (wrong denominator):** demod amp at d_ref DRIVE_END =
  0.442–0.501 vs median *detected-peak* amplitude 0.213–0.224 →
  ratio 2.08–2.26 < 5 → FAIL as computed. Phase drift 13.7–18.9°
  (< 30° ✓, all seeds).
- **Pivot (registered instrument, pre-stated before running):**
  denominator = median of the window *spectrum*. Ratios: **20.5 /
  24.3 / 20.4** in DRIVE_END (3/3 seeds ≥ 5), phase drift unchanged →
  **GI3 PASS as pre-registered.** In SETTLED the ratio is ~3700× —
  after ring-down, the held tone is effectively the only coherent
  content left at the site. Spectral-peak check agrees: in SETTLED
  the peak nearest the tone sits within 0.0024 of it in every G cell,
  and the tone peak is the maximum-amplitude peak at the site
  (0.327–0.474 vs median candidate 0.002). During DRIVE_END the peak
  is pulled 0.10–0.16 off the tone by the driven comb — the demod,
  which is phase-coherent, is the right instrument there, as
  registered.

## 6. Stability (GI2, GI6, transducer carried from iteration 3)

- **GI2 PASS:** ledger closes with the inducer as its own line in
  51/51 runs; max closure fraction 0.055%.
- **Transducer unchanged in character:** sat_frac ≤ 0.236 in every
  run, zero pool clamps in 51 runs, pool bounded, extraction/reinjection
  accounting exact. Inducer proximity raises its duty cycle (0.089
  OFF → 0.180 at d = 0.02) without destabilizing it.
- **GI6 FAIL as literally pre-stated — mechanism named.** The gate
  required negative SETTLED local-energy slopes at A79. Inducer cells
  sit at ≈ −0.003/s (flat, G cells) to +0.017/s (worst: H_d0.10
  seed 2), vs OFF's −3.5/s ring-down. But SETTLED in this iteration
  is *inducer-on by design* — the site energy is **held**, not rung
  down; the criterion assumed a source-free ring-down. The hard
  stability facts all pass: the worst run's *global* E_A slope is
  −0.006/s with stored energy bounded in a narrow band (2.91e3–
  3.78e3, no growth trend), caps bind, clamps are zero, closure is
  0.05%. No runaway anywhere; the gate's letter fails because the
  regime it described no longer exists once a tone is held.

## 7. Gate scorecard

| gate | verdict |
|---|---|
| GI1 drive recovery ≥ 95%/cell | **FAIL as counted** for inducer cells (0.944–0.989; OFF/H/PLACEBO 90/90). Mechanism identified per run: rung 2 (21.02; DETUNE analog 21.23) merges with the tone peak 0.81 away (1.5 FFT bins at DRIVE_END) and is labelled inducer — peak-picker crowding, not a missing drive component. |
| GI2 closure ≤ 10% with inducer line | **PASS** — max 0.055%, 51/51. |
| GI3 tone hold (registered instrument) | **PASS** — demod/spectral-median 20.4–24.3× DRIVE_END, ~3700× SETTLED, drift ≤ 18.9°, 3/3 seeds (pivot-verified; first-pass proxy fail disclosed in §5). |
| GI4 alignment beyond controls | **NULL** — ladder arm not closest; shifts arm-blind, placebo-sized. |
| GI5 no reproducible unexplained clusters | **PASS** — zero in 17/17 cells. |
| GI6 no runaway | **FAIL as stated** (site-slope letter; held-tone regime), all hard indicators green — §6. |
| GI7 proximity state-changes reported with distance | **Reported** — decoupling node at d = 0.20 (x = 0.59); graded table §3. |

## 8. Limitations

- Modal toy, exactly as iterations 1–3: no real geometry, modes are
  indices; "distance" is bump-center separation on a 48-point grid;
  poloidal position is a single phase factor in a reduced pairing.
  PU/PD placements are radius-independent in this pairing (disclosed
  degeneracy, §3).
- A_rin0 coincides with A_in5 (same point, same readout) — the ring's
  *independent* channels are k = 1..7 and the out-of-plane ring; the
  out-of-plane ring is one fixed poloidal phase (1.0 rad), an
  approximation, not an axial scan.
- The scatter comparison's inducer-OFF reference is iteration 3's A79
  configuration (bit-identical system, same seeds) restricted to the
  12 base channels; the rings exist only in iteration 4, so ring
  channels enter within-iteration contrasts only. First analysis pass
  used iteration-3 OFF (no transducer) — caught by the replication
  check, corrected, logged; verdicts here use the corrected reference.
- Tone H's rung-4 proximity (0.3603) is an unresolved confound for
  every H-specific number.
- DETUNE SETTLED candidate populations are small (n = 31 at d_ref);
  its medians are indicative only.
- 3 seeds throughout (grid budget); per-seed SDs are correspondingly
  coarse. GI1's crowding artifact is a property of the peak-picker
  (distance = 3 bins) interacting with the 0.81 tone–rung spacing —
  a readout limitation, stated, not hidden.
- Sim-internal only: nothing here is evidence about physical cavities,
  the Riemann hypothesis, or any laboratory system.

## Files

- `riemann_inducer_sim.py` — engine (28-channel readout, inducer,
  ledger), `--calib` / `--calib-b` / `--grid` modes.
- `riemann_inducer_analyze.py` — audit4 (audit3 + flow + inducer
  targets), scatter, distance/alignment tables, gates.
- `riemann_inducer_pivot.py` + `riemann_inducer_pivot.json` — the GI3
  re-measurement (registered instrument).
- `riemann_inducer_scatter.csv` — full-system scatter, 89,857 peaks,
  audit-labelled, with d_line / d_off / d_offB per peak.
- `riemann_inducer_log.csv` — per-run ledger/stability/purity log
  (51 runs). `riemann_inducer_results.json` — all gate inputs.
- `riemann_inducer_peaks_raw.csv`, `riemann_inducer_series.csv`,
  `riemann_inducer_runs.csv`, `riemann_inducer_raw_meta.json` — raw
  outputs. `ITERATION_LOG.md` — pre-registration, calibration, bug
  disclosures, pivot pre-statement, verdict.
