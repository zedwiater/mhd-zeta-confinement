# RIEMANN_DUAL_CAVITY — Riemann-ladder standing waves in the dual cavity

**Order:** Zak, 2026-10-22 13:21 EDT (standing waves from the Riemann line in the
dual-cavity experiment; extra pickups for standing spectral eddies; log everything
plot-ready against the ladder; report any subspectra not explained by known mixing).
**Follow-the-data amendment:** 13:24 EDT — one focused follow-up pass authorized with
gates pre-stated before it runs. Used once (Pass 2, §5). No third pass.
**Pre-registration:** ITERATION_LOG.md, 2026-10-22 13:26 EDT (gates + grid, before any run).
Pilot adjustments logged 13:31 EDT; Pass-2 pre-statement logged 13:47 EDT.

**HONESTY BOUNDARY: everything below is sim-internal. Nothing here is a claim about
CERN, the LHC, or the physical world. A positive would have meant "candidate for
bench follow-up" — and there is no positive. The result is a null, reported unsoftened.**

---

## 1. What was built, and what was reused

`riemann_dual_cavity_sim.py` — two coupled modal cavities (A driven, B undriven,
coupled per-mode at 2% stiffness), 30 modes each, mode shapes φ_m(x) = sin(mπx)
(loaded cavity: eigenfrequencies need not be harmonic — the ladder supplies them).

Reused from lab stock, not reinvented:
- **audio_cavity_sim.py** conventions: modal oscillators with per-mode damping
  γ_m = f_m/Q (Q = 30), fixed per-mode detune (rng(7) uniform), RK4 stepping.
- **Cubic mixing as the intermodulation source** (the lab's superlinear mixing),
  implemented as a physical-space u³ term projected back onto the modes
  (pseudospectral, 48-point grid).
- **transducer-series dual-cavity picture**: cavity B receives energy only through
  the coupling (the residual second cavity); pickups in both cavities.
- Drive per mode normalized to equal resonant amplitude (A_TARGET), so no mode
  is privileged by the 1/f² response falloff.

The sim is **scale-free**: the first-30 Riemann zeros (t = 14.1347 … 101.3179,
from `~/workspace/riemann-tp/edge_frequencies.json`) are used directly as
dimensionless frequencies. No Hz claim is made anywhere; the Schumann anchoring
in the bench test plan is a free scale choice and is **not** used here.

**Extra pickups (the new ask):** 6 per cavity at x = 0.13 / 0.29 / 0.47 / 0.61 /
0.79 / 0.93 (12 total), each logging its own Hann/Blackman-Harris-windowed FFT
spectrum, so a peak can live at some pickups and not others (eddy localization).

## 2. Pre-stated sweep grid and gates

- **Arms (3):** LADDER (t_n) · DETUNE (1.01·t_n) · SHIFT (t_n + mean_gap/2 = +1.50315).
  Pre-stated caveat: ×1.01 preserves all additive coincidence relations by scaling
  (weak control); SHIFT breaks them (discriminating control).
- **Modulations (5 cells/arm):** NONE · AM depth 0.5 at rate 1.50315 and 3.0063 ·
  FM deviation 0.5 at rate 1.50315 and 3.0063 (rates = mean_gap/2 and mean_gap).
- **Seeds:** 5 per cell (drive phases + pickup noise). **75 runs** main grid.
- **Candidate:** spectral peak > 0.45 (freq units) from every driven frequency.
- **Intermod audit** (tol 0.45): fi+fj, |fi−fj|, 2fi, 3fi, fi/2, 2fi+fj, 2fi−fj,
  sidebands fi±r, fi±2r. First matching rule is the label.
- **Gates:** G1 sanity (≥95% drive recovery, else run void) · G2 audit fraction
  reported exactly · G3 verdict (ladder unexplained+eddy rate must beat BOTH
  controls by >2 pooled SD, else NULL) · G4 eddy reality (≥3/5 seeds, consistent
  freq, lit at 2–7 of 12 pickups).

## 3. Pilot adjustments (instrument, not criteria — logged before the grid)

- Pilot 1 (β=0.30, A=0.35): 360/360 peaks drive, zero candidates — cubic force
  ~10⁻³ of linear restoring force; nothing to find, by arithmetic.
- Pilot 2 (β=3.0, A=1.2, Hann): still zero over the floor. Spectral diagnosis:
  strongest non-drive content 5.3×10⁻³ abs sat inside a broad Hann leakage skirt
  (~5×10⁻³) shed by 30 strong tones. Generic cubic products are
  off-resonance-suppressed by 1/|ω_k²−ω_p²|; only near-mode products can surface,
  and the skirt buried exactly those.
- **Adjustment:** Blackman-Harris window (−92 dB sidelobes), β=10 (Duffing
  self-shift ~0.007 ≪ tol; stability margin ~150×), floor max(8×median, 0.2% of
  max), record window 3.5 s (bin 0.286). Pilot 3: 14 candidates in one run, all
  audit-explained — instrument in the discriminating regime. Gates unchanged.

## 4. Main grid — exact numbers

G1: **450/450 drive recoveries (1.000)** in NONE cells — readout sound, run valid.
Total peaks logged: **40,990** (CSV, plot-ready vs the ladder). Drive peaks 26,943.
**Candidate peaks: 14,047 — unexplained at tol 0.45: 0.**

Audit label counts over all candidates (first-match rule):

| label | count | | label | count |
|---|---|---|---|---|
| fi/2 | 3,390 | | 2fi+fj | 1,957 |
| fi+fj | 2,832 | | 2fi | 867 |
| fi−fj | 2,474 | | 3fi | 246 |
| 2fi−fj | 2,281 | | sidebands | (in above at this tol) |

Reproducible clusters (G4, ≥3/5 seeds): **385, all explained**, by arm:

| arm | clusters | unexplained | unexpl.+eddy | eddy fraction | median amp |
|---|---|---|---|---|---|
| LADDER | 126 | 0 | 0 | 0.34 | 0.0141 |
| DETUNE | 134 | 0 | 0 | 0.35 | 0.0180 |
| SHIFT | 125 | 0 | 0 | 0.34 | 0.0130 |

Structure that IS present (in all arms equally — generic cavity physics, not ladder
signal): near-resonant products (≤1.2 from a drive freq) run ~3.5× stronger than
far products (median 0.034 vs 0.0096); unmodulated cells produce ~0–2 reproducible
candidate clusters per arm vs ~60–70 under AM/FM — **the candidate population is
overwhelmingly modulation-born** (sidebands and their mixing products).

**G3 VERDICT (main grid): NULL.** Unexplained+eddy peaks per run: 0.0 ± 0.0 in all
three arms. There is no ladder enrichment to measure — the difference is 0 vs 0.

## 5. Pass 2 — the one authorized follow-up (pre-stated 13:47 EDT)

**Why:** the audit set holds ~2,700 combination targets over the 0–310 band; at
tol 0.45 their coverage is ~complete, so "explained at 0.45" is close to vacuous.
The discriminating question is whether peaks sit on combination tones at
*measurement precision*, and whether the nearly-clean NONE cells hide anything
resolution-limited.

**(a) Tight re-audit of all 14,047 main-grid candidates:**
- tol 0.30 → unexplained: **0**.
- tol 0.15 (≈ interp precision) → unexplained: **7 peaks (0.05%)**, split
  LADDER 2 / DETUNE 3 / SHIFT 2, all single-seed, **zero clusters** (none
  reproduces in ≥3/5 seeds). The candidate peaks sit on combination tones
  essentially to measurement precision; the stragglers are symmetric noise.

**(b) NONE-cell high-resolution runs** (3 arms × 5 seeds, record 6.5 s, bin 0.154,
audited at tol 0.15 from the start; G1 again 450/450): 702 candidates,
17 unexplained (2.4%) — LADDER **1** / DETUNE 5 / SHIFT **11**. Exactly **one**
unexplained+eddy cluster formed in the whole pass — in the **SHIFT control**:
f = 43.454, median amp 0.0015 (≈10× below typical cluster amps, i.e. noise-floor),
lit pickups A5/A6 only, seeds 0/1/3/4; nearest audited combination 43.171
(dist 0.283 — outside tight tol, inside loose). Most plausibly a higher-order
cascade product (4-wave mixing sits deliberately outside the audit set).

**PASS-2 GATE: NOT MET.** The gate required a LADDER-arm unexplained+eddy cluster
exceeding both controls. The ladder arm was the *cleanest* arm in the pass.

## 6. Final verdict

**NULL, hardened.** In this dual-cavity sim, driving with standing waves at the
Riemann-ladder frequencies produces **no spectral structure beyond intermodulation
of the drive set** — no unique subspectra, no ladder-specific eddies — and the
detuned controls produce statistically indistinguishable spectra by every
measure taken (cluster counts 126/134/125, eddy fractions 0.34/0.35/0.34,
amplitudes within noise, unexplained counts 0 everywhere at tol 0.45, symmetric
7-peak scatter at tol 0.15, and the single tight-tolerance eddy cluster landing in
a control arm). If the ladder has special standing-wave physics, this cavity —
at these Q, coupling, and nonlinearity settings — does not show it.

## 7. Limitations (stated, not buried)

- One cavity family: modal oscillators, Q=30, 2% inter-cavity coupling, cubic
  mixing only. Different Q/coupling/nonlinearity orders are different sims;
  this null covers this instrument, not all conceivable ones.
- The audit set stops at 3rd-order combinations; cascaded (4-wave) products
  exist in the physics and outside the audit — the one SHIFT cluster is the
  visible example, at noise-floor amplitude.
- Scale-free: even a positive would have needed the bench test plan's frozen
  scale + ABBA protocol before meaning anything physical. The null needs no
  such rescue: it is a statement about the sim, and it stands as one.

## 8. Files

- `riemann_dual_cavity_sim.py` — sim (grid + gates in header).
- `riemann_dual_cavity_pickups.csv` — main grid: 40,990 peaks (run, arm,
  modulation, seed, cavity, pickup, freq, amp, drive flag, audit label).
- `riemann_dual_cavity_pickups_none_hi.csv` — Pass-2 NONE-hi peaks (6,102).
- `riemann_dual_cavity_analyze.py` / `riemann_dual_cavity_pass2.py` — verdict code.
- `riemann_dual_cavity_results.json` — all verdict numbers + cluster lists + pass2.
- `ITERATION_LOG.md` — pre-registration (13:26), pilot adjustment (13:31),
  Pass-2 pre-statement + main-grid result (13:47), final verdict entry.
