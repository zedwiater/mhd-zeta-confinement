# RIEMANN_TP_CYLINDER — Riemann drive in a TP + cylinder-flow cavity pair

**Order:** Zak, 2026-10-22 14:35 EDT (TP and cylinder-flow cavities, pickups
inside and outside, Riemann-line source, document exchange / subspectra /
radiatives / transients). **Adjustments:** 14:44 EDT pilot fix (delay
estimator); 14:47 EDT (leakage ratio, pulling watch, fine transients,
asymmetry); 14:55 EDT second delta (leakage phase matrix, audit extension,
mode-splitting line). Pre-registrations: ITERATION_LOG.md 14:38 / 14:47 /
14:56 EDT — all gates stated before the runs they govern.

**HONESTY BOUNDARY: sim-internal only. Nothing here is a claim about CERN,
the LHC, or the physical world. Positives are bench candidates. There is
one positive (ring-down pulling, §7) and it is a property of the cavity,
not of the ladder.**

---

## 1. Geometry and approximations (stated, not buried)

- **Cavity A — TP (toroidal–poloidal).** Toroidal modes n = 1..30, shapes
  sin(nπx) on the unwrapped toroidal coordinate (run-1 basis). Plasma-frame
  eigenfrequencies = the 30-zero Riemann ladder (t = 14.1347 … 101.3179,
  dimensionless; scale-free, no Hz claim). **Poloidal flow** = rigid
  poloidal rotation, first order: lab-frame resonance of mode n shifts by
  Δf = FP·n, FP = 0.025. No flow shear, no flow-induced harmonic coupling
  (neglected, stated).
- **Cavity B — cylinder flow.** Cylindrical modes (m, nr, p), f ∝
  √(j²_{m,nr} + (pπ)²) — Bessel-root modal approximation, a = L = 1,
  lowest 30 (eigen span 20.12 … 85.97 after rescaling ×8.366 so the set
  median equals the ladder median). Shapes J_m(j r) cos(mθ) cos(pπz),
  RMS-normalized on 64 volume sample points. **Axial flow** = co-flow
  Doppler shift Δf/f = 0.02 × axial share (shares run 0 … 0.94);
  counter-flow branch neglected (stated).
- **Arms change the drive set only** (lab antenna): LADDER / DETUNE ×1.01 /
  SHIFT +1.50315. The cavity pair is identical in every arm.
- Coupling A↔B per mode index at 2% stiffness (modal approximation of
  aperture coupling). Cubic mixing per cavity in its own basis, run-1
  formula and constants (β = 10, A_TARGET = 1.2, Q = 30, RK4, DT = 0.0004).
- **Radiation port:** γ_rad = 0.25·γ_m extra damping on every mode; modal
  radiated energy is exact bookkeeping. Exterior pickup = port series
  Σ √（2γ_rad,m)·v_m·φ_m(aperture). 12 channels: 5 interior + 1 exterior
  per cavity.
- Protocol: drive 0→4.0 (0.5 ramp), off, record to 8.0. Windows EARLY
  [0.6,2.0], DRIVE_END [2.2,4.0], SETTLED [5.2,7.6]. Grid 3 arms ×
  {NONE, AM@3.0063, FM@3.0063} × 5 seeds = 45 runs; flow sweep FP ∈
  {0, 0.15} supplementary (15 runs each); fine-window pass (the one
  authorized pivot, 20 runs, windows [0.2,1.0]/[1.0,1.8]/[3.2,4.0]/
  [4.0,4.8]/[4.8,5.6]/[5.6,7.6]).

**Gates:** G1 = 443/450 = **0.984 PASS**. G2 ledger closure: median
**0.045%**, max 0.056%, 45/45 runs ≤ 10% — **PASS** by a wide margin.

## 2. Exchange through the system

NONE-cell medians (energies are in the sim's dimensionless units; see the
scale caveat below):

| arm | E_inj | E_diss A | E_diss B | E_rad A | E_rad B | B share of output |
|---|---|---|---|---|---|---|
| LADDER | 8.59e6 | 6.29e6 | 0.584e6 | 1.571e6 | 0.146e6 | 8.5% |
| DETUNE | 2.91e7 | 2.24e7 | 0.846e6 | 5.60e6 | 0.212e6 | 3.6% |
| SHIFT | 4.91e7 | 3.82e7 | 1.056e6 | 9.55e6 | 0.264e6 | 2.7% |

- Radiated share of total output = **0.200 in every arm** — exactly
  γ_rad/γ_tot = 0.25/1.25. A consistency check on the port bookkeeping,
  not a discovery.
- **Scale caveat (measured, mechanism labelled unproven):** absolute
  injected energy is strongly arm-dependent (SHIFT injects 5.7× LADDER).
  A diagnostic run shows SHIFT's A-modal energies ramp steeply with mode
  index (81 → 5.6e5 across the band) while LADDER's are comparatively
  flat (5.2e3 → 7.8e4) — consistent with the shifted set's cubic products
  (2fᵢ−fⱼ) landing back on its own drive set (the +1.503 offset is
  preserved by that combination), an engine coincidence effect. All
  cross-arm claims below therefore rest on **ratios and shares only**.
- **Which components cross** (transfer ratio = B/A modal energy at
  DRIVE_END, NONE, median over seeds). LADDER: **crossed (> 0.1): modes
  15 (2.23), 17 (2.35), 18 (0.97)** — B holds *more* energy than A at
  those indices (resonant sink; the matching A modes dip). Partial
  (0.01–0.1): 11 modes. Home (< 0.01): 16 modes. The three crossed modes
  are exactly the three smallest detunings to a B eigenfrequency
  (Δ = 0.56 / 0.73 / 0.40). Components above f ≈ 86 (ladder modes 24–30)
  have **no B modes to cross into** — the cylinder's spectrum ends at
  85.97 — and stay home (ratios 0.004–0.03).
- **G5 (resonance-selected crossing): NOT MET as pre-stated.** Spearman ρ
  between transfer ratio and −|detuning to nearest B eigen| = **0.389**
  (p = 0.034) — real but below the 0.5 gate. The extremes obey resonance
  selection cleanly; mid-band partial transfers (index-coupling chains +
  cubic mixing) dilute the rank correlation. Reported as measured.

## 3. Subspectra (DRIVE_END, intermod audit as run 1)

- 33,256 peaks logged (main grid, all windows); DRIVE_END: 14,486 peaks,
  **1,255 candidates, ZERO unexplained at tol 0.45, 0.30, and 0.15 — in
  every arm.** Audit labels: 2fᵢ−fⱼ 501, fᵢ−fⱼ 209, 2fᵢ+fⱼ 145, fᵢ/2 132,
  fᵢ+fⱼ 131, 2fᵢ 72, 3fᵢ 65.
- Reproducible clusters (≥3/5 seeds): 107, all explained — LADDER 44 /
  DETUNE 32 / SHIFT 31; eddy fractions (lit at 2–7 of 12 channels)
  0.43 / 0.28 / 0.39.
- **G3 VERDICT: NULL.** Unexplained+eddy rate 0.0 per run in all three
  arms. Changing the geometry to TP + cylinder did not produce unique
  subspectra; run 1's null survives the geometry change.

## 4. Radiatives

Dirtiness = candidate share of spectral power at DRIVE_END (median per run):

| channel | LADDER | DETUNE | SHIFT |
|---|---|---|---|
| A interior | 3.4e-4 | 3.3e-5 | 4.7e-5 |
| A exterior | 2.7e-4 | 2.4e-5 | 9.5e-6 |
| B interior | 4.9e-4 | 2.6e-4 | 5.0e-5 |
| B exterior | 2.8e-3 | 7.0e-4 | 1.4e-4 |

- **G6 (transducer prior — power clean, state lossy): holds for A, fails
  for B.** A's exterior is cleaner than its interior (2.7e-4 < 3.4e-4);
  B's exterior is ~5.9× *dirtier* than its interior (2.8e-3 vs 4.9e-4).
  Exterior drive-component recovery (NONE): A **0.933** (vs interior
  0.984 — mild state loss), B **0.493** — B's port hears barely half the
  drive set, essentially only what crossed into B. Absolute dirtiness is
  small everywhere (≤ 0.3% of power).
- Modal escape fraction (A, LADDER): median **0.198** per component ≈ the
  constructed 0.20 port share — consistency, not discovery.

**Leakage ratio per component** (outside/inside peak amplitude at the
drive frequency, DRIVE_END; exterior amps corrected to
displacement-equivalent — the raw port series is velocity-weighted and an
uncorrected ratio is a units artifact growing ∝ f; that correction is
applied everywhere in this section):

- **Cavity A (LADDER):** rises monotonically-ish from **0.06** (n=1)
  through ~0.9–1.3 (n=8–14) to **1.1–2.1** (n=15–30) — low components are
  trapped, high components punch through, driven by the port gain
  (√γ_rad ∝ √f) and aperture shape factors. Components **14 and 28 read
  0.0**: present inside, below the exterior floor in ≥ half the runs.
  DETUNE and SHIFT tables agree to ~3 decimals at low n and closely
  throughout — **leakage is a property of the cavity + port, not of the
  ladder** (same physical cavity in every arm).
- **Cavity B (LADDER):** 0.0 for 13 of 30 components (nothing measurable
  outside — the component never built up in B), with spikes exactly at
  the crossed components: n=18 → 1.21, n=21 → 2.58, n=26 → 1.35. What
  crosses into B also leaks out of B.

**Leakage phase matrix** (focused pass, NONE cells; H = X_out /
coherent-mean(X_in) at the exact drive frequency, displacement-corrected):

- **Cavity A: |arg H| ≈ 90° for essentially every component in every
  arm** (band medians −90.0 / −90.0 / −52…−60°). The exterior is a pure
  velocity (derivative) readout of the interior — the sign flips track
  the aperture shape sign sin(nπ·0.97). **No extra boundary phase
  selectivity in A**: phase is flat where the amplitude carried the
  structure (the opposite of the delta's guess, reported as measured).
- **Cavity B: phases scatter mode-to-mode** (LADDER by mode: +89.6,
  +88.9, −11.5, +43.5, −89.8, +89.8, −58.0, −109.4, …) **with the same
  scatter in DETUNE and SHIFT** — component-dependent but
  arm-independent: no ladder specificity in phase. H amplitude medians:
  A 3.67 / 3.79 / 3.75, B 1.81 / 1.82 / 1.82 (LADDER/DETUNE/SHIFT).
- Interior phase spread across pickups ≈ 103° for most modes — the
  standing-wave sign-split signature (pickups sit on both signs of the
  mode shape), not decoherence. Modes 1 and 12: spread 0.0 (all five
  pickups on one sign).

## 5. Transients — between cavities and inside

- **Bulk transfer delay (cavity t50 difference, NONE):** LADDER
  **+0.148 s** (A t50 0.400 → B t50 0.548); DETUNE +0.028 s; SHIFT
  +0.020 s. The LADDER arm has the *largest* bulk lag: resonant A modes
  fill fast and meter energy out through the coupling, while controls
  slosh. **Per-mode t50 delays are ≈ 0** (medians −0.008 / −0.016 /
  −0.012 s, at the 4 ms sampling limit): individual coupled channels
  respond essentially synchronously; the bulk lag is an accumulation
  effect, not a per-channel propagation delay.
- **Ring-down:** τ_A / τ_B (SETTLED energy fits, NONE): LADDER
  0.584 / 0.304 s; DETUNE 0.517 / 0.325; SHIFT 0.358 / 0.207.
- **Per-window candidates (main grid, all channels):** EARLY — LADDER
  377 / DETUNE 234 / SHIFT 294; DRIVE_END — 501 / 361 / 393; SETTLED —
  **1270 / 116 / 1503**. Ring-down is where the candidate population
  lives (free oscillation at cavity eigenfrequencies, which are not
  drive frequencies — most extremely in SHIFT, whose eigen set sits
  1.5 below its drive set). All SETTLED candidates are base-audit
  explained at tol 0.45 (see §8 for the tol-0.15 + flow extension).
- **Fine-window pass:** candidate counts by window — T1 [0.2,1.0] 1004,
  T2 [1.0,1.8] 853, DRIVE_END [3.2,4.0] 929, RING1 [4.0,4.8] **2247**,
  RING2 [4.8,5.6] 1900, SETTLED [5.6,7.6] 1159. The outside share of
  candidates falls monotonically after drive-off: 0.21 (DRIVE_END) →
  0.149 → 0.132 → **0.088** (SETTLED): the radiative channel quiets
  faster than the interior — the port reports the driven state, the
  interior keeps the ring-down.
- **Ring-down attribution (Gemini's topological-filtering prediction,
  modal version): NOT SUPPORTED in this engine.** In the fine-pass
  ring-down windows, 100% of non-drive peaks (inside n = 4619, outside
  n = 687) classify as drive-set members or audited products at tol
  0.45; the eigenmode fraction is 0.0% at *both* inside and outside
  pickups — there is no exterior eigenmode excess to find. Caveat
  stated: attribution at tol 0.45 against a ~2,700-target audit set is
  permissive, and the fine windows' FFT bin is 1.25; the §8 extension
  (tol 0.15 + flow terms) is the sharper instrument and its result is
  the same shape.

## 6. Asymmetry A→B vs B→A

From the coupling-power series (fine pass, NONE medians):

- Attack phase: mean P_AB = **+2.01e5** (into B), mean P_BA = **−2.01e5**
  — near-perfect mirror: during attack the coupling potential is
  quasi-steady and power flows one way, A → B. (SHIFT: +3.66e5 /
  −3.69e5.)
- Ring-down: P_AB +966, P_BA +60 (LADDER) — both positive as the
  coupling potential itself decays; net flow is a trickle next to the
  attack flux.
- **Backflow timing:** the B→A power (P_BA) peaks **0.028 s** after
  drive-off (LADDER), 0.040 s (SHIFT) — backflow onset is essentially
  immediate, against a forward *bulk* filling delay of 0.148 s.
  **The poloidal flow does not create a large directional delay
  asymmetry in this engine**: the asymmetry that exists is
  fill-slow / return-fast, and it is a coupling-accumulation effect
  present with the flow at its (perturbative) pre-registered strength.

## 7. Lock-in / pulling (the one positive)

Flow sweep FP ∈ {0, 0.025, 0.15}; SETTLED-window peak nearest each
LADDER component (A channels, NONE), offset regressed against FP·n:

- **Components 1–7: slopes +0.98 … +1.24** — the ring-down peak tracks
  the flow-shifted eigenfrequency essentially 1:1 as the flow changes.
  At FP = 0 the same peaks sit on the drive set. Median slope over the
  measurable components = **0.92** (pre-stated claim gate > 0.5: met,
  on the measurable subset). At FP = 0.15, 46.7% of components' peaks
  sit nearer the eigenfrequency than the drive frequency.
- Components ≥ 8: **inconclusive, not negative** — at FP = 0.15 their
  eigen shift (up to 4.5) exceeds the ±1.5 search window, so no peak
  was captured to regress.
- DRIVE_END peaks stay on the drive set at all FP (forced response).
  Reading: the *free* oscillation locks to the (flow-shifted) cavity;
  the *driven* response stays on the antenna. This is a property of
  any driven cavity with a flow-shifted spectrum — a bench candidate
  for the flow knob, not ladder physics.

## 8. Audit extension (flow terms) — Gemini's break prediction, tested

Extended set: fᵢ ± k·sᵢ (k = 1,2; sᵢ = FP·i), fᵢ±fⱼ ± s, split-pair
midpoints ((fᵢ+sᵢ)+(fⱼ+sⱼ))/2. Anything explained = **ENGINE-EXPECTED**,
never unique. Applied at tol 0.15 (base audit first):

- **DRIVE_END: zero migration** — there was nothing to migrate (base
  audit already 0 unexplained at 0.15, §3).
- **SETTLED (main grid):** flow-labelled peaks — LADDER **72** (22–25
  per cell), DETUNE 9, SHIFT **0**; **still-unexplained after
  extension: 0 in every cell of every set.** Fine ring-down: LADDER 70,
  SHIFT 41 flow-labelled; flow clusters (≥3/5 seeds): SETTLED LADDER 3 /
  DETUNE 1 / SHIFT 0; RING LADDER 9 / SHIFT 1.
- **Mechanism, stated plainly so this is not misread as a discovery:**
  in the LADDER arm the ring-down content sits at eigenfrequencies =
  drive + sᵢ *by construction* (eigen = ladder + FP·n, drive = ladder),
  so the flow terms fᵢ + sᵢ catch exactly the peaks §7 already measured
  (pulling). The LADDER-weighted counts are that identity re-measured,
  not an independent excess. In SHIFT the same eigen content is
  absorbed by the base combination set.
- **Verdict on the prediction:** the uniqueness null does **not**
  break. Advection sidebands exist precisely where the engine's
  construction puts them, every one is audit-explained, and nothing
  survives as unexplained in any arm, window, or tolerance.

## 9. Mode-splitting line

Amp²-weighted spectral width within ±2.0 of each component (A interior,
NONE), band medians:

| window | FP | LADDER n1–10 / n11–20 / n21–30 | SHIFT n1–10 / n11–20 / n21–30 |
|---|---|---|---|
| DRIVE_END | 0 | 0.000 / 0.504 / 0.380 | 0.001 / 0.698 / 0.429 |
| DRIVE_END | 0.025 | 0.000 / 0.641 / 0.355 | 0.001 / 0.693 / 0.037 |
| DRIVE_END | 0.15 | 0.001 / 0.707 / 0.366 | 0.000 / 0.674 / 0.073 |
| SETTLED | any | ~0.002–0.009 / ~0.02 / unmeasurable | same shape |

Low-band components are single clean peaks (width ≈ 0) at all flows.
Mid-band DRIVE_END width (0.37–0.71) is **neighbour crowding** (drive
spacing 1.2–5 with products in between), not a resolved doublet: it
rises mildly with FP in LADDER (0.504 → 0.641 → 0.707) and is flat in
SHIFT (0.698 / 0.693 / 0.674). SETTLED peaks are single and narrow
(~0.02); the high band is below the floor in ring-down. **No
flow-dependent mode splitting is resolved in this engine** — the
flow's spectral signature is the §7 peak *migration*, not a splitting
*width*. (DETUNE coverage: FP = 0.025 only.)

## 10. Limitations

- Modal approximations throughout: rigid-rotation flow at first order
  (no shear, no harmonic coupling by flow); Bessel-root cylinder set
  with a = L = 1 and a median rescale; index-matched 2% coupling as the
  aperture model; B shapes RMS-normalized on 64 sample points.
- FP = 0.025 is perturbative here: top-mode shift 0.75 vs modal
  linewidth ≈ 4.2. The FP = 0.15 sweep arm exists precisely to probe
  beyond that regime; the main grid's flow effects are small by
  construction and reported at their measured size.
- The audit stops at 3rd-order combinations + flow terms; 4-wave
  cascade products live outside it (run 1 saw one, at noise floor).
- Absolute energies are arm-dependent (§2 caveat); ratios/shares only
  across arms. Ring-down attribution is limited by audit density and,
  in the fine pass, by the 1.25 FFT bin.
- **Gemini prediction 1 (flow noise masking a subspectrum) is
  untestable in this engine**: the modal engine has no turbulence, no
  vortex shedding, no flow noise, and this report does not pretend
  otherwise.
- Scale-free: even the §7 positive needs the bench plan's frozen scale
  + ABBA protocol before it means anything physical.

## 11. Files

- `riemann_tp_cylinder_sim.py` — sim (geometry, gates in header).
- `riemann_tp_cylinder_pickups.csv` — main grid, 33,256 peaks (run,
  arm, modulation, seed, cavity, pickup, kind, window, freq, amp,
  drive flag, audit label). `riemann_tp_cylinder_transients.csv` —
  per-mode window energies, transfer ratios, delays, inj/rad per mode.
  `riemann_tp_cylinder_energy.csv` — stored-energy + coupling-power
  series. `riemann_tp_cylinder_raw_meta.json` — per-run ledger +
  eigen sets.
- `riemann_fp0_pickups.csv`, `riemann_fp015_pickups.csv` — flow sweep.
- `riemann_fine_pickups.csv` / `riemann_fine_energy.csv` — fine-window
  pass (19,254 peaks) + coupling-power series.
- `riemann_tp_cylinder_leakage_phase.csv` — complex leakage matrix
  (focused pass, NONE cells).
- `riemann_tp_cylinder_analyze.py`, `riemann_tp_cylinder_phase.py`,
  `riemann_tp_cylinder_audit2.py` — analysis.
- `riemann_tp_cylinder_results.json` — all verdict numbers.
- `ITERATION_LOG.md` — pre-registrations 14:38 / 14:47 / 14:56, pilot
  14:44, verdict entry.
