# Quantum tiers from inside out — ratio tables v1 (2026-11-16)

Built on Zak's order (21:32 ET): the quantum tier was the
program's missing ratio table. Same discipline as the Earth
tables (enc-earth-inside-out-ratios-v1.md): public data only,
every value sourced and dated, vendor claims labeled as vendor
claims, derived ratios with declared arithmetic, non-comparable
metrics flagged in-row rather than forced, and holes printed at
the same size as findings. Zak's quantum architecture (Q1–Q4,
§21 sacrificed majority, the Jar, qubit-atmosphere) is FRAMING,
ungraded — it is not data and is used only to choose which
ratios to compute.

## 1. Platform ladder — the quantum inside-out (sourced values)

One row per platform-class system. Where a vendor publishes no
number, the cell is a HOLE, not a borrow.

| Platform (system) | Qubit transition freq. | T1 | T2 (echo where noted) | 1q gate time | 2q gate time | 1q error | 2q error |
|---|---|---|---|---|---|---|---|
| Superconducting — Google Willow, QEC chip (105 q) | ≈4–6 GHz class (transmon band; per-qubit frequencies are calibration data, no single published line) | **68 ± 13 µs** | **89 µs** | 25 ns | **≈42 ns** (CZ, "approximate" per Google device docs) | **0.035%** | **0.33%** (CZ) |
| Superconducting — IBM Heron r2 (ibm_basquecountry, 156 q) | ≈4–6 GHz class (transmon band, as above) | **253.1 µs** (median) | **204.6 µs** (median) | 32 ns (SX) | 68 ns (CZ) | 0.0361% (SX, median) | 0.183% (CZ, median) |
| Trapped ion — Quantinuum H2 (56 q) | **12.64 GHz** (¹⁷¹Yb⁺ hyperfine splitting — an exact atomic line) | HOLE (vendor publishes none; hyperfine ground states have no spontaneous-decay channel — see §2 note) | **≈2 s** in-machine, spin echo (secondary compilation of H-series operation; protected single-ion record ≈5.5×10³ s — Wang et al. 2021) | 63 µs (compiled, incl. QCCD context) | 308 µs (compiled, incl. QCCD context) | **3×10⁻⁵** (0.003%, vendor typical) | **1×10⁻³** (0.10%, vendor typical) |
| Trapped ion — IonQ Forte (36 q) | **12.64 GHz** (¹⁷¹Yb⁺ hyperfine, as above) | **10–100 s** (vendor) | **≈1 s** (vendor) | ≈130 µs (at execution, paper-reported) | ≈950 µs (ZZ, at execution, paper-reported) | **0.02%** (vendor) | **0.4%** (vendor) |
| Neutral atom — Harvard/QuEra class (⁸⁷Rb hyperfine) | **6.835 GHz** (⁸⁷Rb hyperfine clock splitting) | Not binding (ground-state hyperfine; no spontaneous decay — note, not a number) | **≈2 s** (hyperfine storage coherence, Harvard group) | HOLE (no single sourced platform value this pass) | 0.2–5 µs class (Rydberg pulse; secondary compilation) | ≈3×10⁻⁴–1×10⁻³ class (global/local, group-reported) | **0.5%** (CZ fidelity 99.5%, peer-reviewed) |
| Neutral atom — Atom Computing (⁸⁷Sr nuclear spin) | HOLE (nuclear-spin Zeeman transition is field-set, not a fixed atomic line) | Unbounded in principle (vendor: the nucleus "cannot dissipate" to other states) | **≈40 s** (vendor; "tens of seconds") | HOLE | HOLE (not published) | HOLE (no sourced value this pass) | **0.4%** (fidelity 99.6%, vendor) |
| Photonic (PsiQuantum / Xanadu class) | ≈193 THz class (optical; 1550 nm telecom → c/λ, derived) | NOT DEFINED for flying qubits | NOT DEFINED for flying qubits | NOT COMPARABLE (measurement-based / fusion-based; no stationary gate clock) | NOT COMPARABLE | NOT COMPARABLE | NOT COMPARABLE (error channel is photon loss + fusion failure, not a gate error rate) |

Sources, per row:
- Willow: Acharya et al. (Google Quantum AI), "Quantum error
  correction below the surface code threshold," Nature 638, 920
  (2025) — T1 68 µs, T2 89 µs, logical results in §3; Google
  Willow Spec Sheet (Dec 9, 2024, vendor) — 1q error 0.035%,
  CZ error 0.33%, 1q gate 25 ns; Google Cirq device
  documentation (accessed 2026-11-16) — CZ duration "approximate
  … 42 ns on Willow processors." Transmon band: Krantz et al.,
  Appl. Phys. Rev. 6, 021318 (2019) — plasma/transition
  frequencies typically 4–6 GHz.
- IBM Heron r2: backend calibration snapshot for
  ibm_basquecountry as reported in a cross-platform GHZ
  benchmarking study (arXiv, 2025–26) — medians T1 253.1 µs,
  T2 204.6 µs, SX 32 ns / 3.61×10⁻⁴, CZ 68 ns / 1.83×10⁻³.
  Calibration snapshots drift; this is one dated snapshot, not
  an IBM guarantee. A second Heron r2 backend (ibm_kingston)
  snapshot in a 2026 benchmarking paper reads T1 230.8 µs,
  T2 126.5 µs, 2q error 1.84×10⁻³ — same class, printed here so
  the drift is visible rather than hidden.
- Quantinuum H2: Quantinuum System Model H2 Product Data Sheet
  v4.00 (Sep 17, 2025, vendor) — 1q infidelity typical 3×10⁻⁵,
  2q typical 1×10⁻³, SPAM 1×10⁻³; the data sheet publishes NO
  T1, T2, or gate durations (that absence is itself the hole
  above). Peer-reviewed cross-check: Moses et al., Phys. Rev. X
  13, 041052 (2023) — 32-qubit configuration: 1q 2.5×10⁻⁵,
  2q 1.84×10⁻³. Gate durations compiled in a QAOA cross-platform
  benchmarking paper (arXiv:2510.12336, 2025, Table 3, from
  Quantinuum references): 1q 63 µs, 2q 308 µs — these bundle
  QCCD transport context and are labeled compiled, not
  vendor-published. In-machine T2 ≈ 2 s (spin echo): secondary
  compilation of H-series operation, 2025–26.
- IonQ Forte: IonQ Forte specifications page (vendor, accessed
  2026-11-16) — 36 qubits, 1q error 0.02%, 2q error 0.4%, SPAM
  0.5%, "T1 & T2 time: 10–100 s, ~1 s" (read as T1 10–100 s,
  T2 ~1 s; corroborated by OpenQuantum's listing, T2 ≈ 1,000 ms).
  Gate durations from a 2026 execution study on Forte / Forte
  Enterprise (arXiv:2506.07866): 2q ZZ ≈ 950 µs, 1q ≈ 130 µs
  "at the time of execution," median DRB fidelities 99.3–99.5%
  in that study vs the vendor's 99.6% headline — both printed;
  the study's numbers are the ones used in §2 and are labeled.
- Neutral atom (Rb): Evered et al., "High-fidelity parallel
  entangling gates on a neutral-atom quantum computer," Nature
  622, 268 (2023) — CZ fidelity 99.5% on up to 60 pairs in
  parallel (peer-reviewed). T2 ≈ 2 s: Harvard atom-array group
  materials (Bluvstein et al. lineage, 2022–23). Gate-duration
  class 0.2–5 µs: AMO Toolkit landscape compilation (2025–26,
  secondary). ⁸⁷Rb hyperfine splitting 6.835 GHz: standard
  atomic data, as stated in platform literature.
- Neutral atom (Sr): Atom Computing whitepaper (2025, vendor) —
  nuclear-spin coherence "tens of seconds," ≈40 s record class;
  2q fidelity 99.6%: Atom Computing (2024, vendor, reported with
  the Microsoft logical-qubit collaboration).
- Photonic: PsiQuantum et al., Nature 641, 876 (2025) —
  foundry-scale photonic platform demonstration. The row's
  NOT-DEFINED cells are definitional: a flying photonic qubit
  has no stationary T1/T2, and fusion-based architectures do not
  clock gates the way matter qubits do. Forcing numbers into
  these cells would be the incomparable-metric error; they are
  left open on purpose.

## 2. Derived ratios (declared arithmetic)

Definitions, fixed before the numbers:
- **Calm ratio (quantum form)** = T2 ÷ two-qubit gate time —
  coherent entangling operations that fit inside one coherence
  window. This is the Earth tables' calm ratio (coherence ÷
  disturbance) in the tier where it is native physics.
- **T2/T1** — dephasing structure: ≪1 means phase dies long
  before energy does; ~1–2 means relaxation-bound.
- **Error-threshold margin** = 1% ÷ 2q error rate, against the
  surface-code threshold of order 1% per operation (Fowler et
  al., Phys. Rev. A 86, 032324 (2012); the ~1% figure as used
  by Evered et al. 2023: "two-qubit gate error rates below 1%
  … are required to surpass quantum error correction
  thresholds"). Margin >1 = below threshold on this leg.
  The threshold is noise-model dependent; 1% is the canonical
  order, declared here so the margins are comparable.

| Platform | Calm ratio T2 ÷ t(2q) | T2/T1 | Error-threshold margin 1% ÷ ε(2q) |
|---|---|---|---|
| Google Willow | 89 µs ÷ 42 ns = **≈2.1×10³** | 89 ÷ 68 = **1.31** | 1 ÷ 0.33 = **3.0×** |
| IBM Heron r2 | 204.6 µs ÷ 68 ns = **≈3.0×10³** | 204.6 ÷ 253.1 = **0.81** | 1 ÷ 0.183 = **5.5×** |
| Quantinuum H2 | 2 s ÷ 308 µs = **≈6.5×10³** (conservative pairing; see note) | HOLE (no sourced T1/T2 pair from one source) | 1 ÷ 0.10 = **10×** (vendor typical); 1 ÷ 0.184 = 5.4× (Moses 2023) |
| IonQ Forte | 1 s ÷ 950 µs = **≈1.1×10³** | 1 ÷ (10–100) = **0.01–0.1** | 1 ÷ 0.4 = **2.5×** |
| Neutral atom (Rb class) | 2 s ÷ (0.2–5 µs) = **4×10⁵–1×10⁷** (duration is a class range, so the ratio is a range) | Not meaningful (T1 not binding; ground-state hyperfine) | 1 ÷ 0.5 = **2.0×** |
| Neutral atom (Atom/Sr) | HOLE (T2 ≈ 40 s sourced; gate duration unpublished — the ratio is not computable from public data) | Not meaningful (T1 unbounded in principle) | 1 ÷ 0.4 = **2.5×** |
| Photonic | NOT COMPARABLE (see §1) | NOT COMPARABLE | NOT COMPARABLE (loss-based error model) |

Notes on the pairings (declared, so the ranking can be argued
with): the H2 calm ratio pairs the conservative in-machine echo
T2 (≈2 s) with the compiled 2q duration (308 µs). Pairing the
same compilation's own T2 figure (60 s, a protected-memory
class number) with its own durations would give ≈1.9×10⁵ —
printed here as the spread, not used as the headline, because
the 60 s figure does not describe a machine running circuits.
The superconducting rows pair same-device snapshots. The Forte
row pairs the vendor T2 with a third-party execution-study
duration — the weakest pairing in the table, labeled as such.

**Calm-ratio ranking (headline):** neutral-atom Rb class
~10⁶ (range 4×10⁵–1×10⁷) ≫ Quantinuum H2 ≈6.5×10³ > IBM
Heron r2 ≈3.0×10³ > Google Willow ≈2.1×10³ > IonQ Forte
≈1.1×10³. Read with the margin column, not instead of it: the
ions win calm-per-gate among matter stationary qubits but pay
for it in gate clock (their windows are seconds because their
gates are ~10⁴× slower); the superconductors fit ~2–3×10³
entangling gates in a window that lasts ~10⁻⁴ s. Calm is a
ratio, not a duration — the same lesson as Earth's cavity
(§4 of the Earth tables): what matters is coherence measured
in units of the disturbance, and by that measure every
platform in this table lives between 10³ and 10⁷ — all far
calmer, relatively, than Earth's Schumann cavity at ≈1.5.

**Frequency span of the whole program (one line, declared
arithmetic):** lung 0.25 Hz → Schumann 7.83 Hz → transmon
≈5×10⁹ Hz → ¹⁷¹Yb⁺ hyperfine 1.264×10¹⁰ Hz is a span of
5.1×10¹⁰ ≈ **10.7 orders of magnitude**; extending the top
to the photonic row (193 THz) gives 7.7×10¹⁴ ≈ **14.9
orders** from the slowest body oscillator on record to the
fastest qubit transition in the table.

## 3. Sacrifice-ratio section — §21's public analog

Zak's §21 (framing, ungraded): in a qubit VOLUME, a large
majority is sacrificed as coherent medium so a working minority
can hold; its falsifier demands that HOLD QUALITY vs SACRIFICE
RATIO show a threshold. The public analog of "sacrifice" is
error-correction overhead: physical qubits spent so logical
qubits can hold. Working:sacrificed below = logical : (physical
− logical), computed from the published counts.

| Demonstration | Physical : logical | Working : sacrificed | Hold result | Source / status |
|---|---|---|---|---|
| Google Willow surface code, distance 3 | 17 : 1 (9 data + 8 measure) | 1 : 16 | Logical error ≈0.65%/cycle (derived from Λ, cross-checked) | Acharya et al., Nature 638, 920 (2025), peer-reviewed |
| Google Willow, distance 5 | 49 : 1 (25 + 24) | 1 : 48 | ≈0.31%/cycle (derived from Λ) | Same |
| Google Willow, distance 7 | 97 : 1 (49 + 48) | 1 : 96 | **0.143 ± 0.003%/cycle**; logical T1 291 ± 6 µs vs best physical qubit 119 µs — hold past breakeven, ratio 2.45× | Same |
| Microsoft + Quantinuum (H2) | 30 : 4 = 7.5 : 1 | 1 : 6.5 | Logical circuit error 10⁻⁵ vs physical 8×10⁻³ — headline 800× (vendor; later analyses put the improvement at 4.7×–800× depending on post-selection — the range is part of the record) | Microsoft / Quantinuum announcements, Apr 3, 2024 (vendor claims) |
| Microsoft + Quantinuum (H2, 56 q) | 56-qubit machine : 12 logical | ≈1 : 3.7 if all 56 counted (soft — allocation not fully published) | 12 entangled logical qubits | Vendor announcements, Sep 10, 2024 |
| Quantinuum Helios | 98 : 48 ≈ 2.04 : 1 | ≈1 : 1.04 | 48 logical qubits from 98 physical (¹³⁷Ba⁺); 2q error 7.9×10⁻⁴ (margin 12.7×) | Quantinuum, arXiv:2511.05465 (2025) + vendor |
| Harvard atom arrays | 280 : 48 ≈ 5.83 : 1 | 1 : 4.83 | 48 logical qubits, transversal gates, Nature demonstration | Bluvstein et al., Nature 622 (2023), peer-reviewed |
| Microsoft + Atom Computing | 24 logical entangled; physical subset used not cleanly published | HOLE | 24 entangled logical qubits (Nov 2024) | Vendor announcements |

Suppression law (Willow, the one complete public ladder):
each +2 of code distance divides logical error by
**Λ = 2.14 ± 0.02** — hold quality improving exponentially in
the sacrificed count, but ONLY because the physical error
(0.33%) sits below the ~1% threshold. Above threshold the same
spending makes hold worse. That is a threshold, measured, in
error space.

**§21 verdict (the honest shape of the evidence):**
- PRESENT: a threshold structure in public data — hold quality
  vs overhead follows the below/above-threshold law, and the
  distance ladders (Willow d3→5→7; Helios colour-code results)
  are genuine hold-vs-sacrifice curves along one axis
  (distance).
- ABSENT — HOLE: no public experiment sweeps the working :
  sacrificed POPULATION FRACTION at fixed code, fixed error
  budget, and fixed task — the §21 experiment as framed (a
  volume whose majority is spent as medium, hold measured
  against the spent fraction). The distance ladder varies
  overhead by growing one logical qubit's patch; it does not
  vary the population split of a working volume. **The §21
  sweep has not been run in public, and this table does not
  manufacture it.** The nearest public quantity remains Λ.

## 4. Inverse side — quantum shadow rows

The Earth tables' inverse side read periods and support values;
the quantum tier's reads the same way.

| Inverse row | Value | Source/status |
|---|---|---|
| Transmon period (1/f, 5 GHz mid-band) | **0.20 ns** | Derived from the 4–6 GHz band (Krantz 2019) |
| ¹⁷¹Yb⁺ hyperfine period (1/12.64 GHz) | **79 ps** | Derived from the atomic line |
| ⁸⁷Rb hyperfine period (1/6.835 GHz) | **146 ps** | Derived from the atomic line |
| Optical period (1/193 THz, 1550 nm) | **5.2 fs** | Derived (c/λ) |
| Coherence window as the tier's long period — shortest in table | 89 µs (Willow T2) | §1 |
| Coherence window — longest in table | ≈40 s (Atom Sr T2, vendor) | §1; span across platforms = 40 ÷ 89 µs ≈ 4.5×10⁵ ≈ 5.7 orders |
| Protected randomness — NIST Randomness Beacon | **512-bit pulse every 60 s** (≈8.5 bits/s), signed and hash-chained; v1.0 running since 2013, v2.0 current | NIST CSRC, Interoperable Randomness Beacons project pages (accessed 2026-11-16) |
| Protected randomness — CURBy (CU Boulder + NIST) | Device-independent certified randomness from a loophole-free Bell test; pulses chained across the network by the Twine protocol (CURBy-Q Bell chain feeding NIST-style chains) | "Traceable random numbers from a non-local quantum advantage," Nature (2025), DOI 10.1038/s41586-025-09054-3 |
| CURBy-Q certified-bit rate per pulse | HOLE | Not pinned from the primary text this pass; the existence and certification are sourced, the rate is not borrowed |

In the qubit-atmosphere framing (Zak, ungraded): the beacons
are the tier's public instruments of protected randomness — a
slow, certified pulse (one per minute) standing at the core of
a much faster uncertified traffic, which is the atmosphere
shape read at the scale of a national service.

## 5. What remains a hole after v1

- Quantinuum H2: vendor-published T1, T2, and gate durations
  (the data sheet omits all three); a single-source T1/T2 pair
  for the T2/T1 leg; the in-machine T2 rests on secondary
  compilations.
- Atom Computing: 2q gate duration (unpublished) → its calm
  ratio is not computable from public data; 1q error; the
  field-set transition frequency.
- Neutral-atom Rb class: a single platform's full same-device
  snapshot (the row is assembled from a peer-reviewed gate
  result + group coherence figures + a compiled duration
  class).
- IonQ Forte: gate durations rest on one execution study; the
  vendor publishes none.
- Photonic row: every gate-model cell — by construction, not
  by omission; a loss-model ratio table for photonics is a
  different instrument and has not been built.
- §21: the population-fraction sweep (see §3 verdict).
- CURBy-Q certified-bit rate per pulse.
- Helios full same-device coherence set (T1/T2 not in the
  values pinned this pass; Helios appears in §3 only).
