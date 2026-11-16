# System Catalog — Tier-Proof Program

**Power and data flow across tiers, and what the language means against
the data, the physical phenomenon, and the assumed phenomenon.**

Added 2026-11-16, on Zak Wiater's order, as the system catalog for the
tier-proof program (source draft: *Sim Revision — Tier-Proof Program*,
Wave-resonance working files, 2026-11-16, DRAFT/ungraded). This catalog
is a map of what would have to be shown, with the receipts named in
advance. It is not a result. Every claim below carries its status:
**observed** (in the data), **assumed** (the presumed substrate
mapping), or **graded by** (the test that decides).

Zak's framing, verbatim (2026-11-16):

> "If everything is a system that is part of a system it's related to
> in physics and we only need the proof to align science grab a subset
> of public data that shows the bifurcation and super fluidic and
> solidic properties across tiers and show the natural and mathematical
> reciepts along the way including how Riemann and NS will function
> across those tiers and in a super fluidic and solidic universe how
> those systems make all the randomness connected through second sound
> particlization and transduction. Especially duality's expression
> across tiers and super states."

> "a qubit processing environment is atmospheric based stability from
> superposition atmosphere out in layers that stablize and protect
> superposition atmosphere where true randomness and uncertainty can
> exist and be calculated/expressed uniformly like our planet full of
> true random expression of photonic transduction and all its possible
> expressions."

---

## 1. The system catalog (public data subset, verified 2026-11-16)

"Natural receipt" = what nature measured. "Mathematical receipt" = the
structure the numbers must show for the mapping to count. Every row
carries its own scrambled control and falsifier.

| Dataset (custodian / access) | Tier · property | Natural receipt | Mathematical receipt | Scrambled control | Falsifier |
|---|---|---|---|---|---|
| ERA5 reanalysis (ECMWF/Copernicus CDS; free, registration) | Climate · bifurcation | Hourly wind/temperature fields, 1940–present, around regime shifts (blocking, monsoon onset) | Pre-shift variance + lag-1 autocorrelation rise (critical slowing), same exponent class across events | Same statistics on phase-randomized surrogate fields | Slowing absent, or identical in surrogates |
| NOAA CPC Oceanic Niño Index (NOAA CPC; direct download) | Climate · bifurcation (coupled oscillator) | ONI series 1950–present: the Pacific's own warm/cold band record, El Niño/La Niña flips | Dwell-time asymmetry + stable spectral peak across decades, consistent with a bifurcating oscillator | Time-shuffled ONI (keeps distribution, destroys phase order) | Flip statistics indistinguishable from shuffled |
| GISP2 / NGRIP ice cores (NOAA NCEI Paleoclimatology; direct download) | Climate, deep time · bifurcation | ~25 Dansgaard–Oeschger abrupt warmings in the last glacial, each within decades | Step shape + waiting-time distribution of a noise-driven bistable system, not a linear response | Smoothed/shuffled core series at same resolution | Steps resolve into gradual ramps at full resolution |
| Planck CMB maps + power spectra (ESA Planck Legacy Archive; public) | Cosmos · solidic (frozen pattern) | The early universe's pattern, frozen at recombination: angular power spectrum Cℓ, acoustic peaks | Peak ratios fixed by one parameter set across ℓ — a mode accounting | Gaussian random-phase skies with the same Cℓ | Any claimed cross-tier statistic that also appears in the random-phase control |
| GWTC catalogs + strain (GWOSC; public, CC BY) | Cosmos · duality | Mergers read as pure wave (strain) in the detector, as discrete events in the catalog | Waveform phase evolution matches GR chirp; event rates Poisson-clean | Time-shifted detector coincidences (standard background) | Claimed structure present in the time-shift background |
| Pulsar glitch tables (ATNF Pulsar Catalogue + Jodrell Bank; public) | Cosmos · superfluidic + solidic crust | A solid crust's spin, kicked by an interior it cannot touch mechanically: Δν/ν ~10⁻⁹–10⁻⁶, post-glitch relaxation | Two-component angular-momentum bookkeeping closes with a ~1% decoupled reservoir across pulsars | Crust-quake-only model populations | Bookkeeping fails without a decoupled (superfluid) reservoir |
| Johns Hopkins Turbulence Database (Johns Hopkins; public web-service/cutout) | Flow tier (NS anchor) · bifurcation of flow regimes | NS integrated honestly at Reλ≈430, 1024³ DNS with full time history | Inertial-range power law; anomalous but stable structure-function exponents | Gaussianized fields with identical spectrum | Cascade statistics reproducible by the Gaussian control |
| Sycamore processor datasets (Google Quantum AI; Dryad, DOI 10.5061/dryad.k6t1rj8) | Quantum · duality + decoherence | A working qubit volume's output, decoherence included (XEB circuits, bitstring samples) | XEB fidelity decay vs circuit depth follows the declared error budget | Classically spoofed / depolarized sample sets | Fidelity decay inconsistent with any single error-budget model |
| NIST Randomness Beacon + CURBy (NIST / CU Boulder; public API) | Quantum · true randomness | Signed 512-bit pulses at 60 s cadence; CURBy Bell-test randomness (Nature 2025) — randomness whose receipt is a Bell violation | Pulse statistics uniform + hash-chained; CHSH violated beyond the classical bound | Classical PRNG streams through the identical test battery | Any statistic separating beacon output from ideal randomness beyond stated bounds |
| Odlyzko zeta-zero tables (Univ. of Minnesota; public download; first 2,001,052 zeros) | Math tier · spectral accounting | Arithmetic's own spectrum — computed, not measured: zero spacings γₙ₊₁−γₙ after unfolding | Spacing distribution matches GUE random-matrix statistics (Montgomery–Odlyzko) | Poisson (uncorrelated) spacing ensemble at the same density | Spacings fit Poisson better than GUE |

**Honest null (kept in the catalog on purpose):** no public,
machine-readable He-II second-sound dataset could be verified on
2026-11-16. Second-sound velocities/attenuation exist as published
literature tables (e.g. Donnelly–Barenghi compilations). The quantum
tier therefore carries He-II as a **literature-digitization task with
page-level provenance per value**, not as downloaded data.

The subset's skeleton, visible before any gate runs: every domain has
one dataset where a system **flips** (bifurcation), one where a pattern
**holds** (solidic), one where a component **flows against another
without friction** (superfluidic), and one where the same thing is read
as **wave in one instrument and event in another** (duality).

---

## 2. Power and data flow table

Two different things flow in this program, and the catalog keeps them
in separate columns because confusing them is the classic failure:
**power** (energy — what NS routes inside a tier) and **data** (pattern
— what Riemann-style accounting tracks across tiers). Walls are where
both change form: transduction converts power across a wall, and only
*projections* of data ever cross one (the interior is read inductively,
from wall modulation — §20 doctrine).

| Flow path | Power flow (energy) | Data flow (pattern) | Transduction point / wall | Receipt that the flow happened |
|---|---|---|---|---|
| Within a tier (any fluid tier) | NS/MHD cascade: energy injected at large scales, routed scale-to-scale, dissipated at small scales; bifurcation routes (pitchfork, Hopf) decide which regime the flow settles into | The cascade's spectrum is the tier's mode census — who holds the energy, at what scale | None — flow stays inside the tier | JHTDB spectra; ERA5 fields; ONI phase record |
| Wave → particle (particlization) | A wave of the medium delivers its energy in discrete transactions: photon arrivals, glitch kicks, beacon pulses | Each transaction is one counted event in a waiting-time series | The detector / wall where the wave is absorbed | GWTC event catalog vs strain; glitch tables; NIST pulse chain |
| Particle → wave (settling) | Discrete events re-enter the medium as waves: post-glitch relaxation, thermal second sound, regime re-formation after a flip | Event statistics relax into a spectrum again | The medium's surface / crust / band boundary | Post-glitch relaxation curves; D–O step recovery shapes |
| Across a wall (transduction proper) | Power converts form across the wall's travel — heat ↔ wave ↔ work — never at a point | State does **not** cross. A projection crosses: wall modulation carrying the pattern of the interior | The wall itself, traveling | Wall-delta readout beats a no-cavity / no-wall control at recovering injected deltas (Q4 protocol) |
| Across tiers (accounting) | No direct power flow is claimed between tiers | Riemann accounting: each tier presents its modes as an unfolded spectrum; spacing law (GUE / Poisson / scrambled) is compared across tiers | The unfolding protocol — the same mathematics applied at every tier | §4 of the source draft: same non-Poisson law in ≥2 tiers under the same protocol, plus a joint scrambled control |
| Qubit atmosphere, outward | Decoherence is the only thing removed and paid, exported at an inverse factor against the stability of the outer layers (atmospheric gradient) | Nothing leaves except the wall delta — the answer *is* a data flow read off the outermost wall | Layer boundaries, innermost (superposition atmosphere) to outermost (sacrificial layer) | Hold quality vs sacrifice ratio × layer count (two-axis sweep); layered-jar shake test vs mixed at equal budget |
| Qubit atmosphere, inward | Light bathing tightens coherence — power in, spectrum-specific, through shadows | Inputs must arrive *as* decoherence inside the stabilized band (Q3): data enters encoded in the wall's own language | The atmosphere's outer surface | Bathing extends hold vs dark matched control, and only in-band (spectrum-specificity check) |
| Processing tier (George / core app — Muse's mapping, labeled) | Compute and thermal budget: each brain's residency is a power allocation; the pinned brain is never evicted | Packets sized to the receiving brain cross between organs; the conduit (core) routes, it does not think | The handoff / dwell boundary between organs | Timed tier receipts: sub-minute simple answers, sub-4-minute tool tasks, check-ins beyond 10 minutes (Zak's human-use tiers) |

---

## 3. The language, applied — against data, physical phenomenon, assumed phenomenon

Each term of Zak's language does a specific job. This table keeps the
three relations separate so the language can never borrow credit from
one column to spend in another.

| Term (application) | In relation to the DATA (catalog) | In relation to the PHYSICAL phenomenon (observed) | In relation to the ASSUMED phenomenon (presumed substrate) |
|---|---|---|---|
| **Bifurcation** — every generative split makes two receipted branches that rejoin | ONI flips, D–O steps, ERA5 regime onsets: flips with dwell times and step shapes on record | Observed: bistable and oscillating systems flip — ENSO, glacial climate, flow regimes at pitchfork/Hopf onset | Assumed: the same split-and-rejoin is the universe's generative operation at every tier (nested daughters from twisted parents). Graded by: flip statistics beating shuffled controls in ≥2 tiers under one protocol |
| **Superfluidic** — a component flows against another without friction, so the pair can store and release against each other | Pulsar glitch bookkeeping closes only with a decoupled ~1% reservoir; He-II literature values (digitization pending) | Observed: neutron-star interiors (inferred from glitches), He-II two-fluid behavior in the lab | Assumed: a superfluid substrate underlies the tiers, which is why decoupled-reservoir bookkeeping keeps appearing. Graded by: bookkeeping failing *without* the reservoir, tier by tier — one clean failure wounds the assumption |
| **Solidic** — a pattern holds its shape across time / a crossing instead of blurring | Planck Cℓ: a pattern frozen at recombination, readable 13.8 Gyr later; pulsar crust as the solid partner in the glitch pair | Observed: frozen-in patterns persist — CMB acoustic peaks, crusts, flux-freezing in MHD plasmas | Assumed: solidity is what a wave looks like when it holds — pattern persistence is a wave state, not a separate substance. Graded by: pattern-survival statistics vs random-phase / Gaussian controls |
| **Second sound** — temperature (state) travels as a wave, not a blur, so pattern can cross a medium coherently | No machine-readable dataset (honest null, §1); the catalog's receipts for it are currently literature tables | Observed: second sound is measured physics in He-II (and second-sound-like modes in other superfluids) | Assumed: second sound is the general state channel of the substrate — the slow, meaning-carrying wave under the fast raw flow, at every tier. Graded by: the digitization task first; then any tier where state demonstrably crosses *without* a wave channel wounds it |
| **Particlization** — a wave of the medium arrives as discrete transactions | GW events, glitch kicks, beacon pulses, photon-counting: waiting-time series in every domain | Observed: detectors register waves as events; the event/wave pair is instrument routine | Assumed: particlization is one operation, tier-invariant — discreteness is how waves *pay*, everywhere. Graded by: waiting-time statistics a scrambled control cannot reproduce, in ≥2 tiers |
| **Transduction** — conversion across a wall, completed over the wall's travel; the wall traveling alone | Wall-delta protocols (Q4) and every detector row: the catalog only ever contains wall readouts | Observed: every measurement in the catalog is a transduction — no instrument reads an interior directly | Assumed: this is a law, not an instrument limitation — interiors are knowable *only* inductively, in physics generally (§20). Graded by: any demonstrated non-wall readout channel of a cavity interior |
| **Duality** — one system, two complete expressions, switchable or instrument-dependent | GW strain vs event catalog; Sycamore bitstrings vs circuit description; ONI phase vs storm events | Observed: wave/particle, warm/cold band, scatter/settle pairs recur in the records themselves | Assumed: duality is a *skill* systems have — an expression they can be switched between (photonic macrophage: patrol in wave mode, act in particle mode). Graded by: the macrophage unit test — one encoded trigger, one measured scatter↔absorb flip |
| **Tier** — a scale-band with its own fluid, its own spectrum, its own floor | The catalog is tier-indexed: climate / cosmos / quantum / math / processing rows never mix custodians | Observed: each band genuinely has distinct dynamics, instruments, and timescales | Assumed: tiers are *floors of one structure* — same operations (flip, hold, flow, pay) recurring per floor, related by ratios, not by direct mechanism. Graded by: the cross-tier spacing-law comparison (§2, accounting row) |
| **Atmosphere** (qubit lens) — layered stability that houses randomness instead of exiling it | Sycamore error budgets; NIST/CURBy certified randomness as the "well-kept atmosphere emits this" reference point | Observed: Earth's atmosphere demonstrably stabilizes a surface where photonic transduction expresses in every form (weather, chemistry, life); layered shielding is standard qubit engineering practice in weaker form | Assumed (NEW in the 2026-11-16 lens): a processing volume should be built as stratigraphy — superposition atmosphere at the core, sacrificial layers outward — and true randomness is a *resident* the layers protect, calculable and expressible uniformly inside. Graded by: layered-jar vs mixed-jar at equal sacrifice budget |
| **Receipt** — no claim without the artifact it rests on | This catalog itself: every row names custodian, access route, control, falsifier | Observed: the A–C ladder lesson — an effect survived while its mechanism attribution failed its scrambled control; both halves are receipts | Assumed: none. The receipt discipline is the one part of the language that assumes nothing — it is how the other nine terms earn or lose their grades |

---

## Provenance

- Framing and language: Zak Wiater, 2026-11-16 (Wave Physics working
  sessions), recorded verbatim in §1 of the source draft.
- Doctrine cited: transduction §20 and qubit Q1–Q4 / §21 / Jar census
  (2026-11-07/08, Wave-resonance paper and census); A–C ladder gate
  discipline (standing brief, Wave Physics).
- Dataset verification: public access and custodians checked
  2026-11-16; the He-II null is stated, not smoothed over.
- Status: catalog and tables are UNGRADED program documents. Gates
  G-T0…G-T5 in the source draft decide what, if anything, graduates.
