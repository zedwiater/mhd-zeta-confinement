# Riemann encasement — the logic note (2026-11-04)

Status: framework statement for the cavity program, written to be attacked.
The Riemann connection in this repository is a **conjecture, stated
formally in Lean 4** — the bridge "all nontrivial zeros are confined
stable configurations" *is* the Riemann hypothesis, and nothing here
claims it as proved. What this note does is pin the physical reading the
cavity iterations have been converging on, with the observations that
would break it.

## The claim

The critical line has to be **encased** by the 0/1 walls of the critical
strip for a daughter spectrum to exist at all. Daughters develop *inside*
the strip. They return to the parent two ways: by bleeding back once a
configuration resonates enough power and data, or by rectification
through a transducer.

## Why the encasement is load-bearing

No walls, no standing waves, no discrete spectrum. A mode set needs a
cavity; the strip is the cavity.

The strip is also exactly the region where both of the parent's ledgers
fail:

- For Re(s) > 1, the Euler product rules — the multiplicative account
  book, fully deterministic.
- For Re(s) ≤ 0, the trivial zeros rule — closed form, also
  deterministic.
- ζ's only pole sits at s = 1, *on* the wall, not inside the room.

So the strip is the gap between two bookkeeping systems: the only region
whose contents are not read off a ledger but have to be earned
dynamically. Inside the encasement, the line Re(s) = 1/2 is the seat at
maximum distance from both walls — the one position that is neither
overfed by the injection side nor starved by the drainage side. In the
ladder program's language it is the marginality axis, λ = 0: modes on it
neither grow nor decay; they persist, and persistence is the return
observable.

## Daughters: span inherited, zeros earned

A daughter tier receives the full spectrum — tiers slide the
distribution's interaction peak without narrowing its width. But the
daughter does not inherit the parent's on-line set. Her zeros are
re-earned by her own transduction, at her own scale. This asymmetry is
the tier filter: spectrum is universal, stability is local.

The zero set's two symmetries read as skin physics: conjugation
(ρ ↔ ρ̄, the helicity pairing) and the functional equation's reflection
(s ↔ 1 − s) about the line itself — parent side and daughter side as
mirror images across the stability axis.

## The two return routes

**Bleed-back** is the slow channel: power and data seeping home once a
configuration resonates enough — in instrument terms, once it meets the
persistence criteria and its throughput clears the gate capacity. The
daughter pays from her own radiative budget; the radiative output *is*
the power return and the payment for the data return in one flow.

**Rectification** is the directed channel, and the rectifier also
**recirculates** — at starburst and quasar scale (a quasar being a
parent-side load adjustment fired by transducer overload) and at the
smallest scale alike, where a photon's encode/decode stripping pass
partly recirculates. One-way-ness lives in the transformation, not in a
wall: what returns is not what left. A photon that passes a transducer
and escapes with **new speed, phase, and twist** carries the receipt
that exchange happened.

Exchange is not confined to the razor edge of the line. It couples
**inductively in a band** marginally outside it — consistent with the
ladder data, where condensates live in a band around λ = 0 (pooled
median |λ| = 0.00299, most inside ε = 0.005, tails to ~0.04) rather than
on a measure-zero seat where nothing could ever couple.

## Falsifiers, named in advance

1. **The receipt needs its control.** Changed speed/phase/twist across a
   transducer pass must beat the changed-triple rate across empty
   propagation of the same length. Dispersion alone changes phase with
   no exchange; if the rates match, the receipt is void.
2. **The band width is unmeasured.** ε = 0.005 in the ladder gates was a
   choice, not a derived coupling width. Booking return events against
   parent |λ| would trace the falloff; until then "marginally" is a
   word, not a number.
3. **The daughter test.** In the stack build now running (Wave Ladder
   iteration 13, Waves C/D), if the daughter tiers' persistent
   configurations scatter off marginality after the full cascade, the
   encasement reading fails in-sim as well as on paper.
4. **The open coordinate.** Which wall the bleed crosses — and whether
   bleed and rectification are one wall-crossing read at two scales —
   is not yet pinned. It is recorded here as an open slot, not smoothed
   over.

## Relation to the cavity iterations

Iterations 1–13 of this directory are the encasement experiments read
one element at a time: walls that generate the waves (bubble geometry),
the skin as the net exporter (helicity is a throughput, not an
inventory), pair nulling as the interior rectification primitive
(counter-phase pairs collapsing net work ~185×), nested rectification
with anti-phased residuals. This note is the frame those results are
interpreted in — offered so the next person can break the frame, not
just the numbers.

*Measurements before metaphors — and proofs before either.*
