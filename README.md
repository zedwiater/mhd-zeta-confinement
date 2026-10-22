# MHD Zeta Confinement — a Lean 4 Formalization

[![Lean Action CI](https://github.com/zedwiater/mhd-zeta-confinement/actions/workflows/lean_action_ci.yml/badge.svg)](https://github.com/zedwiater/mhd-zeta-confinement/actions/workflows/lean_action_ci.yml)
![Lean](https://img.shields.io/badge/Lean-4.36.0--rc1-blue)
![Mathlib](https://img.shields.io/badge/Mathlib-pinned%20d888143-lightgrey)

> *A plasma pinch is a claim about a balance of forces. This repository is where that claim gets proved — or killed — in public.*

## What this is

A [Lean 4](https://lean-lang.org/) workspace for formalizing the physics of **MHD zeta (Z-) pinch confinement**: a column of plasma held together by the magnetic field of its own axial current. Instead of arguing the confinement criteria in prose, the goal is to state them as theorems — with every assumption on the table and a machine checking each step.

This is deliberately adversarial work. A formal proof cannot be charmed, rounded, or hand-waved past a bad step: either the pinch conditions close, or they don't, and the gap shows up as an unprovable goal instead of a plausible paragraph.

## The physics on the table

A Z-pinch lives or dies by a short list of classical results. These are the formalization targets:

| Result | Statement, in words | Status |
|---|---|---|
| **Bennett relation** | A self-pinched column in equilibrium ties current, temperature, and line density into one balance: the magnetic pressure of the pinch carries the plasma pressure. | ⬜ target |
| **Kruskal–Shafranov condition** | The kink-stability limit: wind the field lines too shallowly around the column and the pinch tears itself sideways. | ⬜ target |
| **Taylor relaxation / Woltjer invariants** | A turbulent plasma relaxes toward a minimum-energy state that conserves magnetic helicity — the frozen-in flux does the bookkeeping. | ⬜ target |

The table is the honest part of this README: nothing above is claimed as proved yet. As theorems land, rows flip from ⬜ to ✅ with a link to the file that proves them.

## Repository layout

```
MhdZetaTest.lean        Root module of the library — imports everything that ships
MhdZetaTest/Basic.lean  First definitions module (scaffold)
lakefile.toml           Package definition — Lean library `MhdZetaTest`
lake-manifest.json      Exact dependency pins (Mathlib and friends)
lean-toolchain          Toolchain pin: leanprover/lean4:v4.36.0-rc1
.github/workflows/      Lean Action CI — builds the project on every push
```

## Building it yourself

You need [elan](https://github.com/leanprover/elan) (the Lean toolchain manager). Then:

```bash
git clone https://github.com/zedwiater/mhd-zeta-confinement.git
cd mhd-zeta-confinement
lake exe cache get   # pull the prebuilt Mathlib cache — do not skip this
lake build
```

`lake exe cache get` is the difference between a coffee break and an afternoon: Mathlib is enormous, and the cache gives you its compiled form instead of building it from source.

## Why formalize a pinch?

The pinch is one of the oldest ideas in plasma physics — and one of the most argued about. Equilibrium derivations quietly assume symmetry, steadiness, or boundary behaviour that real columns don't have. Writing the argument in Lean forces each of those assumptions to be named in the theorem statement itself, where it can be attacked directly.

That is the working rule here, the same one used across the wider project: **measurements before metaphors — and proofs before either.**

## Status

🌱 **Scaffold stage.** The project builds against a pinned Mathlib; the library root and module skeleton are in place; the first definitions are next. CI runs on every push, so the badge at the top is the current truth about whether it compiles.

## Related work

- [Wave-resonance](https://github.com/zedwiater/Wave-resonance) — the parent research line: wave-resonance computing, cavity transduction, and energy-accounted inference, where the MHD framing behind this formalization is developed.

---

*Part of an open, standard-setting research program: the code and the proofs are public so anyone can check them, build on them, or prove them wrong.*
