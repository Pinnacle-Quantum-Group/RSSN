# RSSN T3 — Derivation of the Density→Ordinal Map δ (closing L3.4)

**Status of this note:** Derivation. Resolves `LEMMA_DERIVATIONS.md` item
**L3.4** ("δ Derivation — density-to-ordinal mapping not derived from axioms —
**OPEN**") by deriving δ from the Abel-function structure already formalized in
`formal_proofs/RSSN/HierarchyPlacement.lean` (L3.1, L3.2) and the shape-operator
recurrence in `formal_proofs/RSSN/ShapeOperators.lean`.

**Scope note (per author decision, "document both, decide later"):** This note
does **not** modify any existing claim in `README.md`, `LEMMA_DERIVATIONS.md`,
`TOPOLOGY_MAP.md`, `PUBLICATION_ROADMAP.md`, or the Lean sources. It records the
derivation and flags one convention choice (the integer base ordinal α) as an
open decision. The derived **formula for δ is convention-independent**; only the
integer label α shifts by a global +1 between the two conventions (§5).

---

## 0. The problem

`README.md` §5.1 (Theorem 3) *posits*

> S(n) ≈ f_{α + δ(D_k(n))}(n),  with Triangle: α = 3, δ(D_k(n)) = D_k(n) = 1/n.

The functional form δ(D) = D is asserted, with the informal justification that
"density modulates recursion depth and adds a fractional component to the ordinal
index." `LEMMA_DERIVATIONS.md` correctly flags this as **OPEN**: the map from a
density value to a fractional ordinal offset is never derived — it is named.

This note derives it. The result is that **the fractal density is the
leading-order fractional ordinal offset**, forced (not chosen) by the Abel
function of the fast-growing hierarchy.

---

## 1. Ingredients already in the repository

All three ingredients below are already formalized and (for the integer parts)
proven; nothing here is new machinery, only its assembly.

**(I) The fast-growing hierarchy (standard Löb–Wainer), `HierarchyPlacement.lean`:**

```
f_0(n) = n + 1
f_1(n) = 2n
f_2(n) = 2^n · n
f_{k+1}(n) = f_k^[n](n)           -- n-fold iteration  (the successor rule)
```

**(II) The Abel function, `HierarchyPlacement.lean` (L3.1):** for the analytic
representative of f_α, an `AbelFunction` provides a strictly monotone `A_α` with
`A_α(f_α(x)) = A_α(x) + 1`. The proven lemma `L3_1_abel_iteration` states

> **(1)**   A_α( f_α^[m](x) ) = A_α(x) + m.

**(III) The shape operators, `ShapeOperators.lean`:** the RSSN shapes are
*defined by the very same successor rule* as the hierarchy:

```
Triangle(n) = n^n
Square(n)   = Triangle^[n](n)      -- def `square`  = `triangleIter n n`
Circle(n)   = Square^[n](n)        -- def `circle`  = `squareIter n n`
```

This is the crux that makes the derivation clean: **a shape and its successor are
related by exactly the hierarchy successor rule `S_{m+1}(n) = S_m^[n](n)`.**

---

## 2. The integer part of the ordinal is exact, not posited

Because `Square = Triangle^[n]`, `Circle = Square^[n]`, …, the shapes form their
own fast-growing tower with `Triangle` as base generator. Applying (1) with the
shape's own Abel function `A_S`:

> A_S( S^[n](n) ) = A_S(n) + n.

So each shape sits **exactly one** fast-growing level above the previous one —
this is a structural identity (already proven in spirit by
`square_ge_triangle`, `circle_ge_square`, and the iterate-monotonicity lemmas),
not an approximation and not an assumption:

```
if  Triangle ≈ f_β   then   Square ≈ f_{β+1},   Circle ≈ f_{β+2},   …
```

The only thing left to determine is **(a)** the integer base β of `Triangle`, and
**(b)** the *fractional* offset δ of each shape within its level. Both come from
the Abel function.

---

## 3. The fractional offset δ is forced by the Abel function

Define fractional hierarchy levels by **Abel interpolation**, anchored so that
offset 0 is the argument n and offset 1 is the successor f_{α+1}(n):

> **(2)**   f_{α+θ}(n) := A_α^{-1}( A_α(n) + θ·n ),   θ ∈ [0, 1].

Endpoints check out:
- θ = 1 gives `A_α^{-1}(A_α(n) + n) = A_α^{-1}(A_α(f_α^[n](n))) = f_{α+1}(n)` by (1);
- one Abel step, `A_α(n)+1`, lands on `f_α(n)`, i.e. offset 1/n.

Now write a function S as `S(n) = f_{α+δ}(n)`. Equation (2) is **invertible**
because `A_α` is strictly monotone, so δ is *determined*, not free:

> **(★)**   δ_S(n) = ( A_α(S(n)) − A_α(n) ) / n.

This is the derivation of δ: it is the Abel-normalized "how far S(n) has climbed
from n toward f_{α+1}(n)," in units where one full successor step is 1. No
appeal to density yet — δ is read off the hierarchy alone.

The formula (★) is what was missing in `README.md` §5.1. Two immediate sanity
checks:
- **Successor:** δ for `f_α` itself is `(A_α(f_α(n)) − A_α(n))/n = 1/n` by the
  Abel equation. The bottom rung of a level carries offset 1/n.
- **Top of level:** δ for `f_{α+1}` is `(A_α(f_α^[n](n)) − A_α(n))/n = n/n = 1`
  by (1). Reaching the next integer level is offset 1, as it must be.

---

## 4. Evaluating δ for Triangle: the density falls out

Take `Triangle(n) = n^n`. Under the standard convention its base is α = 2
(justified in §5). The relevant Abel function is that of the **repository's** rung
`f_2(x) = 2^x · x` (`formal_proofs/RSSN/HierarchyPlacement.lean:19`), so it obeys
the *exact* equation

> **(A₂)**   A_2(2^x · x) = A_2(x) + 1   — **not** the pure-exponential
> `A_2(2^x) = A_2(x) + 1`. Writing `f_2(x) = 2^{x + log₂ x}`, the `·x` adds only a
> *sub-linear* `log₂ x` to the exponent.

**Bracket (rigorous, convention-exact).** For `n ≥ 3`,
`f_2(n) = 2^n·n ≤ n^n < f_2(f_2(n))` — the left because `n^n/(2^n·n) = n^{n-1}/2^n
→ ∞`, the right because `f_2²(n)` is doubly exponential. With
`A_2(f_2^[k](n)) = A_2(n) + k` (the proven `L3_1_abel_iteration`) and `A_2`
strictly monotone,

```
1 < A_2(n^n) − A_2(n) < 2,    so    1/n < δ_Triangle(n) < 2/n   (all n ≥ 3).
```

**Sharpening to the lower edge.** The Abel *difference* `A_2(u) − A_2(v)` is
invariant under applying `f_2^{-1}` to both arguments (immediate from (A₂)). Since
`f_2^{-1}(y) = log₂ y − log₂log₂ y + o(1)`, one application sends
`n^n ↦ n·log₂ n + o(·)` and `n ↦ log₂ n + o(·)`; iterating `f_2^{-1}`, the two
orbits collapse — their separation is only of order `log₂log₂ n`, crushed by each
further de-exponentiation. So the fractional part above the first iterate vanishes:

> **(†)**   A_2(n^n) − A_2(n) → 1,    so   **δ_Triangle(n) = 1/n + o(1/n)**,

now derived from the repository's `f_2(x) = 2^x·x` itself — the `·x` factor enters
only at the `log₂log₂` level and does not touch the leading order. *(The first
draft used the pure-exponential Abel equation `A_2(2^x)=A_2(x)+1`; thanks to the
Codex review for catching that. The corrected derivation above changes the proof,
not the result: the unconditional bracket already pins `δ ∈ (1/n, 2/n)`, and the
`·x` correction is sub-`log₂log₂`.)*

And the proven density of Triangle is `D(n) = 1/n` (`triangle_ratio_constant`,
`triangle_density_converges`). Therefore

> **Result (L3.4).**   δ_Triangle(n) = D(n) + o(D(n)).

The fractal density **is** the leading-order fractional ordinal offset. The
posited `δ(D) = D` of README §5.1 is *recovered as the leading term* of the
Abel-derived δ — now with a reason: density measures the same Abel-normalized
climb that the fractional ordinal does.

**The other shapes.** By §2, `Square ≈ f_{β+1}`, `Circle ≈ f_{β+2}`. Their
densities converge to 0 (`square_density_vanishes_model`), consistent with their
δ → 0 within their own level: each successor's bottom rung is offset 1/n → 0, and
the density tracks it. The fractional shapes (e.g. `Pentagon(n) = Circle^{D₃(n)}(n)`
in the FractalMega example, README §7) are *defined* as Abel-fractional iterates
by the density — exactly the δ-as-Abel-offset picture, read in the other
direction.

---

## 5. The one open convention choice: the integer base α

The **derived formula (★) is convention-independent**. What is *not* fixed by the
derivation is the integer label α of the base rung, and here the repository is
currently inconsistent — which is the real reason L3.4 reads "OPEN":

| Convention | f at the exponential rung | Triangle = n^n placement | Triangle base α |
|---|---|---|---|
| **Standard (Löb–Wainer)** — used by `HierarchyPlacement.lean` | f_2(n) = 2^n·n | f_2 ≲ Triangle < f_3 | **α = 2** |
| **Shifted** — used by README §5.1 prose ("f_3(n) = 2^n") | f_3(n) = 2^n | f_3 ≲ Triangle < f_4 | **α = 3** |

The committed Lean already chose the **standard** convention and **proves**
`Triangle(n) < f_3(n)` (`L3_3_triangle_below_f3`), which is *incompatible* with the
README's α = 3 (that would put Triangle *above* f_3). So:

- Under the **standard** convention (recommended for consistency with the
  proofs): Triangle α = 2, Square α = 3, Circle α = 4, …
- Under the **shifted** convention: Triangle α = 3, Square α = 4, … — but this
  requires re-defining `f₃` in `HierarchyPlacement.lean` and re-proving L3.3.

Either way, δ is given by (★) and `δ_Triangle = D + o(D)`. **Only the integer α
shifts by a global +1; the density→offset content is unchanged.** Author to
decide which convention becomes canonical; this note deliberately changes neither.

---

## 6. Formalization sketch (for `HierarchyPlacement.lean`, not yet committed)

The derivable core slots directly onto the existing `AbelFunction` structure and
`L3_1_abel_iteration`. Sketch (to be type-checked against Lean 4.5.0 / Mathlib in
CI before committing — no toolchain is available in the authoring environment):

```lean
/-- Derived density→ordinal offset (★): δ_S(n) = (A(S n) − A n) / n. -/
noncomputable def deltaOffset (af : AbelFunction) (S : ℝ → ℝ) (n : ℝ) : ℝ :=
  (af.A (S n) - af.A n) / n

/-- Bottom rung of a level carries offset 1/n. -/
theorem delta_base (af : AbelFunction) (n : ℝ) (hn : n ≠ 0) :
    deltaOffset af af.T n = 1 / n := by
  unfold deltaOffset; rw [af.abel_eq]; field_simp

/-- Reaching the successor f_{α+1}(n) = f_α^[n](n) is offset 1 (uses L3.1). -/
theorem delta_succ (af : AbelFunction) (n : ℕ) (hn : (n:ℝ) ≠ 0) :
    deltaOffset af (fun x => af.T^[n] x) (n : ℝ) = 1 := by
  unfold deltaOffset; rw [L3_1_abel_iteration af n (n:ℝ)]; field_simp
```

The analytic step (†) — `A_2(n^n) − A_2(n) → 1` — is the one genuinely analytic
obligation; it is rigorous on paper (super-log estimate above) and can be
formalized via the bracketing `1 < A_2(n^n) − A_2(n) < 2` (from `f_2 < n^n <
f_2∘f_2`, already provable with the `g₂`-iterate lemmas in the file) plus a
super-log continuity bound for the `o(1)` term.

---

## 7. Summary

- δ is **derived**, not posited: **δ_S(n) = (A_α(S(n)) − A_α(n)) / n**, forced by
  Abel interpolation (2) of the fast-growing levels.
- The integer part of each shape's ordinal is **exact** (`S_{m+1} = S_m^[n]`).
- For Triangle, **δ = D + o(D)** — the fractal density is the leading-order
  fractional ordinal offset, recovering README §5.1's `δ(D) = D` with a reason.
- One convention choice remains open (base α: 2 vs 3); it shifts only the integer
  label and is left to the author. No existing claim is modified by this note.

This reduces L3.4 from **OPEN** to **DERIVED (modulo the documented convention
choice)**, unblocking Paper 3 and Paper 6 in `PUBLICATION_ROADMAP.md` once the
convention is fixed and the §6 lemmas are type-checked in CI.
