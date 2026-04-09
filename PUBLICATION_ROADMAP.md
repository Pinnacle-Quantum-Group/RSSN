# PQG Framework: Publication Roadmap

Ordered by: closest to tight, unlocks most downstream, most novel, TDD-confirmed.

---

## Paper 1: RSF Foundation Results (Ready Now)

**Title:** "Recursive Structure Foundation: Russell's Transcendence and Recursive Foundation"

**Content:** RSF T3 (Russell's Transcendence) + RSF T5 (Recursive Foundation)

**Why first:**
- Both TIGHT from Axiom 7 alone -- no dependencies on other frameworks
- Anchor theorems: everything else builds on these
- Novel contribution: generation/membership distinction eliminates Russell paradox
- Clean, self-contained, shortest path to publication

**Status:** All TDD tests pass. No known gaps.

**Venue suggestions:** Journal of Symbolic Logic, Notre Dame J. Formal Logic

---

## Paper 2: Cardinality Transcendence (Ready Now)

**Title:** "Replacing Cardinality with Recursive Density: Transcendence and Continuum Resolution"

**Content:** RSF T1 + RSF T2 + RSSN T5

**Why second:**
- TIGHT (L2.1 density non-collapse + L2.2 no hidden hierarchy)
- High novelty: replaces entire cardinality framework
- Dissolves Continuum Hypothesis (doesn't decide it -- makes it meaningless)
- Builds on Paper 1's foundation results

**Dependencies:** Paper 1 (RSF axioms)

**Status:** All TDD tests pass. Critical gap (density secretly = cardinality) CLOSED.

**Venue suggestions:** Annals of Pure and Applied Logic, J. Mathematical Logic

---

## Paper 3: Core RSSN Results (Needs RLA Published)

**Title:** "Recursive Shape-Structured Notation: Convergence, Reflection, and Non-Commutativity"

**Content:** RSSN T1, T2, T7 + Bridge B.1

**Why third:**
- T7 is TIGHT (anchor, explicit computation)
- T1 is TIGHT (via L1.1 + L1.2 + L1.3 using RSF + Bridge)
- T2 is TIGHT (L2.1-L2.4)
- Establishes the notation system's mathematical validity

**Dependencies:** Paper 1 (RSF axioms), RLA paper (Bridge B.1)

**Key correction:** T3 hierarchy placement must be corrected (Triangle < f_3, not = f_3)

**Status:** Core results pass. T3 placement corrected. delta derivation remains OPEN.

**Venue suggestions:** Journal of Symbolic Computation, Advances in Mathematics

---

## Paper 4: Singularity Resolution + Bekenstein (Strong Physics)

**Title:** "Fractal Tensor Calculus: Singularity Resolution and the Bekenstein Bound"

**Content:** FTC T4 + FTC T6 + IT connection

**Why fourth:**
- TIGHT for black holes (L5.1-L5.6 complete chain)
- Strongest physical content: R^(n*) -> 0 at BH singularity
- D*_BH = e^{-pi} universal (mass-independent)
- Entropy = 2*pi nats = Bekenstein capacity exactly
- High impact if correct

**Dependencies:** Papers 1-3 (RSF + RSSN + FTC axioms), IT paper

**Status:** All BH tests pass. Naked singularities OPEN.

**Venue suggestions:** Physical Review D, Classical and Quantum Gravity

---

## Paper 5: RLA Geometric Backbone (Independent Track)

**Title:** "Recursive Lie Algebras of Vector Fields" (already written)

**Content:** RLA T2, T5, T7, P8, L3

**Why this track:**
- ALL results TIGHT -- most rigorous document in the framework
- Self-contained: does not require RSSN/RSF/FTC
- Publishable independently as differential geometry contribution
- Enables Bridge B.1 which connects everything else

**Status:** Ready. All proofs complete with appendices.

**Venue suggestions:** Journal of Geometry and Physics, Communications in Mathematical Physics

---

## Paper 6: Corrected Hierarchy + Uncertainty (Future Work)

**Title:** "RSSN in the Fast-Growing Hierarchy: Corrected Placements and Uncertainty Relations"

**Content:** RSSN T3 (corrected) + T8 (Robertson form)

**Why last:**
- T3 has OPEN gap (delta derivation)
- T8 original form fails (n=1 counterexample) -- needs reformulation
- Robertson form is tight but less novel
- Ecalle machinery (L3.1) is established but delta mapping is not

**Status:** Corrections identified. New derivation needed for delta.

---

## Recommended Order

```
Phase 1 (parallel, immediate):
  Paper 1: RSF Foundation        [TIGHT, no dependencies]
  Paper 5: RLA Geometric          [TIGHT, independent]

Phase 2 (after Phase 1):
  Paper 2: Cardinality             [TIGHT, needs Paper 1]
  Paper 3: Core RSSN               [TIGHT, needs Papers 1+5]

Phase 3 (after Phase 2):
  Paper 4: Singularity + Bekenstein [TIGHT for BH, needs 1-3]

Phase 4 (ongoing):
  Paper 6: Hierarchy + Uncertainty  [OPEN gaps, future work]
```

---

## Cross-Framework Coherence Check

Before submitting any paper, verify the core principle instantiation:

| Framework | Value at saturation | Test |
|-----------|--------------------|---------|
| RSF | D in (0,1] finite | `fractal_density_triangle(n) == 1/n` |
| RSSN | D_k(n) converges | `fractal_density_sequence` monotone bounded |
| RLA | D smooth, positive | `rla_alpha_discrete` finite for Triangle |
| FTC | D* = e^{-pi} | `E_NEG_PI ~ 0.0432` |
| IT | I_local = 2*pi | `2 * math.pi ~ 6.283` |

All values consistent. No contradictions found across frameworks.
