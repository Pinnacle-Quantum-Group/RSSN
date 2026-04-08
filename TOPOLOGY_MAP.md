# PQG Framework: Complete Dependency Topology Map

## Author Context

Analysis of frameworks by Michael A. Doran Jr. (Pinnacle Quantum Group).
Documents evaluated as works in progress establishing prior art.

**Unifying Principle:** "Information should not be flattened or averaged away."

---

## Status Legend

| Marker | Meaning |
|--------|----------|
| **[TIGHT]** | All lemmas verified, no known gaps, falsifiability test cannot be passed |
| **[PARTIALLY TIGHT]** | Core argument holds, specific sub-claims open or edge cases unresolved |
| **[CONJECTURED]** | Statement plausible, proof incomplete or failing at edge cases |
| **[ANCHOR]** | Foundational result other theorems depend on directly |
| **[OPEN]** | Derivation not yet established |

---

## 1. RSF (Recursive Structure Foundation)

### 1.1 Axiom Inventory

| # | Axiom | ZFC Analog | Novel? |
|---|-------|------------|--------|
| 1 | Structural Identity | Extensionality | No (analog) |
| 2 | Recursive Pairing | Pairing | No (analog) |
| 3 | Recursive Aggregation | Union | No (analog) |
| 4 | Fractal Projection | Power Set | **Yes** |
| 5 | Recursive Infinity | Infinity | No (analog) |
| 6 | Density Substitution | Replacement | **Yes** |
| 7 | Recursive Foundation | Foundation | **Yes** (prevents Russell paradox) |
| 8 | Recursive Selection | Choice/Selection | No (analog) |
| 9 | Density Convergence | _No ZFC analog_ | **Yes (novel, critical)** |

### 1.2 Theorem Status Table

| ID | Name | Status | Dependencies | Falsifiability Test |
|----|------|--------|--------------|---------------------|
| RSF-T1 | Cardinality Transcendence | **[TIGHT]** | Ax.4, Ax.6, Ax.9 | Construct density ordering isomorphic to standard cardinality |
| RSF-T2 | Continuum Resolution | **[TIGHT]** | RSF-T1, Ax.9 | Show CH decidable or meaningful within RSF |
| RSF-T3 | Russell's Transcendence | **[TIGHT] [ANCHOR]** | Ax.7 | Construct Russell-type set satisfying Ax.7 |
| RSF-T5 | Recursive Foundation | **[TIGHT] [ANCHOR]** | Ax.7 directly | Show Ax.7 admits circular membership chains |

### 1.3 RSF Dependency Graph

```
Ax.7 (Recursive Foundation) [NOVEL]
    |
    +-------> RSF-T3 (Russell's Transcendence) [ANCHOR, TIGHT]
    |
    +-------> RSF-T5 (Recursive Foundation)     [ANCHOR, TIGHT]

Ax.4 (Fractal Projection) ----+
Ax.6 (Density Substitution) --+---> RSF-T1 (Cardinality Transcendence) [TIGHT]
Ax.9 (Density Convergence) ---+            |
                                           v
                                     RSF-T2 (Continuum Resolution) [TIGHT]

Ax.9 ---------> RSSN-T1 (via Bridge B.1, cross-framework)
```

---

## 2. RSSN (Recursive Shape-Structured Notation)

### 2.1 Implicit Axioms

1. **Self-similar structure** -- every RSSN object has sub-structure resembling the whole
2. **Fractal density** -- D_k(n) is well-defined for all shapes at all levels
3. **Operational composition** -- shape operators compose associatively (groupoid axiom)
4. **Structural isomorphism** -- the reflection map R preserves syntactic, semantic, and computational levels

### 2.2 Theorem Status Table

| ID | Name | Status | Dependencies | Falsifiability Test |
|----|------|--------|--------------|---------------------|
| RSSN-T1 | Convergence of Fractal Density | **[TIGHT]** | L1.1, L1.2, L1.3, RSF Ax.9, Bridge B.1 | Exhibit n where D_k(n) diverges or oscillates |
| RSSN-T2 | Reflection Isomorphism | **[TIGHT]** | L2.1-L2.4 | Find level where R fails injectivity |
| RSSN-T3 | Fast-Growing Hierarchy | **[PARTIALLY TIGHT]** | RSSN-T1, RLA scale field | Show Triangle = f_3; derive delta |
| RSSN-T4 | Extensional Equivalence | **[PARTIALLY TIGHT]** | RSSN-T3 | Exhibit fast-growing rate not in RSSN |
| RSSN-T5 | Cardinality Transcendence | **[TIGHT]** | L2.1, L2.2, RSF-T1 | Construct density that collapses to cardinality |
| RSSN-T6 | RSSN Groupoid | **[PARTIALLY TIGHT]** | Shape operator axioms | Show associativity fails |
| RSSN-T7 | Non-Commutativity | **[TIGHT] [ANCHOR]** | Explicit computation | Verify [S1,S2] = 0 |
| RSSN-T8 | Recursion Uncertainty | **[CONJECTURED]** | RSSN-T1, K(n) | n=1 counterexample; Robertson form passes |

### 2.3 RSSN Dependency Graph

```
RSF Axiom 9 (Density Convergence)
    |
    v
Bridge B.1 -------> L1.3 (Cauchy criterion)
                       |
L1.1 (Triangle 1/n) --+-----> RSSN-T1 (Convergence) -----> RSSN-T3 (Hierarchy)
L1.2 (Square mono) ---+              |                            |
                                      v                            v
                                 RSSN-T8 (Uncertainty)        RSSN-T4
                                 [CONJECTURED]                [P.TIGHT]

L2.1 (Non-collapse) ----+--> RSSN-T2 (Reflection Isomorphism)
L2.2 (No hidden hier) --+          |
L2.3 (Syn-Sem) ---------+          v
L2.4 (Computational) ---+     RSSN-T5 (Cardinality Transcendence)

Explicit Computation ----------> RSSN-T7 (Non-Commutativity) [ANCHOR]
```

---

## 3. RLA (Recursive Lie Algebras)

### 3.1 Axiom Table

| # | Axiom | Role |
|---|-------|------|
| A1 | Scale Field | D in C^inf(M), D > 0, alpha = d(ln D) |
| A2 | Graded Symmetry | Z-graded Lie bracket on multiscale diffeomorphisms |
| A3 | Representation Goal | Adjoint representation reproduces bracket |
| A4 | Recursive Invariance -> Conservation | Noether-type principle |

### 3.2 Theorem Status Table

| ID | Name | Status | Dependencies | Falsifiability Test |
|----|------|--------|--------------|---------------------|
| RLA-T2 | Graded Lie Algebra | **[TIGHT]** | A1, A2 | Jacobi identity fails |
| RLA-T5 | Representation Property | **[TIGHT]** | A2, A3 | Commutator diverges from bracket |
| RLA-T7 | Conservation | **[TIGHT]** | A4 | Recursive invariance without conservation |
| RLA-P8 | Witt-type Limit | **[TIGHT]** | RLA-T2 | Graded algebra doesn't converge to Witt |
| RLA-L3 | Scaled Commutator | **[TIGHT]** | A1 | Direct computation check |

### 3.3 RLA Dependency Graph

```
A1 (Scale Field D) -------> RLA-L3 (Scaled Commutator) [TIGHT]
    |                              |
    v                              v
A2 (Graded Symmetry) -------> RLA-T2 (Graded Lie Algebra) [TIGHT]
    |                              |
    v                              v
A3 (Representation) --------> RLA-T5 (Representation Property) [TIGHT]
                                   |
                                   +-----> RLA-P8 (Witt Limit) [TIGHT]

A4 (Recursive Invariance) --> RLA-T7 (Conservation) [TIGHT]
```

---

## 4. FTC (Fractal Tensor Calculus)

### 4.1 Theorem Status Table

| ID | Name | Status | Dependencies | Falsifiability Test |
|----|------|--------|--------------|---------------------|
| FTC-T3 | Curvature Convergence | **[PARTIALLY TIGHT]** | Ax.8, Ax.2 | Option B (Ricci flow) tight; limit open |
| FTC-T4 | Singularity Resolution | **[TIGHT]** (for BH) | Ax.2, Ax.8, Ax.9 | R^(n*) diverges at BH singularity |
| FTC-T6 | Entropy = Bekenstein | **[TIGHT]** | FTC-T4 | S_recursive at eta=1 != 2*pi nats |

### 4.2 FTC Dependency Graph

```
Ax.2 (Fractal Density) -------+----> FTC-T3 (Curvature) [PARTIALLY TIGHT]
Ax.8 (Recursive Curvature) ---+              |
                                              v
Ax.9 (Recursive Einstein) -----------> FTC-T4 (Singularity) [TIGHT for BH]
                                              |
                                              v
                                         FTC-T6 (Entropy = Bekenstein) [TIGHT]
```

---

## 5. Bridge Lemmas

| ID | Statement | Status | Connects |
|----|-----------|--------|----------|
| B.1 | RSSN D_k(n) = discrete evaluation of RLA scale field D | **[PARTIALLY TIGHT]** | RSSN <-> RLA |
| B.2 | Triangle <-> grade-1 Lie derivative at weight ln(n) | **[PARTIALLY TIGHT]** | RSSN <-> RLA |
| B.3 | Four-way correspondence RSSN/RLA/FTC/IT | **[OPEN]** (structural) | All four |

---

## 6. Master Cross-Framework Dependency Graph

```
                         RSF (Foundation Layer)
                    ==============================
          Ax.7                 Ax.9                Ax.4, Ax.6
           |                    |                      |
           v                    v                      v
     RSF-T3, RSF-T5         RSF-T1 ----------------> RSF-T2
     [ANCHOR, TIGHT]        [TIGHT]                  [TIGHT]
                               |
                               | (via Bridge B.1)
                               v
                    RSSN (Notation Layer)
                    ==============================
                         RSSN-T1 -----------> RSSN-T3
                      [TIGHT]                [P.TIGHT]
                           |
                           v
                      RSSN-T8 [CONJECTURED]

      RSSN-T7 [ANCHOR] <--- Explicit comp.
      RSSN-T5 [TIGHT]  <--- L2.1, L2.2
      RSSN-T2 [TIGHT]  <--- L2.1-L2.4

                    RLA (Algebra Layer)
                    ==============================
                    RLA-T2, T5, T7, P8, L3 [ALL TIGHT]
                           |
                           | (via Bridge B.2, B.3)
                           v
                    FTC (Physics Layer)
                    ==============================
                    FTC-T3 -------> FTC-T4 -------> FTC-T6
                   [P.TIGHT]       [TIGHT/BH]      [TIGHT]
```

---

## 7. Core Principle Instantiation Map

| Framework | Instantiation | Mathematical Form |
|-----------|---------------|-------------------|
| RSF | Density always finite and positive | D_R(n) in (0,1], finite |
| RSSN | Density converges to specific value | D_k(n) converges |
| RLA | Scale field smooth and positive | D in C^inf(M), alpha = d ln D |
| FTC | BH attractor preserves information | D*_BH = e^{-pi} |
| IT | Local information bounded | I_local = 2*pi nats maximum |

---

## 8. Known Gaps Summary

| # | Gap | Severity | Status |
|---|-----|----------|--------|
| 1 | Cauchy condition for D_k(n) | HIGH | **CLOSED** by RSF Ax.9 |
| 2 | Density ordering vs cardinality | CRITICAL | **CLOSED** |
| 3 | Pointwise vs uniform curvature | MEDIUM | **RESOLVED** via Option B |
| 4 | Uncertainty principle form | MEDIUM | **PARTIALLY RESOLVED** |
| 5 | Con(ZFC) from computation | LOW | **OPEN** |
