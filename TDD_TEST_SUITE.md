# PQG Framework: TDD Falsifiability Test Suite

One test per theorem, stated precisely enough for independent execution.

**Method:** For each theorem, state what would have to be false for it to fail. If no precise falsifiability condition can be stated, flag for reformulation.

---

## RSSN Tests

| ID | Theorem | Falsifiability Condition | Expected Outcome | Status |
|----|---------|--------------------------|------------------|--------|
| F-RSSN-T1 | Convergence of D_k(n) | Exhibit n,k where D_k(n) oscillates without limit | No such n,k exists for well-defined structures (RSF Ax.9) | **PASS** |
| F-RSSN-T2 | Reflection Isomorphism | Find syntactic expression where R fails injectivity | Cannot: strict growth hierarchy ensures injectivity | **PASS** |
| F-RSSN-T3 | Fast-Growing Hierarchy | Show Triangle(n) >= f_3(n) for large n | Triangle(4)=256 < f_3(4)=65536. CORRECTS original placement. | **PASS (corrected)** |
| F-RSSN-T4 | Extensional Equivalence | Exhibit fast-growing function not expressible in RSSN | Not constructed; claim plausible but unverified for all rates | **OPEN** |
| F-RSSN-T5 | Cardinality Transcendence | Construct density ordering isomorphic to cardinality | Dyadic rationals vs powers-of-2: same |X|, different D(X) | **PASS** |
| F-RSSN-T6 | RSSN Groupoid | Show associativity or identity fails | Not tested for all triples; holds for computed cases | **PARTIAL** |
| F-RSSN-T7 | Non-Commutativity | Verify [Triangle, Square](2) = 0 | [T,S](2) != 0 by explicit computation | **PASS (anchor)** |
| F-RSSN-T8 | Uncertainty Principle | Evaluate Delta_D * Delta_R at n=1 | Delta_D=0 at n=1, so 0 >= K(1)/2 is FALSE. Original form fails. | **FAIL (original), PASS (Robertson)** |

## RSF Tests

| ID | Theorem | Falsifiability Condition | Expected Outcome | Status |
|----|---------|--------------------------|------------------|--------|
| F-RSF-T1 | Cardinality Transcendence | Construct D such that D(X)>D(Y) iff \|X\|>\|Y\| | Impossible: same cardinality, different density demonstrated | **PASS** |
| F-RSF-T2 | Continuum Resolution | Show CH is decidable within RSF | CH becomes trivially YES (density spectrum is dense) | **PASS** |
| F-RSF-T3 | Russell's Transcendence | Construct Russell set satisfying Axiom 7 | Axiom 7 prevents self-membership by construction | **PASS** |
| F-RSF-T5 | Recursive Foundation | Show Axiom 7 admits circular membership | F_1(G_0) intersect G_0 = empty by axiom | **PASS (anchor)** |

## FTC Tests

| ID | Theorem | Falsifiability Condition | Expected Outcome | Status |
|----|---------|--------------------------|------------------|--------|
| F-FTC-T3 | Curvature Convergence | Construct compact K where R^(n) diverges | Option B (fixed scale via Ricci flow) converges | **PASS (Option B)** |
| F-FTC-T4 | Singularity Resolution | Compute R^(n*) at BH singularity; show divergence | R^(n*) -> 0 (super-exponential decay). Singularity replaced. | **PASS (BH)** |
| F-FTC-T6 | Entropy = Bekenstein | Compute S_recursive at eta=1; show != 2*pi | S = e^pi * e^{-pi} * pi + pi = 2*pi nats exactly | **PASS** |

## RLA Tests

| ID | Theorem | Falsifiability Condition | Expected Outcome | Status |
|----|---------|--------------------------|------------------|--------|
| F-RLA-T2 | Graded Lie Algebra | Jacobi identity fails for specific triple | Cancellation via dα=0 ensures Jacobi | **PASS** |
| F-RLA-T5 | Representation | Commutator != graded bracket | Leibniz expansion + weight collection gives exact match | **PASS** |
| F-RLA-T7 | Conservation | Recursive invariance without conservation | Noether machinery + weight terms force conservation | **PASS** |
| F-RLA-P8 | Witt Limit | Discrete grading doesn't limit to continuous | λ in R interpolates; λ in Z recovers discrete | **PASS** |

## Bridge Tests

| ID | Lemma | Falsifiability Condition | Expected Outcome | Status |
|----|-------|--------------------------|------------------|--------|
| F-B1 | RSSN-RLA Correspondence | D_k(n) not recoverable as D on lattice | Triangle: D=1/n is C^inf. Square: D->0 requires regularization. | **PARTIAL** |
| F-B2 | Generator Identification | Triangle action != grade-1 Lie derivative | Structural match; quantitative mismatch (ln(n)^2 vs n*ln(n)) | **PARTIAL** |
| F-B3 | Four-way Correspondence | Framework A predicts value contradicting Framework B | No contradiction found; all give D*=e^{-pi} at saturation | **PASS (structural)** |

## Gap Tests

| ID | Gap | Test | Result |
|----|-----|------|--------|
| G1 | T1 Cauchy | Bounded oscillating sequence | Confirms gap exists; RSF Ax.9 closes it |
| G2 | Cardinality | Density bijection attempt | Density independent of cardinality. CLOSED. |
| G3 | Uniform curvature | x-dependent R on compact set | Rate varies with R; Option B resolves. RESOLVED. |
| G4 | Uncertainty | Scalar commutation test | Scalars commute; Robertson form works. PARTIALLY RESOLVED. |
| G5 | Con(ZFC) | Pathology class mapping | A-C covered; D (large cardinals) OPEN. |

---

## Execution

All tests implemented in Python: `tests/test_*.py`

Run: `cd rssn && python -m pytest tests/ -v`
