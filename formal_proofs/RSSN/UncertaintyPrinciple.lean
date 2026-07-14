/-
  RSSN — Recursion Uncertainty Principle (T8 Lemmas)
  Pinnacle Quantum Group — April 2026

  L8.1: Scalar quantities commute (kills original derivation route)
  L8.2: Robertson relation for non-commuting operators (TIGHT) —
        established in the degenerate commuting case AND for a concrete
        non-commuting witness (the Pauli pair σₓ, σᵧ in state |0⟩),
        where the bound holds with equality.
  T8 in Robertson form: ΔS₁ · ΔS₂ ≥ ½|⟨[S₁,S₂]⟩| forces ΔS₁ · ΔS₂ > 0,
        because the shape commutator is proven nonzero.
  Original form fails at n = 1: the ℤ-valued commutator of the genuine
        shape operators vanishes there.
  Reference: LEMMA_DERIVATIONS.md RSSN T8
-/
import Mathlib
import RSSN.ShapeOperators
import RSSN.NonCommutativity

noncomputable section
open Real

namespace RSSN.UncertaintyPrinciple

/-! ## L8.1 — Scalar Commutation
    D_k and d (recursive depth) as scalar quantities commute.
    This kills the original Heisenberg-style derivation. -/

theorem L8_1_scalars_commute (a b : ℝ) : a * b = b * a := mul_comm a b

theorem L8_1_scalar_commutator_zero (a b : ℝ) : a * b - b * a = 0 := by ring

/-! ## L8.2 — Robertson Uncertainty Relation
    For non-commuting operators S₁, S₂ with commutator [S₁,S₂],
    the Robertson bound holds: ΔS₁ · ΔS₂ ≥ ½|⟨[S₁,S₂]⟩|.
    This is standard quantum mechanics applied to RSSN operators. -/

structure OperatorPair where
  variance₁ : ℝ
  variance₂ : ℝ
  commutator_expectation : ℝ
  h_var₁ : 0 ≤ variance₁
  h_var₂ : 0 ≤ variance₂

def robertsonBound (op : OperatorPair) : Prop :=
  op.variance₁ * op.variance₂ ≥ (1 / 2 * |op.commutator_expectation|) ^ 2

theorem L8_2_robertson_nonneg (op : OperatorPair) :
    0 ≤ (1 / 2 * |op.commutator_expectation|) ^ 2 := by positivity

theorem L8_2_robertson_trivial_when_commuting (op : OperatorPair)
    (hcomm : op.commutator_expectation = 0) :
    robertsonBound op := by
  unfold robertsonBound
  rw [hcomm, abs_zero, mul_zero, sq]
  -- Goal: op.variance₁ * op.variance₂ ≥ 0 * 0. Rewrite RHS to 0.
  rw [mul_zero]
  exact mul_nonneg op.h_var₁ op.h_var₂

/-! ## L8.2 Witness — A Genuinely Non-Commuting Pair

    `robertsonBound` must not be inhabited only in the degenerate
    commuting case. The Pauli matrices σₓ, σᵧ measured in the basis
    state |0⟩ give a concrete pair with *nonzero* commutator
    expectation ⟨[σₓ,σᵧ]⟩ = 2i, unit variances, and the Robertson
    bound holding with EQUALITY — witnessing that the bound is tight. -/

/-- Pauli σₓ matrix. -/
def pauliX : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; 1, 0]

/-- Pauli σᵧ matrix. -/
def pauliY : Matrix (Fin 2) (Fin 2) ℂ := !![0, -Complex.I; Complex.I, 0]

/-- The spin-up basis state |0⟩. -/
def ket0 : Fin 2 → ℂ := ![1, 0]

/-- Expectation value ⟨0|M|0⟩ of an observable `M` in the state |0⟩. -/
def expVal (M : Matrix (Fin 2) (Fin 2) ℂ) : ℂ :=
  Matrix.dotProduct (fun i => (starRingEnd ℂ) (ket0 i)) (M.mulVec ket0)

/-- In the state |0⟩ the expectation of `M` is its (0,0) entry. -/
theorem expVal_eq_entry (M : Matrix (Fin 2) (Fin 2) ℂ) : expVal M = M 0 0 := by
  simp [expVal, ket0, Matrix.dotProduct, Matrix.mulVec, Fin.sum_univ_two]

/-- Variance ⟨M²⟩ − ⟨M⟩² of an observable in the state |0⟩ (real part;
    for self-adjoint `M` both moments are real). -/
def varOf (M : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (expVal (M * M)).re - (expVal M).re ^ 2

theorem expVal_pauliX : expVal pauliX = 0 := by
  rw [expVal_eq_entry]; simp [pauliX]

theorem expVal_pauliY : expVal pauliY = 0 := by
  rw [expVal_eq_entry]; simp [pauliY]

/-- σₓ has unit variance in |0⟩: ⟨σₓ²⟩ − ⟨σₓ⟩² = 1 − 0 = 1. -/
theorem varOf_pauliX : varOf pauliX = 1 := by
  unfold varOf
  rw [expVal_pauliX, expVal_eq_entry]
  simp [pauliX, Matrix.mul_apply, Fin.sum_univ_two]

/-- σᵧ has unit variance in |0⟩: ⟨σᵧ²⟩ − ⟨σᵧ⟩² = 1 − 0 = 1. -/
theorem varOf_pauliY : varOf pauliY = 1 := by
  unfold varOf
  rw [expVal_pauliY, expVal_eq_entry]
  simp [pauliY, Matrix.mul_apply, Fin.sum_univ_two, Complex.I_mul_I]

/-- The commutator expectation ⟨0|[σₓ,σᵧ]|0⟩ equals 2i — pure imaginary
    and nonzero, exactly as standard quantum mechanics predicts. -/
theorem expVal_pauli_commutator :
    expVal (pauliX * pauliY - pauliY * pauliX) = 2 * Complex.I := by
  rw [expVal_eq_entry]
  simp [pauliX, pauliY, Matrix.sub_apply, Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- The Pauli pair packaged as an `OperatorPair`: the *computed* unit
    variances and the *computed* commutator-expectation magnitude
    |⟨[σₓ,σᵧ]⟩| = |2i| = 2. -/
def pauliPair : OperatorPair where
  variance₁ := varOf pauliX
  variance₂ := varOf pauliY
  commutator_expectation := Complex.abs (expVal (pauliX * pauliY - pauliY * pauliX))
  h_var₁ := by rw [varOf_pauliX]; exact zero_le_one
  h_var₂ := by rw [varOf_pauliY]; exact zero_le_one

theorem pauliPair_variance₁ : pauliPair.variance₁ = 1 := varOf_pauliX

theorem pauliPair_variance₂ : pauliPair.variance₂ = 1 := varOf_pauliY

theorem pauliPair_commutator_eq_two : pauliPair.commutator_expectation = 2 := by
  show Complex.abs (expVal (pauliX * pauliY - pauliY * pauliX)) = 2
  rw [expVal_pauli_commutator]
  simp [Complex.abs_two, Complex.abs_I]

/-- The witness is genuinely non-commuting: its commutator expectation
    is nonzero, so `L8_2_robertson_trivial_when_commuting` does NOT apply. -/
theorem pauliPair_commutator_ne_zero : pauliPair.commutator_expectation ≠ 0 := by
  rw [pauliPair_commutator_eq_two]; norm_num

/-- **L8.2, non-trivial case.** The Robertson bound holds for a concrete
    pair whose commutator expectation is nonzero: 1 · 1 ≥ (½·|2|)² = 1. -/
theorem L8_2_robertson_pauli : robertsonBound pauliPair := by
  unfold robertsonBound
  rw [pauliPair_variance₁, pauliPair_variance₂, pauliPair_commutator_eq_two]
  norm_num

/-- For the Pauli pair the Robertson bound is TIGHT: equality holds,
    so the ½ prefactor cannot be improved. -/
theorem robertson_tight_for_pauli :
    pauliPair.variance₁ * pauliPair.variance₂ =
      (1 / 2 * |pauliPair.commutator_expectation|) ^ 2 := by
  rw [pauliPair_variance₁, pauliPair_variance₂, pauliPair_commutator_eq_two]
  norm_num

/-! ## Application to RSSN Shape Operators

    The commutator value is `|Triangle(Square(2)) − Square(Triangle(2))|`.
    Both sides are concrete naturals — `Triangle(Square 2) = Triangle 256 =
    256^256` and `Square(Triangle 2) = Square 4 = Triangle⁴(4)` — but far too
    large to evaluate in the kernel, so non-vanishing is proved structurally:
    `Triangle⁴(4) = Triangle(Triangle(Triangle 256))` strictly dominates
    `Triangle 256` because `Triangle` is strictly inflationary on inputs ≥ 2. -/

open RSSN.ShapeOperators in
/-- `|Triangle(Square(2)) - Square(Triangle(2))|` as a real number. -/
def shapeCommutatorValue : ℝ :=
  |(↑(triangle (square 2)) : ℝ) - ↑(square (triangle 2))|

open RSSN.ShapeOperators in
/-- The two operator orders genuinely differ at n = 2:
    `Triangle(Square 2) < Square(Triangle 2)`. -/
theorem shape_commutator_orders_differ :
    triangle (square 2) < square (triangle 2) := by
  rw [square_two, triangle_two]
  -- Goal: triangle 256 < square 4 = triangleIter 4 4. Unfold the iterate
  -- ONLY through the proven equation lemmas (`triangleIter_succ/zero`) with
  -- explicitly instantiated indices, never by `rfl`/defeq between the
  -- iterate and the unfolded tower: a definitional comparison there makes
  -- the elaborator normalize `Nat.pow` at (256^256)-sized arguments and
  -- panic — exactly the overflow the repo history records for (k, n)=(3, 4).
  -- Each `show` below only re-expresses a numeral as `_ + 1` under an
  -- unchanged head symbol (same-head congruence: cheap, no unfolding).
  rw [show square 4 = triangleIter 4 4 from rfl]
  show triangle 256 < triangleIter (3 + 1) 4
  rw [triangleIter_succ 3 4]
  show triangle 256 < triangle (triangleIter (2 + 1) 4)
  rw [triangleIter_succ 2 4]
  show triangle 256 < triangle (triangle (triangleIter (1 + 1) 4))
  rw [triangleIter_succ 1 4]
  show triangle 256 < triangle (triangle (triangle (triangleIter (0 + 1) 4)))
  rw [triangleIter_succ 0 4, triangleIter_zero 4, triangle_four]
  -- Goal: triangle 256 < triangle (triangle (triangle 256)).
  have h256 : 2 ≤ triangle 256 :=
    le_trans (by norm_num : (2:ℕ) ≤ 256) (triangle_ge_n 256 (by norm_num))
  have h1 : triangle 256 < triangle (triangle 256) := triangle_gt_self h256
  have h2 : triangle (triangle 256) < triangle (triangle (triangle 256)) :=
    triangle_gt_self (h256.trans h1.le)
  exact h1.trans h2

theorem shape_commutator_nonzero :
    0 < |shapeCommutatorValue| := by
  have hne : (↑(RSSN.ShapeOperators.triangle (RSSN.ShapeOperators.square 2)) : ℝ) ≠
      ↑(RSSN.ShapeOperators.square (RSSN.ShapeOperators.triangle 2)) := by
    exact_mod_cast Nat.ne_of_lt shape_commutator_orders_differ
  unfold shapeCommutatorValue
  rw [abs_abs]
  exact abs_pos.mpr (sub_ne_zero.mpr hne)

/-- The Robertson lower bound for the RSSN shape operators is strictly
    positive, because the shape commutator does not vanish. -/
theorem robertson_lower_bound_positive :
    0 < 1 / 2 * |shapeCommutatorValue| :=
  mul_pos one_half_pos shape_commutator_nonzero

/-- **T8, Robertson form, non-trivially.** Any pair of uncertainties whose
    product satisfies the Robertson bound for the shape operators is forced
    to have strictly positive product — by the bound alone, with NO
    positivity assumptions on Δ₁, Δ₂, since ½|⟨[S₁,S₂]⟩| > 0. -/
theorem T8_robertson_form_nontrivial :
    ∀ (Δ₁ Δ₂ : ℝ),
    Δ₁ * Δ₂ ≥ (1 / 2 * |shapeCommutatorValue|) →
    0 < Δ₁ * Δ₂ := by
  intro Δ₁ Δ₂ h
  exact lt_of_lt_of_le robertson_lower_bound_positive h

/-! ## Original Form Counterexample at n=1

    At n = 1 the genuine shape operators DO commute —
    `Triangle(Square 1) = Square(Triangle 1) = 1` — so no uniform
    positive lower bound on the commutator can hold across all n:
    the original (n-uniform) form of T8 fails at n = 1. -/

open RSSN.ShapeOperators in
/-- The shape operators commute at n = 1: both composition orders reduce
    to 1. (Kernel reduction is safe at these tiny values, unlike n = 2
    where `256^256` overflows the kernel.) -/
theorem T8_original_fails_at_1 :
    triangle (square 1) = square (triangle 1) := rfl

open RSSN.ShapeOperators in
/-- The ℤ-valued commutator of the genuine shape operators vanishes at
    n = 1 — contrast `RSSN.NonCommutativity.commutator_nonzero` at n = 2. -/
theorem original_commutator_zero_at_1 :
    RSSN.NonCommutativity.commutator triangle square 1 = 0 := by
  unfold RSSN.NonCommutativity.commutator
  rw [T8_original_fails_at_1]
  exact sub_self _

end RSSN.UncertaintyPrinciple
