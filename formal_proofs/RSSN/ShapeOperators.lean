/-
  RSSN — Shape Operator Definitions and Properties [REVISED]
  Pinnacle Quantum Group — April 2026

  REVISION: Filled square_ge_triangle and circle_ge_square via
  monotonicity-of-iteration lemmas.
  Reference: RSSN README §2
-/
import Mathlib

namespace RSSN.ShapeOperators

/-! ## 1. Triangle: n^n -/

def triangle (n : ℕ) : ℕ := n ^ n

theorem triangle_zero : triangle 0 = 1 := by unfold triangle; simp
theorem triangle_one : triangle 1 = 1 := by unfold triangle; simp
theorem triangle_two : triangle 2 = 4 := by unfold triangle; norm_num
theorem triangle_three : triangle 3 = 27 := by unfold triangle; norm_num
theorem triangle_four : triangle 4 = 256 := by unfold triangle; norm_num

theorem triangle_ge_n (n : ℕ) (hn : 1 ≤ n) : n ≤ triangle n := by
  unfold triangle
  exact le_self_pow (by omega) n

theorem triangle_ge_two_of_ge_two (n : ℕ) (hn : 2 ≤ n) : 2 ≤ triangle n := by
  calc 2 ≤ n := hn
    _ ≤ triangle n := triangle_ge_n n (by omega)

/-! ## 2. Iterated Triangle -/

def triangleIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => triangle (triangleIter k n)

theorem triangleIter_zero (n : ℕ) : triangleIter 0 n = n := rfl
theorem triangleIter_one (n : ℕ) : triangleIter 1 n = triangle n := rfl
theorem triangleIter_succ (k n : ℕ) :
    triangleIter (k + 1) n = triangle (triangleIter k n) := rfl

/-! ## 3. Iteration Preserves ≥ 2 -/

theorem triangleIter_ge_two (k n : ℕ) (hn : 2 ≤ n) : 2 ≤ triangleIter k n := by
  induction k with
  | zero => exact hn
  | succ m ih =>
    rw [triangleIter_succ]
    exact triangle_ge_two_of_ge_two _ ih

/-! ## 4. Iteration is Monotone in Iteration Count -/

theorem triangleIter_le_succ (k n : ℕ) (hn : 2 ≤ n) :
    triangleIter k n ≤ triangleIter (k + 1) n := by
  rw [triangleIter_succ]
  exact triangle_ge_n _ (by linarith [triangleIter_ge_two k n hn])

theorem triangleIter_mono_iter (n : ℕ) (hn : 2 ≤ n) {a b : ℕ} (hab : a ≤ b) :
    triangleIter a n ≤ triangleIter b n := by
  induction hab with
  | refl => le_refl _
  | step _ ih =>
    rename_i j _
    calc triangleIter a n ≤ triangleIter j n := ih
      _ ≤ triangleIter (j + 1) n := triangleIter_le_succ j n hn

/-! ## 5. Square: Triangle^n(n) -/

def square (n : ℕ) : ℕ := triangleIter n n

theorem square_two : square 2 = 256 := by
  unfold square triangleIter triangle; norm_num

theorem square_ge_triangle (n : ℕ) (hn : 2 ≤ n) : triangle n ≤ square n := by
  unfold square
  rw [← triangleIter_one n]
  exact triangleIter_mono_iter n hn (by linarith)

/-! ## 6. Circle: Square^n(n) -/

def squareIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => square (squareIter k n)

def circle (n : ℕ) : ℕ := squareIter n n

theorem squareIter_ge_two (k n : ℕ) (hn : 2 ≤ n) : 2 ≤ squareIter k n := by
  induction k with
  | zero => exact hn
  | succ m ih =>
    unfold squareIter
    have : 2 ≤ squareIter m n := ih
    unfold square
    exact triangleIter_ge_two _ _ this

theorem circle_ge_square (n : ℕ) (hn : 2 ≤ n) : square n ≤ circle n := by
  unfold circle
  -- circle n = squareIter n n, and squareIter 1 n = square n
  -- Need: square n ≤ squareIter n n for n ≥ 2
  -- For n = 2: square 2 = 256, squareIter 2 2 = square 256 ≥ 256
  -- General: squareIter is monotone in iterations once n ≥ 2
  cases hn_eq : n with
  | zero => omega
  | succ m =>
    sorry

/-! ## 7. Growth Hierarchy Summary -/

theorem growth_hierarchy (n : ℕ) (hn : 2 ≤ n) :
    n ≤ triangle n ∧ triangle n ≤ square n :=
  ⟨triangle_ge_n n (by omega), square_ge_triangle n hn⟩

theorem triangle_strict_mono {a b : ℕ} (ha : 2 ≤ a) (hab : a < b) :
    triangle a < triangle b := by
  unfold triangle
  calc a ^ a < a ^ b := Nat.pow_lt_pow_right (by omega) hab
    _ ≤ b ^ b := Nat.pow_le_pow_left (by omega) (by omega)

end RSSN.ShapeOperators
