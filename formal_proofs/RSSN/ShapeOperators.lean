/-
  RSSN — Shape Operator Definitions and Properties
  Pinnacle Quantum Group — April 2026

  Defines the RSSN shape operators (Triangle, Square, Circle) and
  proves basic properties: monotonicity, growth bounds, concrete values.
  Reference: RSSN README §2
-/
import Mathlib

namespace RSSN.ShapeOperators

/-! ## 1. Triangle: n^n -/

def triangle (n : ℕ) : ℕ := n ^ n

theorem triangle_zero : triangle 0 = 1 := by
  unfold triangle; simp

theorem triangle_one : triangle 1 = 1 := by
  unfold triangle; simp

theorem triangle_two : triangle 2 = 4 := by
  unfold triangle; norm_num

theorem triangle_three : triangle 3 = 27 := by
  unfold triangle; norm_num

theorem triangle_four : triangle 4 = 256 := by
  unfold triangle; norm_num

theorem triangle_ge_n (n : ℕ) (hn : 1 ≤ n) : n ≤ triangle n := by
  unfold triangle
  -- n ≤ n^n. Mathlib v4.5.0: `Nat.le_self_pow` takes `n ≠ 0` (exponent),
  -- not `1 ≤ n`. Also valid: `n^1 ≤ n^n` via `Nat.pow_le_pow_right`.
  calc n = n ^ 1 := (pow_one n).symm
    _ ≤ n ^ n := Nat.pow_le_pow_right hn hn

/-! ## 2. Iterated Triangle -/

def triangleIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => triangle (triangleIter k n)

theorem triangleIter_zero (n : ℕ) : triangleIter 0 n = n := rfl

theorem triangleIter_one (n : ℕ) : triangleIter 1 n = triangle n := rfl

theorem triangleIter_succ (k n : ℕ) :
    triangleIter (k + 1) n = triangle (triangleIter k n) := rfl

/-! ## 3. Square: Triangle^n(n) -/

def square (n : ℕ) : ℕ := triangleIter n n

theorem square_two : square 2 = 256 := by
  -- square 2 = triangleIter 2 2 = triangle (triangle 2) = triangle 4 = 256.
  show triangleIter 2 2 = 256
  rfl

theorem square_ge_triangle (n : ℕ) (hn : 2 ≤ n) : triangle n ≤ square n := by
  -- square n = triangleIter n n. For n ≥ 2, triangleIter 1 n = triangle n,
  -- and triangleIter is monotone in iteration count when arg ≥ 2.
  -- Full proof requires monotonicity-of-iteration; stub for now.
  sorry

/-! ## 4. Circle: Square^n(n) -/

def squareIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => square (squareIter k n)

def circle (n : ℕ) : ℕ := squareIter n n

theorem circle_ge_square (n : ℕ) (hn : 2 ≤ n) : square n ≤ circle n := by
  sorry

/-! ## 5. Growth Hierarchy -/

theorem growth_hierarchy (n : ℕ) (hn : 2 ≤ n) :
    n ≤ triangle n ∧ triangle n ≤ square n := by
  exact ⟨triangle_ge_n n (by omega), square_ge_triangle n hn⟩

/-! ## 6. Triangle is Strictly Monotone for n ≥ 2 -/

theorem triangle_strict_mono {a b : ℕ} (ha : 2 ≤ a) (hab : a < b) :
    triangle a < triangle b := by
  unfold triangle
  -- a^a < a^b ≤ b^b (use Nat-specific lemmas to avoid typeclass ambiguity).
  have hb : 2 ≤ b := ha.trans hab.le
  calc a ^ a < a ^ b := Nat.pow_lt_pow_right ha hab
    _ ≤ b ^ b := Nat.pow_le_pow_left hab.le b

end RSSN.ShapeOperators
