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

theorem triangle_gt_self {n : ℕ} (hn : 2 ≤ n) : n < triangle n := by
  unfold triangle
  -- n = n^1 < n^n since the base n ≥ 2 > 1 and the exponent grows 1 < n.
  calc n = n ^ 1 := (pow_one n).symm
    _ < n ^ n := pow_lt_pow_right (show (1:ℕ) < n by omega) (by omega)

/-! ## 2. Iterated Triangle -/

def triangleIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => triangle (triangleIter k n)

theorem triangleIter_zero (n : ℕ) : triangleIter 0 n = n := rfl

theorem triangleIter_one (n : ℕ) : triangleIter 1 n = triangle n := rfl

theorem triangleIter_succ (k n : ℕ) :
    triangleIter (k + 1) n = triangle (triangleIter k n) := rfl

/-- Every iterate keeps the argument as a lower bound: `n ≤ triangle^[k](n)`
    for `n ≥ 1`, since each application of `triangle` is inflationary. -/
theorem triangleIter_ge_arg (n : ℕ) (hn : 1 ≤ n) : ∀ k, n ≤ triangleIter k n := by
  intro k
  induction k with
  | zero => exact le_of_eq (triangleIter_zero n).symm
  | succ k ih =>
    calc n ≤ triangleIter k n := ih
      _ ≤ triangle (triangleIter k n) := triangle_ge_n _ (hn.trans ih)
      _ = triangleIter (k + 1) n := (triangleIter_succ k n).symm

/-- `triangleIter` is monotone in the iteration count (argument ≥ 1). -/
theorem triangleIter_mono_iter (n : ℕ) (hn : 1 ≤ n) :
    ∀ {j k : ℕ}, j ≤ k → triangleIter j n ≤ triangleIter k n := by
  intro j k hjk
  induction hjk with
  | refl => exact le_refl _
  | @step m _ ih =>
    calc triangleIter j n ≤ triangleIter m n := ih
      _ ≤ triangle (triangleIter m n) :=
          triangle_ge_n _ (hn.trans (triangleIter_ge_arg n hn m))
      _ = triangleIter (m + 1) n := (triangleIter_succ m n).symm

/-! ## 3. Square: Triangle^n(n) -/

def square (n : ℕ) : ℕ := triangleIter n n

theorem square_two : square 2 = 256 := by
  -- square 2 = triangleIter 2 2 = triangle (triangle 2) = triangle 4 = 256.
  show triangleIter 2 2 = 256
  rfl

theorem square_ge_triangle (n : ℕ) (hn : 2 ≤ n) : triangle n ≤ square n := by
  -- square n = triangleIter n n; triangle n = triangleIter 1 n; and
  -- triangleIter is monotone in the iteration count for arguments ≥ 1.
  show triangle n ≤ triangleIter n n
  calc triangle n = triangleIter 1 n := (triangleIter_one n).symm
    _ ≤ triangleIter n n := triangleIter_mono_iter n (by omega) (by omega)

/-! ## 4. Circle: Square^n(n) -/

def squareIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => square (squareIter k n)

def circle (n : ℕ) : ℕ := squareIter n n

theorem squareIter_zero (n : ℕ) : squareIter 0 n = n := rfl

theorem squareIter_one (n : ℕ) : squareIter 1 n = square n := rfl

theorem squareIter_succ (k n : ℕ) :
    squareIter (k + 1) n = square (squareIter k n) := rfl

/-- `square` is inflationary for arguments ≥ 1 (inherited from `triangleIter`). -/
theorem square_ge_self (n : ℕ) (hn : 1 ≤ n) : n ≤ square n :=
  triangleIter_ge_arg n hn n

/-- Every `squareIter` iterate keeps the argument as a lower bound (arg ≥ 1). -/
theorem squareIter_ge_arg (n : ℕ) (hn : 1 ≤ n) : ∀ k, n ≤ squareIter k n := by
  intro k
  induction k with
  | zero => exact le_of_eq (squareIter_zero n).symm
  | succ k ih =>
    calc n ≤ squareIter k n := ih
      _ ≤ square (squareIter k n) := square_ge_self _ (hn.trans ih)
      _ = squareIter (k + 1) n := (squareIter_succ k n).symm

/-- `squareIter` is monotone in the iteration count (argument ≥ 1). -/
theorem squareIter_mono_iter (n : ℕ) (hn : 1 ≤ n) :
    ∀ {j k : ℕ}, j ≤ k → squareIter j n ≤ squareIter k n := by
  intro j k hjk
  induction hjk with
  | refl => exact le_refl _
  | @step m _ ih =>
    calc squareIter j n ≤ squareIter m n := ih
      _ ≤ square (squareIter m n) :=
          square_ge_self _ (hn.trans (squareIter_ge_arg n hn m))
      _ = squareIter (m + 1) n := (squareIter_succ m n).symm

theorem circle_ge_square (n : ℕ) (hn : 2 ≤ n) : square n ≤ circle n := by
  -- circle n = squareIter n n; square n = squareIter 1 n; monotone in count.
  show square n ≤ squareIter n n
  calc square n = squareIter 1 n := (squareIter_one n).symm
    _ ≤ squareIter n n := squareIter_mono_iter n (by omega) (by omega)

/-! ## 5. Growth Hierarchy -/

theorem growth_hierarchy (n : ℕ) (hn : 2 ≤ n) :
    n ≤ triangle n ∧ triangle n ≤ square n := by
  exact ⟨triangle_ge_n n (by omega), square_ge_triangle n hn⟩

/-! ## 6. Triangle is Strictly Monotone for n ≥ 2 -/

theorem triangle_strict_mono {a b : ℕ} (ha : 2 ≤ a) (hab : a < b) :
    triangle a < triangle b := by
  unfold triangle
  -- a^a < a^b ≤ b^b. v4.5.0: `Nat.pow_lt_pow_right` doesn't exist; the
  -- generic `pow_lt_pow_right` resolves via the ℕ ordered semiring instance.
  have hb : 2 ≤ b := ha.trans hab.le
  calc a ^ a < a ^ b := pow_lt_pow_right (show (1:ℕ) < a by omega) hab
    _ ≤ b ^ b := Nat.pow_le_pow_left hab.le b

end RSSN.ShapeOperators
