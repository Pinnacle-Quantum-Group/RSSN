/-
  RSSN — Non-Commutativity (T7) [REVISED v3]
  Pinnacle Quantum Group — April 2026

  REVISION: Cleaner proof using triangleIter monotonicity.
  Proves [Triangle, Square](2) ≠ 0 via explicit chain.
-/
import Mathlib

namespace RSSN.NonCommutativity

/-! ## 1. Operators -/

def triangle (n : ℕ) : ℕ := n ^ n

def triangleIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => triangle (triangleIter k n)

def square (n : ℕ) : ℕ := triangleIter n n

/-! ## 2. Concrete Values -/

theorem triangle_2 : triangle 2 = 4 := by unfold triangle; norm_num
theorem triangle_4 : triangle 4 = 256 := by unfold triangle; norm_num

theorem square_2 : square 2 = 256 := by
  unfold square triangleIter triangle; norm_num

theorem triangleIter_4_4_step1 : triangleIter 1 4 = 256 := by
  unfold triangleIter triangle; norm_num

theorem triangleIter_4_4_step2 : triangleIter 2 4 = 256 ^ 256 := by
  unfold triangleIter; rw [triangleIter_4_4_step1]; unfold triangle

/-! ## 3. triangle(square 2) computation -/

theorem triangle_of_square_2 : triangle (square 2) = 256 ^ 256 := by
  rw [square_2]; unfold triangle

/-! ## 4. square(triangle 2) computation -/

theorem square_of_triangle_2 : square (triangle 2) = triangleIter 4 4 := by
  rw [triangle_2]; unfold square

/-! ## 5. triangleIter is monotone in iterations (for arg ≥ 2) -/

theorem triangle_ge_n {n : ℕ} (hn : 1 ≤ n) : n ≤ triangle n := by
  unfold triangle; exact le_self_pow (by omega) n

theorem triangleIter_ge_two (k : ℕ) {n : ℕ} (hn : 2 ≤ n) : 2 ≤ triangleIter k n := by
  induction k with
  | zero => exact hn
  | succ m ih =>
    show 2 ≤ triangle (triangleIter m n)
    calc 2 ≤ triangleIter m n := ih
      _ ≤ triangle (triangleIter m n) := triangle_ge_n (by linarith)

theorem triangleIter_strict_mono_succ {k : ℕ} {n : ℕ} (hn : 2 ≤ n)
    (h_iter_ge : 2 ≤ triangleIter k n) :
    triangleIter k n < triangleIter (k + 1) n := by
  show triangleIter k n < triangle (triangleIter k n)
  unfold triangle
  exact Nat.lt_pow_self (by linarith) (by linarith)

/-! ## 6. triangleIter 4 4 > triangleIter 2 4 -/

theorem step3_gt_step2 : triangleIter 2 4 < triangleIter 3 4 :=
  triangleIter_strict_mono_succ (by norm_num) (triangleIter_ge_two 2 (by norm_num))

theorem step4_gt_step3 : triangleIter 3 4 < triangleIter 4 4 :=
  triangleIter_strict_mono_succ (by norm_num) (triangleIter_ge_two 3 (by norm_num))

theorem step4_gt_step2 : triangleIter 2 4 < triangleIter 4 4 :=
  step3_gt_step2.trans step4_gt_step3

/-! ## 7. Non-Commutativity -/

theorem noncommutative_at_2 : triangle (square 2) ≠ square (triangle 2) := by
  rw [triangle_of_square_2, square_of_triangle_2]
  rw [← triangleIter_4_4_step2]
  exact ne_of_lt step4_gt_step2

/-! ## 8. Commutator -/

def commutator (f g : ℕ → ℕ) (n : ℕ) : ℤ :=
  ↑(f (g n)) - ↑(g (f n))

theorem commutator_sign : commutator triangle square 2 < 0 := by
  unfold commutator
  rw [triangle_of_square_2, square_of_triangle_2]
  have h : 256 ^ 256 < triangleIter 4 4 := by
    rw [← triangleIter_4_4_step2]; exact step4_gt_step2
  omega

theorem commutator_nonzero : commutator triangle square 2 ≠ 0 :=
  ne_of_lt commutator_sign

end RSSN.NonCommutativity
