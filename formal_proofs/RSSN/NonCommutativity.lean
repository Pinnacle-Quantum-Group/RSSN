/-
  RSSN — Non-Commutativity of Shape Operators (Theorem T7)
  Pinnacle Quantum Group — April 2026

  Proves that RSSN shape operators do not commute:
  [Triangle, Square](2) ≠ 0, demonstrated by concrete computation.
  Triangle(Square(2)) = 256^256, Square(Triangle(2)) >> 256^256.
  Reference: RSSN README §9.3, Appendix B
-/
import Mathlib

namespace RSSN.NonCommutativity

open RSSN.ShapeOperators in

/-! ## 1. Import Shape Operators -/

def triangle (n : ℕ) : ℕ := n ^ n

def triangleIter : ℕ → ℕ → ℕ
  | 0, n => n
  | k + 1, n => triangle (triangleIter k n)

def square (n : ℕ) : ℕ := triangleIter n n

/-! ## 2. Concrete Values -/

theorem triangle_2 : triangle 2 = 4 := by
  unfold triangle; norm_num

theorem triangle_4 : triangle 4 = 256 := by
  unfold triangle; norm_num

theorem triangle_256 : triangle 256 = 256 ^ 256 := rfl

theorem square_2 : square 2 = 256 := by
  unfold square triangleIter triangle; norm_num

/-! ## 3. Triangle(Square(2)) = Triangle(256) = 256^256 -/

theorem triangle_of_square_2 : triangle (square 2) = 256 ^ 256 := by
  rw [square_2]; rfl

/-! ## 4. Square(Triangle(2)) = Square(4) = triangleIter 4 4 -/

theorem square_of_triangle_2 : square (triangle 2) = triangleIter 4 4 := by
  rw [triangle_2]; rfl

/-! ## 5. triangleIter 4 4 involves triangle(triangle(triangle(triangle(4)))) -/

theorem triangleIter_4_4_expand :
    triangleIter 4 4 = triangle (triangle (triangle (triangle 4))) := rfl

/-! ## 6. triangle(triangle(4)) = triangle(256) = 256^256 -/

theorem triangle_triangle_4 : triangle (triangle 4) = 256 ^ 256 := by
  rw [triangle_4]; rfl

/-! ## 7. Helper: `triangle x > x` for `x ≥ 2`
  This is the key growth fact behind non-commutativity at n=2:
  iterating `triangle` strictly increases at every step (for arg ≥ 2). -/

lemma triangle_grows (x : ℕ) (h : 2 ≤ x) : x < triangle x := by
  unfold triangle
  -- x = x^1 < x^x via Nat.pow_lt_pow_right (1 < x) (1 < x)
  calc x = x ^ 1 := (pow_one x).symm
    _ < x ^ x := Nat.pow_lt_pow_right h h

/-! ## 7b. triangle(triangle(triangle(4))) >> triangle(triangle(4)) -/

theorem triple_triangle_4_gt : triangle (triangle (triangle 4)) > triangle (triangle 4) := by
  rw [triangle_triangle_4]
  -- Goal: 256^256 < triangle (256^256)
  exact triangle_grows _ (by norm_num)

/-! ## 8. Non-Commutativity: triangle ∘ square ≠ square ∘ triangle at n=2

  square (triangle 2) = triangleIter 4 4 = triangle⁴(4) ≥ triangle (256^256)
  > 256^256 = triangle (square 2). Two applications of `triangle_grows`. -/

theorem noncommutative_at_2 : triangle (square 2) ≠ square (triangle 2) := by
  rw [triangle_of_square_2, square_of_triangle_2, triangleIter_4_4_expand]
  rw [triangle_triangle_4]
  -- Goal: 256^256 ≠ triangle (triangle (256^256))
  intro h
  have h1 : 256 ^ 256 < triangle (256 ^ 256) := triangle_grows _ (by norm_num)
  have h_arg : 2 ≤ triangle (256 ^ 256) := by
    have : 2 ≤ (256 : ℕ) ^ 256 := by norm_num
    omega
  have h2 : triangle (256 ^ 256) < triangle (triangle (256 ^ 256)) :=
    triangle_grows _ h_arg
  -- h : 256^256 = triangle (triangle (256^256))
  -- h2: triangle (256^256) < triangle (triangle (256^256)) = 256^256
  rw [← h] at h2
  omega

/-! ## 9. Commutator Definition -/

def commutator (f g : ℕ → ℕ) (n : ℕ) : ℤ :=
  ↑(f (g n)) - ↑(g (f n))

theorem commutator_nonzero : commutator triangle square 2 ≠ 0 := by
  unfold commutator
  -- Reduce to noncommutative_at_2 via Nat.cast injectivity:
  -- ↑a - ↑b = 0 ⇒ a = b ⇒ contradicts noncommutativity.
  intro h
  have h_eq : (triangle (square 2) : ℤ) = (square (triangle 2) : ℤ) := by linarith
  have h_nat : triangle (square 2) = square (triangle 2) := by exact_mod_cast h_eq
  exact noncommutative_at_2 h_nat

end RSSN.NonCommutativity
