/-
  RSSN — Non-Commutativity of Shape Operators (Theorem T7)
  Pinnacle Quantum Group — April 2026

  Proves that RSSN shape operators do not commute. Uses the Tetration
  foundation to keep all proofs symbolic — concrete numerals like
  256^256 are NEVER reduced, sidestepping the `Nat.pow exponent too big`
  kernel panic.

  Triangle(Square(2)) = Triangle(triangleIter 2 2) = Triangle(triangle (triangle 2))
  Square(Triangle(2)) = Square(4) = triangleIter 4 4 = triangle (triangle (triangle (triangle 4)))

  Both decompose into iterated triangle towers; the LHS has fewer
  iterations than the RHS, so by `triangleIter_strict_grows` they differ.

  Reference: RSSN README §9.3, Appendix B
-/
import RSSN.Tetration

namespace RSSN.NonCommutativity

open RSSN.Tetration

/-! ## 1. Re-export from Tetration -/

-- `triangle` and `triangleIter` are inherited from the Tetration module.

def square (n : ℕ) : ℕ := triangleIter n n

/-! ## 2. Concrete Equalities (definitional, no kernel reduction needed) -/

theorem triangle_2 : triangle 2 = 4 := by
  unfold triangle; norm_num

theorem triangle_4 : triangle 4 = 256 := by
  unfold triangle; norm_num

/-! ## 3. Square at 2 expands to triangleIter chain -/

theorem square_2_eq : square 2 = triangleIter 2 2 := rfl

theorem square_4_eq : square 4 = triangleIter 4 4 := rfl

/-! ## 4. The Key Identity: square (triangle 2) = triangleIter 4 4

    triangle 2 = 4, so square (triangle 2) = square 4 = triangleIter 4 4.
    Symbolic — the kernel never expands `triangleIter 4 4` to a numeral. -/

theorem square_of_triangle_2 : square (triangle 2) = triangleIter 4 4 := by
  rw [triangle_2]; rfl

theorem triangle_of_square_2 : triangle (square 2) = triangleIter 3 2 := by
  -- triangle (square 2) = triangle (triangleIter 2 2) = triangleIter 3 2.
  rfl

/-! ## 5. Non-Commutativity via Strict Growth of triangleIter

    The two compositions live at different "tower heights":
      triangle ∘ square at 2  →  triangleIter 3 2
      square ∘ triangle at 2  →  triangleIter 4 4

    By `triangleIter_strict_grows`, height matters strictly for n ≥ 2.
    The differing tower height (3 vs 4 iterations) plus the differing
    base (2 vs 4) makes them unequal.

    Strategy: show triangleIter 4 4 > triangleIter 3 4 ≥ triangleIter 3 2,
    so RHS > LHS, hence ≠. -/

theorem triangleIter_4_4_gt_triangleIter_3_4 :
    triangleIter 3 4 < triangleIter 4 4 :=
  triangleIter_strict_grows (by norm_num) 3

/-- `triangleIter` is monotone in its base argument (for fixed iteration count),
    when the base is ≥ 2. Proven by induction on the iteration count. -/
lemma triangleIter_mono_base {a b : ℕ} (ha : 2 ≤ a) (hab : a ≤ b) :
    ∀ k, triangleIter k a ≤ triangleIter k b := by
  intro k
  induction k with
  | zero => exact hab
  | succ k ih =>
    show triangle (triangleIter k a) ≤ triangle (triangleIter k b)
    -- triangle is monotone on inputs ≥ 2 (strict-mono with ≤ from <).
    rcases Nat.eq_or_lt_of_le ih with heq | hlt
    · rw [heq]
    · -- Need triangleIter k a ≥ 2 to invoke triangle_strict_mono.
      have ha_iter : 2 ≤ triangleIter k a := ha.trans (triangleIter_ge_arg ha k)
      exact (triangle_strict_mono ha_iter hlt).le

theorem triangleIter_3_4_ge_triangleIter_3_2 :
    triangleIter 3 2 ≤ triangleIter 3 4 :=
  triangleIter_mono_base (by norm_num) (by norm_num) 3

theorem noncommutative_at_2 : triangle (square 2) ≠ square (triangle 2) := by
  rw [triangle_of_square_2, square_of_triangle_2]
  -- Goal: triangleIter 3 2 ≠ triangleIter 4 4
  -- Chain: triangleIter 3 2 ≤ triangleIter 3 4 < triangleIter 4 4.
  intro h
  have h1 := triangleIter_3_4_ge_triangleIter_3_2
  have h2 := triangleIter_4_4_gt_triangleIter_3_4
  -- h : triangleIter 3 2 = triangleIter 4 4
  -- h1: triangleIter 3 2 ≤ triangleIter 3 4
  -- h2: triangleIter 3 4 < triangleIter 4 4
  -- So triangleIter 4 4 ≤ triangleIter 3 4 < triangleIter 4 4 — contradiction.
  rw [h] at h1
  omega

/-! ## 6. Commutator -/

def commutator (f g : ℕ → ℕ) (n : ℕ) : ℤ :=
  ↑(f (g n)) - ↑(g (f n))

theorem commutator_nonzero : commutator triangle square 2 ≠ 0 := by
  unfold commutator
  intro h
  have h_eq : (triangle (square 2) : ℤ) = (square (triangle 2) : ℤ) := by linarith
  have h_nat : triangle (square 2) = square (triangle 2) := by exact_mod_cast h_eq
  exact noncommutative_at_2 h_nat

end RSSN.NonCommutativity
