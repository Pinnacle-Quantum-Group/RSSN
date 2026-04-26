/-
  RSSN — Tetration Foundation (Knuth Up-Arrow ↑↑)
  Pinnacle Quantum Group — April 2026

  Knuth's tetration `tet a k = a ↑↑ k` defined by:
    tet a 0     = 1
    tet a (k+1) = a ^ tet a k

  Identifies the RSSN shape ladder via tetration:
    Triangle(n) = n^n              = tet n 2
    Square(n)   = Triangle^[n](n)  ≥ tet n (n+1)
    Circle(n)   = Square^[n](n)    (pentation height)

  This module provides the core monotonicity / growth lemmas that
  collapse per-shape arguments to one-liners. The classical
  non-commutativity of tetration (Goodstein 1947, Knuth 1976) becomes
  the single algebraic input, and the architectural argument for
  Triangle/Square/Circle non-commutativity follows uniformly.
-/
import Mathlib

namespace RSSN.Tetration

/-! ## 1. Tetration Definition -/

def tet : ℕ → ℕ → ℕ
  | _, 0     => 1
  | a, k + 1 => a ^ tet a k

@[simp] lemma tet_zero (a : ℕ) : tet a 0 = 1 := rfl

@[simp] lemma tet_succ (a k : ℕ) : tet a (k + 1) = a ^ tet a k := rfl

lemma tet_one (a : ℕ) : tet a 1 = a := by
  show a ^ tet a 0 = a
  simp

lemma tet_two (a : ℕ) : tet a 2 = a ^ a := by
  show a ^ tet a 1 = a ^ a
  rw [tet_one]

/-! ## 2. Positivity / Lower Bounds -/

lemma tet_pos {a : ℕ} (ha : 1 ≤ a) : ∀ k, 0 < tet a k := by
  intro k
  induction k with
  | zero => simp
  | succ k _ =>
    show 0 < a ^ tet a k
    exact Nat.pos_pow_of_pos _ (by omega)

lemma tet_ge_one {a : ℕ} (ha : 1 ≤ a) (k : ℕ) : 1 ≤ tet a k := tet_pos ha k

/-! ## 3. The Workhorse: `n < a^n` for a ≥ 2, n ≥ 1

    Used by every height-monotonicity argument below. -/

lemma nat_lt_pow_self {a : ℕ} (ha : 2 ≤ a) :
    ∀ {n : ℕ}, 1 ≤ n → n < a ^ n := by
  -- Intro `n` explicitly (it's an implicit binder in the goal).
  -- Don't `intro hn` first: the eliminator reverts dependent hypotheses,
  -- so `hn` ends up in the goal as an implication. Intro per branch.
  intro n
  induction n with
  | zero =>
    intro hn
    -- hn : 1 ≤ 0 is decidably False.
    exact absurd hn (by decide)
  | succ k ih =>
    intro _
    rcases Nat.eq_zero_or_pos k with hk | hk
    · subst hk
      -- Goal: 1 < a^1 = a; ha : 2 ≤ a closes it.
      simp; omega
    · have ihk := ih hk
      -- Goal: k+1 < a^(k+1) = a * a^k
      rw [pow_succ]
      have hapk : 1 ≤ a ^ k := Nat.one_le_iff_ne_zero.mpr (pow_ne_zero _ (by omega))
      -- ihk : k < a^k. Then a^k ≥ k+1. Then a^k * a ≥ (k+1) * 2 = 2k + 2 > k+1.
      have h1 : k + 1 ≤ a ^ k := ihk
      nlinarith [h1, hapk, ha]

/-! ## 4. Height Monotonicity (Strict, for base ≥ 2) -/

lemma tet_lt_tet_succ {a : ℕ} (ha : 2 ≤ a) {k : ℕ} (_hk : 1 ≤ k) :
    tet a k < tet a (k + 1) := by
  show tet a k < a ^ tet a k
  -- The premise `1 ≤ k` isn't strictly needed (the bound holds for k=0 too
  -- since tet a 0 = 1 < a), but we keep it to match the callsite's
  -- semantic intent. `nat_lt_pow_self` infers n from the goal type.
  exact nat_lt_pow_self ha (tet_ge_one (by omega) k)

lemma tet_strict_mono_height {a : ℕ} (ha : 2 ≤ a) :
    ∀ {m n : ℕ}, 1 ≤ m → m ≤ n → tet a m ≤ tet a n := by
  intro m n hm hmn
  induction hmn with
  | refl => exact le_refl _
  | step h ih =>
    -- Bind the underscored hyp `h : m ≤ n'` so omega can see it.
    exact ih.trans (tet_lt_tet_succ ha (hm.trans h)).le

/-! ## 5. tet a k > a for k ≥ 2, a ≥ 2 -/

lemma a_lt_tet_two {a : ℕ} (ha : 2 ≤ a) : a < tet a 2 := by
  rw [tet_two]
  -- a < a^a — direct via `nat_lt_pow_self`, which is exactly this fact.
  exact nat_lt_pow_self ha (by omega)

lemma a_lt_tet {a : ℕ} (ha : 2 ≤ a) {k : ℕ} (hk : 2 ≤ k) : a < tet a k := by
  -- a < tet a 2 ≤ tet a k by height monotonicity.
  exact (a_lt_tet_two ha).trans_le
    (tet_strict_mono_height ha (by omega) hk)

/-! ## 6. Connection to RSSN Shape Operators -/

def triangle (n : ℕ) : ℕ := n ^ n

lemma triangle_eq_tet (n : ℕ) : triangle n = tet n 2 := by
  rw [triangle, tet_two]

/-- Single-line corollary: `triangle x > x` for `x ≥ 2`. -/
lemma triangle_grows {x : ℕ} (hx : 2 ≤ x) : x < triangle x := by
  rw [triangle_eq_tet]; exact a_lt_tet_two hx

/-- `triangleIter k n = triangle^[k] n`. -/
def triangleIter : ℕ → ℕ → ℕ
  | 0, n     => n
  | k + 1, n => triangle (triangleIter k n)

@[simp] lemma triangleIter_zero (n : ℕ) : triangleIter 0 n = n := rfl
@[simp] lemma triangleIter_succ (k n : ℕ) :
    triangleIter (k + 1) n = triangle (triangleIter k n) := rfl

/-- All iterates stay ≥ n (and hence ≥ 2 if n ≥ 2). -/
lemma triangleIter_ge_arg {n : ℕ} (hn : 2 ≤ n) :
    ∀ k, n ≤ triangleIter k n := by
  intro k
  induction k with
  | zero => exact le_refl n
  | succ k ih =>
    -- triangleIter (k+1) n = triangle (triangleIter k n) ≥ triangleIter k n ≥ n
    show n ≤ triangle (triangleIter k n)
    have hge : 2 ≤ triangleIter k n := hn.trans ih
    exact ih.trans (triangle_grows hge).le

/-- `triangle` is strictly monotone on inputs `≥ 2`. -/
lemma triangle_strict_mono {x y : ℕ} (hx : 2 ≤ x) (hxy : x < y) :
    triangle x < triangle y := by
  unfold triangle
  -- x^x < y^y. Use x^x ≤ y^x < y^y.
  have hy : 2 ≤ y := hx.trans hxy.le
  -- v4.5.0: `Nat.pow_lt_pow_right` is unknown; use generic `pow_lt_pow_right`
  -- (resolves via the StrictOrderedSemiring instance on ℕ).
  calc x ^ x
      ≤ y ^ x := Nat.pow_le_pow_left hxy.le x
    _ < y ^ y := pow_lt_pow_right (show (1:ℕ) < y by omega) hxy

/-- For `n ≥ 2`, iterating `triangle` strictly grows.
    Source of non-commutativity for the Triangle/Square ladder. -/
lemma triangleIter_strict_grows {n : ℕ} (hn : 2 ≤ n) (k : ℕ) :
    triangleIter k n < triangleIter (k + 1) n := by
  induction k with
  | zero =>
    show n < triangle n
    exact triangle_grows hn
  | succ k ih =>
    show triangle (triangleIter k n) < triangle (triangleIter (k + 1) n)
    have h_arg_ge : 2 ≤ triangleIter k n := hn.trans (triangleIter_ge_arg hn k)
    exact triangle_strict_mono h_arg_ge ih

end RSSN.Tetration
