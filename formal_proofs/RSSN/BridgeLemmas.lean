/-
  RSSN — Bridge Lemmas (B.1, B.2)
  Pinnacle Quantum Group — April 2026

  B.1: RSSN D_k(n) = discrete evaluation of RLA scale field D
  B.2: Triangle(n) corresponds to grade-1 Lie derivative at weight w = ln(n)
  These connect the RSSN notation layer to the RLA algebraic backbone.
  Reference: LEMMA_DERIVATIONS.md Bridge Lemmas
-/
import Mathlib

noncomputable section
open Real Filter Topology BigOperators  -- BigOperators needed for `∑`

namespace RSSN.BridgeLemmas

/-! ## B.1 — RSSN-RLA Correspondence
    D_k(n) in RSSN is the discrete evaluation of RLA scale field D
    at recursion level k. The twisted bracket α-twist encodes the
    density gradient. -/

structure DiscreteScaleField where
  D : ℕ → ℝ
  D_pos : ∀ k, 0 < D k
  D_bounded : ∀ k, D k ≤ 1

def discreteAlpha (sf : DiscreteScaleField) (k : ℕ) : ℝ :=
  log (sf.D (k + 1)) - log (sf.D k)

theorem B1_alpha_from_density (sf : DiscreteScaleField) (k : ℕ) :
    discreteAlpha sf k = log (sf.D (k + 1) / sf.D k) := by
  unfold discreteAlpha
  rw [log_div (ne_of_gt (sf.D_pos (k + 1))) (ne_of_gt (sf.D_pos k))]

theorem B1_alpha_telescopes (sf : DiscreteScaleField) (n : ℕ) :
    ∑ k in Finset.range n, discreteAlpha sf k =
    log (sf.D n) - log (sf.D 0) := by
  induction n with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ, ih]
    unfold discreteAlpha; ring

theorem B1_constant_field_zero_alpha (c : ℝ) (hc : 0 < c) (hc1 : c ≤ 1) :
    ∀ k, discreteAlpha ⟨fun _ => c, fun _ => hc, fun _ => hc1⟩ k = 0 := by
  intro k; unfold discreteAlpha; simp

/-! ## B.2 — Generator Identification
    Triangle(n) = n^n = e^{n ln n}.
    Setting D = n, ι_X(α) = ln(n), weight w = ln(n). -/

theorem B2_triangle_as_exp (n : ℕ) (hn : 0 < n) :
    (↑n : ℝ) ^ (n : ℕ) = exp (↑n * log ↑n) := by
  rw [← rpow_nat_cast (↑n : ℝ) n]
  rw [rpow_def_of_pos (Nat.cast_pos.mpr hn)]
  congr 1; push_cast; ring

def triangleWeight (n : ℕ) : ℝ := log ↑n

theorem B2_weight_positive (n : ℕ) (hn : 2 ≤ n) :
    0 < triangleWeight n := by
  unfold triangleWeight
  exact log_pos (by exact_mod_cast hn)

theorem B2_grade1_scaling (n : ℕ) (hn : 0 < n) :
    exp (1 * triangleWeight n * ↑n) = (↑n : ℝ) ^ (n : ℕ) := by
  unfold triangleWeight
  rw [one_mul, ← rpow_nat_cast (↑n : ℝ) n, rpow_def_of_pos (Nat.cast_pos.mpr hn)]
  -- `congr 1` already closes by reflexivity after the rewrite chain;
  -- the trailing `ring` would error with "no goals to be solved".
  congr 1

/-! ## B.1 + B.2 Combined: Density Gradient = Bracket Twist -/

structure BracketTwistData where
  density : ℕ → ℝ
  weight : ℝ
  alpha : ℕ → ℝ
  h_alpha : ∀ k, alpha k = log (density (k + 1)) - log (density k)

theorem bridge_twist_from_density (bt : BracketTwistData) (k : ℕ) :
    bt.alpha k = log (bt.density (k + 1)) - log (bt.density k) :=
  bt.h_alpha k

theorem bridge_zero_twist_constant_density (bt : BracketTwistData)
    (hconst : ∀ k, bt.density k = bt.density 0) :
    ∀ k, bt.alpha k = 0 := by
  intro k; rw [bt.h_alpha, hconst (k + 1), hconst k, sub_self]

end RSSN.BridgeLemmas
