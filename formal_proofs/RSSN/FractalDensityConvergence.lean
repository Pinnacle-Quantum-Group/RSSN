/-
  RSSN — Fractal Density Convergence (Theorem T1)
  Pinnacle Quantum Group — April 2026

  Proves convergence of fractal density D_k(n) = lim_{i→∞} F_i(n)/G_i
  for each shape operator. Triangle density = 1/n (constant ratio).
  Square density converges to 0 (exponentially shrinking ratio).
  Reference: RSSN README §3.3
-/
import Mathlib

noncomputable section
open Filter Topology BigOperators

namespace RSSN.FractalDensityConvergence

/-! ## 1. Density Sequence Structure -/

structure FractalDensitySeq where
  F : ℕ → ℝ
  G : ℕ → ℝ
  hF_nonneg : ∀ i, 0 ≤ F i
  hG_pos : ∀ i, 0 < G i

def FractalDensitySeq.ratio (s : FractalDensitySeq) (i : ℕ) : ℝ :=
  s.F i / s.G i

/-! ## 2. Triangle Density: Constant Ratio 1/n -/

def triangleDensitySeq (n : ℕ) (hn : 1 ≤ n) : FractalDensitySeq where
  F := fun i => (↑n : ℝ) ^ i
  G := fun i => (↑n : ℝ) ^ (i + 1)
  hF_nonneg := fun i => by positivity
  hG_pos := fun i => by positivity

theorem triangle_ratio_constant (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    (triangleDensitySeq n (by omega)).ratio i = 1 / (↑n : ℝ) := by
  unfold FractalDensitySeq.ratio triangleDensitySeq
  simp
  rw [pow_succ]
  field_simp
  ring

theorem triangle_density_converges (n : ℕ) (hn : 2 ≤ n) :
    Tendsto (triangleDensitySeq n (by omega)).ratio atTop (nhds (1 / (↑n : ℝ))) := by
  -- The ratio is the constant function `1/n`; rewrite via funext so simp can
  -- replace `(triangleDensitySeq n _).ratio` (a function) with `fun _ => 1/n`.
  have h : (triangleDensitySeq n (by omega)).ratio = fun _ => 1 / (↑n : ℝ) := by
    funext i; exact triangle_ratio_constant n hn i
  rw [h]
  exact tendsto_const_nhds

theorem triangle_density_positive (n : ℕ) (hn : 2 ≤ n) :
    0 < 1 / (↑n : ℝ) := by positivity

/-! ## 3. General Convergence: Monotone Bounded Ratio -/

theorem monotone_ratio_converges (s : FractalDensitySeq)
    (hmono : Antitone s.ratio) :
    ∃ L, Tendsto s.ratio atTop (nhds L) ∧ 0 ≤ L := by
  have hbdd : BddBelow (Set.range s.ratio) := by
    exact ⟨0, by rintro _ ⟨i, rfl⟩; exact div_nonneg (s.hF_nonneg i) (le_of_lt (s.hG_pos i))⟩
  -- ℝ is ConditionallyCompleteLattice — c-version of tendsto_atTop_iInf.
  exact ⟨iInf s.ratio, tendsto_atTop_ciInf hmono hbdd,
    le_ciInf fun i => div_nonneg (s.hF_nonneg i) (le_of_lt (s.hG_pos i))⟩

/-! ## 4. Geometric Decay to Zero -/

theorem geometric_density_vanishes (r : ℝ) (hr0 : 0 ≤ r) (hr1 : r < 1) :
    Tendsto (fun n => r ^ n) atTop (nhds 0) :=
  tendsto_pow_atTop_nhds_zero_of_lt_one hr0 hr1

/-! ## 5. Square Density: Converges to 0 -/

theorem square_density_vanishes_model (n : ℕ) (hn : 2 ≤ n) :
    let r := 1 / (↑n : ℝ)
    Tendsto (fun i => r ^ i) atTop (nhds 0) := by
  apply tendsto_pow_atTop_nhds_zero_of_lt_one
  · positivity
  · rw [div_lt_one (by positivity : (0 : ℝ) < ↑n)]; exact_mod_cast hn

/-! ## 6. Density Bounded in [0, 1] -/

theorem ratio_bounded (s : FractalDensitySeq)
    (hFG : ∀ i, s.F i ≤ s.G i) (i : ℕ) :
    0 ≤ s.ratio i ∧ s.ratio i ≤ 1 := by
  unfold FractalDensitySeq.ratio
  exact ⟨div_nonneg (s.hF_nonneg i) (le_of_lt (s.hG_pos i)),
         div_le_one_of_le (hFG i) (le_of_lt (s.hG_pos i))⟩

end RSSN.FractalDensityConvergence
