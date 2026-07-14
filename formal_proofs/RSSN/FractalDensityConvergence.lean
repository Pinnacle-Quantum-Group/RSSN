/-
  RSSN — Fractal Density Convergence (Theorem T1)
  Pinnacle Quantum Group — April 2026

  Proves convergence of fractal density D_k(n) = lim_{i→∞} F_i(n)/G_i
  for the shape operators, using explicit density sequences:
  * Triangle: F_i = n^i, G_i = n^(i+1) — constant ratio, density = 1/n.
  * Square:   F_i = n^i, G_i = n^(n^i).  G follows the README §3.2 recursion
    G_0 = n, G_{i+1} = G_i^n exactly; the n-term sum defining F_i is
    abstracted to growth by a factor of n per level.  The ratio is antitone,
    dominated by the geometric sequence (1/n)^i, and converges to 0.
  * Circle:   F_i = n^i, G_i given by the README §3.2 recursion
    G_0 = n, G_{i+1} = G_i^{G_i} (same abstraction for F).  The ratio is
    antitone, dominated by the Square ratio, and converges to 0.
  The general antitone-bounded lemma `monotone_ratio_converges` is
  instantiated by both the Square and Circle sequences.
  Reference: RSSN README §3.2–3.3
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

/-! ## 5. Square Density: Converges to 0

The README §3.2 recursion for Square is `G_0 = n`, `G_{i+1} = G_i^n`, whose
closed form is `G_i = n^(n^i)`.  The substructure count `F_i` (an `n`-term
sum over Triangle iterates in the README) is abstracted here to growth by a
factor of `n` per level, i.e. `F_i = n^i`. -/

/-- `i + i ≤ n ^ i` for `n ≥ 2`: the exponent gap that makes the Square ratio
dominated by `(1/n)^i`. -/
private lemma two_mul_le_pow {n : ℕ} (hn : 2 ≤ n) (i : ℕ) : i + i ≤ n ^ i := by
  induction i with
  | zero => simp
  | succ k ih =>
    rcases Nat.eq_zero_or_pos k with hk | hk
    · subst hk; rw [pow_one]; omega
    · have h2 : 2 ≤ n ^ k := le_trans hn (by
        calc n = n ^ 1 := (pow_one n).symm
          _ ≤ n ^ k := Nat.pow_le_pow_right (by omega) hk)
      calc (k + 1) + (k + 1) = (k + k) + 2 := by omega
        _ ≤ n ^ k + n ^ k := Nat.add_le_add ih h2
        _ = 2 * n ^ k := by ring
        _ ≤ n * n ^ k := Nat.mul_le_mul_right _ hn
        _ = n ^ (k + 1) := by rw [pow_succ]

/-- `1 + n ^ i ≤ n ^ (i + 1)` for `n ≥ 2`: the exponent step that makes the
Square ratio antitone. -/
private lemma succ_pow_bound {n : ℕ} (hn : 2 ≤ n) (i : ℕ) : 1 + n ^ i ≤ n ^ (i + 1) := by
  have h1 : 1 ≤ n ^ i := Nat.one_le_pow i n (by omega)
  calc 1 + n ^ i ≤ n ^ i + n ^ i := by omega
    _ = 2 * n ^ i := by ring
    _ ≤ n * n ^ i := Nat.mul_le_mul_right _ hn
    _ = n ^ (i + 1) := by rw [pow_succ]

/-- The Square density sequence: `F_i = n^i` and `G_i = n^(n^i)`.
`G` satisfies the README recursion exactly: `G_0 = n^(n^0) = n` and
`G_{i+1} = n^(n^i · n) = (G_i)^n`. -/
def squareDensitySeq (n : ℕ) (hn : 2 ≤ n) : FractalDensitySeq where
  F := fun i => (↑n : ℝ) ^ i
  G := fun i => (↑n : ℝ) ^ (n ^ i)
  hF_nonneg := fun i => pow_nonneg (Nat.cast_nonneg n) i
  hG_pos := fun _ => pow_pos (by exact_mod_cast (by omega : 0 < n)) _

/-- Sanity check: `squareDensitySeq` really satisfies the README recursion
`G_{i+1} = G_i ^ n`. -/
theorem squareDensitySeq_G_recursion (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    (squareDensitySeq n hn).G (i + 1) = (squareDensitySeq n hn).G i ^ n := by
  show (↑n : ℝ) ^ (n ^ (i + 1)) = ((↑n : ℝ) ^ (n ^ i)) ^ n
  rw [← pow_mul, pow_succ']

/-- The Square ratio is dominated by the geometric sequence `(1/n)^i`. -/
theorem square_ratio_le_geo (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    (squareDensitySeq n hn).ratio i ≤ (1 / (↑n : ℝ)) ^ i := by
  have hn0 : (0 : ℝ) < ↑n := by exact_mod_cast (by omega : 0 < n)
  have hn1 : (1 : ℝ) ≤ ↑n := by exact_mod_cast (by omega : 1 ≤ n)
  show (↑n : ℝ) ^ i / (↑n : ℝ) ^ (n ^ i) ≤ (1 / (↑n : ℝ)) ^ i
  rw [div_pow, one_pow, div_le_div_iff (pow_pos hn0 _) (pow_pos hn0 _), one_mul,
    ← pow_add]
  exact pow_le_pow_right hn1 (two_mul_le_pow hn i)

/-- The Square density ratio `F_i/G_i = n^i / n^(n^i)` converges to 0:
Theorem T1 for the Square operator. -/
theorem square_density_vanishes (n : ℕ) (hn : 2 ≤ n) :
    Tendsto (squareDensitySeq n hn).ratio atTop (nhds 0) := by
  have hn0 : (0 : ℝ) < ↑n := by exact_mod_cast (by omega : 0 < n)
  refine squeeze_zero (fun i => ?_) (square_ratio_le_geo n hn) ?_
  · exact div_nonneg ((squareDensitySeq n hn).hF_nonneg i)
      (le_of_lt ((squareDensitySeq n hn).hG_pos i))
  · exact geometric_density_vanishes (1 / ↑n) (by positivity)
      (by rw [div_lt_one hn0]; exact_mod_cast (by omega : 1 < n))

/-- The Square ratio is antitone, so it also falls under the general
convergence lemma `monotone_ratio_converges`. -/
theorem square_ratio_antitone (n : ℕ) (hn : 2 ≤ n) :
    Antitone (squareDensitySeq n hn).ratio := by
  have hn0 : (0 : ℝ) < ↑n := by exact_mod_cast (by omega : 0 < n)
  have hn1 : (1 : ℝ) ≤ ↑n := by exact_mod_cast (by omega : 1 ≤ n)
  refine antitone_nat_of_succ_le fun i => ?_
  show (↑n : ℝ) ^ (i + 1) / (↑n : ℝ) ^ (n ^ (i + 1)) ≤ (↑n : ℝ) ^ i / (↑n : ℝ) ^ (n ^ i)
  rw [div_le_div_iff (pow_pos hn0 _) (pow_pos hn0 _), ← pow_add, ← pow_add]
  refine pow_le_pow_right hn1 ?_
  have h := succ_pow_bound hn i
  omega

/-- `monotone_ratio_converges` is inhabited: applied to the Square sequence it
yields a nonnegative limit. -/
theorem square_density_limit_exists (n : ℕ) (hn : 2 ≤ n) :
    ∃ L, Tendsto (squareDensitySeq n hn).ratio atTop (nhds L) ∧ 0 ≤ L :=
  monotone_ratio_converges _ (square_ratio_antitone n hn)

/-- Any limit of the Square ratio is 0 — identifies the limit produced by
`square_density_limit_exists`. -/
theorem square_density_limit_eq_zero (n : ℕ) (hn : 2 ≤ n) (L : ℝ)
    (hL : Tendsto (squareDensitySeq n hn).ratio atTop (nhds L)) : L = 0 :=
  tendsto_nhds_unique hL (square_density_vanishes n hn)

/-- Geometric envelope only (kept for reference): the dominating sequence
`(1/n)^i` itself tends to 0.  Superseded by `square_density_vanishes`, which
proves the statement for the actual Square ratio `F_i/G_i`. -/
theorem square_density_vanishes_model (n : ℕ) (hn : 2 ≤ n) :
    let r := 1 / (↑n : ℝ)
    Tendsto (fun i => r ^ i) atTop (nhds 0) := by
  apply tendsto_pow_atTop_nhds_zero_of_lt_one
  · positivity
  · rw [div_lt_one (by exact_mod_cast (by omega : 0 < n) : (0 : ℝ) < ↑n)]
    exact_mod_cast hn

/-! ## 6. Circle Density: Converges to 0

The README §3.2 recursion for Circle is `G_0 = n`, `G_{i+1} = G_i^{G_i}`
(doubly exponential growth); `F_i` is abstracted to `n^i` as for Square.
The Circle configuration space dominates the Square one
(`n^(n^i) ≤ G_i`), so the Circle ratio is squeezed by the Square ratio. -/

/-- Total configuration space for Circle: `G_0 = n`, `G_{i+1} = G_i ^ G_i`
(README §3.2). -/
def circleG (n : ℕ) : ℕ → ℕ
  | 0 => n
  | i + 1 => circleG n i ^ circleG n i

/-- The Circle space dominates the Square space: `n^(n^i) ≤ circleG n i`. -/
private lemma circleG_ge_square {n : ℕ} (hn : 2 ≤ n) : ∀ i, n ^ n ^ i ≤ circleG n i := by
  intro i
  induction i with
  | zero => simp [circleG]
  | succ i ih =>
    have hpos : 0 < n := by omega
    have hni : 1 ≤ n ^ i := Nat.one_le_pow i n hpos
    have hgn : n ≤ circleG n i := by
      calc n = n ^ 1 := (pow_one n).symm
        _ ≤ n ^ n ^ i := Nat.pow_le_pow_right hpos hni
        _ ≤ circleG n i := ih
    simp only [circleG]
    calc n ^ n ^ (i + 1) = (n ^ n ^ i) ^ n := by rw [← pow_mul, pow_succ']
      _ ≤ circleG n i ^ n := Nat.pow_le_pow_left ih n
      _ ≤ circleG n i ^ circleG n i :=
          Nat.pow_le_pow_right (lt_of_lt_of_le hpos hgn) hgn

private lemma circleG_pos {n : ℕ} (hn : 2 ≤ n) (i : ℕ) : 0 < circleG n i :=
  lt_of_lt_of_le Nat.zero_lt_one
    (le_trans (Nat.one_le_pow _ n (by omega)) (circleG_ge_square hn i))

private lemma circleG_ge_base {n : ℕ} (hn : 2 ≤ n) (i : ℕ) : n ≤ circleG n i :=
  calc n = n ^ 1 := (pow_one n).symm
    _ ≤ n ^ n ^ i := Nat.pow_le_pow_right (by omega) (Nat.one_le_pow i n (by omega))
    _ ≤ circleG n i := circleG_ge_square hn i

/-- One recursion step gains at least a factor of `n`:
`n * circleG n i ≤ circleG n (i+1)`.  Drives antitonicity of the ratio. -/
private lemma circleG_mul_le_succ {n : ℕ} (hn : 2 ≤ n) (i : ℕ) :
    n * circleG n i ≤ circleG n (i + 1) := by
  have hgn : n ≤ circleG n i := circleG_ge_base hn i
  have hg2 : 2 ≤ circleG n i := le_trans hn hgn
  simp only [circleG]
  calc n * circleG n i ≤ circleG n i * circleG n i :=
        Nat.mul_le_mul_right _ hgn
    _ = circleG n i ^ 2 := (pow_two _).symm
    _ ≤ circleG n i ^ circleG n i := Nat.pow_le_pow_right (by omega) hg2

/-- The Circle density sequence: `F_i = n^i` and `G_i` the README
`G_{i+1} = G_i^{G_i}` recursion. -/
def circleDensitySeq (n : ℕ) (hn : 2 ≤ n) : FractalDensitySeq where
  F := fun i => (↑n : ℝ) ^ i
  G := fun i => (↑(circleG n i) : ℝ)
  hF_nonneg := fun i => pow_nonneg (Nat.cast_nonneg n) i
  hG_pos := fun i => by
    show (0 : ℝ) < ↑(circleG n i)
    exact_mod_cast circleG_pos hn i

/-- The Circle ratio is dominated by the Square ratio (its configuration
space grows at least as fast). -/
theorem circle_ratio_le_square (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    (circleDensitySeq n hn).ratio i ≤ (squareDensitySeq n hn).ratio i := by
  have hn0 : (0 : ℝ) < ↑n := by exact_mod_cast (by omega : 0 < n)
  show (↑n : ℝ) ^ i / (↑(circleG n i) : ℝ) ≤ (↑n : ℝ) ^ i / (↑n : ℝ) ^ (n ^ i)
  refine div_le_div_of_le_left (pow_nonneg hn0.le i) (pow_pos hn0 _) ?_
  exact_mod_cast circleG_ge_square hn i

/-- The Circle ratio is dominated by the geometric sequence `(1/n)^i`. -/
theorem circle_ratio_le_geo (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    (circleDensitySeq n hn).ratio i ≤ (1 / (↑n : ℝ)) ^ i :=
  le_trans (circle_ratio_le_square n hn i) (square_ratio_le_geo n hn i)

/-- The Circle density ratio `F_i/G_i` converges to 0: Theorem T1 for the
Circle operator. -/
theorem circle_density_vanishes (n : ℕ) (hn : 2 ≤ n) :
    Tendsto (circleDensitySeq n hn).ratio atTop (nhds 0) := by
  have hn0 : (0 : ℝ) < ↑n := by exact_mod_cast (by omega : 0 < n)
  refine squeeze_zero (fun i => ?_) (circle_ratio_le_geo n hn) ?_
  · exact div_nonneg ((circleDensitySeq n hn).hF_nonneg i)
      (le_of_lt ((circleDensitySeq n hn).hG_pos i))
  · exact geometric_density_vanishes (1 / ↑n) (by positivity)
      (by rw [div_lt_one hn0]; exact_mod_cast (by omega : 1 < n))

/-- The Circle ratio is antitone: a second instantiation of the hypothesis of
`monotone_ratio_converges`. -/
theorem circle_ratio_antitone (n : ℕ) (hn : 2 ≤ n) :
    Antitone (circleDensitySeq n hn).ratio := by
  refine antitone_nat_of_succ_le fun i => ?_
  have hG : ∀ j, (0 : ℝ) < ↑(circleG n j) := fun j => by exact_mod_cast circleG_pos hn j
  show (↑n : ℝ) ^ (i + 1) / (↑(circleG n (i + 1)) : ℝ) ≤
    (↑n : ℝ) ^ i / (↑(circleG n i) : ℝ)
  rw [div_le_div_iff (hG (i + 1)) (hG i)]
  have key : n ^ (i + 1) * circleG n i ≤ n ^ i * circleG n (i + 1) := by
    calc n ^ (i + 1) * circleG n i = n ^ i * (n * circleG n i) := by
          rw [pow_succ]; ring
      _ ≤ n ^ i * circleG n (i + 1) :=
          Nat.mul_le_mul_left _ (circleG_mul_le_succ hn i)
  exact_mod_cast key

/-- Existence of the Circle density limit via `monotone_ratio_converges`. -/
theorem circle_density_limit_exists (n : ℕ) (hn : 2 ≤ n) :
    ∃ L, Tendsto (circleDensitySeq n hn).ratio atTop (nhds L) ∧ 0 ≤ L :=
  monotone_ratio_converges _ (circle_ratio_antitone n hn)

/-- Any limit of the Circle ratio is 0. -/
theorem circle_density_limit_eq_zero (n : ℕ) (hn : 2 ≤ n) (L : ℝ)
    (hL : Tendsto (circleDensitySeq n hn).ratio atTop (nhds L)) : L = 0 :=
  tendsto_nhds_unique hL (circle_density_vanishes n hn)

/-! ## 7. Density Bounded in [0, 1] -/

theorem ratio_bounded (s : FractalDensitySeq)
    (hFG : ∀ i, s.F i ≤ s.G i) (i : ℕ) :
    0 ≤ s.ratio i ∧ s.ratio i ≤ 1 := by
  unfold FractalDensitySeq.ratio
  exact ⟨div_nonneg (s.hF_nonneg i) (le_of_lt (s.hG_pos i)),
         div_le_one_of_le (hFG i) (le_of_lt (s.hG_pos i))⟩

/-- `ratio_bounded` is inhabited: the Square ratio lies in `[0, 1]`. -/
theorem square_ratio_mem_unit_interval (n : ℕ) (hn : 2 ≤ n) (i : ℕ) :
    0 ≤ (squareDensitySeq n hn).ratio i ∧ (squareDensitySeq n hn).ratio i ≤ 1 := by
  refine ratio_bounded _ (fun j => ?_) i
  show (↑n : ℝ) ^ j ≤ (↑n : ℝ) ^ (n ^ j)
  exact pow_le_pow_right (by exact_mod_cast (by omega : 1 ≤ n))
    (le_of_lt (Nat.lt_pow_self (by omega) j))

end RSSN.FractalDensityConvergence
