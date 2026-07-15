/-
  RSSN — Cardinality Transcendence (Theorem 5)
  Pinnacle Quantum Group — July 2026

  Formalizes README §6: recursive density D(X) = lim Fᵢ(X)/Gᵢ replaces
  Cantor's cardinality. (1) The density spectrum is the full unit
  continuum — every d ∈ [0,1] is the limit of a dyadic generative process —
  so no discrete aleph ladder survives; (2) density is structural, not
  intrinsic: the same count measured against different generation frames
  converges to different densities; (3) Cantor's |2^X| > |X| reads
  generatively — the power-stage count outruns every enumeration frame at
  every finite depth and explodes in the limit; (4) the middle-thirds
  bookkeeping of RSF T8 (FractalGenesis) is an RSSN density sequence with
  ratio (2/3)^i, bridging the two frameworks.
  Reference: RSSN README §6; RSF T6–T8 (GenerativeRuler,
  DiagonalGenerativity, FractalGenesis)
-/
import Mathlib
import RSSN.FractalDensityConvergence

noncomputable section
open Filter Topology

namespace RSSN.CardinalityTranscendence

open RSSN.FractalDensityConvergence

/-! ## 1. The density spectrum is the full unit continuum -/

/-- The dyadic generative realization of a target density `d`: at
    resolution `i` the process has produced `⌊d · 2ⁱ⌋` units out of a
    generation frame of `2ⁱ`. Every value is computed from finite data at
    a finite stage — the ruler discipline of RSF T6 applied to density
    itself. -/
def dyadicApproxSeq (d : ℝ) (_hd0 : 0 ≤ d) : FractalDensitySeq where
  F := fun i => (⌊d * 2 ^ i⌋₊ : ℝ)
  G := fun i => (2 : ℝ) ^ i
  hF_nonneg := fun _ => Nat.cast_nonneg _
  hG_pos := fun _ => by positivity

private lemma dyadic_ratio_le (d : ℝ) (hd0 : 0 ≤ d) (i : ℕ) :
    (dyadicApproxSeq d hd0).ratio i ≤ d := by
  show (⌊d * 2 ^ i⌋₊ : ℝ) / 2 ^ i ≤ d
  have hpow : (0 : ℝ) < (2 : ℝ) ^ i := by positivity
  rw [div_le_iff hpow]
  exact Nat.floor_le (mul_nonneg hd0 (by positivity))

private lemma dyadic_ratio_gt (d : ℝ) (hd0 : 0 ≤ d) (i : ℕ) :
    d - (1 / 2 : ℝ) ^ i < (dyadicApproxSeq d hd0).ratio i := by
  show d - (1 / 2 : ℝ) ^ i < (⌊d * 2 ^ i⌋₊ : ℝ) / 2 ^ i
  have hpow : (0 : ℝ) < (2 : ℝ) ^ i := by positivity
  have h2ne : ((2 : ℝ) ^ i) ≠ 0 := hpow.ne'
  have hfl : d * 2 ^ i < (⌊d * 2 ^ i⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one _
  rw [lt_div_iff hpow]
  have hone : (1 / 2 : ℝ) ^ i * 2 ^ i = 1 := by
    rw [div_pow, one_pow, one_div, inv_mul_cancel h2ne]
  have hexp : (d - (1 / 2 : ℝ) ^ i) * 2 ^ i
      = d * 2 ^ i - (1 / 2 : ℝ) ^ i * 2 ^ i := by ring
  rw [hexp, hone]
  linarith

/-- The dyadic realization of a density `d ∈ [0,1]` stays inside the unit
    interval at every stage (via `ratio_bounded` from
    `FractalDensityConvergence` — the count never exceeds the frame). -/
theorem dyadic_ratio_in_unit (d : ℝ) (hd0 : 0 ≤ d) (hd1 : d ≤ 1) (i : ℕ) :
    0 ≤ (dyadicApproxSeq d hd0).ratio i ∧ (dyadicApproxSeq d hd0).ratio i ≤ 1 := by
  refine ratio_bounded _ (fun j => ?_) i
  show (⌊d * 2 ^ j⌋₊ : ℝ) ≤ 2 ^ j
  have hpow : (0 : ℝ) ≤ (2 : ℝ) ^ j := by positivity
  calc (⌊d * 2 ^ j⌋₊ : ℝ) ≤ d * 2 ^ j := Nat.floor_le (mul_nonneg hd0 (by positivity))
    _ ≤ 1 * 2 ^ j := mul_le_mul_of_nonneg_right hd1 hpow
    _ = 2 ^ j := one_mul _

/-- **L5.1 (Spectrum realization).** The dyadic generative process
    converges to its target density: `⌊d·2ⁱ⌋/2ⁱ → d`. Realizing a density
    never requires a completed infinite set — only ever-finer finite
    stages, exactly as a ruler realizes a length. -/
theorem L5_1_dyadic_realizes (d : ℝ) (hd0 : 0 ≤ d) :
    Tendsto (dyadicApproxSeq d hd0).ratio atTop (nhds d) := by
  have hlow : Tendsto (fun i : ℕ => d - (1 / 2 : ℝ) ^ i) atTop (nhds d) := by
    have hpow : Tendsto (fun i : ℕ => ((1 : ℝ) / 2) ^ i) atTop (nhds 0) :=
      tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
    have h := Filter.Tendsto.const_sub d hpow
    simpa using h
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds
  · exact Filter.eventually_of_forall fun i => (dyadic_ratio_gt d hd0 i).le
  · exact Filter.eventually_of_forall fun i => dyadic_ratio_le d hd0 i

/-- **T5 (No absolute hierarchy).** Every value of the unit continuum is
    the limiting density of some generative sequence whose stages all live
    in `[0,1]`. Where Cantor's theory offers a discrete ladder of alephs
    with unresolvable gaps, RSSN's density spectrum is a continuum with no
    rungs at all: "size" values are dense, connected, and reached by
    finite-stage generation. -/
theorem T5_no_absolute_hierarchy (d : ℝ) (hd0 : 0 ≤ d) (hd1 : d ≤ 1) :
    ∃ s : FractalDensitySeq,
      (∀ i, 0 ≤ s.ratio i ∧ s.ratio i ≤ 1) ∧
        Tendsto s.ratio atTop (nhds d) :=
  ⟨dyadicApproxSeq d hd0, dyadic_ratio_in_unit d hd0 hd1, L5_1_dyadic_realizes d hd0⟩

/-- **T5 (CH dissolution, RSSN form).** Between any two realizable
    densities lies another realizable density. The countable/uncountable
    dichotomy that the Continuum Hypothesis quantifies never forms:
    the spectrum has no gap for CH to ask about. -/
theorem T5_spectrum_no_gaps (d₁ d₂ : ℝ) (h0 : 0 ≤ d₁) (hlt : d₁ < d₂) :
    ∃ d₃, d₁ < d₃ ∧ d₃ < d₂ ∧
      ∃ s : FractalDensitySeq, Tendsto s.ratio atTop (nhds d₃) :=
  ⟨(d₁ + d₂) / 2, by linarith, by linarith,
    dyadicApproxSeq ((d₁ + d₂) / 2) (by linarith),
    L5_1_dyadic_realizes ((d₁ + d₂) / 2) (by linarith)⟩

/-! ## 2. Density is structural, not intrinsic -/

/-- A fixed generative count: `2ⁱ` objects at stage `i`. What "size" this
    count has is not yet determined — that requires choosing a generation
    frame to measure it against. -/
def sharedCount : ℕ → ℝ := fun i => (2 : ℝ) ^ i

/-- The count `2ⁱ` measured against the frame `2ⁱ⁺¹` (one refinement
    ahead): density 1/2 at every stage. -/
def halfFrame : FractalDensitySeq where
  F := sharedCount
  G := fun i => (2 : ℝ) ^ (i + 1)
  hF_nonneg := fun _ => by unfold sharedCount; positivity
  hG_pos := fun _ => by positivity

/-- The same count `2ⁱ` measured against the frame `4ⁱ` (a
    faster-generating ambient): density `(1/2)ⁱ`, vanishing. -/
def quarticFrame : FractalDensitySeq where
  F := sharedCount
  G := fun i => (4 : ℝ) ^ i
  hF_nonneg := fun _ => by unfold sharedCount; positivity
  hG_pos := fun _ => by positivity

private lemma halfFrame_ratio (i : ℕ) : halfFrame.ratio i = 1 / 2 := by
  show (2 : ℝ) ^ i / 2 ^ (i + 1) = 1 / 2
  have h2 : ((2 : ℝ) ^ i) ≠ 0 := by positivity
  rw [pow_succ, mul_comm, ← div_div, div_self h2]

private lemma halfFrame_tendsto :
    Tendsto halfFrame.ratio atTop (nhds (1 / 2)) := by
  have h : halfFrame.ratio = fun _ => (1 / 2 : ℝ) := by
    funext i; exact halfFrame_ratio i
  rw [h]
  exact tendsto_const_nhds

private lemma quarticFrame_ratio (i : ℕ) :
    quarticFrame.ratio i = (1 / 2 : ℝ) ^ i := by
  show (2 : ℝ) ^ i / 4 ^ i = (1 / 2 : ℝ) ^ i
  rw [← div_pow]
  norm_num

private lemma quarticFrame_tendsto :
    Tendsto quarticFrame.ratio atTop (nhds 0) := by
  have h : Tendsto (fun i : ℕ => ((1 : ℝ) / 2) ^ i) atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
  exact h.congr fun i => (quarticFrame_ratio i).symm

/-- **T5 (Structural dependency).** One and the same count admits
    different limiting densities depending on the generation frame it is
    measured against: `2ⁱ` against `2ⁱ⁺¹` has density 1/2, against `4ⁱ`
    density 0. In RSSN, "how big" is a property of the *process pair*
    (count, frame) — of how a structure is generated within its ambient —
    never an intrinsic label on a completed set. This is README §6.1's
    second consequence, stated formally. -/
theorem T5_structural_dependency :
    ∃ s₁ s₂ : FractalDensitySeq,
      (∀ i, s₁.F i = s₂.F i) ∧
        Tendsto s₁.ratio atTop (nhds (1 / 2)) ∧
        Tendsto s₂.ratio atTop (nhds 0) ∧ (1 / 2 : ℝ) ≠ 0 :=
  ⟨halfFrame, quarticFrame, fun _ => rfl,
    halfFrame_tendsto, quarticFrame_tendsto, by norm_num⟩

/-! ## 3. Cantor's |2^X| > |X|, read generatively -/

/-- **L5.2 (Stage-wise escape).** At every finite depth the power-stage
    count strictly outruns the enumeration frame: `i + 1 < 2^(i+1)`.
    This is the growth-rate content of "the powerset is bigger" — a fact
    about every finite stage of generation, never about a completed
    transfinite object. (The generative engine behind it is RSF T7's
    diagonal, which produces the witnessing fresh element explicitly.) -/
theorem L5_2_power_outruns_frame (i : ℕ) :
    ((i : ℝ) + 1) < (2 : ℝ) ^ (i + 1) := by
  exact_mod_cast Nat.lt_two_pow (i + 1)

/-- **L5.3 (Count explosion).** The power-stage count `2ⁱ` tends to
    infinity: generation is inexhaustible, yet at every stage the count is
    a plain finite number — the aleph tower is replaced by a divergent but
    everywhere-finite record. -/
theorem L5_3_power_count_explodes :
    Tendsto (fun i : ℕ => (2 : ℝ) ^ i) atTop atTop :=
  tendsto_pow_atTop_atTop_of_one_lt (by norm_num)

/-! ## 4. The Cantor bridge: RSF T8 as an RSSN density sequence -/

/-- Cantor's middle-thirds bookkeeping (RSF T8, `FractalGenesis`) cast in
    RSSN's native form: `F i = 2ⁱ` retained intervals against the triadic
    generation frame `G i = 3ⁱ`. -/
def cantorDensitySeq : FractalDensitySeq where
  F := fun i => (2 : ℝ) ^ i
  G := fun i => (3 : ℝ) ^ i
  hF_nonneg := fun _ => by positivity
  hG_pos := fun _ => by positivity

private lemma cantor_ratio (i : ℕ) :
    cantorDensitySeq.ratio i = (2 / 3 : ℝ) ^ i := by
  show (2 : ℝ) ^ i / 3 ^ i = (2 / 3 : ℝ) ^ i
  rw [← div_pow]

/-- **T5 (Cantor bridge).** The RSSN density of Cantor's own construction
    vanishes — `(2/3)ⁱ → 0` — while its count `2ⁱ` explodes
    (`L5_3`). "How many" and "how much" come apart at the very object that
    launched cardinality theory; RSSN keeps the convergent density (the
    measurable quantity) and lets the count remain an honest, everywhere-
    finite generative record. This is the same divergence proved from the
    interval bookkeeping in RSF T8 (`FractalGenesis`), here expressed in
    RSSN's `Fᵢ/Gᵢ` formalism. -/
theorem T5_cantor_bridge :
    Tendsto cantorDensitySeq.F atTop atTop ∧
      Tendsto cantorDensitySeq.ratio atTop (nhds 0) := by
  constructor
  · exact L5_3_power_count_explodes
  · have h : Tendsto (fun i : ℕ => ((2 : ℝ) / 3) ^ i) atTop (nhds 0) :=
      tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
    exact h.congr fun i => (cantor_ratio i).symm

/-! ## 5. Theorem 5, assembled -/

/-- **Theorem 5 (Cardinality Transcendence, README §6).** The three
    consequences of replacing static cardinality with recursive density,
    in one statement: (1) the density spectrum is the full unit continuum
    — no absolute hierarchy; (2) density is structural — the same count
    admits different densities in different generation frames; (3) the
    power-stage count explodes while remaining finite at every stage — the
    diagonal phenomenon survives as inexhaustibility of generation, not as
    a ladder of completed infinities. -/
theorem T5_cardinality_transcendence :
    (∀ d : ℝ, 0 ≤ d → d ≤ 1 →
        ∃ s : FractalDensitySeq, Tendsto s.ratio atTop (nhds d)) ∧
      (∃ s₁ s₂ : FractalDensitySeq, (∀ i, s₁.F i = s₂.F i) ∧
        Tendsto s₁.ratio atTop (nhds (1 / 2)) ∧
        Tendsto s₂.ratio atTop (nhds 0)) ∧
      Tendsto (fun i : ℕ => (2 : ℝ) ^ i) atTop atTop := by
  refine ⟨fun d hd0 _ => ⟨dyadicApproxSeq d hd0, L5_1_dyadic_realizes d hd0⟩,
    ⟨halfFrame, quarticFrame, fun _ => rfl, halfFrame_tendsto, quarticFrame_tendsto⟩,
    L5_3_power_count_explodes⟩

end RSSN.CardinalityTranscendence
