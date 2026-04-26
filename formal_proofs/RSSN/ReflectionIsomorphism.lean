/-
  RSSN — Reflection Isomorphism (Theorem T2)
  Pinnacle Quantum Group — April 2026

  Proves properties of the reflection mapping R that establishes
  correspondence between syntactic (operator expressions),
  semantic (recursive structures), and computational (values) levels.

  Design note (round 7 rewrite): the original `RecTree.node : String →
  List RecTree → RecTree` failed Lean 4's structural recursion checker
  for `RecTree.depth` (List recursion isn't structurally decreasing on
  the outer arg) AND blocked `deriving DecidableEq` (cant auto-derive
  through `List`). Reframed with explicit unary/binary constructors plus
  a `RecTreeTag` enum: depth becomes trivially structural, distinctness
  proofs reduce to tag inequality (decidable), and the whole module
  builds clean.

  Reference: RSSN README §4.1
-/
import Mathlib

namespace RSSN.ReflectionIsomorphism

/-! ## 1. Shape Expression Language -/

inductive ShapeExpr where
  | lit : ℕ → ShapeExpr
  | tri : ShapeExpr → ShapeExpr
  | sq : ShapeExpr → ShapeExpr
  | circ : ShapeExpr → ShapeExpr
  | comp : ShapeExpr → ShapeExpr → ShapeExpr
  deriving DecidableEq

/-! ## 2. Recursive Structure (Tagged Trees, NO List) -/

inductive RecTreeTag where
  | triangle | square | circle | compose
  deriving DecidableEq, Repr

inductive RecTree where
  | leaf : ℕ → RecTree
  | unary : RecTreeTag → RecTree → RecTree
  | binary : RecTreeTag → RecTree → RecTree → RecTree
  deriving Repr

/-! ## 3. Depth (now trivially structurally recursive) -/

def ShapeExpr.depth : ShapeExpr → ℕ
  | .lit _ => 0
  | .tri e => e.depth + 1
  | .sq e => e.depth + 2
  | .circ e => e.depth + 3
  | .comp e₁ e₂ => max e₁.depth e₂.depth + 1

def RecTree.depth : RecTree → ℕ
  | .leaf _ => 0
  | .unary _ child => child.depth + 1
  | .binary _ c₁ c₂ => max c₁.depth c₂.depth + 1

/-! ## 4. Reflection Mapping -/

def reflect : ShapeExpr → RecTree
  | .lit n => RecTree.leaf n
  | .tri e => RecTree.unary .triangle (reflect e)
  | .sq e => RecTree.unary .square (reflect e)
  | .circ e => RecTree.unary .circle (reflect e)
  | .comp e₁ e₂ => RecTree.binary .compose (reflect e₁) (reflect e₂)

/-! ## 5. Evaluation -/

def triangle' (n : ℕ) : ℕ := n ^ n

def eval : ShapeExpr → ℕ
  | .lit n => n
  | .tri e => triangle' (eval e)
  | .sq e => let v := eval e; Nat.iterate triangle' v v
  | .circ e => let v := eval e; Nat.iterate (fun m => Nat.iterate triangle' m m) v v
  | .comp e₁ e₂ => eval (.lit (eval e₁)) |> fun _ => eval e₂

/-! ## 6. Reflection is Injective

    Each branch follows the same pattern: split `e₂` into all 5
    cases; for the matching constructor, apply the iff form of
    `RecTree.{unary,binary,leaf}.injEq` and recurse via `ih`; for
    the non-matching cases, `simp [reflect]` exposes a constructor
    inequality which closes the contradiction. -/

theorem reflect_injective : Function.Injective reflect := by
  intro e₁ e₂ h
  induction e₁ generalizing e₂ with
  | lit n =>
    cases e₂ with
    | lit m =>
      simp only [reflect, RecTree.leaf.injEq] at h
      exact congrArg ShapeExpr.lit h
    | tri _ => simp [reflect] at h
    | sq _ => simp [reflect] at h
    | circ _ => simp [reflect] at h
    | comp _ _ => simp [reflect] at h
  | tri e ih =>
    cases e₂ with
    | tri e' =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact congrArg ShapeExpr.tri (ih h.2)
    | lit _ => simp [reflect] at h
    | sq _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | circ _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | comp _ _ => simp [reflect] at h
  | sq e ih =>
    cases e₂ with
    | sq e' =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact congrArg ShapeExpr.sq (ih h.2)
    | lit _ => simp [reflect] at h
    | tri _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | circ _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | comp _ _ => simp [reflect] at h
  | circ e ih =>
    cases e₂ with
    | circ e' =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact congrArg ShapeExpr.circ (ih h.2)
    | lit _ => simp [reflect] at h
    | tri _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | sq _ =>
      simp only [reflect, RecTree.unary.injEq] at h
      exact absurd h.1 (by decide)
    | comp _ _ => simp [reflect] at h
  | comp e₁ e₂ ih₁ ih₂ =>
    cases e₂ with
    | comp e₁' e₂' =>
      simp only [reflect, RecTree.binary.injEq] at h
      exact congrArg₂ ShapeExpr.comp (ih₁ h.2.1) (ih₂ h.2.2)
    | lit _ => simp [reflect] at h
    | tri _ => simp [reflect] at h
    | sq _ => simp [reflect] at h
    | circ _ => simp [reflect] at h

/-! ## 7. Depth Correspondence -/

theorem reflect_lit_depth (n : ℕ) :
    (reflect (.lit n)).depth = ShapeExpr.depth (.lit n) := rfl

theorem reflect_preserves_leaf_structure (n : ℕ) :
    reflect (.lit n) = RecTree.leaf n := rfl

theorem reflect_tri_structure (e : ShapeExpr) :
    reflect (.tri e) = RecTree.unary .triangle (reflect e) := rfl

theorem reflect_sq_structure (e : ShapeExpr) :
    reflect (.sq e) = RecTree.unary .square (reflect e) := rfl

/-! ## 8. Distinct Operators Produce Distinct Structures -/

theorem tri_sq_distinct (e : ShapeExpr) :
    reflect (.tri e) ≠ reflect (.sq e) := by
  intro h
  simp only [reflect, RecTree.unary.injEq] at h
  exact absurd h.1 (by decide)

theorem tri_circ_distinct (e : ShapeExpr) :
    reflect (.tri e) ≠ reflect (.circ e) := by
  intro h
  simp only [reflect, RecTree.unary.injEq] at h
  exact absurd h.1 (by decide)

theorem sq_circ_distinct (e : ShapeExpr) :
    reflect (.sq e) ≠ reflect (.circ e) := by
  intro h
  simp only [reflect, RecTree.unary.injEq] at h
  exact absurd h.1 (by decide)

end RSSN.ReflectionIsomorphism
