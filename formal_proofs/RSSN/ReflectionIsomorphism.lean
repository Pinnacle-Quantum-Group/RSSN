/-
  RSSN — Reflection Isomorphism (Theorem T2)
  Pinnacle Quantum Group — April 2026

  Proves properties of the reflection mapping R that establishes
  correspondence between syntactic (operator expressions),
  semantic (recursive structures), and computational (values) levels:
  `reflect` is injective (syntactic level), and evaluation commutes
  with reflection (`eval_reflect`, `applyTo_reflect`), giving the
  commuting triangle of README Theorem 2 — including composition
  preservation R(S₁(S₂(n))) = R(S₁)(R(S₂(n))) (`reflect_comp_spec`).

  Design note (round 7 rewrite): the original `RecTree.node : String →
  List RecTree → RecTree` failed Lean 4's structural recursion checker
  for `RecTree.depth` (List recursion isn't structurally decreasing on
  the outer arg) AND blocked `deriving DecidableEq` (cant auto-derive
  through `List`). Reframed with explicit unary/binary constructors plus
  a `RecTreeTag` enum: depth becomes trivially structural, distinctness
  proofs reduce to tag inequality (decidable), and the whole module
  builds clean.

  Fix note (round 8): round 7's `eval` had a silent semantic bug — its
  `comp` clause computed `eval e₁` and then discarded it, so composition
  was a no-op on its first argument, and `eval` appeared in no theorem.
  The old definition is kept as `evalOld` with `evalOld_comp_discards`
  proving the defect; the corrected `eval` composes via substitution
  semantics (`evalWith`), and §6–§7 prove the tree-level evaluators
  agree with it along `reflect`.

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

/-! ## 5. Evaluation (Computational Level)

    `triangle'`, `square'`, `circle'` interpret the three operators,
    each level iterating the previous one diagonally (the `square'` and
    `circle'` bodies are exactly round 7's inlined `Nat.iterate`
    expressions, factored out and named). `evalWith e v` reads `e` as a
    unary operator applied to `v` — every literal is the hole — so
    `comp` genuinely composes; `eval` closes expressions off. -/

def triangle' (n : ℕ) : ℕ := n ^ n

/-- Square operator: diagonal iteration of `triangle'`. -/
def square' (n : ℕ) : ℕ := Nat.iterate triangle' n n

/-- Circle operator: diagonal iteration of `square'`. -/
def circle' (n : ℕ) : ℕ := Nat.iterate square' n n

/-- Evaluate `e` with every literal replaced by `v` — i.e. read `e` as
    a unary operator and apply it to `v`. -/
def evalWith : ShapeExpr → ℕ → ℕ
  | .lit _, v => v
  | .tri e, v => triangle' (evalWith e v)
  | .sq e, v => square' (evalWith e v)
  | .circ e, v => circle' (evalWith e v)
  | .comp e₁ e₂, v => evalWith e₁ (evalWith e₂ v)

def eval : ShapeExpr → ℕ
  | .lit n => n
  | .tri e => triangle' (eval e)
  | .sq e => square' (eval e)
  | .circ e => circle' (eval e)
  | .comp e₁ e₂ => evalWith e₁ (eval e₂)

/-- Composition now actually composes: the first operand acts on the
    value of the second. Contrast with `evalOld_comp_discards`. -/
theorem eval_comp (e₁ e₂ : ShapeExpr) :
    eval (.comp e₁ e₂) = evalWith e₁ (eval e₂) := rfl

/-! ### 5a. The round-7 defect, kept on record

    Round 7's `comp` clause read
    `eval (.lit (eval e₁)) |> fun _ => eval e₂`: it computed `eval e₁`
    and then threw the result away, so composition was semantically a
    no-op on its first argument — contradicting T2's composition-
    preservation claim. The old definition is preserved verbatim as
    `evalOld`, and `evalOld_comp_discards` proves the defect, so the
    correction is visible rather than silent. -/

def evalOld : ShapeExpr → ℕ
  | .lit n => n
  | .tri e => triangle' (evalOld e)
  | .sq e => let v := evalOld e; Nat.iterate triangle' v v
  | .circ e => let v := evalOld e; Nat.iterate (fun m => Nat.iterate triangle' m m) v v
  | .comp e₁ e₂ => evalOld (.lit (evalOld e₁)) |> fun _ => evalOld e₂

/-- Witness of the old defect: under `evalOld`, composition discarded
    its first argument entirely. -/
theorem evalOld_comp_discards (e₁ e₂ : ShapeExpr) :
    evalOld (.comp e₁ e₂) = evalOld e₂ := by
  simp [evalOld]

/-! ## 6. Semantic Evaluation on Trees

    The tree-level mirrors of `evalWith` and `eval`. Tag combinations
    never produced by `reflect` (a `unary .compose`, or a `binary`
    node with a non-`compose` tag) default to acting via the first
    child. -/

/-- Apply the operator denoted by a tree to a value (mirror of
    `evalWith`). -/
def RecTree.applyTo : RecTree → ℕ → ℕ
  | .leaf _, v => v
  | .unary .triangle c, v => triangle' (RecTree.applyTo c v)
  | .unary .square c, v => square' (RecTree.applyTo c v)
  | .unary .circle c, v => circle' (RecTree.applyTo c v)
  | .unary .compose c, v => RecTree.applyTo c v
  | .binary .compose c₁ c₂, v => RecTree.applyTo c₁ (RecTree.applyTo c₂ v)
  | .binary _ c₁ _, v => RecTree.applyTo c₁ v

/-- Evaluate a closed tree (mirror of `eval`). -/
def RecTree.eval : RecTree → ℕ
  | .leaf n => n
  | .unary .triangle c => triangle' (RecTree.eval c)
  | .unary .square c => square' (RecTree.eval c)
  | .unary .circle c => circle' (RecTree.eval c)
  | .unary .compose c => RecTree.eval c
  | .binary .compose c₁ c₂ => RecTree.applyTo c₁ (RecTree.eval c₂)
  | .binary _ c₁ _ => RecTree.eval c₁

/-! ## 7. The Commuting Triangle (Theorem T2)

    Reflection commutes with evaluation, both in operator form
    (`applyTo_reflect`) and in closed form (`eval_reflect`); the
    composition-preservation instance R(S₁(S₂(n))) = R(S₁)(R(S₂(n)))
    is `reflect_comp_spec`. Together with `reflect_injective` below,
    this is the syntactic/semantic/computational correspondence of
    README Theorem 2. -/

/-- Operator-level commutation: reflecting an expression and applying
    the resulting tree to `v` agrees with substitution semantics. -/
theorem applyTo_reflect (e : ShapeExpr) (v : ℕ) :
    RecTree.applyTo (reflect e) v = evalWith e v := by
  induction e generalizing v with
  | lit n => rfl
  | tri e ih => simp [reflect, RecTree.applyTo, evalWith, ih]
  | sq e ih => simp [reflect, RecTree.applyTo, evalWith, ih]
  | circ e ih => simp [reflect, RecTree.applyTo, evalWith, ih]
  | comp e₁ e₂ ih₁ ih₂ => simp [reflect, RecTree.applyTo, evalWith, ih₁, ih₂]

/-- Closed-form commutation: semantic evaluation after reflection
    equals computational evaluation. -/
theorem eval_reflect (e : ShapeExpr) : RecTree.eval (reflect e) = eval e := by
  induction e with
  | lit n => rfl
  | tri e ih => simp [reflect, RecTree.eval, eval, ih]
  | sq e ih => simp [reflect, RecTree.eval, eval, ih]
  | circ e ih => simp [reflect, RecTree.eval, eval, ih]
  | comp e₁ e₂ _ih₁ ih₂ =>
    simp [reflect, RecTree.eval, eval, ih₂, applyTo_reflect]

/-- README T2 composition preservation: the reflection of a composite
    acts as the composition of the reflections,
    R(S₁ ∘ S₂)(v) = R(S₁)(R(S₂)(v)). -/
theorem reflect_comp_spec (e₁ e₂ : ShapeExpr) (v : ℕ) :
    RecTree.applyTo (reflect (.comp e₁ e₂)) v
      = RecTree.applyTo (reflect e₁) (RecTree.applyTo (reflect e₂) v) := rfl

/-! ## 8. Reflection is Injective

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

/-! ## 9. Depth Correspondence -/

theorem reflect_lit_depth (n : ℕ) :
    (reflect (.lit n)).depth = ShapeExpr.depth (.lit n) := rfl

theorem reflect_preserves_leaf_structure (n : ℕ) :
    reflect (.lit n) = RecTree.leaf n := rfl

theorem reflect_tri_structure (e : ShapeExpr) :
    reflect (.tri e) = RecTree.unary .triangle (reflect e) := rfl

theorem reflect_sq_structure (e : ShapeExpr) :
    reflect (.sq e) = RecTree.unary .square (reflect e) := rfl

/-! ## 10. Distinct Operators Produce Distinct Structures -/

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
