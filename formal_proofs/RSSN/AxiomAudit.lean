/-
  Axiom audit — the machine-checked claim → theorem map for RSSN.

  Each `#print axioms` below reports the complete axiom footprint of one
  headline result. CI (.github/workflows/lean.yml) runs this file and fails
  if any reported axiom falls outside the allowed trust base:

    propext, Classical.choice, Quot.sound   (standard Lean/Mathlib axioms)
    Lean.ofReduceBool, Lean.trustCompiler   (compiler reflection — used only
                                             by `native_decide` in
                                             HierarchyPlacement, where the
                                             kernel would OOM on the literals)

  This makes the results falsifiable in two concrete ways:
  * an admitted proof anywhere beneath a listed theorem surfaces in the
    reported axioms as `sorryAx` and turns CI red (a plain `lake build`
    merely warns about admitted proofs);
  * any custom `axiom` smuggled into the development is listed by name and
    rejected by the CI allowlist.

  If a theorem is renamed or deleted, this file fails to elaborate, so the
  claim map cannot silently drift out of sync with the proofs.
-/
import RSSN.ShapeOperators
import RSSN.Tetration
import RSSN.FractalDensityConvergence
import RSSN.NonCommutativity
import RSSN.ReflectionIsomorphism
import RSSN.HierarchyPlacement
import RSSN.BridgeLemmas
import RSSN.UncertaintyPrinciple

-- T1: triangle fractal density converges (to 1/n)
#print axioms RSSN.FractalDensityConvergence.triangle_density_converges
-- T2: reflection is injective (structure-preserving embedding)
#print axioms RSSN.ReflectionIsomorphism.reflect_injective
-- T3: triangle sits below f_3 in the fast-growing hierarchy
--     (expected to use Lean.ofReduceBool via native_decide)
#print axioms RSSN.HierarchyPlacement.L3_3_triangle_below_f3
-- T3: hierarchy placement summary
#print axioms RSSN.HierarchyPlacement.hierarchy_summary
-- T7: shape operators do not commute
#print axioms RSSN.NonCommutativity.noncommutative_at_2
-- T7: the commutator is nonzero
#print axioms RSSN.NonCommutativity.commutator_nonzero
-- T8: Robertson uncertainty form is nontrivial for shape operators
#print axioms RSSN.UncertaintyPrinciple.T8_robertson_form_nontrivial
-- Bridge: twist field derived from fractal density
#print axioms RSSN.BridgeLemmas.bridge_twist_from_density
-- Triangle equals tetration at height 2 (links to Knuth arrows)
#print axioms RSSN.Tetration.triangle_eq_tet
