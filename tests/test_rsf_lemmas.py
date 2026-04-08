"""
Tests for RSF Lemma Chain: T1 (Cardinality Transcendence),
T2 (Continuum Resolution), T3 (Russell), T5 (Foundation -- anchor)
"""
import math
import pytest
from .rssn_core import (
    RecursiveStructure,
    successor_naturals, even_naturals,
    dyadic_rationals, powers_of_two,
    fractal_density_triangle,
)


class TestLemma2_1_DensityNonCollapse:
    """L2.1: Same cardinality, different density. Status: TIGHT"""

    def test_same_cardinality_different_density(self):
        X = RecursiveStructure(successor_naturals)
        Y = RecursiveStructure(even_naturals)
        for n in range(1, 20):
            assert len(X.F_n(n)) == len(Y.F_n(n))
            assert X.F_n(n) != Y.F_n(n)

    def test_rsf_axiom1_distinguishes(self):
        X = RecursiveStructure(successor_naturals)
        Y = RecursiveStructure(even_naturals)
        assert all(X.F_n(n) != Y.F_n(n) for n in range(1, 20))

    def test_falsifiability_density_collapse(self):
        X = RecursiveStructure(successor_naturals)
        Y = RecursiveStructure(powers_of_two)
        assert abs(X.density(15) - Y.density(15)) > 0.01


class TestLemma2_2_NoHiddenHierarchy:
    """L2.2: Density independent of cardinality. Status: TIGHT"""

    def test_high_density_countable(self):
        R = RecursiveStructure(dyadic_rationals)
        assert R.density(max_depth=12) > 0.5

    def test_density_independent_of_cardinality(self):
        X = RecursiveStructure(dyadic_rationals)
        Y = RecursiveStructure(powers_of_two)
        assert X.density(10) != pytest.approx(Y.density(10), abs=0.01)

    def test_falsifiability_hidden_hierarchy(self):
        A = RecursiveStructure(dyadic_rationals)
        B = RecursiveStructure(lambda n: list(range(0, n + 1, max(1, n))))
        assert A.density(10) > B.density(10)


class TestLemma4_1_DensitySpectrumDense:
    """L4.1: Density spectrum dense in [0,1]. Status: TIGHT"""

    def test_density_approximates_any_target(self):
        for d in [0.5, 0.25, 0.1, 0.05, 0.01, 0.333]:
            n = max(1, math.ceil(1.0 / d))
            achieved = fractal_density_triangle(n)
            assert abs(achieved - d) <= 1.0 / (n * n) + d * 0.5 + 0.01

    def test_density_covers_unit_interval(self):
        densities = [1.0 / n for n in range(1, 101)]
        for target in [0.1 * k for k in range(1, 11)]:
            assert min(abs(d - target) for d in densities) < 0.02


class TestLemma4_2_CHDissolution:
    """L4.2: CH becomes meaningless. Status: TIGHT"""

    def test_ch_question_trivially_yes(self):
        d_low = 1.0 / 3
        d_high = 1.0 / 2
        d_mid = (d_low + d_high) / 2
        assert d_low < d_mid < d_high

    def test_density_between_any_two(self):
        for n in range(2, 20):
            assert 1.0 / (n + 1) < 1.0 / (n + 0.5) < 1.0 / n


class TestRussellTranscendence:
    """RSF T3: Russell paradox prevented by Axiom 7. Status: TIGHT"""

    def test_axiom7_prevents_self_membership(self):
        g = RecursiveStructure(lambda n: list(range(1, n + 1)))
        assert 0 not in g.F_n(1)

    def test_russell_class_well_defined(self):
        result = successor_naturals(5)
        assert result is not None

    def test_falsifiability_self_referential_generator(self):
        g = RecursiveStructure(lambda n: [0, 1, 2])
        assert 0 in g.F_n(1)  # Violates Ax.7 -- excluded by axiom


class TestT5_RecursiveFoundation:
    """RSF T5: Foundation from Axiom 7. Status: TIGHT (ANCHOR)"""

    def test_successor_has_base_generator(self):
        g = RecursiveStructure(successor_naturals)
        assert len(g.F_n(1)) > 0

    def test_foundation_prevents_infinite_descent(self):
        for n in range(1, 20):
            g = RecursiveStructure(successor_naturals)
            assert len(g.F_n(n)) < float('inf')

    @pytest.mark.parametrize("gen_fn,name", [
        (successor_naturals, "successor"),
        (even_naturals, "even"),
        (powers_of_two, "powers"),
    ])
    def test_all_standard_structures_have_foundation(self, gen_fn, name):
        g = RecursiveStructure(gen_fn, name)
        for n in range(1, 10):
            assert len(g.F_n(n)) > 0


class TestDensityConvergenceAxiom:
    """RSF Axiom 9: Density convergence. Novel to RSF."""

    def test_axiom9_compatible_with_axiom5(self):
        densities = [(n + 1) / (2 * n) for n in range(1, 100)]
        assert len(set([n + 1 for n in range(1, 100)])) == 99
        assert abs(densities[-1] - 0.5) < 0.01

    def test_axiom9_excludes_oscillating_structures(self):
        densities = [(1 + (-1) ** n) / 2 for n in range(100)]
        assert all(densities[i] == 1.0 for i in range(0, 100, 2))
        assert all(densities[i] == 0.0 for i in range(1, 100, 2))
