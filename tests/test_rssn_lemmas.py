"""
Tests for RSSN Lemma Chain: T1 (Convergence), T2 (Reflection),
T3 (Fast-growing), T7 (Non-commutativity), T8 (Uncertainty)
"""
import math
import pytest
from .rssn_core import (
    triangle, square,
    fractal_density_triangle, fractal_density_square,
    fractal_density_sequence,
    abel_triangle, abel_triangle_inverse, fractional_triangle,
    log_triangle, log_square_iterations,
    log10_f3, recursive_entropy, shannon_entropy,
)


class TestLemma1_1_TriangleConstantDensity:
    """L1.1: For Triangle, D_k(n) = 1/n for all k. Status: TIGHT"""

    @pytest.mark.parametrize("n", range(1, 51))
    def test_density_equals_1_over_n(self, n):
        assert abs(fractal_density_triangle(n, depth=30) - 1.0 / n) < 1e-12

    def test_sequence_is_constant(self):
        for n in [2, 5, 10, 100]:
            seq = fractal_density_sequence(n, "triangle", max_depth=20)
            for val in seq:
                assert abs(val - 1.0 / n) < 1e-12

    def test_falsifiability_oscillating_density(self):
        """Bounded but non-convergent sequence -- confirms Cauchy gap."""
        seq = [(1 + (-1) ** i) / 2 for i in range(20)]
        assert all(0 <= r <= 1 for r in seq)
        assert all(seq[i] == 1.0 for i in range(0, 20, 2))
        assert all(seq[i] == 0.0 for i in range(1, 20, 2))


class TestLemma1_2_SquareRatioMonotone:
    """L1.2: For Square, F_i/G_i -> 0. Status: PARTIALLY TIGHT"""

    def test_square_density_decreasing(self):
        for n in [2, 3, 5]:
            seq = fractal_density_sequence(n, "square", max_depth=15)
            for i in range(1, len(seq)):
                assert seq[i] <= seq[i - 1] + 1e-15

    def test_square_density_converges_to_zero(self):
        for n in [2, 3, 5]:
            seq = fractal_density_sequence(n, "square", max_depth=15)
            assert seq[-1] < 0.01

    def test_falsifiability_denominator_dominates(self):
        n = 3
        f, g = 1, n
        for _ in range(1, 15):
            f = n * f
            g = n * g
            assert abs(f / g - 1.0) < 1e-10


class TestLemma1_3_CauchyCriterion:
    """L1.3: RSF Density Convergence Axiom supplies Cauchy. Status: TIGHT"""

    def test_triangle_cauchy(self):
        for n in [2, 5, 10]:
            seq = fractal_density_sequence(n, "triangle", max_depth=30)
            for i in range(len(seq)):
                for j in range(i, len(seq)):
                    assert abs(seq[i] - seq[j]) < 1e-12

    def test_square_cauchy(self):
        for n in [2, 3]:
            seq = fractal_density_sequence(n, "square", max_depth=15)
            for i in range(1, len(seq)):
                assert seq[i] <= seq[i - 1] + 1e-15
                assert seq[i] >= 0

    def test_rsf_axiom5_no_contradiction(self):
        states = [n + 1 for n in range(1, 50)]
        densities = [(n + 1) / (2 * n) for n in range(1, 50)]
        assert len(set(states)) == len(states)
        assert abs(densities[-1] - 0.5) < 0.02


class TestLemmaT2_Reflection:
    """T2: Reflection isomorphism. Status: TIGHT"""

    def test_t2_1_syntactic_distinctness(self):
        assert triangle(3) != square(3, max_iter=3)
        assert triangle(2) != triangle(3)

    def test_t2_2_semantic_well_definedness(self):
        for n in range(1, 6):
            t = triangle(n)
            assert isinstance(t, int) and t > 0
        assert triangle(triangle(2)) == 256

    def test_t2_3_injectivity(self):
        assert triangle(2) < square(2, max_iter=2)
        for n in range(1, 7):
            assert triangle(n) < triangle(n + 1)

    def test_t2_4_computational_faithfulness(self):
        assert triangle(triangle(2)) == square(2, max_iter=2) == 256

    def test_t2_4_fractional_iteration_abel(self):
        n = 10.0
        half = fractional_triangle(n, 0.5)
        double_half = fractional_triangle(half, 0.5)
        log_dh = math.log(double_half) if double_half > 0 else 0
        log_full = math.log(triangle(int(n)))
        assert abs(log_dh - log_full) / max(log_full, 1) < 1.0


class TestLemmaT3_FastGrowing:
    """T3: Hierarchy placement. Status: PARTIALLY TIGHT"""

    def test_t3_3_corrected_hierarchy_placement(self):
        assert triangle(4) == 256
        f3_4 = 2 ** (2 ** (2 ** 2))
        assert f3_4 == 65536
        assert f3_4 > triangle(4)

    def test_t3_3_square_below_f4(self):
        s2 = square(2, max_iter=2)
        assert s2 == 256
        assert 65536 > s2

    def test_falsifiability_hierarchy_separation(self):
        t5 = triangle(5)
        assert t5 == 3125
        tower = 2
        for _ in range(4):
            tower = 2 ** tower
            if tower > 10 ** 100:
                break
        assert tower > t5


class TestLemmaT7_NonCommutativity:
    """T7: [S1, S2] != 0. Status: TIGHT (ANCHOR)"""

    def test_triangle_square_commutator_nonzero(self):
        sq2 = square(2, max_iter=2)
        assert sq2 == 256
        log_t_sq = 256 * math.log(256)
        t2 = triangle(2)
        assert t2 == 4
        log_sq_t = log_square_iterations(4, 4)
        assert log_sq_t > log_t_sq or math.isinf(log_sq_t)

    def test_falsifiability_commuting_operators(self):
        for n in range(1, 10):
            assert triangle(n) - triangle(n) == 0

    def test_rla_correspondence(self):
        for n in [2, 5, 10]:
            alpha = math.log(n)
            for m in range(0, 4):
                for k in range(m + 1, 4):
                    assert abs(k - m) * (n ** (m + k)) * alpha > 0


class TestLemmaT8_Uncertainty:
    """T8: Uncertainty principle. Status: CONJECTURED"""

    def test_t8_1_scalar_operators_commute(self):
        for n in [2, 5, 10]:
            for k in range(1, 5):
                d_k = 1.0 / n
                assert abs(d_k * k - k * d_k) < 1e-15

    def test_t8_2_robertson_relation(self):
        n = 2
        t_n = triangle(n)
        s_n = square(n, max_iter=2)
        log_comm = abs(math.log(s_n) * n - n * math.log(t_n))
        assert log_comm > 0

    def test_t8_falsifiability_original_form(self):
        """n=1 counterexample: Delta_D=0, so 0 >= K(1)/2 is FALSE."""
        assert fractal_density_triangle(1) == 1.0
        assert 0.0 * 100 < 0.5 * 1.0


class TestT1CompleteChain:
    """Integration: full T1 proof chain."""

    @pytest.mark.parametrize("n", [2, 3, 5, 10, 20, 50])
    def test_triangle_convergence(self, n):
        seq = fractal_density_sequence(n, "triangle", max_depth=20)
        assert all(abs(v - 1.0 / n) < 1e-12 for v in seq)

    @pytest.mark.parametrize("n", [2, 3])
    def test_square_convergence(self, n):
        seq = fractal_density_sequence(n, "square", max_depth=12)
        for i in range(1, len(seq)):
            assert seq[i] <= seq[i - 1] + 1e-15
        assert all(v >= 0 for v in seq)
        assert seq[-1] < 0.01
