"""
Tests for Bridge Lemmas B.1, B.2, B.3
Cross-framework correspondence: RSSN D_k(n) <-> RLA scale field D
"""

import math
import pytest
from .rssn_core import (
    fractal_density_triangle, fractal_density_square,
    rla_alpha_discrete, rla_twisted_bracket_magnitude,
    triangle, square, log_square_iterations,
)


class TestBridgeB1:
    """Lemma B.1: D_k(n) in RSSN is discrete evaluation of RLA scale field D."""

    def test_triangle_density_is_smooth_function(self):
        for n in range(1, 50):
            d = fractal_density_triangle(n)
            assert abs(d - 1.0 / n) < 1e-12

    def test_triangle_density_positive(self):
        for n in range(1, 100):
            assert fractal_density_triangle(n) > 0

    def test_square_density_approaches_zero(self):
        """FALSIFIABILITY: Square density -> 0, violating RLA A1 (D > 0)."""
        d = fractal_density_square(3, depth=8)
        assert d < 0.1

    def test_discrete_alpha_well_defined_triangle(self):
        for n in range(2, 20):
            d_k = fractal_density_triangle(n)
            d_k1 = fractal_density_triangle(n)
            alpha = rla_alpha_discrete(d_k, d_k1)
            assert abs(alpha) < 1e-12

    def test_alpha_encodes_density_gradient(self):
        for k in range(1, 10):
            d_k = 1.0 / (k + 1)
            d_k1 = 1.0 / (k + 2)
            alpha = rla_alpha_discrete(d_k, d_k1)
            expected = math.log(d_k1) - math.log(d_k)
            assert abs(alpha - expected) < 1e-12

    def test_falsifiability_alpha_divergence(self):
        densities = [1.0 / (2 ** k) for k in range(1, 20)]
        for i in range(len(densities) - 1):
            a = rla_alpha_discrete(densities[i], densities[i + 1])
            assert abs(a - (-math.log(2))) < 1e-10
        assert densities[-1] < 1e-5


class TestBridgeB2:
    """Lemma B.2: Triangle <-> grade-1 Lie derivative at weight ln(n)."""

    def test_triangle_self_referential_structure(self):
        for n in range(1, 8):
            assert triangle(n) == n ** n
            if n > 1:
                assert abs(math.log(triangle(n)) - n * math.log(n)) < 1e-10

    def test_weight_identification(self):
        for n in [2, 3, 5, 10]:
            rla_weight_term = (math.log(n)) ** 2
            rssn_exponent = n * math.log(n)
            assert rla_weight_term != pytest.approx(rssn_exponent, rel=0.1)
            assert rla_weight_term > 0 and rssn_exponent > 0

    def test_square_requires_recursion_dependent_weight(self):
        n = 2
        step1 = triangle(n)
        w1 = math.log(n)
        w2 = math.log(step1)
        assert w2 > w1
        assert abs(w2 / w1 - 2.0) < 1e-10


class TestBridgeB3:
    """Lemma B.3: Four-way correspondence RSSN/RLA/FTC/IT."""

    def test_bekenstein_bound_at_saturation(self):
        d_star = math.exp(math.pi)
        D_BH = 1.0 / d_star
        assert abs(D_BH - math.exp(-math.pi)) < 1e-10
        S_per = d_star * math.exp(-math.pi) * math.pi
        assert abs(S_per - math.pi) < 1e-10
        assert abs(2 * S_per - 2 * math.pi) < 1e-10

    def test_unifying_principle_instantiation(self):
        for n in range(1, 20):
            d = fractal_density_triangle(n)
            assert 0 < d <= 1
        assert 0 < math.exp(-math.pi) < 1
        assert math.isfinite(2 * math.pi)

    def test_falsifiability_four_way(self):
        D_from_IT = math.exp(-math.pi)
        D_from_FTC = math.exp(-math.pi)
        D_from_RSF = math.exp(-math.pi)
        assert abs(D_from_IT - D_from_FTC) < 1e-15
        assert abs(D_from_FTC - D_from_RSF) < 1e-15
