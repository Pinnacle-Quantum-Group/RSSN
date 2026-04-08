"""
Core RSSN mathematical primitives used across all test modules.

Implements:
- Shape operators (Triangle, Square, Circle)
- Fractal density functions D_k(n)
- Abel functions for fractional iteration
- Fast-growing hierarchy functions f_alpha
- RSF recursive structure primitives
- RLA scale field discretization
"""

import math
from functools import lru_cache
from typing import Callable, Optional, Tuple, List


# ---------------------------------------------------------------------------
# Shape Operators
# ---------------------------------------------------------------------------

def triangle(n: int) -> int:
    """Triangle(n) = n^n"""
    if n <= 0:
        raise ValueError("RSSN operators require positive integer input")
    return n ** n


def square(n: int, max_iter: Optional[int] = None) -> int:
    """Square(n) = Triangle^n(n). Capped to prevent overflow."""
    if n <= 0:
        raise ValueError("RSSN operators require positive integer input")
    iterations = min(n, max_iter) if max_iter else n
    result = n
    for _ in range(iterations):
        result = triangle(result)
    return result


def circle(n: int, max_iter: Optional[int] = None) -> int:
    """Circle(n) = Square^n(n). Capped to prevent overflow."""
    if n <= 0:
        raise ValueError("RSSN operators require positive integer input")
    iterations = min(n, max_iter) if max_iter else n
    result = n
    for _ in range(iterations):
        result = square(result, max_iter=max_iter)
    return result


# ---------------------------------------------------------------------------
# Log-scale Shape Operators (for large values)
# ---------------------------------------------------------------------------

def log_triangle(log_n: float) -> float:
    """log(Triangle(n)) = n * log(n), given log(n) as input."""
    n = math.exp(log_n)
    return n * log_n


def log_square_iterations(n: int, depth: int) -> float:
    """Compute log(Triangle^depth(n)) iteratively."""
    log_val = math.log(n)
    for _ in range(depth):
        val = math.exp(log_val) if log_val < 700 else float('inf')
        if math.isinf(val):
            return float('inf')
        log_val = val * log_val
    return log_val


# ---------------------------------------------------------------------------
# Fractal Density Functions
# ---------------------------------------------------------------------------

def fractal_density_triangle(n: int, depth: int = 20) -> float:
    """
    D_k(n) for Triangle operator.
    F_i(n) = n * F_{i-1}(n), G_i = n * G_{i-1}
    => F_i/G_i = F_0/G_0 = 1/n for all i.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    f = [0] * (depth + 1)
    g = [0] * (depth + 1)
    f[0] = 1
    g[0] = n
    for i in range(1, depth + 1):
        f[i] = n * f[i - 1]
        g[i] = n * g[i - 1]
    return f[depth] / g[depth]


def fractal_density_square(n: int, depth: int = 10) -> float:
    """
    D_k(n) for Square operator.
    Uses log-scale to avoid overflow.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    log_g = [0.0] * (depth + 1)
    log_g[0] = math.log(n)
    for i in range(1, depth + 1):
        log_g[i] = n * log_g[i - 1]
    log_f = [0.0] * (depth + 1)
    log_f[0] = 0.0
    for i in range(1, depth + 1):
        log_f[i] = log_f[i - 1] + math.log(n)
    ratios = [log_f[i] - log_g[i] for i in range(depth + 1)]
    final_log_ratio = ratios[depth]
    if final_log_ratio < -700:
        return 0.0
    return math.exp(final_log_ratio)


def fractal_density_sequence(n: int, operator: str = "triangle",
                              max_depth: int = 20) -> List[float]:
    """Return the sequence {F_i(n)/G_i} for convergence analysis."""
    if operator == "triangle":
        f, g = 1, n
        seq = [f / g]
        for _ in range(max_depth):
            f = n * f
            g = n * g
            seq.append(f / g)
        return seq
    elif operator == "square":
        ratio = 1.0 / n
        seq = [ratio]
        for i in range(max_depth):
            ratio = ratio / n
            seq.append(max(ratio, 0.0))
        return seq
    else:
        raise ValueError(f"Unknown operator: {operator}")


# ---------------------------------------------------------------------------
# Fast-Growing Hierarchy
# ---------------------------------------------------------------------------

def f_alpha(alpha: int, n: int, max_depth: int = 50) -> float:
    """Standard fast-growing hierarchy."""
    if alpha == 0:
        return n + 1
    if alpha == 1:
        return 2 * n
    if alpha == 2:
        return n * (2 ** n)
    if alpha == 3:
        return float('inf')
    return float('inf')


def log10_f3(n: int) -> float:
    """log10(f_3(n)) -- tower of 2s of height n."""
    if n <= 0:
        return 0
    val = 1.0
    for _ in range(n):
        val = 2 ** val
        if val > 1e308:
            return float('inf')
    return math.log10(val) if val > 0 and not math.isinf(val) else float('inf')


# ---------------------------------------------------------------------------
# Abel Functions (for fractional iteration)
# ---------------------------------------------------------------------------

def abel_triangle(n: float) -> float:
    """Abel function for Triangle: A(Triangle(n)) = A(n) + 1."""
    if n <= 1:
        return 0.0
    log_n = math.log(n)
    if log_n <= 1:
        return log_n
    log_log_n = math.log(log_n)
    if log_log_n <= 0:
        return log_n
    return log_n / log_log_n


def abel_triangle_inverse(a: float) -> float:
    """Inverse Abel function (approximate via Newton's method)."""
    n = math.exp(a)
    for _ in range(20):
        val = abel_triangle(n)
        if abs(val - a) < 1e-10:
            break
        h = max(n * 1e-8, 1e-10)
        deriv = (abel_triangle(n + h) - abel_triangle(n - h)) / (2 * h)
        if abs(deriv) < 1e-15:
            break
        n = n - (val - a) / deriv
        n = max(n, 1.01)
    return n


def fractional_triangle(n: float, alpha: float) -> float:
    """Triangle^alpha(n) via Abel function."""
    a = abel_triangle(n)
    return abel_triangle_inverse(a + alpha)


# ---------------------------------------------------------------------------
# RLA Scale Field Discretization
# ---------------------------------------------------------------------------

def rla_alpha_discrete(d_k: float, d_k_plus_1: float) -> float:
    """Discrete 1-form alpha: alpha|_{k->k+1} = ln D_{k+1} - ln D_k"""
    if d_k <= 0 or d_k_plus_1 <= 0:
        return float('nan')
    return math.log(d_k_plus_1) - math.log(d_k)


def rla_twisted_bracket_magnitude(alpha: float, m: int, n: int) -> float:
    """Magnitude of the alpha-twist in the RLA bracket."""
    return abs(n - m) * abs(alpha)


# ---------------------------------------------------------------------------
# RSF Recursive Structure Primitives
# ---------------------------------------------------------------------------

class RecursiveStructure:
    """A recursive structure in RSF, defined by its generator and depth."""

    def __init__(self, generator: Callable[[int], list], name: str = ""):
        self.generator = generator
        self.name = name

    def F_n(self, n: int) -> list:
        return self.generator(n)

    def density(self, max_depth: int = 20) -> float:
        sizes = []
        for n in range(1, max_depth + 1):
            fn = len(self.F_n(n))
            gn = n
            sizes.append(fn / gn if gn > 0 else 0)
        return sizes[-1] if sizes else 0.0

    def density_sequence(self, max_depth: int = 20) -> List[float]:
        seq = []
        for n in range(1, max_depth + 1):
            fn = len(self.F_n(n))
            gn = n
            seq.append(fn / gn if gn > 0 else 0)
        return seq


def successor_naturals(n: int) -> list:
    return list(range(n + 1))

def even_naturals(n: int) -> list:
    return list(range(0, 2 * n + 1, 2))

def dyadic_rationals(n: int) -> list:
    denom = 2 ** min(n, 20)
    return [k / denom for k in range(denom + 1)]

def powers_of_two(n: int) -> list:
    return [2 ** k for k in range(n + 1)]


# ---------------------------------------------------------------------------
# FTC Ricci Flow Construction
# ---------------------------------------------------------------------------

def ricci_flow_metric(g0: float, R_exact: float, t: float) -> float:
    """Ricci flow: g(t) = g0 - 2*R*t + O(t^2)"""
    return g0 - 2 * R_exact * t


def ftc_recursive_density(g0: float, R_exact: float, n: int) -> float:
    """D_n(g) under Ricci flow construction. F_n(g) = g(1/n^2), G_n = g0."""
    if n == 0:
        return 1.0
    t = 1.0 / (n * n)
    g_at_t = ricci_flow_metric(g0, R_exact, t)
    return g_at_t / g0


def ftc_recursive_ricci(g0: float, R_exact: float, n: int) -> float:
    """Recursive Ricci tensor under Ricci flow construction."""
    if n <= 0:
        n = 1
    d_n = ftc_recursive_density(g0, R_exact, n)
    d_n1 = ftc_recursive_density(g0, R_exact, n + 1)
    return (d_n1 - d_n) * g0 * n * n


# ---------------------------------------------------------------------------
# Information-Theoretic Quantities
# ---------------------------------------------------------------------------

def recursive_entropy(density_sequence: List[float]) -> float:
    """S_recursive = sum_n D_n(x) * log(1/D_n(x))"""
    s = 0.0
    for d in density_sequence:
        if d > 0 and d <= 1:
            s += d * math.log(1.0 / d)
    return s


def shannon_entropy(probabilities: List[float]) -> float:
    """Standard Shannon entropy: -sum p_i log p_i"""
    s = 0.0
    for p in probabilities:
        if p > 0:
            s -= p * math.log(p)
    return s


E_PI = math.exp(math.pi)
E_NEG_PI = math.exp(-math.pi)
I_LOCAL = 2 * math.pi
