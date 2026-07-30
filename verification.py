"""Numerical transcription checks for the corrected GOE-factor density.

The proof is in CORRECTION.md.  This script compares the corrected density
with a finite-difference derivative of the conditional CDF.  It intentionally
uses only the Python standard library.
"""

from __future__ import annotations

import math


def adaptive_simpson(
    f,
    a: float,
    b: float,
    *,
    tolerance: float = 2e-10,
    max_depth: int = 24,
) -> float:
    """Integrate a continuous function by recursive adaptive Simpson quadrature."""

    def simpson(left: float, right: float, f_left: float, f_mid: float, f_right: float) -> float:
        return (right - left) * (f_left + 4.0 * f_mid + f_right) / 6.0

    mid = (a + b) / 2.0
    fa, fm, fb = f(a), f(mid), f(b)
    whole = simpson(a, b, fa, fm, fb)

    def recurse(
        left: float,
        right: float,
        f_left: float,
        f_mid: float,
        f_right: float,
        estimate: float,
        tol: float,
        depth: int,
    ) -> float:
        center = (left + right) / 2.0
        left_mid = (left + center) / 2.0
        right_mid = (center + right) / 2.0
        f_left_mid = f(left_mid)
        f_right_mid = f(right_mid)
        left_estimate = simpson(left, center, f_left, f_left_mid, f_mid)
        right_estimate = simpson(center, right, f_mid, f_right_mid, f_right)
        refined = left_estimate + right_estimate
        error = refined - estimate
        if depth == 0 or abs(error) <= 15.0 * tol:
            return refined + error / 15.0
        return recurse(
            left,
            center,
            f_left,
            f_left_mid,
            f_mid,
            left_estimate,
            tol / 2.0,
            depth - 1,
        ) + recurse(
            center,
            right,
            f_mid,
            f_right_mid,
            f_right,
            right_estimate,
            tol / 2.0,
            depth - 1,
        )

    return recurse(a, b, fa, fm, fb, whole, tolerance, max_depth)


def even_chi_cdf(y: float, m: int) -> float:
    """CDF of chi_(2m), evaluated by the integer-shape gamma formula."""

    if y <= 0.0:
        return 0.0
    z = y * y / 2.0
    if z > 745.0:
        return 1.0
    term = 1.0
    partial_sum = 1.0
    for j in range(1, m):
        term *= z / j
        partial_sum += term
    return -math.expm1(-z) - math.exp(-z) * (partial_sum - 1.0)


def conditional_cdf(t: float, m: int) -> float:
    """Compute P(R_m <= t) from the conditioning argument."""

    chi_one_constant = math.sqrt(2.0 / math.pi)

    def integrand(x: float) -> float:
        if x == 0.0:
            return t * chi_one_constant
        if x == 1.0:
            return 0.0
        s = t * x
        g_squared = 0.5 * t * t * (1.0 / (x * x) - x * x)
        g = math.sqrt(g_squared)
        return (
            t
            * chi_one_constant
            * math.exp(-s * s / 2.0)
            * even_chi_cdf(g, m)
        )

    return adaptive_simpson(integrand, 0.0, 1.0)


def corrected_density(t: float, m: int) -> float:
    """Evaluate the corrected density formula."""

    constant = 2.0 ** (2 - m) / (math.gamma(m) * math.sqrt(2.0 * math.pi))

    def integrand(x: float) -> float:
        if x == 0.0:
            return 0.0
        g_squared = 0.5 * t * t * (1.0 / (x * x) - x * x)
        if g_squared == 0.0:
            g_power = 1.0 if m == 1 else 0.0
        else:
            g_power = g_squared ** (m - 1)
        exponent = -0.5 * (t * t * x * x + g_squared)
        return g_power * math.exp(exponent) / (x * x)

    return constant * t * t * adaptive_simpson(integrand, 0.0, 1.0)


def finite_difference_density(t: float, m: int) -> float:
    step = 2e-4 * max(1.0, t)
    return (conditional_cdf(t + step, m) - conditional_cdf(t - step, m)) / (
        2.0 * step
    )


def main() -> None:
    largest_relative_error = 0.0
    for m in (1, 2, 3, 5):
        for t in (0.6, 1.0, 1.8, 3.0):
            exact_formula = corrected_density(t, m)
            numerical_derivative = finite_difference_density(t, m)
            scale = max(1e-10, abs(exact_formula), abs(numerical_derivative))
            relative_error = abs(exact_formula - numerical_derivative) / scale
            largest_relative_error = max(largest_relative_error, relative_error)
            assert relative_error < 2e-5, (
                m,
                t,
                exact_formula,
                numerical_derivative,
                relative_error,
            )

    print(f"All density/CDF checks passed; max relative error={largest_relative_error:.3e}")


if __name__ == "__main__":
    main()

