"""Numerical checks for the representation-theoretic formulas.

The proofs are in DYSON_REPRESENTATION_NOTE.md.  This checker independently
compares the beta-hypergeometric angular moments with quadrature and compares
the fixed-domain density kernel with the original corrected integral.
"""

from __future__ import annotations

import math

from verification import adaptive_simpson, corrected_density


def beta(a: float, b: float) -> float:
    return math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))


def hypergeometric_2f1(
    a: float,
    b: float,
    c: float,
    z: float,
    *,
    tolerance: float = 2e-15,
    max_terms: int = 100_000,
) -> float:
    """Evaluate 2F1 by its absolutely convergent power series at |z| < 1."""

    term = 1.0
    total = 1.0
    for n in range(max_terms):
        term *= (a + n) * (b + n) * z / ((c + n) * (n + 1))
        new_total = total + term
        if abs(term) <= tolerance * max(1.0, abs(new_total)):
            return new_total
        total = new_total
    raise RuntimeError("hypergeometric series did not converge")


def angular_moment_formula(m: int, r: float, k: int) -> float:
    """Return E[B^k (B(2-B))^r] for B ~ Beta(1/2,m)."""

    return (
        2.0**r
        * beta(0.5 + r + k, m)
        / beta(0.5, m)
        * hypergeometric_2f1(
            -r,
            0.5 + r + k,
            m + 0.5 + r + k,
            0.5,
        )
    )


def angular_moment_quadrature(m: int, r: float, k: int) -> float:
    """Independently integrate after B=sin(theta)^2 removes the singularity."""

    normalization = beta(0.5, m)

    def integrand(theta: float) -> float:
        if theta == 0.0:
            return 0.0
        if theta == math.pi / 2.0:
            return 0.0
        sine = math.sin(theta)
        cosine = math.cos(theta)
        b_value = sine * sine
        return (
            2.0
            * cosine ** (2 * m - 1)
            * b_value ** (r + k)
            * (2.0 - b_value) ** r
            / normalization
        )

    return adaptive_simpson(integrand, 0.0, math.pi / 2.0)


def fixed_domain_density(t: float, m: int) -> float:
    """Evaluate the x-integral obtained from x=s^2/t^2."""

    constant = (
        2.0 ** (2 - 2 * m)
        * t ** (2 * m)
        / (math.gamma(m) * math.sqrt(2.0 * math.pi))
    )

    def integrand(x: float) -> float:
        if x == 0.0:
            return 0.0
        if x == 1.0:
            return math.exp(-t * t / 2.0) if m == 1 else 0.0
        log_value = (
            (m - 1) * math.log1p(-x * x)
            - (m + 0.5) * math.log(x)
            - 0.25 * t * t * (x + 1.0 / x)
        )
        return math.exp(log_value)

    return constant * adaptive_simpson(integrand, 0.0, 1.0)


def main() -> None:
    largest_moment_error = 0.0
    for m in (1, 2, 4, 7):
        for r in (0.25, 0.75, 1.5):
            for k in (0, 1, 2):
                closed_form = angular_moment_formula(m, r, k)
                quadrature = angular_moment_quadrature(m, r, k)
                scale = max(1e-14, abs(closed_form), abs(quadrature))
                relative_error = abs(closed_form - quadrature) / scale
                largest_moment_error = max(largest_moment_error, relative_error)
                assert relative_error < 2e-9, (m, r, k, closed_form, quadrature)

    largest_density_error = 0.0
    for m in (1, 2, 3, 5):
        for t in (0.6, 1.0, 1.8, 3.0):
            original = corrected_density(t, m)
            transformed = fixed_domain_density(t, m)
            scale = max(1e-14, abs(original), abs(transformed))
            relative_error = abs(original - transformed) / scale
            largest_density_error = max(largest_density_error, relative_error)
            assert relative_error < 2e-8, (m, t, original, transformed)

    # The first Jacobi mode is proportional to B-E[B].  Check its exact norm.
    for m in (1, 2, 5, 10):
        mean = 1.0 / (2 * m + 1)
        variance = 4.0 * m / ((2 * m + 1) ** 2 * (2 * m + 3))
        direct_second_moment = beta(2.5, m) / beta(0.5, m)
        assert abs((direct_second_moment - mean * mean) - variance) < 2e-14

    print(
        "All representation checks passed; "
        f"max moment error={largest_moment_error:.3e}, "
        f"max density error={largest_density_error:.3e}"
    )


if __name__ == "__main__":
    main()
