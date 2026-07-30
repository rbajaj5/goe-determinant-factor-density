"""Exact checks for DEFORMATION_AND_STABILITY.md.

The script uses Fraction arithmetic.  It checks the deformed mixed-moment
formula against both Ward recurrences, the algebraic closure, the exact
mean-square coupling formula, and the first two coefficients of the
large-m moment expansion.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb


def rising(x: Fraction, n: int) -> Fraction:
    result = Fraction(1)
    for k in range(n):
        result *= x + k
    return result


def mixed_moment(
    m: int,
    c: Fraction,
    a: int,
    b: int,
    n: int,
) -> Fraction:
    total = Fraction(0)
    for j in range(n + 1):
        total += (
            comb(n, j)
            * c**j
            * rising(Fraction(1, 2), a + 2 * n - j)
            * rising(Fraction(m), b + j)
        )
    return 2 ** (a + b + 2 * n) * total


def rising_m_polynomial(n: int) -> list[int]:
    """Coefficients of (m)_n in ascending powers of m."""
    coefficients = [1]
    for shift in range(n):
        next_coefficients = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            next_coefficients[degree] += shift * coefficient
            next_coefficients[degree + 1] += coefficient
        coefficients = next_coefficients
    return coefficients


def scaled_moment_coefficient(n: int, c: Fraction, power: int) -> Fraction:
    """Coefficient of m^(-power) in E[(T/(2cm))^n]."""
    target_degree = n - power
    if target_degree < 0:
        return Fraction(0)

    total = Fraction(0)
    for j in range(n + 1):
        polynomial = rising_m_polynomial(j)
        if target_degree < len(polynomial):
            total += (
                comb(n, j)
                * c**j
                * rising(Fraction(1, 2), 2 * n - j)
                * polynomial[target_degree]
            )
    return Fraction(2**n, 1) * total / c**n


def check_deformed_hierarchy() -> tuple[int, int]:
    differential_checks = 0
    algebraic_checks = 0
    parameters = (
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
        Fraction(3),
        Fraction(5),
    )
    for c in parameters:
        for m in range(1, 7):
            for a in range(4):
                for b in range(4):
                    for n in range(6):
                        base = mixed_moment(m, c, a, b, n)

                        first = (2 * a + 1) * base
                        first -= mixed_moment(m, c, a + 1, b, n)
                        if n:
                            first += n * (
                                4 * mixed_moment(m, c, a + 2, b, n - 1)
                                + 2
                                * c
                                * mixed_moment(m, c, a + 1, b + 1, n - 1)
                            )
                        assert first == 0
                        differential_checks += 1

                        second = (2 * m + 2 * b) * base
                        second -= mixed_moment(m, c, a, b + 1, n)
                        if n:
                            second += (
                                2
                                * c
                                * n
                                * mixed_moment(m, c, a + 1, b + 1, n - 1)
                            )
                        assert second == 0
                        differential_checks += 1

                        algebraic = mixed_moment(m, c, a + 2, b, n)
                        algebraic += c * mixed_moment(m, c, a + 1, b + 1, n)
                        algebraic -= mixed_moment(m, c, a, b, n + 1)
                        assert algebraic == 0
                        algebraic_checks += 1
    return differential_checks, algebraic_checks


def check_quantitative_limit() -> int:
    checks = 0
    parameters = (
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
        Fraction(3),
        Fraction(5),
    )
    for c in parameters:
        for m in range(1, 21):
            mean_w = mixed_moment(m, c, 0, 0, 1) / (2 * c * m)
            mean_difference = mean_w - 1
            assert mean_difference == Fraction(3, 2 * c * m)

            second_w = mixed_moment(m, c, 0, 0, 2) / (2 * c * m) ** 2
            cross_w_u = mixed_moment(m, c, 1, 0, 1) / (2 * c * m)
            mean_square_difference = second_w - 2 * cross_w_u + 3
            assert mean_square_difference == (
                Fraction(3, m) + Fraction(105, 4 * c * c * m * m)
            )
            checks += 2

    for c in parameters:
        for n in range(1, 11):
            leading = 2**n * rising(Fraction(1, 2), n)
            first_correction = (
                2 ** (n - 1)
                * n
                * rising(Fraction(1, 2), n)
                * (
                    Fraction(n - 1)
                    + Fraction(2, 1) * Fraction(2 * n + 1, 2) / c
                )
            )
            assert scaled_moment_coefficient(n, c, 0) == leading
            assert scaled_moment_coefficient(n, c, 1) == first_correction
            checks += 2
    return checks


def main() -> None:
    differential_checks, algebraic_checks = check_deformed_hierarchy()
    limit_checks = check_quantitative_limit()
    total = differential_checks + algebraic_checks + limit_checks
    print(
        f"Verified {total} exact deformation identities: "
        f"{differential_checks} differential, "
        f"{algebraic_checks} algebraic, and "
        f"{limit_checks} quantitative-limit identities."
    )


if __name__ == "__main__":
    main()
