"""Exact checks for ORTHOGONALITY_NOTE.md.

The script uses Fraction arithmetic.  It checks the moment formula by an
independent radial-angular expansion, verifies the Wiener-chaos variance
split, and reconstructs the first monic recurrence coefficients.
"""

from __future__ import annotations

import math
from fractions import Fraction


def rising(a: Fraction, n: int) -> Fraction:
    result = Fraction(1)
    for j in range(n):
        result *= a + j
    return result


def moment_factorial_expansion(m: int, n: int) -> Fraction:
    """E[(U^2+2UV)^n] from independent gamma moments."""

    return 2 ** (2 * n) * sum(
        Fraction(math.comb(n, j) * 2**j)
        * rising(Fraction(1, 2), 2 * n - j)
        * rising(Fraction(m), j)
        for j in range(n + 1)
    )


def moment_radial_angular(m: int, n: int) -> Fraction:
    """E[rho^(4n)] E[(B(2-B))^n] from beta-gamma separation."""

    radial = 2 ** (2 * n) * rising(Fraction(2 * m + 1, 2), 2 * n)
    angular = sum(
        Fraction((-1) ** j * math.comb(n, j) * 2 ** (n - j))
        * rising(Fraction(1, 2), n + j)
        / rising(Fraction(2 * m + 1, 2), n + j)
        for j in range(n + 1)
    )
    return radial * angular


def polynomial_inner_product(
    left: list[Fraction],
    right: list[Fraction],
    moments: list[Fraction],
) -> Fraction:
    result = Fraction(0)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result += left_coefficient * right_coefficient * moments[i + j]
    return result


def first_recurrence_coefficients(m: int) -> tuple[Fraction, Fraction, Fraction]:
    moments = [moment_factorial_expansion(m, n) for n in range(4)]
    alpha_zero = moments[1]
    q_one = [-alpha_zero, Fraction(1)]
    norm_one = polynomial_inner_product(q_one, q_one, moments)
    t_q_one = [Fraction(0), *q_one]
    alpha_one = polynomial_inner_product(t_q_one, q_one, moments) / norm_one
    return alpha_zero, norm_one, alpha_one


def main() -> None:
    for m in range(1, 13):
        for n in range(9):
            assert moment_factorial_expansion(m, n) == moment_radial_angular(m, n)

        mean = Fraction(4 * m + 3)
        second_chaos = Fraction(32 * m * m + 112 * m + 72)
        fourth_chaos = Fraction(32 * m + 24)
        variance = (
            moment_factorial_expansion(m, 2)
            - moment_factorial_expansion(m, 1) ** 2
        )
        assert variance == second_chaos + fourth_chaos

        alpha_zero, beta_one, alpha_one = first_recurrence_coefficients(m)
        assert alpha_zero == mean
        assert beta_one == 16 * (2 * m * m + 9 * m + 6)
        assert alpha_one == Fraction(
            40 * m**3 + 402 * m**2 + 1035 * m + 612,
            2 * m**2 + 9 * m + 6,
        )

    # Fixed moments of T_m/(4m) converge to chi-square_1 moments.
    for n in range(1, 7):
        target = 2**n * rising(Fraction(1, 2), n)
        errors = []
        for m in (20, 80, 320):
            scaled = moment_factorial_expansion(m, n) / (4 * m) ** n
            errors.append(abs(float(scaled - target)))
        assert errors[2] < errors[1] < errors[0]

    print(
        "All orthogonality checks passed: "
        "108 moment identities, 12 chaos splits, "
        "and 36 exact recurrence coefficients."
    )


if __name__ == "__main__":
    main()
