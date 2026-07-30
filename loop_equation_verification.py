"""Exact checks for the loop equations of the exceptional GOE factor.

The script uses rational arithmetic only.  It verifies finite ranges of the
two mixed-moment recurrences against the independent closed formula.
"""

from fractions import Fraction
from math import comb


def rising(x: Fraction, n: int) -> Fraction:
    """Return the rising factorial (x)_n."""
    out = Fraction(1)
    for k in range(n):
        out *= x + k
    return out


def mixed_moment(m: int, a: int, b: int, n: int) -> Fraction:
    """Return E[U^a V^b (U^2 + 2UV)^n] exactly."""
    total = Fraction(0)
    for j in range(n + 1):
        total += (
            comb(n, j)
            * 2**j
            * rising(Fraction(1, 2), a + 2 * n - j)
            * rising(Fraction(m), b + j)
        )
    return 2 ** (a + b + 2 * n) * total


def check_displayed_moments() -> None:
    for m in range(1, 13):
        assert mixed_moment(m, 0, 0, 0) == 1
        assert mixed_moment(m, 0, 0, 1) == 4 * m + 3
        assert mixed_moment(m, 0, 0, 2) == 48 * m * m + 168 * m + 105


def check_loop_hierarchy() -> int:
    checks = 0
    for m in range(1, 9):
        for a in range(5):
            for b in range(5):
                for n in range(7):
                    base = mixed_moment(m, a, b, n)

                    first = (2 * a + 1) * base
                    first -= mixed_moment(m, a + 1, b, n)
                    if n:
                        first += 4 * n * (
                            mixed_moment(m, a + 2, b, n - 1)
                            + mixed_moment(m, a + 1, b + 1, n - 1)
                        )
                    assert first == 0
                    checks += 1

                    second = (2 * m + 2 * b) * base
                    second -= mixed_moment(m, a, b + 1, n)
                    if n:
                        second += (
                            4
                            * n
                            * mixed_moment(m, a + 1, b + 1, n - 1)
                        )
                    assert second == 0
                    checks += 1
    return checks


def main() -> None:
    check_displayed_moments()
    checks = check_loop_hierarchy()
    print(f"Verified {checks} exact mixed-moment loop identities.")


if __name__ == "__main__":
    main()
