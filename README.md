# Density of an exceptional GOE determinant factor

This repository records a correction and completion of the density calculation
in the appendix of Ravi Andrew Bajaj's 2016 MIT course report,
[log-Determinant of a Wigner matrix](https://web.mit.edu/18.338/www/2016s/projects/bajaj_report.pdf).

Bornemann and La Croix discovered an independent-factor representation of the
GOE determinant. One exceptional factor in that representation leads, after
taking a square root, to

```math
R_m=\left(\xi_1^4+2\xi_1^2\xi_{2m}^2\right)^{1/4},
```

where $\xi_1$ and $\xi_{2m}$ are independent chi random variables.
The 2016 appendix derived a direct one-dimensional integral for the density of
$R_m$. The principal GOE source instead records a hypergeometric Mellin
transform and notes that its inverse Mellin transform is not readily written
down.

## Corrected result

For $t>0$, define

```math
g_t(s)=
\sqrt{\frac12\left(\frac{t^4}{s^2}-s^2\right)}
\qquad (0\lt s\lt t).
```

Then

```math
\boxed{
f_{R_m}(t)=
\frac{2^{2-m}t^3}{\Gamma(m)\sqrt{2\pi}}
\int_0^t
g_t(s)^{2m-2}
\frac{\exp\!\left(-\frac{s^2+g_t(s)^2}{2}\right)}{s^2}
\mathrm{d}s .
}
```

The typeset 2016 formula has $2^{1-m}$ where $2^{2-m}$ is
required and $m-2$ where $2m-2$ is required. The conditioning
argument itself gives the corrected formula directly.

See [CORRECTION.md](CORRECTION.md) for the proof, an independent
beta-gamma check, and a precise statement of the priority boundary.

## Representation-theoretic continuation

[DYSON_REPRESENTATION_NOTE.md](DYSON_REPRESENTATION_NOTE.md) treats the
corrected integral as an algebraic orbital integral. It identifies the compact
rank-one quotient, derives the exact Jacobi spherical expansion and its first
anisotropy coefficient, gives the Mellin transform, and rewrites the density
as a hyperbolic Bessel-type integral. It also states precisely how this
rank-one calculation relates to—but does not reprove—Dyson's type-$A$
integral.

[ORTHOGONALITY_NOTE.md](ORTHOGONALITY_NOTE.md) continues this analysis through
the oscillator representation. It gives the complete Jacobi--Laguerre basis,
an exact second/fourth Wiener-chaos decomposition, all integer moments, the
first scalar recurrence coefficients, and a coefficientwise Laguerre limit
for the resulting nonclassical orthogonal-polynomial family.

[LOOP_EQUATION_NOTE.md](LOOP_EQUATION_NOTE.md) lifts the scalar factor to its
two independent gamma coordinates and proves an exact Ward/loop hierarchy,
mixed-moment recurrences, a resolvent hierarchy, and two characterizations
of the lifted law. The note also gives the precise boundary with the
random-matrix universality theorem of Bourgade and Huang: the proof
architecture transfers, but their point-process conclusion does not apply
to a single determinant factor.

[DEFORMATION_AND_STABILITY.md](DEFORMATION_AND_STABILITY.md) places the GOE
factor in the one-parameter family \(T_{m,c}=U^2+cUV\). It proves the exact
deformed loop hierarchy and mixed moments, then gives a nonasymptotic
quadratic-Wasserstein bound for convergence of \(T_{m,c}/(2cm)\) to
\(\chi_1^2\), including the first correction to every fixed moment.

## Verification

The checkers use only the Python standard library. The first compares the
corrected density against a numerical derivative of the independently
computed conditional CDF. The second compares the beta-hypergeometric angular
moments with direct quadrature and independently checks the transformed
density kernel. The fourth checker verifies the two loop recurrences against
the closed mixed-moment formula using exact rational arithmetic:

```text
python verification.py
python representation_verification.py
python orthogonality_verification.py
python loop_equation_verification.py
python deformation_verification.py
```

The computation is a transcription check, not a substitute for the proof.

## Status

- Corrected density: proved in this repository.
- Direct density method: present in the 2016 report and apparently absent from
  the principal GOE source checked here.
- Global priority: not claimed without a broader literature review.
- Peer review: not yet undertaken.

## References

1. Ravi Andrew Bajaj, *log-Determinant of a Wigner matrix*, MIT 18.338
   project report, 2016.
2. Folkmar Bornemann and Michael La Croix,
   [*The Singular Values of the GOE*](https://arxiv.org/abs/1502.05946),
   2015.
3. Paul Bourgade and Jiaoyang Huang,
   [*Loop Equations Characterize Random Matrix Statistics*](https://arxiv.org/abs/2607.07617),
   2026.
