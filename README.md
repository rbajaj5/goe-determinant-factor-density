# Density of an exceptional GOE determinant factor

This repository records a correction and completion of the density calculation
in the appendix of Ravi Andrew Bajaj's 2016 MIT course report,
[log-Determinant of a Wigner matrix](https://web.mit.edu/18.338/www/2016s/projects/bajaj_report.pdf).

Bornemann and La Croix discovered an independent-factor representation of the
GOE determinant. One exceptional factor in that representation leads, after
taking a square root, to

$$
R_m=\left(\xi_1^4+2\xi_1^2\xi_{2m}^2\right)^{1/4},
$$

where \(\xi_1\) and \(\xi_{2m}\) are independent chi random variables.
The 2016 appendix derived a direct one-dimensional integral for the density of
\(R_m\). The principal GOE source instead records a hypergeometric Mellin
transform and notes that its inverse Mellin transform is not readily written
down.

## Corrected result

For \(t>0\), define

$$
g_t(s)=
\sqrt{\frac12\left(\frac{t^4}{s^2}-s^2\right)}
\qquad (0<s<t).
$$

Then

$$
\boxed{
f_{R_m}(t)=
\frac{2^{\,2-m}t^3}{\Gamma(m)\sqrt{2\pi}}
\int_0^t
g_t(s)^{\,2m-2}
\frac{\exp\!\left(-\frac{s^2+g_t(s)^2}{2}\right)}{s^2}
\,ds .
}
$$

The typeset 2016 formula has \(2^{1-m}\) where \(2^{2-m}\) is
required and \(m-2\) where \(2m-2\) is required. The conditioning
argument itself gives the corrected formula directly.

See [CORRECTION.md](CORRECTION.md) for the proof, an independent
beta-gamma check, and a precise statement of the priority boundary.

## Verification

The checker uses only the Python standard library. It compares the corrected
density against a numerical derivative of the independently computed
conditional CDF:

```text
python verification.py
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

