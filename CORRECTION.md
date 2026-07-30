# Correction of the appendix density formula

## 1. Setup

Let \(m\geq 1\), and let \(X\) and \(Y\) be independent with

$$
X\sim\chi_1,
\qquad
Y\sim\chi_{2m}.
$$

Define

$$
R_m=(X^4+2X^2Y^2)^{1/4}.
$$

The densities of the two chi variables are

$$
f_X(s)=\sqrt{\frac{2}{\pi}}e^{-s^2/2},
\qquad s>0,
$$

and

$$
f_Y(y)=
\frac{2^{1-m}}{\Gamma(m)}
y^{2m-1}e^{-y^2/2},
\qquad y>0.
$$

## 2. Conditional distribution

Fix \(t>0\). Conditional on \(X=s\), the event \(R_m\leq t\) is

$$
s^4+2s^2Y^2\leq t^4.
$$

It is empty for \(s>t\). For \(0<s<t\), it is equivalent to

$$
Y\leq
g_t(s):=
\sqrt{\frac12\left(\frac{t^4}{s^2}-s^2\right)}.
$$

Consequently,

$$
\Pr(R_m\leq t)
=
\int_0^t f_X(s)F_Y(g_t(s))\,ds.
$$

At the moving endpoint \(s=t\), we have \(g_t(t)=0\) and hence
\(F_Y(g_t(t))=0\). There is therefore no boundary contribution when
the expression is differentiated with respect to \(t\).

## 3. Differentiation

Direct calculation gives

$$
\frac{\partial g_t(s)}{\partial t}
=
\frac{t^3}{s^2g_t(s)}.
$$

Thus

$$
\begin{aligned}
f_{R_m}(t)
&=
\int_0^t
f_X(s)f_Y(g_t(s))
\frac{\partial g_t(s)}{\partial t}\,ds \\
&=
\frac{2^{\,2-m}t^3}{\Gamma(m)\sqrt{2\pi}}
\int_0^t
g_t(s)^{\,2m-2}
\frac{\exp\!\left(-\frac{s^2+g_t(s)^2}{2}\right)}{s^2}
\,ds.
\end{aligned}
$$

This proves the corrected formula.

## 4. What changes from the 2016 typesetting

Equation (2) in the report displays

- the prefactor \(2^{1-m}/\sqrt{2\pi}\), and
- the power \(g_t(s)^{m-2}\).

The calculation above gives, respectively,

- \(2^{2-m}/\sqrt{2\pi}\), because the \(\chi_1\) density contributes
  \(2/\sqrt{2\pi}\), and
- \(g_t(s)^{2m-2}\), because the \(\chi_{2m}\) density contributes
  \(g_t(s)^{2m-1}\) and the derivative contributes \(g_t(s)^{-1}\).

These are corrections to the final typeset expression. The report's central
conditioning idea remains the proof.

## 5. Independent radial-angular check

There is a second way to organize the same distribution. Put

$$
U=X^2,\qquad V=Y^2.
$$

Then \(U\) and \(V\) are independent gamma variables with the same scale:

$$
U\sim\mathrm{Gamma}\left(\frac12,2\right),
\qquad
V\sim\mathrm{Gamma}(m,2).
$$

The beta-gamma decomposition gives independent variables

$$
S=U+V\sim\chi^2_{2m+1},
\qquad
B=\frac{U}{U+V}
\sim\mathrm{Beta}\left(\frac12,m\right).
$$

Since

$$
R_m^4=U^2+2UV=S^2B(2-B),
$$

we obtain the exact factorization

$$
\boxed{
R_m=\chi_{2m+1}[B(2-B)]^{1/4},
}
$$

where the two factors on the right are independent. This supplies an
independent normalization check and isolates the directional dependence of
the exceptional determinant factor.

## 6. Priority and scope

Bornemann and La Croix identify the exceptional determinant factor and give
its hypergeometric Mellin transform. They state that the corresponding inverse
Mellin transform cannot readily be written down. The 2016 report instead
derives a one-dimensional density integral by conditioning on \(\chi_1\).

Accordingly, the defensible claim is narrow:

> The direct conditional-density derivation appears original to the 2016
> report relative to its principal GOE source.

This repository does not claim an exhaustive global priority search, a new
GOE factorization, or a new random-matrix universality theorem.

