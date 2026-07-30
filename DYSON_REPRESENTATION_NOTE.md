# A rank-one representation-theoretic form of the GOE factor

## 1. The algebraic object behind the integral

Let

```math
Z=(Z_0,Z_1,\ldots,Z_{2m})
```

be a standard Gaussian vector in $\mathbb R^{2m+1}$. Put

```math
X=|Z_0|,
\qquad
Y=(Z_1^2+\cdots+Z_{2m}^2)^{1/2}.
```

The exceptional GOE determinant factor studied in this repository is

```math
R_m=(X^4+2X^2Y^2)^{1/4}.
```

This is the Gaussian pushforward of the fourth root of the polynomial

```math
q(z)=z_0^4+2z_0^2(z_1^2+\cdots+z_{2m}^2).
```

The polynomial is invariant under

```math
K=O(1)\mathbin{\times}O(2m).
```

Thus the integral is not merely a two-variable probability calculation. It is
an orbital integral for the compact rank-one symmetric space

```math
\mathbb {RP}^{2m}
=
O(2m+1)/(O(1)\mathbin{\times}O(2m)).
```

## 2. Radial-angular separation

Write $Z=\rho\omega$, where $\rho=\lVert Z\rVert$ and
$\omega$ is uniform on the unit sphere. Define

```math
B=\omega_0^2.
```

Then $\rho$ and $B$ are independent, with

```math
\rho\sim\chi_{2m+1},
\qquad
B\sim\mathrm {Beta}\left(\frac12,m\right).
```

The polynomial becomes

```math
q(\rho\omega)=\rho^4 B(2-B),
```

and hence

```math
\boxed{
R_m=\rho[B(2-B)]^{1/4}.
}
```

This is the beta-gamma factorization from the correction note, now interpreted
as the polar decomposition of a $K$-invariant polynomial.

## 3. The exact spherical transform

For a complex parameter $s$ with $\mathrm {Re}(s)\gt-2$, independence gives

```math
\mathbb E[R_m^s]
=
\mathbb E[\rho^s]\,
\mathbb E\!\left([B(2-B)]^{s/4}\right).
```

Put $r=s/4$. The radial factor is

```math
\mathbb E[\rho^s]
=
2^{s/2}
\frac{\Gamma(m+\frac12+\frac{s}{2})}
     {\Gamma(m+\frac12)}.
```

Euler's beta integral gives the angular factor

```math
\mathbb E\!\left([B(2-B)]^r\right)
=
2^r
\frac{\mathrm B(\frac12+r,m)}
     {\mathrm B(\frac12,m)}
\,{}_2F_1\!\left(
 -r,\frac12+r;
 m+\frac12+r;
 \frac12
\right).
```

Consequently,

```math
\boxed{
\begin{aligned}
\mathbb E[R_m^s]
={}&
2^{3s/4}
\frac{\Gamma(m+\frac12+\frac{s}{2})}
     {\Gamma(m+\frac12)}
\frac{\mathrm B(\frac12+\frac{s}{4},m)}
     {\mathrm B(\frac12,m)}
\\
&\times{}_2F_1\!\left(
 -\frac{s}{4},\frac12+\frac{s}{4};
 m+\frac12+\frac{s}{4};
 \frac12
\right).
\end{aligned}
}
```

The hypergeometric function is therefore the rank-one angular transform of
the polynomial $q$, rather than an accidental special function.

## 4. Jacobi modes and exact anisotropy

Under $x=2B-1$, the angular probability measure is proportional to

```math
(1-x)^{m-1}(1+x)^{-1/2}\,\mathrm dx
\qquad (-1\lt x\lt1).
```

Its orthogonal polynomials are

```math
P_k^{(m-1,-1/2)}(x).
```

These are the spherical functions for the compact rank-one quotient above.
They correspond to the even spherical harmonics of degree $2k$ on
$S^{2m}$, with Laplace eigenvalue

```math
-2k(2k+2m-1).
```

Thus the angular observable

```math
h_r(B)=[B(2-B)]^r
```

has a canonical representation-theoretic expansion

```math
h_r(B)
=
\sum_{k\geq0}
a_k(r)P_k^{(m-1,-1/2)}(2B-1).
```

The constant coefficient is the angular factor in Section 3. The first
nonconstant mode is already explicit. Let

```math
\mu=\mathbb E[B]=\frac{1}{2m+1},
\qquad
\sigma^2=\mathrm {Var}(B)
=
\frac{4m}{(2m+1)^2(2m+3)}.
```

For $k\geq0$, define

```math
M_k(r)
=
2^r
\frac{\mathrm B(\frac12+r+k,m)}
     {\mathrm B(\frac12,m)}
\,{}_2F_1\!\left(
 -r,\frac12+r+k;
 m+\frac12+r+k;
 \frac12
\right).
```

Then $M_k(r)=\mathbb E[B^kh_r(B)]$, and the orthogonal projection onto
the affine functions is

```math
h_r(B)
\longmapsto
M_0(r)+
\frac{M_1(r)-\mu M_0(r)}{\sigma^2}(B-\mu).
```

The numerator $M_1(r)-\mu M_0(r)$ is the exact first anisotropy coefficient.
Higher coefficients follow from the finite polynomial expansion of each
Jacobi polynomial and the same family $M_k(r)$.

## 5. The density as a hyperbolic orbital integral

The corrected density also has a useful fixed-domain form. Substituting
$x=s^2/t^2$ into the density integral gives

```math
f_{R_m}(t)
=
\frac{2^{2-2m}t^{2m}}{\Gamma(m)\sqrt{2\pi}}
\int_0^1
(1-x^2)^{m-1}x^{-m-1/2}
\exp\!\left[-\frac{t^2}{4}(x+x^{-1})\right]
\,\mathrm dx.
```

With $x=e^{-u}$ this becomes

```math
\boxed{
f_{R_m}(t)
=
\frac{2^{1-m}t^{2m}}{\Gamma(m)\sqrt{2\pi}}
\int_0^\infty
(\sinh u)^{m-1}e^{u/2}
\exp\!\left[-\frac{t^2}{2}\cosh u\right]
\,\mathrm du.
}
```

This is a finite-dimensional Bessel-type orbital integral. Expanding
$(\sinh u)^{m-1}$ expresses it as a finite linear combination of kernels

```math
\mathcal J_\nu(z)
=
\int_0^\infty e^{-z\cosh u+\nu u}\,\mathrm du,
```

whose symmetric part satisfies

```math
\mathcal J_\nu(z)+\mathcal J_{-\nu}(z)=2K_\nu(z).
```

This explains why modified-Bessel behavior is present while the density need
not collapse to one ordinary $K$-Bessel function.

## 6. Relation to Dyson's integral

Dyson's circular integral is the type-$A$ normalization

```math
\frac{1}{(2\pi)^n}
\int_{[0,2\pi]^n}
\prod_{i\lt j}
|e^{i\theta_i}-e^{i\theta_j}|^{2a}
\,\mathrm d\theta_1\cdots\mathrm d\theta_n
=
\frac{\Gamma(1+an)}{\Gamma(1+a)^n}.
```

For the orthogonal symmetry class, the Dyson index is $1$, so $a=1/2$.
Macdonald's root-system program generalizes this type-$A$ identity.

The exceptional factor here is not a new proof of Dyson's identity. Its
connection is structural:

- Dyson's integral is a higher-rank type-$A$ Weyl integral.
- The present factor is a rank-one projective-space orbital integral.
- Both reduce invariant integration to a power-product Jacobian.
- Here the resulting orthogonal functions are the rank-one Jacobi spherical
  functions.

The corrected integral therefore supplies a tractable rank-one test case for
the same representation-theoretic method, while retaining the exact finite
GOE determinant factor that motivated the calculation.

## 7. What is new here

The following are proved consequences of the corrected integral:

1. the explicit projective-space orbital interpretation;
2. the exact radial times angular factorization;
3. the closed Mellin transform derived directly from that factorization;
4. the complete Jacobi-mode algorithm and explicit first anisotropy
   coefficient;
5. the equivalent hyperbolic Bessel-type density integral.

No new random-matrix universality theorem or new proof of the Dyson or
Macdonald constant-term identities is claimed.

## References

1. Folkmar Bornemann and Michael La Croix,
   [*The Singular Values of the GOE*](https://arxiv.org/abs/1502.05946).
2. Freeman Dyson,
   [*Statistical Theory of the Energy Levels of Complex Systems. I*](https://sites.math.rutgers.edu/~zeilberg/akherim/dyson1962.pdf).
3. I. G. Macdonald,
   [*Some Conjectures for Root Systems*](https://doi.org/10.1137/0513070).
4. Lars Vretare,
   [*Elementary Spherical Functions on Symmetric Spaces*](https://doi.org/10.7146/math.scand.a-11667).
