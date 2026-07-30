# Orthogonality of the exceptional GOE factor

## 1. Why the ordinary chi-square analogy is incomplete

Let

```math
T_m=R_m^4=X^4+2X^2Y^2,
```

where $X\sim\chi_1$ and $Y\sim\chi_{2m}$ are independent. A chi-square
variable is a quadratic Gaussian invariant. Its orthogonal polynomials are
Laguerre polynomials because the radial Gaussian measure is a gamma measure.

The variable $T_m$ is instead a quartic invariant under

```math
K=O(1)\mathbin{\times}O(2m).
```

It is not radial under the full group $O(2m+1)$. That is the source of both
its anomalous density and its richer orthogonality: the radial Laguerre ladder
survives, but it is coupled to a projective-space Jacobi ladder.

## 2. The complete separated orthogonal basis

Let $Z\in\mathbb R^{2m+1}$ be standard Gaussian, write
$Z=\rho\omega$, and put $B=\omega_0^2$. For $n,k\geq0$, set

```math
\Psi_{n,k}(\rho,B)
=
\rho^{2k}
L_n^{(2k+m-1/2)}\!\left(\frac{\rho^2}{2}\right)
P_k^{(m-1,-1/2)}(2B-1).
```

Up to normalization, the functions $\Psi_{n,k}$ form an orthogonal basis of
the $K$-invariant subspace of Gaussian $L^2(\mathbb R^{2m+1})$.

The proof separates:

1. $P_k^{(m-1,-1/2)}(2B-1)$ is the unique $K$-fixed spherical harmonic
   type of degree $2k$;
2. after its solid-harmonic factor $\rho^{2k}$ is included, the radial
   measure is

   ```math
   x^{2k+m-1/2}e^{-x}\,\mathrm dx,
   \qquad x=\rho^2/2,
   ```

   whose orthogonal polynomials are
   $L_n^{(2k+m-1/2)}(x)$;
3. spherical-harmonic completeness followed by Laguerre completeness gives
   completeness of the product basis.

This is the concrete special-function model of the
$O(2m+1)$--$\mathfrak {sl}_2$ oscillator decomposition: the angular label
$k$ records the compact representation type and the radial label $n$ records
the Laguerre ladder.

## 3. A finite-band selection rule

The angular part of $T_m$ is

```math
h(B)=B(2-B).
```

In the coordinate $x=2B-1$,

```math
h(B)=\frac{-x^2+2x+3}{4}.
```

Multiplication by $x$ is tridiagonal in a Jacobi basis. Multiplication by
$h$ is therefore five-diagonal:

```math
hP_k^{(m-1,-1/2)}
\in
\mathrm {span}\{
P_{k-2}^{(m-1,-1/2)},\ldots,
P_{k+2}^{(m-1,-1/2)}
\}.
```

Grouping consecutive Jacobi modes into pairs converts this five-term
recurrence into a $2$-by-$2$ block-Jacobi recurrence. This is exactly the
higher-recurrence mechanism treated in the theory of matrix-valued
orthogonal polynomials.

The full operator is multiplication by

```math
T_m=\rho^4h(B).
```

Since it is a degree-four Gaussian polynomial, it changes total Hermite
degree by at most four. In the separated basis, it is consequently a
finite-band self-adjoint operator. The law of $T_m$ is its scalar spectral
measure at the constant vector.

This gives a practical interpretation:

> Orthogonal polynomials for the anomalous factor are obtained by applying
> the Lanczos or Gram--Schmidt procedure to a sparse representation-theoretic
> multiplication operator.

## 4. Exact Wiener-chaos decomposition

Let

```math
H_2(x)=x^2-1,
\qquad
H_4(x)=x^4-6x^2+3,
```

and write

```math
C_m=\sum_{j=1}^{2m}H_2(Z_j).
```

Then direct expansion gives the orthogonal decomposition

```math
\boxed{
\begin{aligned}
T_m={}&4m+3\\
&+(4m+6)H_2(Z_0)+2C_m\\
&+H_4(Z_0)+2H_2(Z_0)C_m.
\end{aligned}
}
```

The second line is second Gaussian chaos and the third line is fourth
Gaussian chaos. They are orthogonal. Their squared norms are

```math
\left\|(T_m-\mathbb ET_m)_{[2]}\right\|_2^2
=
32m^2+112m+72
```

and

```math
\left\|(T_m-\mathbb ET_m)_{[4]}\right\|_2^2
=
32m+24.
```

Therefore

```math
\mathbb E[T_m]=4m+3,
\qquad
\mathrm {Var}(T_m)
=
32m^2+144m+96.
```

This quantifies the anomaly. The quartic chaos is real, but its variance is
only order $m$, whereas the distinguished second-chaos mode has variance of
order $m^2$.

## 5. Exact moments and the new scalar orthogonal family

Put

```math
U=X^2\sim\mathrm {Gamma}\left(\frac12,2\right),
\qquad
V=Y^2\sim\mathrm {Gamma}(m,2).
```

Since $T_m=U^2+2UV$, every integer moment is the finite sum

```math
\boxed{
\mu_n(m):=\mathbb E[T_m^n]
=
2^{2n}
\sum_{j=0}^n
\binom nj 2^j
\left(\frac12\right)_{2n-j}(m)_j.
}
```

The first moments are

```math
\mu_0=1,
\qquad
\mu_1=4m+3,
```

and

```math
\mu_2=48m^2+168m+105.
```

Let $Q_n^{(m)}$ be the monic polynomials orthogonal for the law of $T_m$:

```math
tQ_n^{(m)}(t)
=
Q_{n+1}^{(m)}(t)
+\alpha_n(m)Q_n^{(m)}(t)
+\beta_n(m)Q_{n-1}^{(m)}(t).
```

The moment formula gives, exactly,

```math
\alpha_0(m)=4m+3,
```

```math
\beta_1(m)=16(2m^2+9m+6),
```

and

```math
\alpha_1(m)
=
\frac{40m^3+402m^2+1035m+612}
     {2m^2+9m+6}.
```

Thus

```math
Q_0^{(m)}(t)=1,
\qquad
Q_1^{(m)}(t)=t-(4m+3),
```

and

```math
Q_2^{(m)}(t)
=
(t-\alpha_1(m))(t-\alpha_0(m))-\beta_1(m).
```

These coefficients already differ from the affine-in-$n$ and
quadratic-in-$n$ patterns of the classical gamma/Laguerre family. The moment
formula nevertheless determines the entire family without numerical
integration.

## 6. Laguerre reappears as the large-dimension limit

The anomaly has a clean limiting form:

```math
\frac{T_m}{4m}
=
\frac{U^2}{4m}
+U\frac{V}{2m}.
```

Because $V/(2m)\to1$ in every fixed $L^p$ and $U$ is independent of $V$,

```math
\boxed{
\frac{T_m}{4m}\longrightarrow U=X^2
}
```

in every fixed $L^p$. Equivalently, for each fixed $n$,

```math
\frac{\mu_n(m)}{(4m)^n}
\longrightarrow
2^n\left(\frac12\right)_n,
```

the $n$th moment of $\chi_1^2$.

Let $\widehat Q_{n,m}$ be the monic orthogonal polynomial for
$T_m/(4m)$. Moment-matrix convergence then gives, coefficientwise for every
fixed $n$,

```math
\boxed{
\widehat Q_{n,m}(t)
\longrightarrow
(-2)^n n!L_n^{-1/2}(t/2).
}
```

So the finite-dimensional family is a nonclassical deformation of the
Laguerre system, but it returns exactly to the $\chi_1^2$ Laguerre family as
$m$ grows.

For the determinant factor itself this says

```math
\frac{R_m^2}{2\sqrt m}\longrightarrow |X|,
```

whereas

```math
\frac{R_m}{(4m)^{1/4}}\longrightarrow |X|^{1/2}.
```

## 7. The compact angular pushforward

The angular variable

```math
W=B(2-B)
```

has an explicit density on $(0,1)$:

```math
\mathrm d\nu_m(w)
=
\frac{
(1-w)^{(m-2)/2}
\sqrt{1+\sqrt{1-w}}
}{
2\,\mathrm B(m,\frac12)\sqrt w
}\,\mathrm dw.
```

This is a quadratic pushforward of the projective Jacobi measure. It is
positive and compactly supported, but is not itself a classical Jacobi
weight. Its orthogonal polynomials are equivalently the scalar spectral
polynomials of the five-diagonal operator $h(J_m)$, where $J_m$ is the Jacobi
matrix for the projective-space measure.

This is the precise point at which matrix-valued methods become natural:
$h(J_m)$ becomes block tridiagonal after pairing consecutive Jacobi modes.

## 8. What the literature adds

The surrounding literature supplies four established pieces:

1. rank-one compact spherical functions are Jacobi polynomials;
2. Gaussian radial-angular bases are spherical-harmonic times Laguerre bases;
3. higher odd-term recurrences can be recoded as matrix-valued three-term
   recurrences;
4. Bessel or Macdonald-function weights lead to related multiple-orthogonal
   systems.

The exact quartic $T_m$, its chaos split, its moment family, and its Laguerre
limit above are specific to the exceptional GOE factor. No source located in
the searches below states this combined orthogonality description.

## 9. Next questions

The most concrete open continuations are:

1. find a closed formula for all $\alpha_n(m)$ and $\beta_n(m)$;
2. determine whether the $Q_n^{(m)}$ satisfy a finite-order differential
   equation with polynomial coefficients;
3. compute the block-Jacobi matrix of $h(J_m)$ explicitly;
4. determine whether the hyperbolic density kernel yields a useful
   multiple-orthogonality characterization;
5. prove Plancherel or approximation estimates adapted directly to the
   $T_m$ spectral variable.

## References

1. Lars Vretare,
   [*Elementary Spherical Functions on Symmetric Spaces*](https://doi.org/10.7146/math.scand.a-11667).
2. Sundaram Thangavelu,
   [*Lectures on Hermite and Laguerre Expansions*](https://press.princeton.edu/books/paperback/9780691604213/lectures-on-hermite-and-laguerre-expansions).
3. Jürgen Prestin and Christian Wülker,
   [*Fast Fourier Transforms for Spherical Gauss--Laguerre Basis Functions*](https://arxiv.org/abs/1604.05140).
4. A. J. Durán and Walter Van Assche,
   [*Orthogonal Matrix Polynomials and Higher Order Recurrence Relations*](https://arxiv.org/abs/math/9310220).
5. Gert Heckman and Maarten van Pruijssen,
   [*Matrix Valued Orthogonal Polynomials for Gelfand Pairs of Rank One*](https://arxiv.org/abs/1310.5134).
6. Walter Van Assche and S. B. Yakubovich,
   [*Multiple Orthogonal Polynomials Associated with Macdonald Functions*](https://arxiv.org/abs/math/0101188).
