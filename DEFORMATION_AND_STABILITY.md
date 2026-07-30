# A gamma deformation and a quantitative stability theorem

## 1. The deformation

Let

```math
U\sim\mathrm{Gamma}\left(\frac12,2\right),
\qquad
V\sim\mathrm{Gamma}(m,2),
\qquad
U\ \text{and}\ V\ \text{independent},
```

where \(m>0\). For a parameter \(c>0\), define

```math
T_{m,c}=U^2+cUV,
\qquad
W_{m,c}=\frac{T_{m,c}}{2cm}.
```

The exceptional GOE determinant factor in this repository is the case
\(c=2\):

```math
T_{m,2}=R_m^4,
\qquad
W_{m,2}=\frac{R_m^4}{4m}.
```

This one-parameter family separates the features that come from the two
gamma coordinates from the numerical coefficient \(2\) in the GOE factor.

## 2. Exact mixed moments

For nonnegative integers \(a,b,n\), put

```math
M^{(c)}_{a,b,n}
=
\mathbb E\left[U^aV^bT_{m,c}^n\right].
```

### Theorem 1

```math
\boxed{
M^{(c)}_{a,b,n}
=
2^{a+b+2n}
\sum_{j=0}^{n}
\binom nj c^j
\left(\frac12\right)_{a+2n-j}
(m)_{b+j}.
}
```

Here \((x)_k=x(x+1)\cdots(x+k-1)\) is the rising factorial.

#### Proof

The binomial theorem gives

```math
T_{m,c}^n
=
\sum_{j=0}^{n}
\binom nj c^j U^{2n-j}V^j.
```

Independence and the gamma moment formulas

```math
\mathbb E[U^r]
=
2^r\left(\frac12\right)_r,
\qquad
\mathbb E[V^s]
=
2^s(m)_s
```

give the result term by term. \(\square\)

## 3. The deformed loop hierarchy

The two gamma integration-by-parts identities from
[LOOP_EQUATION_NOTE.md](LOOP_EQUATION_NOTE.md) do not depend on \(c\).
Only the derivatives of the pushforward map change:

```math
\partial_u T_{m,c}=2u+cv,
\qquad
\partial_v T_{m,c}=cu.
```

Consequently, for every differentiable \(f\) for which the expectations
exist,

```math
\boxed{
\begin{aligned}
0={}&
\mathbb E\left[
\big((2a+1)U^aV^b-U^{a+1}V^b\big)f(T_{m,c})
\right.\\
&\left.\hspace{18mm}
+\big(4U^{a+2}V^b+2cU^{a+1}V^{b+1}\big)
f'(T_{m,c})
\right],
\end{aligned}
}
```

and

```math
\boxed{
\begin{aligned}
0={}&
\mathbb E\left[
\big((2m+2b)U^aV^b-U^aV^{b+1}\big)f(T_{m,c})
\right.\\
&\left.\hspace{30mm}
+2cU^{a+1}V^{b+1}f'(T_{m,c})
\right].
\end{aligned}
}
```

Taking \(f(t)=t^n\) yields the exact recurrences

```math
\boxed{
\begin{aligned}
M^{(c)}_{a+1,b,n}
={}&(2a+1)M^{(c)}_{a,b,n}\\
&+n\left(
4M^{(c)}_{a+2,b,n-1}
+2cM^{(c)}_{a+1,b+1,n-1}
\right),
\end{aligned}
}
```

and

```math
\boxed{
M^{(c)}_{a,b+1,n}
=
(2m+2b)M^{(c)}_{a,b,n}
+2cnM^{(c)}_{a+1,b+1,n-1}.
}
```

Terms multiplied by \(n\) are omitted when \(n=0\). The defining polynomial
also gives the algebraic closure

```math
\boxed{
M^{(c)}_{a+2,b,n}
+cM^{(c)}_{a+1,b+1,n}
=
M^{(c)}_{a,b,n+1}.
}
```

The resolvent and Laplace-observable systems deform in the same way. For

```math
G^{(c)}_{a,b}(z)
=
\mathbb E\left[\frac{U^aV^b}{z-T_{m,c}}\right],
```

the multiplication closure is

```math
\boxed{
G^{(c)}_{a+2,b}(z)
+cG^{(c)}_{a+1,b+1}(z)
=
zG^{(c)}_{a,b}(z)-M^{(c)}_{a,b,0}.
}
```

Thus the exact hierarchy in the GOE case is one member of a continuous
family of gamma-pushforward hierarchies.

## 4. Quantitative convergence to the chi-square law

The qualitative limit

```math
\frac{R_m^4}{4m}\ \Longrightarrow\ \chi_1^2
```

was recorded in [ORTHOGONALITY_NOTE.md](ORTHOGONALITY_NOTE.md). The same
coupling gives a nonasymptotic rate for the full deformation.

### Theorem 2

Let \(\mathcal L(X)\) denote the law of \(X\), and let \(W_2\) be the
quadratic Wasserstein distance. Then

```math
\boxed{
W_2\left(
\mathcal L(W_{m,c}),
\chi_1^2
\right)
\leq
\sqrt{
\frac{3}{m}
+\frac{105}{4c^2m^2}
}.
}
```

More precisely, on the defining coupling with \(U\sim\chi_1^2\),

```math
\boxed{
\mathbb E[W_{m,c}-U]
=
\frac{3}{2cm},
}
```

and

```math
\boxed{
\mathbb E\left[(W_{m,c}-U)^2\right]
=
\frac{3}{m}
+\frac{105}{4c^2m^2}.
}
```

#### Proof

Write

```math
W_{m,c}-U
=
\frac{U(V-2m)}{2m}
+\frac{U^2}{2cm}.
```

The gamma moments needed below are

```math
\mathbb E[U^2]=3,
\qquad
\mathbb E[U^4]=105,
\qquad
\mathbb E[V-2m]=0,
\qquad
\mathbb E[(V-2m)^2]=4m.
```

The mean formula follows immediately. After squaring, the cross term has
expectation zero by independence and \(\mathbb E[V-2m]=0\). Therefore

```math
\begin{aligned}
\mathbb E\left[(W_{m,c}-U)^2\right]
&=
\frac{\mathbb E[U^2]\mathbb E[(V-2m)^2]}{4m^2}
+\frac{\mathbb E[U^4]}{4c^2m^2}\\
&=
\frac3m+\frac{105}{4c^2m^2}.
\end{aligned}
```

The Wasserstein bound follows because \(W_2^2\) is the infimum of the
mean-square cost over all couplings, and the displayed construction is one
admissible coupling. \(\square\)

For the exceptional GOE factor \(c=2\), this becomes

```math
\boxed{
W_2\left(
\mathcal L\left(\frac{R_m^4}{4m}\right),
\chi_1^2
\right)
\leq
\sqrt{
\frac3m+\frac{105}{16m^2}
}.
}
```

## 5. First correction to every fixed moment

The exact formula also determines the first correction to the limiting
chi-square moments.

### Theorem 3

For each fixed positive integer \(n\),

```math
\boxed{
\begin{aligned}
\mathbb E[W_{m,c}^n]
={}&
2^n\left(\frac12\right)_n\\
&+
\frac{
2^{n-1}n\left(\frac12\right)_n
\left[
(n-1)+\frac{2}{c}\left(n+\frac12\right)
\right]
}{m}
+O_{n,c}\left(m^{-2}\right).
\end{aligned}
}
```

In the GOE case \(c=2\),

```math
\boxed{
\mathbb E\left[
\left(\frac{R_m^4}{4m}\right)^n
\right]
=
2^n\left(\frac12\right)_n
+\frac{
2^{n-2}n(4n-1)\left(\frac12\right)_n
}{m}
+O_n\left(m^{-2}\right).
}
```

#### Proof

Divide the formula of Theorem 1, with \(a=b=0\), by \((2cm)^n\).
The summand \(j=n\) contributes

```math
2^n\left(\frac12\right)_n
\frac{(m)_n}{m^n}
=
2^n\left(\frac12\right)_n
\left(
1+\frac{n(n-1)}{2m}+O_n(m^{-2})
\right).
```

The summand \(j=n-1\) contributes

```math
\frac{
2^n n\left(\frac12\right)_{n+1}
}{cm}
+O_{n,c}(m^{-2}).
```

All remaining summands are \(O_{n,c}(m^{-2})\). Using

```math
\left(\frac12\right)_{n+1}
=
\left(n+\frac12\right)\left(\frac12\right)_n
```

and collecting the two first-order contributions proves the formula.
\(\square\)

For \(n=0\), the normalized moment is identically \(1\).

## 6. Scope

The deformation theorem shows that the loop hierarchy and the chi-square
limit are structural consequences of the two gamma coordinates, rather than
isolated identities caused by the coefficient \(2\). The Wasserstein rate
and the first moment correction are exact consequences of this model.

These results do not assert that every \(c\) arises from a random-matrix
factorization, nor do they claim global literature priority for the
elementary gamma calculations.
