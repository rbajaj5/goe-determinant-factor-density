# Exact loop equations for the exceptional GOE factor

## 1. Setup

Let

```math
U\sim\mathrm{Gamma}\left(\frac12,2\right),
\qquad
V\sim\mathrm{Gamma}(m,2),
\qquad
U\text{ and }V\text{ independent},
```

where \(m>0\), and set

```math
T=U^2+2UV.
```

For the determinant factor studied in this repository, \(T=R_m^4\).
The gamma lift \((U,V)\) makes it possible to derive an exact hierarchy of
integration-by-parts identities. These identities are a finite-dimensional
analogue of the loop-equation strategy used in random-matrix theory.

## 2. The underlying Ward identities

For every continuously differentiable test function \(g\) for which the
expectations below exist,

```math
\boxed{
\mathbb E\left[
2U\,\partial_u g(U,V)+(1-U)g(U,V)
\right]=0
}
```

and

```math
\boxed{
\mathbb E\left[
2V\,\partial_v g(U,V)+(2m-V)g(U,V)
\right]=0.
}
```

For compactly supported \(g\), these follow immediately by integration by
parts against the joint density

```math
\frac{
u^{-1/2}v^{m-1}e^{-(u+v)/2}
}{
2^{m+1/2}\Gamma(1/2)\Gamma(m)
}.
```

The identities extend to the test functions used below by approximation,
because the gamma tails dominate every polynomial factor.

### Exponential linearization

Define the joint Laplace transform

```math
\Phi(s,t)=\mathbb E\left[e^{-sU-tV}\right],
\qquad s,t\geq0.
```

Taking \(g(u,v)=e^{-su-tv}\) in the two Ward identities gives the closed
linear system

```math
\boxed{
(1+2s)\partial_s\Phi+\Phi=0,
\qquad
(1+2t)\partial_t\Phi+2m\Phi=0,
\qquad
\Phi(0,0)=1.
}
```

Its unique solution is

```math
\boxed{
\Phi(s,t)=(1+2s)^{-1/2}(1+2t)^{-m}.
}
```

This supplies another characterization: any probability law on the
nonnegative quadrant whose Laplace transform satisfies this system is the
product of the two gamma laws above. Unlike a moment-generating-function
argument, this Laplace-transform formulation needs no positive exponential
moment.

## 3. Exact polynomial loop hierarchy

### Theorem 1

Let \(a,b\) be nonnegative integers. For every continuously differentiable
function \(f\) for which the displayed expectations are finite,

```math
\boxed{
\begin{aligned}
0={}&
\mathbb E\left[
\big((2a+1)U^aV^b-U^{a+1}V^b\big)f(T)
\right.\\
&\left.\hspace{31mm}
+4U^{a+1}V^b(U+V)f'(T)
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
\big((2m+2b)U^aV^b-U^aV^{b+1}\big)f(T)
\right.\\
&\left.\hspace{31mm}
+4U^{a+1}V^{b+1}f'(T)
\right].
\end{aligned}
}
```

#### Proof

Apply the first Ward identity to

```math
g(u,v)=u^av^bf(u^2+2uv).
```

Since

```math
\partial_u T=2(u+v),
```

the product and chain rules give the first formula. Applying the second Ward
identity to the same \(g\), using

```math
\partial_v T=2u,
```

gives the second formula. This proves both identities. \(\square\)

## 4. An exact resolvent hierarchy

For \(z\) outside the nonnegative real axis, define the mixed resolvents

```math
G_{a,b}(z)
=
\mathbb E\left[\frac{U^aV^b}{z-T}\right].
```

Taking \(f(t)=(z-t)^{-1}\) in Theorem 1 and observing that

```math
\mathbb E\left[\frac{U^aV^b}{(z-T)^2}\right]
=-\partial_zG_{a,b}(z)
```

gives

```math
\boxed{
(2a+1)G_{a,b}
-G_{a+1,b}
-4\partial_z\left(G_{a+2,b}+G_{a+1,b+1}\right)
=0
}
```

and

```math
\boxed{
(2m+2b)G_{a,b}
-G_{a,b+1}
-4\partial_zG_{a+1,b+1}
=0.
}
```

The polynomial definition of \(T\) supplies a second, algebraic part of the
system:

```math
\boxed{
G_{a+2,b}(z)+2G_{a+1,b+1}(z)
=
zG_{a,b}(z)-M_{a,b,0}.
}
```

Indeed,

```math
\frac{T}{z-T}=\frac{z}{z-T}-1,
```

and multiplication by \(U^aV^b\) followed by expectation proves the
identity. Differentiating it also gives

```math
\partial_zG_{a+2,b}
+2\partial_zG_{a+1,b+1}
=
G_{a,b}+z\partial_zG_{a,b}.
```

Thus the mixed resolvents satisfy an exact differential-algebraic loop
system. It is not yet a closed scalar differential equation for \(G_{0,0}\);
eliminating all mixed resolvents would be an additional result.

### Laplace-observable hierarchy for \(T\)

Define

```math
F_{a,b}(\lambda)
=
\mathbb E\left[U^aV^b e^{-\lambda T}\right],
\qquad \lambda\geq0.
```

Taking \(f(t)=e^{-\lambda t}\) in Theorem 1 gives the linear relations

```math
\boxed{
\begin{aligned}
F_{a+1,b}
={}&(2a+1)F_{a,b}\\
&-4\lambda\left(F_{a+2,b}+F_{a+1,b+1}\right),
\end{aligned}
}
```

and

```math
\boxed{
F_{a,b+1}
=(2m+2b)F_{a,b}
-4\lambda F_{a+1,b+1}.
}
```

Differentiation with respect to \(\lambda\) supplies the algebraic part:

```math
\boxed{
\partial_\lambda F_{a,b}
=
-F_{a+2,b}-2F_{a+1,b+1}.
}
```

Thus an exponential observable turns the exact Ward hierarchy into a linear
differential-algebraic system. The initial values at \(\lambda=0\) are

```math
F_{a,b}(0)
=
2^{a+b}\left(\frac12\right)_a(m)_b.
```

## 5. Exact mixed-moment recurrences

Put

```math
M_{a,b,n}
=
\mathbb E\left[U^aV^bT^n\right].
```

Taking \(f(t)=t^n\) in Theorem 1 gives, for \(n\geq 1\),

```math
\boxed{
\begin{aligned}
M_{a+1,b,n}
={}&(2a+1)M_{a,b,n}\\
&+4n\left(M_{a+2,b,n-1}+M_{a+1,b+1,n-1}\right),
\end{aligned}
}
```

and

```math
\boxed{
M_{a,b+1,n}
=(2m+2b)M_{a,b,n}
+4nM_{a+1,b+1,n-1}.
}
```

The same equations hold for \(n=0\) after the terms multiplied by \(n\) are
deleted.

These recurrences are independently closed by the exact formula

```math
\boxed{
M_{a,b,n}
=
2^{a+b+2n}
\sum_{j=0}^{n}
\binom nj 2^j
\left(\frac12\right)_{a+2n-j}
(m)_{b+j}.
}
```

Indeed, expand

```math
T^n=(U^2+2UV)^n
```

and use the gamma moments

```math
\mathbb E[U^r]=2^r\left(\frac12\right)_r,
\qquad
\mathbb E[V^s]=2^s(m)_s.
```

The case \(a=b=0\) recovers the moment formula in
[ORTHOGONALITY_NOTE.md](ORTHOGONALITY_NOTE.md).

## 6. The polynomial hierarchy characterizes the law

### Theorem 2

Let \((\widetilde U,\widetilde V)\) be a random vector in the nonnegative
quadrant, and suppose that

```math
\mathbb E\left[e^{c(\widetilde U+\widetilde V)}\right]<\infty
```

for some \(c>0\). Assume that, for every \(a,b\geq0\), its mixed moments
\(\widetilde M_{a,b}=\mathbb E[\widetilde U^a\widetilde V^b]\) satisfy the
base level of the two loop recurrences:

```math
\widetilde M_{a+1,b}=(2a+1)\widetilde M_{a,b},
\qquad
\widetilde M_{a,b+1}=(2m+2b)\widetilde M_{a,b}.
```

Then \(\widetilde U\) and \(\widetilde V\) are independent and

```math
\widetilde U\sim\mathrm{Gamma}\left(\frac12,2\right),
\qquad
\widetilde V\sim\mathrm{Gamma}(m,2).
```

In particular, the polynomial loop hierarchy, together with normalization
and an exponential-moment bound, characterizes the lifted law.

#### Proof

Starting from \(\widetilde M_{0,0}=1\), the two recurrences give

```math
\widetilde M_{a,b}
=
2^{a+b}
\left(\frac12\right)_a(m)_b.
```

The exponential-moment assumption makes the joint moment-generating function
analytic in a neighborhood of the origin. Therefore

```math
\begin{aligned}
\mathbb E\left[e^{s\widetilde U+t\widetilde V}\right]
&=
\sum_{a,b\geq0}
\frac{\widetilde M_{a,b}}{a!\,b!}s^at^b\\
&=
(1-2s)^{-1/2}(1-2t)^{-m}.
\end{aligned}
```

This is the product of the two stated gamma moment-generating functions, so
it determines both marginals and their independence. \(\square\)

The point of the hypothesis is analytic determinacy, not a hidden
random-matrix assumption. It can be replaced by any standard condition that
makes the joint moment problem determinate.

## 7. A density-level characterization

### Theorem 3

Let \(q\) be a probability density on the positive quadrant that is locally
absolutely continuous. Suppose that, for every compactly supported smooth
test function \(g\),

```math
\int
\left(2u\,\partial_u g+(1-u)g\right)q\,\mathrm du\,\mathrm dv=0
```

and

```math
\int
\left(2v\,\partial_v g+(2m-v)g\right)q\,\mathrm du\,\mathrm dv=0.
```

Then

```math
\boxed{
q(u,v)
=
\frac{
u^{-1/2}v^{m-1}e^{-(u+v)/2}
}{
2^{m+1/2}\Gamma(1/2)\Gamma(m)
}.
}
```

Consequently, the two Ward identities characterize the joint law of
\((U,V)\), and hence characterize the pushforward law of \(T=U^2+2UV\).

#### Proof

Distributional integration by parts in the first identity gives

```math
-\partial_u(2uq)+(1-u)q=0.
```

For almost every fixed \(v\), this first-order equation has the form

```math
q(u,v)=c(v)u^{-1/2}e^{-u/2}.
```

The second identity similarly forces

```math
q(u,v)=d(u)v^{m-1}e^{-v/2}.
```

The two factorizations together determine \(q\) up to a constant.
Normalization gives the displayed density. \(\square\)

## 8. Relation to random-matrix loop equations

Bourgade and Huang prove that the universal local point processes
\(\mathrm{Sine}_\beta\) and \(\mathrm{Airy}_\beta\), for rational
\(\beta>0\), are characterized by their loop-equation hierarchies. Their
universality argument has three essential ingredients:

1. a random point process and its resolvent hierarchy;
2. approximate loop equations obtained from integration by parts or a
   discrete switching identity;
3. analytic input such as local laws and resolvent stability that controls
   the approximation error.

The result above shares the proof architecture

```text
integration by parts -> hierarchy of identities -> characterization of a law,
```

but not the asymptotic random-matrix conclusion. The exceptional determinant
factor is a single scalar random variable, not a microscopic eigenvalue point
process. It has no bulk or edge scaling limit to which the
\(\mathrm{Sine}_\beta\) or \(\mathrm{Airy}_\beta\) uniqueness theorem can be
applied.

Their analytic uniqueness proof also has to select admissible solution
branches. They first exponentiate integrated Stieltjes transforms to obtain
characteristic-ratio observables satisfying deformed
Calogero--Moser--Sutherland equations. In the edge problem, the subsequent
uniqueness analysis involves square-free derivative reduction, a radial block
system with a spectral gap, and a Volterra contraction.

None of that machinery is silently imported here. The exponentials
\(e^{-sU-tV}\) and \(e^{-\lambda T}\) play an analogous linearizing role, but
their equations are the first-order gamma system and the mixed
differential-algebraic system displayed above, not deformed
Calogero--Moser--Sutherland equations. After lifting the scalar variable to
two independent gamma coordinates, the loop equations are exact, and their
base polynomial level is triangular. Theorem 2 therefore selects the law
directly from its mixed moments, without a local law, cumulant truncation, or
asymptotic error estimate.

There is also an important algebraic difference. In the matrix calculation,
differentiating the matrix resolvent and summing over matrix indices produces
quadratic trace terms such as the square of the normalized Stieltjes
transform. For this scalar factor, differentiating \((z-T)^{-1}\) introduces
mixed \(U,V\) weights instead. The multiplication identity in Section 4
controls those weights, but it does not by itself produce a closed quadratic
equation for \(G_{0,0}\).

For random regular graphs, Bourgade and Huang additionally split into a
high-probability local-law event and its complement, linearize a discrete
switching identity, and prove that the Taylor remainder is negligible. There
is no corresponding bad-event or switching remainder here: the gamma
integration-by-parts identities and the exponential-observable system are
finite-dimensional equalities. A perturbative theorem for approximately
gamma inputs would require a new quantitative error hypothesis, such as a
Stein discrepancy; it is not implied by the exact calculation.

## 9. Status

- The Ward identities, loop hierarchy, moment recurrences, and
  two characterization theorems are proved above.
- The identities are specializations of standard gamma integration by parts;
  global novelty is not claimed.
- The exponential-observable linearization is an exact gamma analogue of the
  organizational move in Bourgade--Huang, not an identification with their
  deformed Calogero--Moser--Sutherland system.
- The exact mixed hierarchy is useful for systematic moment calculations and
  for checking future differential equations for the density or Stieltjes
  transform.
- No random-matrix universality theorem is claimed for this scalar factor.

## Reference

Paul Bourgade and Jiaoyang Huang,
[*Loop Equations Characterize Random Matrix Statistics*](https://arxiv.org/abs/2607.07617),
arXiv:2607.07617, 2026.
