# Optimal Power Flow

AC optimal power flow in Cartesian coordinates. It finds the least-cost
operating point within the component limits.

## Formulation

Minimize

```math
\begin{aligned}
f &= \sum_g u_g (c_{0g} + c_{1g} P_g + c_{2g} P_g^2)
\end{aligned}
```

subject to

```math
\begin{aligned}
0 &= \sum_{c \to k} P_{ck} && \text{each bus } k \\
0 &= \sum_{c \to k} Q_{ck} && \text{each bus } k \\
(V_k^{\min})^2 &\le V_{\mathrm{r}k}^2 + V_{\mathrm{i}k}^2 \le (V_k^{\max})^2 && \text{each bus } k \\
0 &= V_{\mathrm{i}k} && \text{each reference bus } k \\
P_{\ell k}^2 + Q_{\ell k}^2 &\le (S_\ell^{\max})^2 && \text{each terminal } k \text{ of each branch } \ell \\
P_g^{\min} &\le P_g \le P_g^{\max} && \text{each in-service generator } g \\
Q_g^{\min} &\le Q_g \le Q_g^{\max} && \text{each in-service generator } g \\
0 &= P_g && \text{each out-of-service generator } g \\
0 &= Q_g && \text{each out-of-service generator } g
\end{aligned}
```

## Derivatives

Each component supplies `objective(x)` and `constraints(x, g)`, which
`ComponentModel` differentiates with Enzyme:

Driver   | Derivative                                                | Mode
---------|-----------------------------------------------------------|-----
`DfDx`   | Objective gradient                                        | Reverse
`DgDx`   | Constraint Jacobian                                       | Forward
`D2LDx2` | Lower triangle of the Hessian of $\sigma f + \lambda^T g$ | Forward over reverse

## Data

`parseMatpowerCase` reads a
[MATPOWER case](https://matpower.org/docs/ref/matpower/caseformat.html):

MATPOWER                                               | Component   | Parameters
-------------------------------------------------------|-------------|-----------
bus `BUS_TYPE` 3                                       | `Bus`       | `reference`
bus `VMIN`, `VMAX`                                     | `Bus`       | $V^{\min}$, $V^{\max}$
bus `PD`, `QD` when either is nonzero                  | `Load`      | $P$, $Q$
bus `GS`, `BS` when either is nonzero                  | `Shunt`     | $G$, $B$
branch `BR_R`, `BR_X`, `BR_B`                          | `Branch`    | $R$, $X$, $B$
branch `TAP`, `SHIFT`                                  | `Branch`    | $\tau$, $\theta$
branch `RATE_A` when positive                          | `Branch`    | $S^{\max}$
branch `BR_STATUS`                                     | `Branch`    | $u$
gen `PMIN`, `PMAX`, `QMIN`, `QMAX`                     | `Generator` | $P^{\min}$, $P^{\max}$, $Q^{\min}$, $Q^{\max}$
gen `GEN_STATUS`                                       | `Generator` | $u$
gencost `COST` of a polynomial of at most second order | `Generator` | $c_2$, $c_1$, $c_0$

A zero `TAP` gives $\tau = 1$ for a line. The branch also gets
$B_\mathrm{mag} = B (1 / \tau^2 - 1) / 2$ so that the line charging at bus 1
scales by $1 / \tau^2$ as in MATPOWER. Values in MW, Mvar, and MVA and costs in
P [MW] become system base values, and `SHIFT` in degrees becomes radians. Every
island needs a reference bus. The model has no isolated buses, angle-difference
limits, or other cost models, so `parseMatpowerCase` throws for a case that uses
them.
