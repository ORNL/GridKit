# AC Optimal Power Flow

GridKit implements a steady-state polar AC OPF with fixed network topology. All
electrical quantities are per unit on `va_base`; angles and phase shifts are in
radians.

## Formulation

The decision vector contains bus voltage magnitude and angle, followed by
generator active and reactive power:

$$
x = \{v_i,\theta_i\}_{i\in\mathcal N}
    \cup \{p_g,q_g\}_{g\in\mathcal G},
\qquad V_i=v_i e^{\mathrm j\theta_i}.
$$

The objective is the sum of the per-unit generator costs:

$$
f(x)=\sum_{g\in\mathcal G}a_g
\left(c_{0g}+c_{1g}p_g+c_{2g}p_g^2\right),
$$

where $a_g$ is one for an online generator and zero otherwise.

For branch $\ell=(i,j)$, with the transformer on the from side,

$$
y_\ell=\frac{1}{R_\ell+\mathrm jX_\ell},\qquad
y_\ell^c=G_\ell+\mathrm jB_\ell,\qquad
T_\ell=\tau_\ell e^{\mathrm j\phi_\ell},
$$

$$
\begin{bmatrix}I_{\ell i}\\I_{\ell j}\end{bmatrix}
=a_\ell
\begin{bmatrix}
(y_\ell+y_\ell^c/2)/|T_\ell|^2 & -y_\ell/T_\ell^*\\
-y_\ell/T_\ell & y_\ell+y_\ell^c/2
\end{bmatrix}
\begin{bmatrix}V_i\\V_j\end{bmatrix},
\qquad
S_{\ell k}=V_k I_{\ell k}^*.
$$

$S_{\ell k}$ is power leaving bus $k$ and entering the branch. The two bus
balance equations are the real and imaginary parts of

$$
\sum_{g\in\mathcal G_i}S_g
+\sum_{d\in\mathcal D_i}S_d^{\mathrm{inj}}
-\sum_{h\in\mathcal H_i}a_h(G_h-\mathrm jB_h)v_i^2
-\sum_{\ell\in\delta(i)}S_{\ell i}=0.
$$

Load state uses a signed-injection convention: ordinary demand has negative
`p` and `q`. A rated branch adds

$$
0\le |S_{\ell i}|^2\le(s_\ell^{max})^2,
\qquad
0\le |S_{\ell j}|^2\le(s_\ell^{max})^2.
$$

Bus voltage and generator power use their supplied bounds. The angle of the
lowest-numbered bus is fixed to its input-state angle to remove rotational
freedom; it is not a slack-bus model.

The objective gradient, constraint Jacobian, and lower-triangular Hessian of

$$
\mathcal L(x,\lambda,\sigma)=\sigma f(x)+\lambda^\mathsf T g(x)
$$

are generated exactly with Enzyme and assembled into fixed sparse structures.

## State

The case structure is defined by `.opf.json`; operating values come from the
shared `.state.json` format.

| Element | Input state |
|---|---|
| Bus | Required finite `vr`, `vi` |
| Generator | `online` defaults to true; `p`, `q` default to zero |
| Load | Required signed `p`, `q`; `online` defaults to true |
| Branch | `open` defaults to false, `tap` to one, `phase` to zero |
| Shunt | `online` defaults to true |

Solved state updates bus voltage, generator power, and bus current injections
while preserving unrelated state. The generic `active` field is not used by
OPF. See [INPUT_FORMAT.md](INPUT_FORMAT.md) for the JSON fields.

## References

- [MATPOWER User's Manual, AC Optimal Power Flow](https://matpower.org/docs/MATPOWER-manual-8.1.pdf)
- [MATPOWER's Extensible Optimal Power Flow Architecture](https://matpower.org/docs/MATPOWER-OPF.pdf)
- [PowerModels.jl: An Open-Source Framework for Exploring Power Flow Formulations](https://arxiv.org/abs/1711.01728)
