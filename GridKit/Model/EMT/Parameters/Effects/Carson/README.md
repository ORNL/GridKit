# Carson Model

`Carson` computes frequency-dependent earth-return resistance and inductance
matrices for overhead conductors using the Deri-Semlyen approximation to
Carson's homogeneous-earth return term.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
`tower` | [-] | - | Static conductor tower | [`Tower`](../../Geometry/Tower/README.md)
$\sigma_g$ | [S/m] | `earth.conductivity` | Earth conductivity | $\sigma_g\ge0$
$\varepsilon_g$ | [F/m] | `earth.permittivity` | Earth permittivity | $\varepsilon_g\ge\varepsilon_0$
$\mu_0$ | [H/m] | - | Vacuum permeability | constant
$\varepsilon_0$ | [F/m] | - | Vacuum permittivity | constant

### Parameter Validation

```math
\sigma_g\ge0,\qquad
\varepsilon_g\ge\varepsilon_0,\qquad
\omega > 0,\qquad
(\gamma_g^{\mathrm r})^2+(\gamma_g^{\mathrm i})^2 > 0
```

The connected tower model owns its geometry validation.

### Model Derived Parameters

For conductor pair $(i,j)$,

```math
\begin{aligned}
d_{ij} &= |x_i-x_j| \\
s_{ij} &= h_i+h_j \\
\Delta_g &= (\gamma_g^{\mathrm r})^2+(\gamma_g^{\mathrm i})^2 \\
p^{\mathrm r} &= \frac{\gamma_g^{\mathrm r}}{\Delta_g} \\
p^{\mathrm i} &= -\frac{\gamma_g^{\mathrm i}}{\Delta_g} \\
H_{ij}^{\mathrm r} &= s_{ij}+2p^{\mathrm r} \\
H_{ij}^{\mathrm i} &= 2p^{\mathrm i} \\
A_{ij} &= d_{ij}^2+(H_{ij}^{\mathrm r})^2-(H_{ij}^{\mathrm i})^2 \\
B_{ij} &= 2H_{ij}^{\mathrm r}H_{ij}^{\mathrm i} \\
M_{ij} &= A_{ij}^2+B_{ij}^2 \\
\theta_{ij} &= \operatorname{atan2}(B_{ij},A_{ij}) \\
\Lambda_{ij} &= \frac{1}{4}\ln M_{ij}-\ln D'_{ij}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\eta_g$ | [1/m$^2$] | Magnitude of earth-return propagation radicand | real
$\gamma_g^{\mathrm r}$ | [1/m] | Real earth-return propagation constant component | real
$\gamma_g^{\mathrm i}$ | [1/m] | Imaginary earth-return propagation constant component | real
$\mathbf{R}^{\mathrm{carson}}$ | [$\Omega$/m] | Carson earth-return resistance | $\mathbb{R}^{K\times K}$
$\mathbf{L}^{\mathrm{carson}}$ | [H/m] | Carson earth-return inductance | $\mathbb{R}^{K\times K}$

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
0 &= -\eta_g^2
  + (\omega^2\mu_0\varepsilon_g)^2
  + (\omega\mu_0\sigma_g)^2 \\
0 &= -2(\gamma_g^{\mathrm r})^2
  + \eta_g - \omega^2\mu_0\varepsilon_g \\
0 &= -2(\gamma_g^{\mathrm i})^2
  + \eta_g + \omega^2\mu_0\varepsilon_g \\
0 &= -4\pi R_{ij}^{\mathrm{carson}}
  - \omega\mu_0\theta_{ij} \\
0 &= -2\pi L_{ij}^{\mathrm{carson}}
  + \mu_0\Lambda_{ij}
\end{aligned}
```

## Initialization

Initialize $\eta_g$, $\gamma_g^{\mathrm r}$, $\gamma_g^{\mathrm i}$,
$\mathbf{R}^{\mathrm{carson}}$, and $\mathbf{L}^{\mathrm{carson}}$ from the
algebraic equations at the current $\omega$.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\sigma_g$ | [S/m] | Earth conductivity | real
$\varepsilon_g$ | [F/m] | Earth permittivity | real
$\eta_g$ | [1/m$^2$] | Magnitude of earth-return propagation radicand | real
$\gamma_g^{\mathrm r}$ | [1/m] | Real earth-return propagation constant component | real
$\gamma_g^{\mathrm i}$ | [1/m] | Imaginary earth-return propagation constant component | real
$\mathbf{R}^{\mathrm{carson}}$ | [$\Omega$/m] | Carson earth-return resistance | $\mathbb{R}^{K\times K}$
$\mathbf{L}^{\mathrm{carson}}$ | [H/m] | Carson earth-return inductance | $\mathbb{R}^{K\times K}$
