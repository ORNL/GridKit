# H Model

`H` computes the modal-domain finite-length propagation function from modal
propagation constants.

The model is real-valued in implementation. The real and imaginary parts of
$h_m$ are owned as separate differential states.

`H` uses derivative signals from `Gamma`. In logarithmic frequency sweeps,
`LogEvaluator` converts solver derivatives to derivatives with respect to
angular frequency before `H` is evaluated.

```math
h_m = \exp(-\gamma_m \ell),\qquad
\gamma_m = a_m + j b_m
```

Equivalently,

```math
\begin{aligned}
h^{\mathrm r}_m &=
e^{-a_m\ell}\cos(b_m\ell) \\
h^{\mathrm i}_m &=
-e^{-a_m\ell}\sin(b_m\ell)
\end{aligned}
```

The differential form used during the sweep is

```math
\dot h_m = -\ell \dot\gamma_m h_m,\qquad
\dot\gamma_m = \dot a_m + j\dot b_m .
```

## Model Parameters

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\ell$ | [m] | Path length | from [`Path`](../../Geometry/Path/README.md)

### Parameter Validation

The path length must satisfy

```math
\ell > 0 .
```

The connected `Gamma` model owns modal-vector validation.

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$h^{\mathrm r}_m$ | [-] | Real part of modal propagation function for mode $m$ | real
$h^{\mathrm i}_m$ | [-] | Imaginary part of modal propagation function for mode $m$ | real

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\dot a_m$ | [s/m] | Modal attenuation derivative for mode $m$ | from [`Gamma`](../Gamma/README.md)
$\dot b_m$ | [s/m] | Modal phase derivative for mode $m$ | from [`Gamma`](../Gamma/README.md)

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$a_m$ | [1/m] | Modal attenuation constant for mode $m$ | from [`Gamma`](../Gamma/README.md)
$b_m$ | [1/m] | Modal phase constant for mode $m$ | from [`Gamma`](../Gamma/README.md)

## Model Equations

### Differential Equations

For each mode $m=1,\dots,K$:

```math
\begin{aligned}
0 &= -\dot h^{\mathrm r}_m
  - \ell\left(\dot a_m h^{\mathrm r}_m
  - \dot b_m h^{\mathrm i}_m\right) \\
0 &= -\dot h^{\mathrm i}_m
  - \ell\left(\dot a_m h^{\mathrm i}_m
  + \dot b_m h^{\mathrm r}_m\right)
\end{aligned}
```

### Algebraic Equations

None.

## Initialization

Initialize $h^{\mathrm r}_m$ and $h^{\mathrm i}_m$ from the current $a_m$ and
$b_m$ inputs using the closed-form expression above. Initialize the derivatives
from $\dot h_m=-\ell\dot\gamma_m h_m$.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{h}^{\mathrm r}$ | [-] | Real part of modal propagation function | $\mathbb{R}^K$
$\mathbf{h}^{\mathrm i}$ | [-] | Imaginary part of modal propagation function | $\mathbb{R}^K$

## Monitors

Monitor | Symbol | Units | Shape | Description
------- | ------ | ----- | ----- | -----------
`H` | $\mathbf{h}^{\mathrm r}+j\mathbf{h}^{\mathrm i}$ | [-] | $K$ | Modal propagation function

CSV output uses explicit real and imaginary columns:

```text
Overhead_H_real_m
Overhead_H_imag_m
```

No complex-valued states, residuals, monitor storage, or model API are used.
