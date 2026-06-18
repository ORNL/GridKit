# SkinEffect Model

`SkinEffect` computes frequency-dependent internal resistance and inductance
vectors for overhead conductors using the simplified skin-effect formulation of
Monteiro et al. [1].

## Model Parameters

For $K$ physical conductors:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$r_i$ | [m] | `radius` | Outer radius | from [`Conductor`](../../Geometry/Conductor/README.md)
$\sigma_i$ | [S/m] | `conductivity` | Conductivity | from [`Conductor`](../../Geometry/Conductor/README.md)
$\mu_i$ | [H/m] | `permeability` | Permeability | from [`Conductor`](../../Geometry/Conductor/README.md)

### Parameter Validation

The connected conductor model owns vector validation. This model requires

```math
\omega>0 .
```

The formulation uses the solid cylindrical-conductor approximation in [1].

### Model Derived Parameters

The number of retained root terms is fixed by the model:

```math
N_s=8 .
```

For conductor $i$ and retained root $k=1,\dots,N_s$,

```math
\begin{aligned}
\xi_k &= \left(k-\frac{1}{4}\right)\pi \\
R_{ik}^{\mathrm{br}} &=
  \frac{\xi_k^2}{4\pi\sigma_i r_i^2} \\
L_i^{\mathrm{br}} &=
  \frac{\mu_i}{4\pi}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{Q}$ | [$(\Omega/\mathrm{m})^2$] | Branch denominator | $\mathbb{R}^{K\times N_s}$
$\mathbf{G}$ | [m/$\Omega$] | Real internal-admittance sum | $\mathbb{R}^{K}$
$\mathbf{H}$ | [s m/$\Omega$] | Inductive internal-admittance sum | $\mathbb{R}^{K}$
$\mathbf{W}$ | [$(\mathrm{m}/\Omega)^2$] | Internal-impedance denominator | $\mathbb{R}^{K}$
$\mathbf{r}^{\mathrm{skin}}$ | [$\Omega$/m] | Internal conductor resistance | $\mathbb{R}^{K}$
$\mathbf{l}^{\mathrm{skin}}$ | [H/m] | Internal conductor inductance | $\mathbb{R}^{K}$

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
0 &= -Q_{ik}
  + (R_{ik}^{\mathrm{br}})^2
  + \omega^2(L_i^{\mathrm{br}})^2 \\
0 &= -G_i
  + \sum_{k=1}^{N_s}
  \frac{R_{ik}^{\mathrm{br}}}{Q_{ik}} \\
0 &= -H_i
  + \sum_{k=1}^{N_s}
  \frac{L_i^{\mathrm{br}}}{Q_{ik}} \\
0 &= -W_i
  + G_i^2
  + \omega^2H_i^2 \\
0 &= -W_i r_i^{\mathrm{skin}}
  + G_i \\
0 &= -W_i l_i^{\mathrm{skin}}
  + H_i
\end{aligned}
```

## Initialization

Initialize all algebraic variables from the algebraic equations at the current
$\omega$.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{r}^{\mathrm{skin}}$ | [$\Omega$/m] | Internal conductor resistance | $\mathbb{R}^{K}$
$\mathbf{l}^{\mathrm{skin}}$ | [H/m] | Internal conductor inductance | $\mathbb{R}^{K}$

## References

[1] J. H. A. Monteiro, E. C. M. Costa, A. J. G. Pinto, S. Kurokawa,
O. M. O. Gatous, and J. Pissolato, "Simplified skin-effect formulation
for power transmission lines," IET Science, Measurement & Technology,
vol. 8, no. 2, pp. 47-53, 2014. doi:10.1049/iet-smt.2013.0072.
