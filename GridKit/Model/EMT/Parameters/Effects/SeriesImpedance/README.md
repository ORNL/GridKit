# SeriesImpedance Model

`SeriesImpedance` builds the full-conductor per-unit-length series resistance
and inductance matrices for one angular frequency.

## Model Parameters

None. `SeriesImpedance` combines effect-model outputs for $K$ physical
conductors.

### Parameter Validation

None.

### Model Derived Parameters

For a vector $\mathbf a\in\mathbb{R}^K$, let
$\operatorname{diag}(\mathbf a)$ denote the diagonal matrix with $\mathbf a$ on
the main diagonal.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}^{\mathrm{int}}$ | [$\Omega$/m] | Internal conductor resistance | diagonal
$\mathbf{L}^{\mathrm{int}}$ | [H/m] | Internal conductor inductance | diagonal
$\mathbf{R}^{\mathrm{ext}}$ | [$\Omega$/m] | External series resistance | full matrix
$\mathbf{L}^{\mathrm{ext}}$ | [H/m] | External series inductance | full matrix
$\mathbf{R}'$ | [$\Omega$/m] | Full-conductor series resistance | $\mathbb{R}^{K\times K}$
$\mathbf{L}'$ | [H/m] | Full-conductor series inductance | $\mathbb{R}^{K\times K}$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{r}^{\mathrm{skin}}$ | [$\Omega$/m] | Internal conductor resistance vector | from connected `SkinEffect` models
$\mathbf{l}^{\mathrm{skin}}$ | [H/m] | Internal conductor inductance vector | from connected `SkinEffect` models
$\mathbf{L}^{\mathrm{geo}}$ | [H/m] | External geometric inductance | from connected `GeometricInductance` model
$\mathbf{R}^{\mathrm{carson}}$ | [$\Omega$/m] | Carson earth-return resistance | from connected `Carson` models
$\mathbf{L}^{\mathrm{carson}}$ | [H/m] | Carson earth-return inductance | from connected `Carson` models

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
\mathbf{0} &= -\mathbf{R}^{\mathrm{int}}
  + \operatorname{diag}(\mathbf{r}^{\mathrm{skin}}) \\
\mathbf{0} &= -\mathbf{L}^{\mathrm{int}}
  + \operatorname{diag}(\mathbf{l}^{\mathrm{skin}}) \\
\mathbf{0} &= -\mathbf{R}^{\mathrm{ext}}
  + \mathbf{R}^{\mathrm{carson}} \\
\mathbf{0} &= -\mathbf{L}^{\mathrm{ext}}
  + \mathbf{L}^{\mathrm{geo}}
  + \mathbf{L}^{\mathrm{carson}} \\
\mathbf{0} &= -\mathbf{R}'
  + \mathbf{R}^{\mathrm{int}}
  + \mathbf{R}^{\mathrm{ext}} \\
\mathbf{0} &= -\mathbf{L}'
  + \mathbf{L}^{\mathrm{int}}
  + \mathbf{L}^{\mathrm{ext}}
\end{aligned}
```

## Initialization

None.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}'$ | [$\Omega$/m] | Full-conductor series resistance | $\mathbb{R}^{K\times K}$
$\mathbf{L}'$ | [H/m] | Full-conductor series inductance | $\mathbb{R}^{K\times K}$
