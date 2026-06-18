# InsulatorLeakage Model

`InsulatorLeakage` computes an empirical direct shunt conductance effect for
the physical conductors.

## Model Parameters

For $K$ physical conductors:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$g_i^{\mathrm{leak}}$ | [S/m] | `leakage_conductance` | Insulator leakage conductance to ground | $g_i^{\mathrm{leak}}\ge0$

### Parameter Validation

```math
g_i^{\mathrm{leak}}\ge0,\qquad i=1,\dots,K
```

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
$\mathbf{G}^{\mathrm{leak}}$ | [S/m] | Insulator leakage shunt conductance | diagonal

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
\mathbf{0} = -\mathbf{G}^{\mathrm{leak}}
  + \operatorname{diag}(\mathbf{g}^{\mathrm{leak}})
```

## Initialization

Validate the leakage conductance vector.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}^{\mathrm{leak}}$ | [S/m] | Insulator leakage shunt conductance | $\mathbb{R}^{K\times K}$
