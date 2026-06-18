# ShuntAdmittance Model

`ShuntAdmittance` builds the full-conductor per-unit-length shunt conductance
and capacitance matrices from connected shunt-effect outputs. Potential-derived
conductance is zero; direct shunt losses are represented explicitly and set to
zero until a leakage model is wired.

## Model Parameters

None. `ShuntAdmittance` combines effect-model outputs for $K$ physical
conductors.

### Parameter Validation

None.

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}^{\mathrm{ext}}$ | [S/m] | Direct shunt conductance effects | $\mathbb{R}^{K\times K}$
$\mathbf{C}^{\mathrm{ext}}$ | [F/m] | Direct shunt capacitance effects | $\mathbb{R}^{K\times K}$
$\mathbf{G}'$ | [S/m] | Full-conductor shunt conductance | $\mathbb{R}^{K\times K}$
$\mathbf{C}'$ | [F/m] | Full-conductor shunt capacitance | $\mathbb{R}^{K\times K}$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}^{\mathrm{pot}}$ | [S/m] | Potential-derived shunt conductance | zero matrix from connected `ShuntPotential` model
$\mathbf{C}^{\mathrm{pot}}$ | [F/m] | Potential-derived shunt capacitance | from connected `ShuntPotential` model

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
\mathbf{0} &= -\mathbf{G}^{\mathrm{ext}} \\
\mathbf{0} &= -\mathbf{C}^{\mathrm{ext}} \\
\mathbf{0} &= -\mathbf{G}'
  + \mathbf{G}^{\mathrm{pot}}
  + \mathbf{G}^{\mathrm{ext}} \\
\mathbf{0} &= -\mathbf{C}'
  + \mathbf{C}^{\mathrm{pot}}
  + \mathbf{C}^{\mathrm{ext}}
\end{aligned}
```

## Initialization

The direct shunt-effect terms are initialized to zero. The assembled matrices
are initialized from the connected `ShuntPotential` output:

```math
\begin{aligned}
\mathbf{G}^{\mathrm{ext}} &= \mathbf{0} \\
\mathbf{C}^{\mathrm{ext}} &= \mathbf{0} \\
\mathbf{G}' &= \mathbf{G}^{\mathrm{pot}} \\
\mathbf{C}' &= \mathbf{C}^{\mathrm{pot}}
\end{aligned}
```

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}'$ | [S/m] | Full-conductor shunt conductance | $\mathbb{R}^{K\times K}$
$\mathbf{C}'$ | [F/m] | Full-conductor shunt capacitance | $\mathbb{R}^{K\times K}$
