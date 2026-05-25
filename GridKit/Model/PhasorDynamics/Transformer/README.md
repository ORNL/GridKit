# **Transformer Model**

The transformer model represents an ideal off-nominal transformer ratio in
series with a branch equivalent circuit. Positive terminal current
contributions to bus residuals are oriented leaving the transformer and
entering the adjacent buses.

## Block Diagram

## Model Parameters

Symbol      | Units   | Description                     | Typical Value | Note
------------|---------|---------------------------------|---------------|------
$R$         | [p.u.]  | Series resistance               |               |
$X$         | [p.u.]  | Series reactance                |               |
$G$         | [p.u.]  | Total shunt conductance         |               |
$B$         | [p.u.]  | Total shunt susceptance         |               |
$a$         | [p.u.]  | Off-nominal transformer ratio   | 1             |
$\theta$    | [rad]   | LTC phase shift                 | 0             |

### Parameter Validation

### Model Derived Parameters

```math
\begin{aligned}
  m &= \dfrac{1}{a} \\
  G_{br} &= \dfrac{R}{R^2 + X^2} \\
  B_{br} &= -\dfrac{X}{R^2 + X^2} \\
  \mathbf{M} &=
    \begin{bmatrix}
      m & 0 \\
      0 & e^{j\theta}
    \end{bmatrix} \\
  \mathbf{Y} &=
    \begin{bmatrix}
      -Y_{br} & Y_{br} \\
      Y_{br} & -Y_{br}
    \end{bmatrix}
    +
    \dfrac{1}{2}
    \begin{bmatrix}
      -Y_{sh} & 0 \\
      0 & -Y_{sh}
    \end{bmatrix} \\
  \mathbf{Y}' &= \mathbf{M}^{\dagger}\mathbf{Y}\mathbf{M}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
-------|-------|-------------|-----
$I_r^{\mathrm{from}}$ | [p.u.] | Current injection, real component, from bus | Positive into the from bus
$I_i^{\mathrm{from}}$ | [p.u.] | Current injection, imaginary component, from bus | Positive into the from bus
$I_r^{\mathrm{to}}$ | [p.u.] | Current injection, real component, to bus | Positive into the to bus
$I_i^{\mathrm{to}}$ | [p.u.] | Current injection, imaginary component, to bus | Positive into the to bus

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
-------|-------|-------------|-----
$V_r^{\mathrm{from}}$ | [p.u.] | Terminal voltage, real component, from bus | Owned by bus object
$V_i^{\mathrm{from}}$ | [p.u.] | Terminal voltage, imaginary component, from bus | Owned by bus object
$V_r^{\mathrm{to}}$ | [p.u.] | Terminal voltage, real component, to bus | Owned by bus object
$V_i^{\mathrm{to}}$ | [p.u.] | Terminal voltage, imaginary component, to bus | Owned by bus object

## Model Equations

### Differential Equations

None.

### Algebraic Equations

The transformer model can be written in phasor form as
$\mathbf{0} = -\mathbf{I} + \mathbf{Y}'\mathbf{V}$, or expanded as:

```math
\begin{aligned}
  0 &= -I_r^{\mathrm{from}}
       &+ G'_{ff}V_r^{\mathrm{from}} &- B'_{ff}V_i^{\mathrm{from}}
       &+ G'_{ft}V_r^{\mathrm{to}}   &- B'_{ft}V_i^{\mathrm{to}} \\
  0 &= -I_i^{\mathrm{from}}
       &+ B'_{ff}V_r^{\mathrm{from}} &+ G'_{ff}V_i^{\mathrm{from}}
       &+ B'_{ft}V_r^{\mathrm{to}}   &+ G'_{ft}V_i^{\mathrm{to}} \\
  0 &= -I_r^{\mathrm{to}}
       &+ G'_{tf}V_r^{\mathrm{from}} &- B'_{tf}V_i^{\mathrm{from}}
       &+ G'_{tt}V_r^{\mathrm{to}}   &- B'_{tt}V_i^{\mathrm{to}} \\
  0 &= -I_i^{\mathrm{to}}
       &+ B'_{tf}V_r^{\mathrm{from}} &+ G'_{tf}V_i^{\mathrm{from}}
       &+ B'_{tt}V_r^{\mathrm{to}}   &+ G'_{tt}V_i^{\mathrm{to}}
\end{aligned}
```

## Network Interface

## Initialization

## Model Outputs

Output | Units | Description | Note
-------|-------|-------------|-----
