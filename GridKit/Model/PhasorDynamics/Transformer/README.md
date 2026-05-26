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
  Y_{\text{br}} &= \dfrac{1}{R + jX} \\
  Y_{\text{sh}} &= G + jB \\
  \mathbf{Y} &=
    \begin{bmatrix}
      -Y_{\text{br}} & Y_{\text{br}} \\
      Y_{\text{br}} & -Y_{\text{br}}
    \end{bmatrix}
    +
    \dfrac{1}{2}
    \begin{bmatrix}
      -Y_{\text{sh}} & 0 \\
      0 & -Y_{\text{sh}}
    \end{bmatrix} \\
  \mathbf{M} &=
    \begin{bmatrix}
      a^{-1} & 0 \\
      0 & e^{j\theta}
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
$I_{\text{r}}^{\mathrm{from}}$ | [p.u.] | Current injection, real component, from bus | Positive into the from bus
$I_{\text{i}}^{\mathrm{from}}$ | [p.u.] | Current injection, imaginary component, from bus | Positive into the from bus
$I_{\text{r}}^{\mathrm{to}}$ | [p.u.] | Current injection, real component, to bus | Positive into the to bus
$I_{\text{i}}^{\mathrm{to}}$ | [p.u.] | Current injection, imaginary component, to bus | Positive into the to bus

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
-------|-------|-------------|-----
$V_{\text{r}}^{\mathrm{from}}$ | [p.u.] | Terminal voltage, real component, from bus | Owned by bus object
$V_{\text{i}}^{\mathrm{from}}$ | [p.u.] | Terminal voltage, imaginary component, from bus | Owned by bus object
$V_{\text{r}}^{\mathrm{to}}$ | [p.u.] | Terminal voltage, real component, to bus | Owned by bus object
$V_{\text{i}}^{\mathrm{to}}$ | [p.u.] | Terminal voltage, imaginary component, to bus | Owned by bus object

## Model Equations

### Differential Equations

None.

### Algebraic Equations

The transformer model can be written in phasor form as
$\mathbf{0} = -\mathbf{I} + \mathbf{Y}'\mathbf{V}$, or expanded as:

```math
\begin{aligned}
  0 &= -I_{\text{r}}^{\mathrm{from}}
       &+ G'_{\text{ff}}V_r^{\mathrm{from}} &- B'_{\text{ff}}V_{\text{i}}^{\mathrm{from}}
       &+ G'_{\text{ft}}V_r^{\mathrm{to}}   &- B'_{\text{ft}}V_{\text{i}}^{\mathrm{to}} \\
  0 &= -I_{\text{i}}^{\mathrm{from}}
       &+ B'_{\text{ff}}V_r^{\mathrm{from}} &+ G'_{\text{ff}}V_{\text{i}}^{\mathrm{from}}
       &+ B'_{\text{ft}}V_r^{\mathrm{to}}   &+ G'_{\text{ft}}V_{\text{i}}^{\mathrm{to}} \\
  0 &= -I_{\text{r}}^{\mathrm{to}}
       &+ G'_{\text{tf}}V_r^{\mathrm{from}} &- B'_{\text{tf}}V_{\text{i}}^{\mathrm{from}}
       &+ G'_{\text{tt}}V_r^{\mathrm{to}}   &- B'_{\text{tt}}V_{\text{i}}^{\mathrm{to}} \\
  0 &= -I_{\text{i}}^{\mathrm{to}}
       &+ B'_{\text{tf}}V_r^{\mathrm{from}} &+ G'_{\text{tf}}V_{\text{i}}^{\mathrm{from}}
       &+ B'_{\text{tt}}V_r^{\mathrm{to}}   &+ G'_{\text{tt}}V_{\text{i}}^{\mathrm{to}}
\end{aligned}
```

## Initialization

## Model Outputs

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$I_{\text{r}}^{\mathrm{from}}$ | [p.u.] | Current injection, real component, from bus | Oriented entering the from bus
$I_{\text{i}}^{\mathrm{from}}$ | [p.u.] | Current injection, imaginary component, from bus | Oriented entering the from bus
$I_{\text{r}}^{\mathrm{to}}$   | [p.u.] | Current injection, real component, to bus | Oriented entering the to bus
$I_{\text{i}}^{\mathrm{to}}$   | [p.u.] | Current injection, imaginary component, to bus | Oriented entering the to bus
