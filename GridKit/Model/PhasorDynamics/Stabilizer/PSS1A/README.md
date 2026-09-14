# PSS1A

Specification draft; not implemented. The equations and initialization remain unverified.

## Block Diagram

![](../../../../../docs/Figures/PSS1A.JPG)

Figure 1: Power system stabilizer PSS1A model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol | Units | JSON | Description | Typical Value | Note
-------|-------|------|-------------|---------------|-----
$I_{\mathrm{cs}}$ | [-] | TBD | Stabilizer input code | 2 |
$A_1$ | [s] | TBD | Notch denominator coefficient | 0 |
$A_2$ | [s²] | TBD | Notch denominator coefficient | 0 |
$T_1$ | [s] | TBD | Lead–lag 1 numerator time constant | 0.25 |
$T_2$ | [s] | TBD | Lead–lag 1 denominator time constant | 0.03 |
$T_3$ | [s] | TBD | Lead–lag 2 numerator time constant | 0.25 |
$T_4$ | [s] | TBD | Lead–lag 2 denominator time constant | 0.03 |
$T_5$ | [s] | TBD | Washout numerator time constant | 20 |
$T_6$ | [s] | TBD | Transducer time constant | 0.02 |
$K_\mathrm{s}$ | [p.u.] | TBD | Stabilizer gain | 10 |
$L_\mathrm{s}^{\max}$ | [p.u.] | TBD | Maximum stabilizer output | 0.1 |
$L_\mathrm{s}^{\min}$ | [p.u.] | TBD | Minimum stabilizer output | -0.1 |
$V_\mathrm{cu}$ | [p.u.] | TBD | Upper cutout threshold | 0 |
$V_\mathrm{cl}$ | [p.u.] | TBD | Lower cutout threshold | 0 |

### Parameter Validation

TBD.

### Model Derived Parameters

TBD.

## Model Ports

Name     | Port   | Init | Description
---------|--------|------|------------------------------------------
`input`  | Input  | TBD  | Stabilizer input $u$ selected by $I_{\mathrm{cs}}$
`vct`    | Input  | TBD  | Cutout signal $V_{\mathrm{ct}}$
`output` | Output | TBD  | Limited stabilizer output $V_{\mathrm{ss}}$

## Model Variables

These were the variables listed in the old documentation.

1. rotor speed deviation (p.u.)
2. bus frequency deviation (p.u.) - default
3. generator electrical power in Gen MVA Base (p.u.)
4. generator accelerating power (p.u.)
5. bus voltage (p.u.)
6. derivative of p.u. bus voltage

### Internal Variables

#### Differential

TBD.
#### Algebraic

TBD.
### External Variables

#### Differential

None.

#### Algebraic

Symbol   | Units  | Description                                 | Note
---------|--------|---------------------------------------------|-----------------------
$u$      | [p.u.] | Stabilizer input signal                     |
$V_{\mathrm{ct}}$ | [p.u.] | Cutout signal (compared to $V_\mathrm{cl},V_\mathrm{cu}$) | from the block diagram

## Model Equations

### Internal Equations

#### Differential

```math
\begin{aligned}
\dot{V_{1}} &= \dfrac{1}{ T_{6} }( V_{SI}-V_{1} ) \\
\dot{x_{1}} &= -\dfrac{ V_{2} }{ T_{5} } \\
\dfrac{d^{2}V_{3}}{dt^{2}}+\dfrac{A_{1}}{A_{2}}\dfrac{dV_{3}}{dt}&=\dfrac{1}{A_{2}}(V_{2}-1) \\
\dfrac{dx_{2}}{dt}&=\dfrac{1}{T_{2}}(V_{3}-V_{4}) \\
\dfrac{dx_{3}}{dt}&=\dfrac{1}{T_{4}}(V_{4}-V_{5})
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
V_{2} &= x_{1} + K_\mathrm{s} V_{1} \\
V_{4}&=x_{2}+\dfrac{T_{1}}{T_{2}}V_{3} \\
V_{5}&=x_{3}+\dfrac{T_{3}}{T_{4}}V_{4} \\
V_{\mathrm{ss}} &= \begin{cases}
   L_\mathrm{s}^{\max} &\text{if } V_{5}>L_\mathrm{s}^{\max} \\
   L_\mathrm{s}^{\min} &\text{if } V_{5}<L_\mathrm{s}^{\min} \\
   V_{5}
\end{cases}
\end{aligned}
```

### External Equations

None.

## Initialization

TBD.

## Monitors

TBD.
