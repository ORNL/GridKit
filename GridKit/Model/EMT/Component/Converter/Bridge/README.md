# Bridge Model

`Bridge` represents a two-level converter with $N$ legs, two DC terminals, and
an LCL filter to the EMT bus. Switching functions set the pole voltages and DC
current injections. Grid-side currents $\mathbf{i}_g$ are injected from the
bridge into the EMT bus.

## Block Diagram

![Bridge model diagram](../../../../../../docs/Figures/EMT/Bridge/diagram.png)

Figure 1: Representative three-phase bridge with binary contact positions;
the equations apply to $N$ legs and define continuous interpolation

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of legs | Required, positive integer
$R_\mathrm{s}$ | $[\Omega]$ | `Rs` | Converter-side resistance | Nonnegative
$L_\mathrm{s}$ | [H] | `Ls` | Converter-side inductance | Required, positive
$C$ | [F] | `C` | Filter capacitance | Required, positive
$R_g$ | $[\Omega]$ | `Rg` | Grid-side resistance | Nonnegative
$L_g$ | [H] | `Lg` | Grid-side inductance | Required, positive

### Parameter Validation

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
R_\mathrm{s} &\ge 0 \\
L_\mathrm{s} &> 0 \\
C &> 0 \\
R_g &\ge 0 \\
L_g &> 0
\end{aligned}
```

### Model Derived Parameters

Define the leg-index set

```math
\mathcal{N} = \{1,\ldots,N\}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{s}$ | `s` | Input | [-] | Switching functions | $\mathbf{s} \in \mathbb{R}^N$
$\mathbf{v}$, $\mathbf{i}_g$ | `bus` | Bus | [V], [A] | Reads the bus voltage and injects the grid-side current | $\mathbf{v}, \mathbf{i}_g \in \mathbb{R}^N$
$v_+$, $i_+$ | `vpos` | Bus | [V], [A] | Positive DC terminal voltage and current injection | Scalar
$v_-$, $i_-$ | `vneg` | Bus | [V], [A] | Negative DC terminal voltage and current injection | Scalar
$\mathbf{i}$ | `i` | Output | [A] | Converter-side currents | $\mathbf{i} \in \mathbb{R}^N$
$\mathbf{v}_\mathrm{o}$ | `vo` | Output | [V] | Capacitor voltages | $\mathbf{v}_\mathrm{o} \in \mathbb{R}^N$
$\mathbf{i}_g$ | `ig` | Output | [A] | Grid-side currents | $\mathbf{i}_g \in \mathbb{R}^N$

All currents are positive from the bridge toward the connected buses. The
capacitor voltages share the bus voltage reference. Each switching function
satisfies $0 \le s_n \le 1$ and may come from
[PWM](../../../Operators/Modulation/PWM/README.md).

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Converter-side currents | $\mathbf{i} \in \mathbb{R}^N$
$\mathbf{v}_\mathrm{o}$ | [V] | Capacitor voltages | $\mathbf{v}_\mathrm{o} \in \mathbb{R}^N$
$\mathbf{i}_g$ | [A] | Grid-side currents | $\mathbf{i}_g \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^N$
$v_+$ | [V] | Positive DC terminal voltage owned by EMT bus |
$v_-$ | [V] | Negative DC terminal voltage owned by EMT bus |

#### Algebraic

None.

## Model Equations

The pole voltages are

```math
e_n = (1-s_n)v_- + s_n v_+,
\quad n \in \mathcal{N}
```

### Internal Equations

#### Differential

For each leg $n \in \mathcal{N}$,

```math
\begin{aligned}
0 &= L_\mathrm{s}\dfrac{\mathrm{d}i_n}{\mathrm{d}t}
     + R_\mathrm{s}i_n + v_{\mathrm{o},n} - e_n \\
0 &= C\dfrac{\mathrm{d}v_{\mathrm{o},n}}{\mathrm{d}t}
     + i_{g,n} - i_n \\
0 &= L_g\dfrac{\mathrm{d}i_{g,n}}{\mathrm{d}t}
     + R_g i_{g,n} + v_n - v_{\mathrm{o},n}
\end{aligned}
```

#### Algebraic

None.

### External Equations

The bridge supplies input $d$ of the [bus current balance](../../../Bus/README.md#model-equations):

```math
\mathbf{i}_d \leftarrow \mathbf{i}_g
```

The DC current injections are

```math
\begin{aligned}
i_+ &\leftarrow -\sum_{n \in \mathcal{N}} s_n i_n \\
i_- &\leftarrow -\sum_{n \in \mathcal{N}} (1-s_n)i_n
\end{aligned}
```

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Converter-side currents | $\mathbf{i} \in \mathbb{R}^N$
`vo` | [V] | Capacitor voltages | $\mathbf{v}_\mathrm{o} \in \mathbb{R}^N$
`ig` | [A] | Grid-side currents | $\mathbf{i}_g \in \mathbb{R}^N$
`ipos` | [A] | Positive DC current injection | Scalar
`ineg` | [A] | Negative DC current injection | Scalar
`e` | [V] | Pole voltages | $\mathbf{e} \in \mathbb{R}^N$
