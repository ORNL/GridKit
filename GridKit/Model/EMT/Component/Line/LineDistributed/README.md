# LineDistributed Model

`LineDistributed` represents an $N$-phase, $K$-conductor distributed EMT line.
Each Bus owns its characteristic-admittance current and states. The line
owns the reflected outputs, propagation states, and histories.

## Block Diagram

![LineDistributed model block diagram](../../../../../../docs/Figures/EMT/LineDistributed/diagram.png)

Figure 1: LineDistributed terminal interconnection

The conductor-to-phase mappings are shown in the equations and omitted from
the diagram for clarity.

## Model Parameters

Define the phase- and conductor-index sets

```math
\mathcal{N} = \{1,\ldots,N\},
\qquad
\mathcal{K} = \{1,\ldots,K\}.
```

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{c}$ | [-] | `conductors` | Conductor phase-index list | $\mathbf{c} \in \mathcal{N}^K$

### Parameter Validation

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
K &\in \mathbb{Z}_{>0} \\
\mathbf{c} &\in \mathcal{N}^K \\
\{c_k \mid k \in \mathcal{K}\} &= \mathcal{N}
\end{aligned}
```

### Derived Parameters

```math
P_{\phi,nk} =
\begin{cases}
1, & n = c_k \\
0, & n \ne c_k
\end{cases},
\quad n \in \mathcal{N},\quad k \in \mathcal{K}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$ | `v1` | Input | [V] | Terminal 1 voltage connection used to assemble its Norton source | $\mathbb{R}^N$
$\mathbf{v}_2$ | `v2` | Input | [V] | Terminal 2 voltage connection used to assemble its Norton source | $\mathbb{R}^N$
$\mathbf{i}_1^\mathrm{c}$ | `Ish1` | Input | [A] | Characteristic-admittance current from terminal 1 Bus | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{c}$ | `Ish2` | Input | [A] | Characteristic-admittance current from terminal 2 Bus | $\mathbb{R}^K$
$\mathbf{i}_1^\mathrm{ref}$ | `i_ref1` | Output | [A] | Reflected current at terminal 1 | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{ref}$ | `i_ref2` | Output | [A] | Reflected current at terminal 2 | $\mathbb{R}^K$
$\mathbf{i}_1^\mathrm{inc}$ | `i_inc1` | Output | [A] | Incident current supplied to terminal 1 Bus | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{inc}$ | `i_inc2` | Output | [A] | Incident current supplied to terminal 2 Bus | $\mathbb{R}^K$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{y}_1^\mathrm{c}$ | Characteristic admittance at terminal 1 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}^\mathrm{c}}$ | `Yc` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{y}_2^\mathrm{c}$ | Characteristic admittance at terminal 2 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}^\mathrm{c}}$ | `Yc` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{h}_{21}$ | Propagation from terminal 2 to terminal 1 | [Propagation](../../../Operators/Shift/Propagation/README.md) | Composite | `H` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{h}_{12}$ | Propagation from terminal 1 to terminal 2 | [Propagation](../../../Operators/Shift/Propagation/README.md) | Composite | `H` | $\mathbb{R}^K$ | $\mathbb{R}^K$

`Yc` and `H` each provide one coefficient set. The model assumes a reciprocal
uniform line, so

```math
\mathbf{Y}_1^{\mathrm{c}}(s)=\mathbf{Y}_2^{\mathrm{c}}(s),
\qquad
\mathbf{H}_{12}(s)=\mathbf{H}_{21}(s).
```

The terminal admittances have independent states in their buses;
the propagation instances have independent states and histories in the line.
In case files, `bus1` and `bus2` select the terminal buses using the usual
scalar voltage aliases. Assembly creates their Norton sources and connects
the characteristic and incident currents shown above.

### Submodel Validation

The characteristic-admittance fits must be stable, proper, and positive real.
Together with the propagation fits, they must produce a passive line model.
The runtime validates dimensions, properness, and stable poles. Positive
realness and passivity of the combined fitted line are offline fitting
requirements; pole stability alone does not establish either property.

### Submodel Wiring

```math
\begin{aligned}
\mathbf{i}_1^\mathrm{inc} &\leftarrow
  \mathbf{h}_{21}[\mathbf{i}_2^\mathrm{ref}] \\
\mathbf{i}_2^\mathrm{inc} &\leftarrow
  \mathbf{h}_{12}[\mathbf{i}_1^\mathrm{ref}]
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_1^\mathrm{ref}$ | [A] | Reflected current at terminal 1 | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{ref}$ | [A] | Reflected current at terminal 2 | $\mathbb{R}^K$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_1^\mathrm{c}$ | [A] | Characteristic-admittance current owned by terminal 1 Norton source | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{c}$ | [A] | Characteristic-admittance current owned by terminal 2 Norton source | $\mathbb{R}^K$
$\mathbf{i}_1^\mathrm{inc}$ | [A] | Output of the terminal 2-to-1 propagation operator | $\mathbb{R}^K$
$\mathbf{i}_2^\mathrm{inc}$ | [A] | Output of the terminal 1-to-2 propagation operator | $\mathbb{R}^K$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = -\mathbf{i}_e^\mathrm{ref}+2\mathbf{i}_e^\mathrm{c}-\mathbf{i}_e^\mathrm{inc},
\qquad e\in\{1,2\}
```

### External Equations

None. Each bus registers its incident-current output
$\mathbf{i}_e^\mathrm{inc}$ with positive sign and its characteristic-admittance
current $\mathbf{i}_e^\mathrm{c}$ with negative sign. The net injection shown
in the diagram is

```math
\Delta\mathbf{i}_e = \mathbf{P}_\phi
  (\mathbf{i}_e^\mathrm{inc}-\mathbf{i}_e^\mathrm{c}),\quad e\in\{1,2\}.
```

## Initialization

Supply each terminal's reflected-current prehistory in the state file's
[`history`](../../../STATE.md#history) section. The history contains an
angular frequency, instantaneous current values, and their time derivatives.
It is constant when $\omega=0$ and harmonic otherwise:

```math
\mathbf{i}_e^{\mathrm{ref}}(t)
  =\mathbf{i}_e^{\mathrm{ref}}(t_0)\cos\!\left(\omega(t-t_0)\right)
   +\dfrac{1}{\omega}\left.\dfrac{\mathrm{d}\mathbf{i}_e^{\mathrm{ref}}}{\mathrm{d}t}\right|_{t_0}
      \sin\!\left(\omega(t-t_0)\right).
```

Each propagation instance initializes its rational states and filtered
history from its local reflected-current prehistory. The opposite terminal
receives its delayed output. A zero history represents energization from an
unenergized line. Initial reflected outputs, if also specified in `devices`,
must agree with the history endpoint. IDA then reconciles the network's
algebraic variables and differential-state derivatives.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i_c1a`, `i_c1b`, `i_c1c` | [A] | Characteristic-admittance current at terminal 1 | $\mathbf{i}_1^\mathrm{c} \in \mathbb{R}^K$
`i_c2a`, `i_c2b`, `i_c2c` | [A] | Characteristic-admittance current at terminal 2 | $\mathbf{i}_2^\mathrm{c} \in \mathbb{R}^K$
`i_inc1a`, `i_inc1b`, `i_inc1c` | [A] | Incident current at terminal 1 | $\mathbf{i}_1^\mathrm{inc} \in \mathbb{R}^K$
`i_inc2a`, `i_inc2b`, `i_inc2c` | [A] | Incident current at terminal 2 | $\mathbf{i}_2^\mathrm{inc} \in \mathbb{R}^K$
`i_ref1a`, `i_ref1b`, `i_ref1c` | [A] | Reflected current at terminal 1 | $\mathbf{i}_1^\mathrm{ref} \in \mathbb{R}^K$
`i_ref2a`, `i_ref2b`, `i_ref2c` | [A] | Reflected current at terminal 2 | $\mathbf{i}_2^\mathrm{ref} \in \mathbb{R}^K$

The three-phase implementation exposes scalar monitor names by appending
`a`, `b`, or `c`, for example `i_inc1a`. Reflected and incident current
outputs use the same phase suffixes.

## Development

The implemented three-phase formulation requires $N=K=3$ and
$\mathbf{c}=[1,2,3]^\mathsf{T}$, hence $\mathbf{P}_\phi=\mathbf{I}_3$.
The generalized conductor mapping above documents the extension; other
dimensions are rejected by this implementation, as for `LineLumped`.
