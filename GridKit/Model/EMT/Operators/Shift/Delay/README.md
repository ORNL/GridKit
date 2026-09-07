# Delay Model

`Delay` applies a constant transport delay to each of $M$ input channels:

```math
\mathbf{D}_{\boldsymbol{\tau}}(s)
  = \mathrm{diag}\left(\exp(-s\tau_1),\ldots,\exp(-s\tau_M)\right).
```

At runtime, accepted-step input samples are reconstructed with cubic Hermite
interpolation. Each channel owns an algebraic delayed-output variable. When a
step extends beyond its delay, the delayed value depends implicitly on the
current trial input. A scalar delay is the $M=1$ case.

## Block Diagram

![Delay operator block diagram](../../../../../../docs/Figures/EMT/Delay/diagram.png)

Figure 1: Delay model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$M$ | [-] | `M` | Channel count | Required, positive integer
$\boldsymbol{\tau}$ | [s] | `tau` | Channel delays | Required, each positive

### Parameter Validation

```math
\begin{aligned}
M &\in \mathbb{Z}_{>0} \\
\tau_m &> 0
\end{aligned}
```

### Derived Parameters

```math
\tau_{\min} = \min(\boldsymbol{\tau}),
\qquad
\tau_{\max} = \max(\boldsymbol{\tau})
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^M$
$\mathbf{y}$ | `out` | Output | $[u]$ | Delayed output port | $\mathbf{y} \in \mathbb{R}^M$

The output is a bound DAE signal, with its derivative supplied by the solver.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

History samples are implementation data. Only the delayed outputs occupy DAE
variables and residual rows.

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | $[u]$ | Delayed output | $\mathbf{y} \in \mathbb{R}^M$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | $\mathbf{u} \in \mathbb{R}^M$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0=-y_m(t)+u_m(t-\tau_m),\qquad m\in\{1,\ldots,M\}
```

### External Equations

None.

## History Realization

### History Record

Accepted input history is stored as the knot sequence

```math
(t_j,\ \mathbf{u}_j,\ \mathbf{u}'_j),
\qquad t_0 \le t_1 \le \cdots \le t_n,
```

where $\mathbf{u}_j$ and $\mathbf{u}'_j$ are the input value and derivative at
$t_j$. The first knot holds the initialized values. Channels share knot times
and read the record at $t-\tau_m$ independently.

![Delay history record and per-channel taps](../../../../../../docs/Figures/EMT/DelayHistory/diagram.png)

Figure 2: Delay history record and channel taps

Later knots are appended only at accepted steps. Rejected trials and monitor
samples never enter the record. At a discontinuity, two knots retain the left
and right limits at the same time. Old knots are removed after retaining the
bracket needed by the longest delay.

### Interpolation

Suppressing the channel index, a lookup at $\xi < t_0$ uses the analytic
prehistory. For $\xi \in (t_j,t_{j+1}]$, the bracketing knots define the cubic
Hermite interpolant

```math
\begin{aligned}
h_j &= t_{j+1}-t_j,
\qquad
\theta = \dfrac{\xi-t_j}{h_j} \\
u(\xi)
  &= (1-\theta)^2(1+2\theta)\,u_j
   + \theta(1-\theta)^2\,h_j\,u'_j \\
  &\quad
   + \theta^2(3-2\theta)\,u_{j+1}
   - \theta^2(1-\theta)\,h_j\,u'_{j+1},
\end{aligned}
```

For smooth input and exact knot data, the interpolant is $C^1$ with nominal
fourth-order value accuracy in the knot spacing. The first segment after a
restart uses linear interpolation because consistent-condition calculation
does not update algebraic input derivatives.

### Steps Longer Than a Delay

Let $t_n$ be the accepted frontier and $t=t_n+h$ a trial endpoint. If
$\xi=t-\tau_m>t_n$, define $\theta=(\xi-t_n)/h$. The overlap polynomial uses
the accepted value and slope and the current trial value:

```math
u_m(\xi)
  \approx (1-\theta^2)u_{m,n}
    +h\theta(1-\theta)u'_{m,n}
    +\theta^2u_m(t).
```

The $\theta^2$ coefficient participates in the sparse Jacobian, including
the chain rule for computed input signals. It does not introduce a trial
input-derivative dependency. Immediately after a restart the overlap uses
linear interpolation, with coefficient $\theta$.

The solver remains adaptive. `limit_step: true` requests the reference
method-of-steps bound $h\le\tau_{\min}$; the default is `false`. This option
is also accepted by `Propagation` and applies to its delay bank. Unrestricted
overlap has a third-order local value error for smooth exact input data;
IDA's error test alone is not an independent bound on history reconstruction
error. Compare refined tolerances and the bounded-step realization when
assessing transient accuracy.

### Discontinuities

An accepted input jump schedules its arrival at $t+\tau_m$. IDA stops at the
next arrival, retains the left limit, then calculates consistent right-limit
conditions and restarts. Repeated arrivals are propagated through subsequent
reflections. A relative jump threshold of $10^{-10}$ suppresses roundoff-level
restart cascades. Smooth dispersive tails continue under the usual adaptive
error control.

## Initialization

A prehistory for $\mathbf{u}(t)$ and
$\mathrm{d}\mathbf{u}/\mathrm{d}t$ must be specified over

```math
t \in [t_0-\tau_{\max},t_0].
```

The delay does not synthesize prehistory. A supplied endpoint differing from
the consistent initial input represents a switching event at $t_0$; both
limits are retained. At $t_0$,

```math
\begin{aligned}
\mathbf{y}(t_0)
  &\leftarrow \mathbf{u}(t_0-\boldsymbol{\tau}) \\
\left.\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}\right|_{t_0}
  &\leftarrow
  \left.\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  \right|_{t_0-\boldsymbol{\tau}}.
\end{aligned}
```

## Monitors

None.
