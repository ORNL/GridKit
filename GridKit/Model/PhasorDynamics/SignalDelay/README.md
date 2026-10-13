# Delayed signal

A delayed signal publishes a scalar source at an earlier time. Multiple delays
of the same source share one history. Smooth source signals and constant delays
are assumed.

## Model Parameters

| Name | Description | Unit |
| --- | --- | --- |
| $\tau$ | Delay | s |

### Parameter Validation

The component requires finite $\tau>0$. A signal declaration permits $\tau=0$,
which the composer implements as an alias without an additional variable.
Source IDs must exist. Delayed signals cannot also have device producers.
Declaration cycles are invalid; feedback through devices is supported.

### Model Derived Parameters

Chains of delay declarations resolve to their original source and the sum of
their delays. A source history retains the largest delay requested by its taps.

## Model Ports

| Port | Direction | Description |
| --- | --- | --- |
| $u$ | Input | Source signal |
| $z$ | Output | Delayed source signal |

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

$z$: delayed output, in the source's units.

### External Variables

#### Differential

$u$, when the connected source is a differential variable.

#### Algebraic

$u$, when the connected source is algebraic. A prescribed time function or
constant source has no corresponding system variable.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

$$F_z(t)=-z(t)+u(t-\tau)=0.$$

This is a history-dependent equation. Let $t_n$ be the latest accepted time and
$q=t-\tau$. Queries before the initial time use prehistory $\phi(q)$. Queries
through $t_n$ use accepted history. Queries inside the current trial step use
the current source iterate:

$$\widehat u(q)=a(q,t)+b(q,t)u(t).$$

For up to three recent accepted anchors $(s_i,u_i)$, ending at $s_{m-1}=t_n$,

$$b(q,t)=\prod_{i=0}^{m-1}\frac{q-s_i}{t-s_i},$$

$$a(q,t)=\sum_{i=0}^{m-1}u_i\frac{q-t}{s_i-t}
         \prod_{j\ne i}\frac{q-s_j}{s_i-s_j}.$$

Startup uses the available anchors. The Jacobian entries are
$\partial F_z/\partial z=-1$ and $\partial F_z/\partial u=b$.
For a query in accepted history, $b=0$. The source column stays in the sparse
pattern as the query crosses the accepted-history boundary.

With one anchor and $h=t-t_n>\tau$, this reduces to

$$z(t)=\frac{\tau}{h}u_n+\left(1-\frac{\tau}{h}\right)u(t).$$

Thus steps may exceed a delay. Integration error, interpolation error and
nonlinear convergence determine the usable step size.

### External Equations

None. Consumers read the ordinary output variable $z$.

## Initialization

### Input Initialization

By default, delayed connections behave as ordinary connections during component
initialization, including initialization writes. The consistency solve uses
$z=u$. Its final source value defines constant prehistory $\phi(s)=u(t_0)$.

For a nonconstant past, set a source-level callback before configuring the solver:

```cpp
system.setSignalPrehistory(source_id,
    [](double t) { return std::cos(t); });
```

The callback must cover all requested times before $t_0$ and agree with the
initial source value within history tolerance. Its smoothness is the caller's
responsibility.

### Internal Initialization

$$z(t_0)=\phi(t_0-\tau).$$

A fresh model initialization clears stored samples. Solver retries and
continuation at the accepted time preserve history. Restarting from a different
time requires a fresh model initialization and a corresponding prehistory.

### Output Initialization

Link $z$ to the declared delayed signal node.

## Monitors

The optional signal declaration field `"mon": ["z"]` selects the delayed output.

## Numerical History

Each accepted interval stores the same interpolation polynomial used in its
solve, represented by four cubic Bernstein coefficients. Recent interpolation
anchors are stored separately. History outside the largest requested delay is
discarded after output sampling.

Before accepting a step, compare its history polynomial with a lower-degree
interpolant and source samples from IDA dense output at quarter, half and
three-quarter points and at delayed query locations. Also check dense delayed
outputs against history when they use the current source iterate or span finer
stored intervals. If a step spans more than two stored intervals at a delayed
input, also probe their shifted midpoints to resolve
previously captured features as they arrive. IDA controls integration
error separately. These comparisons are error estimates, not rigorous bounds on
an arbitrary smooth input. Rejected candidates are retried from the accepted
state, with no history commit or output. Restarting IDA discards its internal
multistep history, so retries have a performance cost. The next step cap is

$$h\min\!\left(2,\max\!\left(0.5,\frac{0.8}{\sqrt[3]{\max(\eta,10^{-12})}}\right)\right),$$

where $\eta$ is the normalized history error estimate.
This predicts refinement before a rejection; IDA still
chooses its step under that cap. No cap depends on the shortest delay.
A temporary stop enforces a reduced cap when IDA has already selected a larger
step. Otherwise the stop remains at the simulation end, allowing step growth.

Adjacent stored intervals may merge into a cubic matching the outer values and
slopes. Bernstein coefficients on a common interval give the bound

$$\|P-Q\|_\infty\le\max_k|p_k-q_k|.$$

After subdivision of a candidate to each child interval, its new bound is the
maximum of each child's accumulated compression error plus this difference.
Merges spend at most one quarter of the source absolute tolerance. Interpolation
uses the remaining absolute tolerance and the relative tolerance. Without
compression, interpolation uses the full tolerance. Failed merge checks are
reused until one of the two intervals changes. Intervals
containing the current delayed queries and recent anchors are protected.

Compression is disabled by default. Setting `compress` to `false` retains every
accepted polynomial within the required window. This removes compression error;
integration and interpolation remain numerical approximations to the true delay.

## Usage

```json
"signals": [
  { "signal_id": 10, "name": "omega" },
  { "signal_id": 11, "name": "omega_delayed",
    "source": 10, "delay": 0.02, "mon": ["z"] }
]
```

Connect device inputs to signal 11 normally. `DynamicSimulation` and
`ContingencyAnalysis` attach the history controller automatically. Programmatic
forward simulations attach it before solver configuration:

```cpp
configureSignalDelays(system, ida, {
    .rel_tol = relative_tolerance,
    .abs_tol = absolute_tolerance,
    .compress = false});
ida.configureSimulation();
ida.initializeSimulation(t0);
ida.runSimulation(tf);
```

Include `GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp`.
Quadrature and adjoint stepping do not support this controller. A requested fixed
step fails explicitly if history accuracy requires refinement; adaptive stepping
can retry. Delayed inputs must remain smooth across any simulation events.

## Testing

Unit tests cover nonuniform interpolation, residuals and Jacobians, initialization,
declaration validation, history commits, retention and cumulative compression
bounds. Integration tests use the exact feedback solution $x(t)=e^{\lambda t}$
of $\dot x=-\gamma x+kx(t-\tau)$ with
$\gamma=ke^{-\lambda\tau}-\lambda$. They check convergence, steps longer than the
delay, compression, rejection rollback and continuation. A narrow prescribed
pulse with a long delay checks that resolved source features reach the output.
