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

Let $h=t-t_n$, $x=(q-t_n)/h$, and let $d_n$ be the stored slope at $t_n$.
The trial segment is cubic Hermite interpolation:

$$\widehat u(q)=H_{00}(x)u_n+hH_{10}(x)d_n
                  +H_{01}(x)u(t)+hH_{11}(x)d(t),$$

$$H_{00}=1-3x^2+2x^3,\quad H_{10}=x-2x^2+x^3,\quad
  H_{01}=3x^2-2x^3,\quad H_{11}=x^3-x^2.$$

The right slope comes from the interpolant through up to three recent accepted
values $(s_i,u_i)$ and the live endpoint $u(t)$:

$$d(t)=a_d+b_du(t),\qquad b_d=\sum_i\frac{1}{t-s_i},$$

$$a_d=\sum_i\frac{u_i}{s_i-t}
                   \prod_{j\ne i}\frac{t-s_j}{s_i-s_j}.$$

Thus $a=H_{00}u_n+hH_{10}d_n+hH_{11}a_d$ and
$b=H_{01}+hH_{11}b_d$. On acceptance, this endpoint slope becomes the next
segment's left slope. Adjacent segments share both value and slope.
Startup uses the available anchors.

The Jacobian entries remain $\partial F_z/\partial z=-1$ and
$\partial F_z/\partial u=b$. For a query in accepted history, $b=0$.
The source column stays in the sparse pattern across the history boundary.
Steps may exceed any delay: the current endpoint participates in the same
implicit solve as the delayed output. Accuracy and nonlinear convergence
control the usable step size; there is no restriction $h\le\min\tau$.

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
responsibility. The initial slope is zero for constant prehistory and a
second-order backward difference for a supplied function.

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

The IDAS step-check extension runs after native nonlinear/error tests and before
committing a step. It compares the trial history with source samples at six
Chebyshev-Lobatto nodes and the midpoint. This measures the reconstruction error
directly, without a lower-order interpolant that can overestimate it.

IDAS dense output has degree at most five. The same six samples recover its
Bernstein coefficients. On every intersected history interval, subdivide this
polynomial and elevate the stored cubic to degree five. Then

$$\|z_{\rm IDAS}(\cdot)-\widehat u(\cdot-\tau)\|_\infty
  \le \max_i|c_i-b_i|.$$

The relative error scale uses a lower bound on the absolute history value on
that interval (zero if its Bernstein range contains zero). Consequently this
check also detects a stored feature between the dense-output sample times.
It covers accepted and current-trial history without a segment-count heuristic.
Arbitrary prescribed source functions and prehistory are sampled; their errors
are estimates, not rigorous continuous bounds.

For normalized error $\eta$, predict a history cap with

$$h_{\rm history}=h\,\operatorname{clip}
 \left(\frac{0.8}{\max(\eta,10^{-12})^{1/4}},\;0.1,\;2\right).$$

On rejection, limit the factor to at most $0.5$. IDAS restores its trial state,
retains multistep history and retries at the reduced step. An accepted step uses
the smaller of IDAS's proposed step and the history cap. Only accepted steps
commit source history and publish output. No restart or temporary stop is needed.

The optional solver setting `h_min` specifies the minimum adaptive step magnitude.
Both native and history error tests fail if accuracy cannot be attained there.
As in IDAS, landing on the final stop time may require a shorter last step.
A step floor is not a source sampling guarantee: an arbitrarily narrow smooth
feature can evade finitely many probes. It does not impose a shortest-delay cap.

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
ida.setMinStep(1e-6); // Optional user-selected floor, in seconds.
ida.configureSimulation();
ida.initializeSimulation(t0);
ida.runSimulation(tf);
```

Include `GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp`.
Signal delays require the local
[IDAS step-check patch](../../../../buildsystem/spack_repo/gridkit/packages/sundials/README.md).
The GridKit Spack recipe applies it to SUNDIALS 7.8.0. Unpatched SUNDIALS remains
usable for simulations without delays; requesting history control fails explicitly.
Quadrature and adjoint stepping do not support this controller. A requested fixed
step fails explicitly if history accuracy requires refinement; adaptive stepping
can retry. Delayed inputs must remain smooth across any simulation events.

## Testing

Unit tests cover nonuniform C1 interpolation, residuals and Jacobians, initialization,
declaration validation, history commits, retention and cumulative compression
bounds. Integration tests use the exact feedback solution $x(t)=e^{\lambda t}$
of $\dot x=-\gamma x+kx(t-\tau)$ with
$\gamma=ke^{-\lambda\tau}-\lambda$. They check convergence, steps longer than the
delay, compression, rejection rollback and continuation. A narrow prescribed
pulse with a long delay checks that resolved source features reach the output.
A direct IDAS test checks candidate versus committed dense polynomials and their
derivatives, retained order on retry, and decreasing maximum steps. Solver tests
also cover the step floor and exceptions raised inside a trial check.
