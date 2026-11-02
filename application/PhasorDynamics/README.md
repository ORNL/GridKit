# Input file for GridKit phasor dynamics application

## Root elements

- `system_model_file`: Path to the system model file.[^1]
- `dt_monitor`: Monitor output time interval for recorded simulation results.
  Default: 0, with no intermediate monitoring.
- `tmax`: A floating-point value for max time.
- `rel_tol`: Relative solver tolerance. Default: 1.0e-7.
- `abs_tol`: Absolute solver tolerance override. Default: 1.0e-9.
- `cosim_rel_tol`: Relative co-simulation interface residual tolerance.
  Default: 1.0e-8.
- `cosim_step_tol`: KINSOL scaled step tolerance for the co-simulation
  interface. Default: 1.0e-12.
- `cosim_function_error`: Estimated relative error of co-simulation residual
  evaluations, used to choose finite difference perturbations. Default: 1.0e-10.
- `dt_fixed`: Fixed solver time step size, or 0 for adaptive stepping.
  Default: 0.
- `max_steps`: Maximum number of solver time steps, 0 for the IDA default,
  or a negative number for unlimited steps. Default: 0.
- `max_order`: Maximum IDA integration method order from 1 to 5. Default: 5;
  fixed stepping is capped at 2.
- `consistent_ic_type`: IDA consistent initial condition calculation type;
  one of `"y"`, `"ya_ydp"`. Default: `"ya_ydp"`. Co-simulation requires
  `"ya_ydp"`.
- `events`: An array of event groups (see [Events](#events) below).
- `output_file`: Path to output file. Optional.
- `output_layout`: Layout of `output_file`. This may be one of `"wide"`
  (default) or `"long"`. The long layout writes one row per variable and
  sample. Event times are written twice, before and after the event, so
  use `sample` rather than `t` to identify a sample. Not supported with
  `reference_file`.
- `reference_file`: A string containing the name of the case. Optional.
- `error_type`: One of `"relative"` (default), `"absolute"`.
- `error_tolerance`: A floating-point value for greatest allowable total error.
  Default: 1.0e-4.
- `abs_err_threshold`: A floating-point value for the smallest value at which
  to scale relative error. Default: machine epsilon for double-precision.

[^1]: See system model [case format].

[case format](../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md)

## Co-simulation interface solve

KINSOL is used to solve the coupling conditions at each simulation step. The
client coordinates this solve at initialization and after both partitions
have advanced to a co-simulation step boundary.

Let $\mathbf{u}$ contain the client's input signals, supplied by the server.
At time $t_n$, let $\mathcal{C}_n(\mathbf{u})$ denote the client output signals
after applying those inputs and restoring algebraic consistency. Similarly,
$\mathcal{S}_n$ maps the server's input signals to its consistent outputs.
The interface must satisfy

```math
\begin{aligned}
  \mathbf{v}_n &= \mathcal{C}_n(\mathbf{u}_n), \\
  \mathbf{u}_n &= \mathcal{S}_n(\mathbf{v}_n), \\
  F_n(\mathbf{u}) &= \mathcal{S}_n\!\left(\mathcal{C}_n(\mathbf{u})\right).
\end{aligned}
```

For each trial $\mathbf{u}$, the client applies the signals and calls
IDA's consistent IC calculation. It sends the resulting outputs to the
server in an `UPDATE` message. The server applies them, restores its own
algebraic consistency, and returns its outputs. Both partitions keep their
differential states fixed during these evaluations.

The first trial is the server's most recent output vector
$\mathbf{u}^{(0)}$. The client chooses a scale for each signal and supplies
the following residual to KINSOL:

```math
\begin{aligned}
  s_i &= \max\!\left(1, |u_i^{(0)}|\right), \\
  r_i(\mathbf{u}) &= \frac{F_{n,i}(\mathbf{u}) - u_i}{s_i}, \\
  \mathbf{r}(\mathbf{u}_n) &= \mathbf{0}.
\end{aligned}
```

If the initial trial already satisfies the interface tolerance, a KINSOL
iteration is not performed. Otherwise, KINSOL performs Newton iterations
without line search, with a finite difference Jacobian:

```math
\begin{aligned}
  J_{i,j}(\mathbf{u}) &= \frac{\partial r_i}{\partial u_j} \approx \frac{r_i(\mathbf{u} + h_j\mathbf{e}_j) - r_i(\mathbf{u})}{h_j}, \\
  \mathbf{J}(\mathbf{u}^{(k)})\,\Delta\mathbf{u}^{(k)} &= -\mathbf{r}(\mathbf{u}^{(k)}), \\
  \mathbf{u}^{(k+1)} &= \mathbf{u}^{(k)} + \Delta\mathbf{u}^{(k)}.
\end{aligned}
```

Here, $k$ indexes Newton iterations, $\mathbf{e}_j$ is the $j$-th coordinate
unit vector, and $h_j$ is the finite difference perturbation selected
by KINSOL. After each solve, each signal is compared against `cosim_rel_tol`:

```math
\max_i \frac{|F_{n,i}(\mathbf{u}) - u_i|}{\max(1, |u_i|)} < \texttt{cosim\_rel\_tol}.
```

If this check fails after KINSOL returns, the client retries once with a
function tolerance 100 times smaller, clamped to the smallest possible
subnormal `double`. If the check still fails after retrying, a runtime
error is thrown with a message specifying that the interface did not
converge. Otherwise, the client sends `COMMIT`, then both partitions record
any pending monitor output and the simulation proceeds.

## Events

Each event group describes a system event that occurs at a given time point.

- `time`: A floating-point value for the time the event occurs.
- `type`: Event type; one of `"fault_on"`, `"fault_off"`, `"in_service"`,
  `"out_of_service"`.
- `bus`: Number of the faulted bus. Omitted for `ContingencyAnalysis`, which
  faults every bus in turn using the same `R` and `X`.
- `R`: Fault resistance in per unit. Required for `"fault_on"`.
- `X`: Fault reactance in per unit. Required for `"fault_on"`.
- `device`: `id` of the branch to put in service (`"in_service"`) or take out of
  service (`"out_of_service"`), in the top-level system or any subsystem.
