# Optimal Power Flow

AC optimal power flow in Cartesian coordinates: the least-cost operating point
within the component limits.

## Formulation

```math
\begin{aligned}
\min_x \quad & f(x) \\
\text{subject to} \quad & g_l \le g(x) \le g_u \\
& x_l \le x \le x_u
\end{aligned}
```

The objective and each constraint sum component contributions, and an equality
has $g_l = g_u = 0$. Every island needs a reference bus, which the system does
not check. Cartesian voltages make a solution a starting point for phasor
dynamics.

## Components

A component owns blocks of the variables and constraints. At each port it reads
the variables of another component and adds to its constraints, so a device
adds its power to the balance of its bus. `tagInequality()` marks each
constraint as an inequality or an equality, as `tagDifferentiable()` marks the
variables in PhasorDynamics. `allocate()` composes the system and builds the
sparsity patterns, and `initialize()` sets the bounds and the flat start.
Enzyme differentiates the components for the gradient, the constraint Jacobian,
and the Hessian of $\sigma f + \lambda^T g$ at the multipliers `lambda()`.

## Data

`parseSystemModelData` reads a JSON case in the
[case format](../PhasorDynamics/INPUT_FORMAT.md) with the classes
[`Bus`](Bus/README.md), [`Branch`](Branch/README.md),
[`Generator`](Generator/README.md), and [`Load`](Load/README.md), in per unit
on the system base `va_base`. It reads the `number` and `params` of each bus and
the `class`, `id`, `ports`, and `params` of each device, and ignores other keys.
An unknown class, parameter, or port is an error. Real parameters are finite
and accept integers. `reference` and `in_service` are booleans, and $u$ is 1 in
service and 0 out of service. `Bus` and `Branch` share their class names with
PhasorDynamics but not their parameters, so a case file belongs to one family.
