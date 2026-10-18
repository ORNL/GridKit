# Optimal Power Flow

AC optimal power flow in Cartesian coordinates. It finds the least-cost
operating point within the component limits.

## Formulation

Minimize $f(x)$ subject to the constraints $g_l \le g(x) \le g_u$ and the
variable bounds $x_l \le x \le x_u$. The objective and each constraint sum
component contributions. Each constraint and variable bound belongs to one
component, which tags the constraint as an inequality or leaves it an equality
with $g_l = g_u = 0$. An inequality without a finite limit has infinite bounds.
Every island needs a reference bus.

Cartesian voltages make a solution directly a starting point for a phasor
dynamics simulation. Branch power is quadratic in these voltages; ZIP
constant-current demand also depends on their magnitude.

## Components

`SystemModel` implements `OptimizationEvaluator` and owns the global vectors
and sparse matrices. Every component, including `Bus`, views its own blocks
of those vectors. Components do not inherit the complete nonlinear program.
Service status is a parameter of each component that uses it.

Each model enumerates its internal variables, its external variables, its
internal constraints, and its external constraint contributions as
`InternalVariablesT`, `ExternalVariablesT`, `InternalConstraintsT`, and
`ExternalConstraintsT`. `NUM_PORT_VARIABLES` and `NUM_PORT_CONSTRAINTS` split
the external entries among the ports. `Component` derives every size from
these. Each kernel receives the internal variables `x` and the external
variables `external`, each indexed by its enumeration. A bus defines what a
port to it shares: its voltage and its power balance.
`connect(port, component, offset)` selects contiguous blocks owned by any
component; variable and constraint offsets are independent. Connections must
target components owned by the same system. `allocate()` binds the views,
which prevents further connection changes, tags the inequalities, and builds
the derivative patterns. `initialize()` sets the bounds and the flat start.
Every variable and multiplier starts at zero, every equality at zero, and every
inequality unbounded, and a model overrides `initialize()` for its limits and
any nonzero initial value.

### Internal Constraints

`internalConstraints(x, external, g)` contributes to the constraints owned by
the component. These may depend on both internal and external variables. Each
constraint is the sum of its owner's contribution and those of every component
with a port to it. `tagInequality()` tags each of them as an inequality or an
equality, as `tagDifferentiable()` tags each variable in PhasorDynamics.
Variable bounds remain separate from these constraint bounds.

### External Constraints

`externalConstraints(x, external, g)` contributes to constraints owned by the
connected components, equalities and inequalities alike. These contributions
are added to the owning constraints; they are not separate constraints. For
example, a generator contributes power to a bus's balance without reading its
voltage.

Case bus numbers are resolved when constructing the system. Electrical islands
follow in-service branches, independently of component dependencies.

`SystemModel` also accepts a vector of owned components for direct composition.
That constructor validates component parameters and connections. Electrical
island and reference checks use the network metadata supplied by the case-data
constructor.

## Derivatives

Each component supplies `objective(x, external)` and the constraint kernels above,
which `Component` differentiates with Enzyme. Internal and external
contributions share one differentiation pass:

Driver   | Derivative                          | Mode
---------|-------------------------------------|-----
`DfDx`   | Objective gradient                  | Reverse
`DgDx`   | Constraint Jacobian                 | Forward
`D2LDx2` | Hessian of $\sigma f + \lambda^T g$ | Forward over reverse

`Component` retains the lower triangle in global indices for assembly. The
Hessian is at the multipliers `lambda()`. A solver sets them, and leaves them at
the solution with the bound multipliers `zLower()` and `zUpper()`.

## Data

`parseSystemModelData` reads a JSON input file in the
[case format](../PhasorDynamics/INPUT_FORMAT.md), with the classes
[`Bus`](Bus/README.md), [`Branch`](Branch/README.md),
[`Generator`](Generator/README.md), and [`Load`](Load/README.md). Values are
per unit on the system base `va_base`.

The parser reads `params.va_base`, the `number` and `params` of each bus, and
the `class`, `id`, `ports`, and `params` of each device, and ignores every
other key. An unknown device class, parameter, or port is an error. Each model
reads its parameters with `Model::ParameterReader`, so a real accepts an
integer, and `reference` and `in_service` must be booleans. `Bus` and `Branch`
share their class names with PhasorDynamics but not their parameters, so a
case file belongs to one family.

A bus carries its shunt admittance to ground in `G` and `B`. Loads use `P`,
`Q`, `Pi`, `Qi`, `Pz`, and `Qz`: the constant-power, constant-current, and
constant-impedance demand at 1 p.u. voltage.
