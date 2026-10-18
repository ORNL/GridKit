# Optimal Power Flow

AC optimal power flow in Cartesian coordinates. It finds the least-cost
operating point within the component limits.

## Formulation

Minimize $f(x)$ subject to the equality constraints $c(x) = 0$, the inequality
constraints $d_l \le d(x) \le d_u$, and the variable bounds
$x_l \le x \le x_u$. The objective and equality rows sum component contributions.
Each inequality row and variable bound belongs to one component. An inequality
without a finite limit has infinite bounds. Every island needs a reference bus.

Cartesian voltages make a solution directly a starting point for a phasor
dynamics simulation. Branch power is quadratic in these voltages; ZIP
constant-current demand also depends on their magnitude.

## Components

`SystemModel` implements `OptimizationEvaluator` and owns the global vectors
and sparse matrices. Every component, including `Bus`, views its own blocks
of those vectors. Components do not inherit the complete nonlinear program.
Service status is a parameter of each component that uses it.

`NUM_INTERNAL_VARIABLES`, `NUM_INTERNAL_EQUALITIES`, and
`NUM_INTERNAL_INEQUALITIES` declare the owned dimensions. `NUM_PORT_VARIABLES`
and `NUM_PORT_EQUALITIES` declare, at each port, the variables read and the
equality rows contributed to. `Component` derives every other size from these
counts. Local arrays hold internal entries first, then external entries in port
order. A bus defines what a port to it shares: its voltage and its power
balance. `connect(port, component, offset)` selects contiguous blocks owned by
any component; variable and equality offsets are independent. Connections must
target components owned by the same system. Allocation binds the views and
prevents further connection changes.

### Internal Constraints

#### Equality

`internalEqualities(x, c)` contributes to rows owned by the component. Each
equality row is owned by one component and is the sum of its owner's
contribution and those of every component with a port to it; the assembled row
equals zero.

#### Inequality

`internalInequalities(x, d)` evaluates rows owned by the component. These may
depend on both internal and external variables. Variable bounds remain separate
from these function bounds.

### External Constraints

`externalEqualities(x, c)` contributes to equality rows owned by the connected
components. These contributions are added to the owning rows; they are not
separate constraints equated to zero. For example, a generator contributes
power to a bus's balance rows without reading its voltage.

Case bus numbers are resolved when constructing the system. Electrical islands
follow in-service branches, independently of component dependencies.

`SystemModel` also accepts a vector of owned components for direct composition.
That constructor validates component parameters and connections. Electrical
island and reference checks use the network metadata supplied by the case-data
constructor.

## Derivatives

Each component supplies `objective(x)` and the constraint kernels above,
which `Component` differentiates with Enzyme. Internal and external equality
contributions share one differentiation pass:

Driver   | Derivative                                            | Mode
---------|-------------------------------------------------------|-----
`DfDx`   | Objective gradient                                    | Reverse
`DgDx`   | Equality or inequality constraint Jacobian            | Forward
`D2LDx2` | Hessian of $\sigma f + \lambda_c^T c + \lambda_d^T d$ | Forward over reverse

`Component` retains the lower triangle in global indices for assembly.

## Data

`parseSystemModelData` reads a JSON input file in the
[case format](../PhasorDynamics/INPUT_FORMAT.md), with the classes
[`Bus`](Bus/README.md), [`Branch`](Branch/README.md),
[`Generator`](Generator/README.md), and [`Load`](Load/README.md). Values are
per unit on the system base `va_base`.

A bus carries its shunt admittance to ground in `G` and `B`. Loads use `P`,
`Q`, `Pi`, `Qi`, `Pz`, and `Qz`: the constant-power, constant-current, and
constant-impedance demand at 1 p.u. voltage.
