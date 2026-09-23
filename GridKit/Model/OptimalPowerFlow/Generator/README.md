# Generator

Dispatchable active and reactive power injection with a quadratic cost. The
`OptimalDispatch` application makes one for every GENROU, GENSAL,
GenClassical, and REGCA device in the case, with limits and costs from the
MATPOWER case on the system base.

## Model Parameters

Symbol               | Units  | JSON   | Description                   | Typical Value | Note
---------------------|--------|--------|-------------------------------|---------------|-----
$P^{\min}$           | [p.u.] | `Pmin` | Active power lower limit      | 0             | Default unbounded
$P^{\max}$           | [p.u.] | `Pmax` | Active power upper limit      | 1             | Default unbounded
$Q^{\min}$           | [p.u.] | `Qmin` | Reactive power lower limit    | -0.5          | Default unbounded
$Q^{\max}$           | [p.u.] | `Qmax` | Reactive power upper limit    | 0.5           | Default unbounded
$c_0$                |        | `c0`   | Constant cost coefficient     | 0             | Default 0
$c_1$                |        | `c1`   | Linear cost coefficient       | 1000          | Default 0
$c_2$                |        | `c2`   | Quadratic cost coefficient    | 100           | Default 0

Values are on the system base.

### Parameter Validation

- $P^{\min} \le P^{\max}$
- $Q^{\min} \le Q^{\max}$
- $c_0$, $c_1$, and $c_2$ are finite

### Model Derived Parameters

None.

## Model Ports

Name  | Port | Init  | Description
------|------|-------|------------
`bus` | Bus  | Known | Terminal bus

## Model Variables

### Internal Variables

Symbol | Units  | Description                    | Note
-------|--------|--------------------------------|-----
$P$    | [p.u.] | Active power into the bus      | System base
$Q$    | [p.u.] | Reactive power into the bus    | System base

### External Variables

Symbol           | Units  | Description                  | Note
-----------------|--------|------------------------------|-----
$V_\mathrm{r}$   | [p.u.] | Terminal real voltage        | Owned by the bus
$V_\mathrm{i}$   | [p.u.] | Terminal imaginary voltage   | Owned by the bus

## Model Equations

### Objective

```math
a \left(c_0 + c_1 P + c_2 P^2\right)
```

### External Constraints

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= P \\
\Delta Q^\mathrm{bus} &\mathrel{+}= Q
\end{aligned}
```

### Bounds

```math
\begin{aligned}
P^{\min} &\le P \le P^{\max} \\
Q^{\min} &\le Q \le Q^{\max}
\end{aligned}
```

An offline generator has $P = Q = 0$.

## Initialization

The state terminal current and bus voltage give the starting point, which is
zero without them:

```math
P + jQ \leftarrow V_\mathrm{bus} I^*
```

$a \leftarrow 0$ if the state marks the device `online` as false, and $a
\leftarrow 1$ otherwise.
