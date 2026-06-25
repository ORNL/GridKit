# PowerFlowData Translation

## Purpose

Translate only values available from
`GridKit::PhasorDynamics::SystemModelData`. PowerFlow input supplies study
fields not present in PhasorDynamics: bus type, reference bus, statuses,
limits, ratings, area, zone, and standalone bus shunts.

## Buses

| PowerFlow | PhasorDynamics |
| --------- | -------------- |
| `baseMVA` | `va_base / 1e6` |
| `bus_i` | `bus_id` |
| `Vm` | `Vr0`, `Vi0` |
| `Va` | `Vr0`, `Vi0` |
| `baseKV` | `v_base / 1e3` |

```math
\begin{aligned}
\left|V_n\right| &= \sqrt{\left(V_n^r\right)^2 + \left(V_n^i\right)^2} \\
\theta_n &= \operatorname{atan2}\left(V_n^i, V_n^r\right)\frac{180}{\pi}
\end{aligned}
```

PowerFlow bus type and reference bus are not translated from PhasorDynamics
`bus` or `infinite_bus`.

## Loads

TBD.

```math
\begin{aligned}
R(S,V) &= ... \\
X(S,V) &= ...
\end{aligned}
```

Problem is that we need an option to intialize parameters as a function of the power flow solution, which is how we do this now.

## Generators

Each supported PhasorDynamics machine translates to one
`PowerFlowData::GenData` row. `SystemSteadyStateModel` then uses
`GeneratorFactory` to instantiate the PowerFlow generator model selected by
the PowerFlow bus type.

| PowerFlow | PhasorDynamics |
| --------- | -------------- |
| `gen.bus` | `Genrou`, `Gensal`, or `GenClassical` port `bus` |
| `Pg` | `Genrou`, `Gensal`, or `GenClassical` parameter `p0` |
| `Qg` | `Genrou`, `Gensal`, or `GenClassical` parameter `q0` |
| `mBase` | `Genrou`, `Gensal`, or `GenClassical` parameter `mva`, when present |

Generator limits, status, voltage setpoint, ramp data, and cost data are not
translated from PhasorDynamics data. Do not combine multiple generators at the
same bus during translation.

## Branches

| PowerFlow | PhasorDynamics |
| --------- | -------------- |
| `fbus`, `tbus` | `bus1`, `bus2` |
| `g`, `b` | derived `g`, `b` |
| `G`, `B` | `G`, `B` |
| $\tau$ | `tap`; default 1 |
| $\theta$ | `phase` $180 / \pi$; default 0 |

Branch shunt admittance is total line charging and is split equally between
terminals.
