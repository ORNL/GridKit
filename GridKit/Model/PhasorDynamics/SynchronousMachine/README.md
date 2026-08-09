# General Synchronous Machine Model

> [!WARNING]
> Generator voltage-compensation parameters $R_\mathrm{Comp}$ and $X_\mathrm{Comp}$ are not currently modeled.

## Convention

![](../../../../docs/Figures/SM1.JPG)

Figure 1: Synchronous Machine. Figure courtesy of
[PowerWorld](https://www.powerworld.com/files/Synchronous-Machines.pdf/)

The following conventions are used for the d-q reference frame.
- The q-axis leads the d-axis
- The Rotor angle is w.r.t. to q-axis

## Types

- Classical Generator (See [GenClassical](GenClassical/README.md))
- Round Rotor (See [GENROU](GENROU/README.md))
- Salient Rotor/Pole (See [GENSAL](GENSAL/README.md))
- GENPWS
- GENTPF
- GENTPJ
- GENQEC

### Per-Unit Basis

In relevant models, the terminal impedances are on the generator impedance base.
To convert to the network base, the following must be performed.

```math
\begin{aligned}
  Z_{term} &
  \mapsto Z_{term}\dfrac{S_{base,sys}}{S_{base,machine}}
\end{aligned}
```

For example, say the terminal impedence is $Z=0.05$ in per-unit on the
machine's base of $S_{base,machine}=50$  MW, and the system base is
$S_{base,sys}=100$ MW. Then the terminal impedance on the system
base is calculated as follows.

``` math
\begin{aligned}
  Z_{sys} = 0.05\dfrac{100 \text{MW}}{50 \text{MW}} = 0.1
\end{aligned}
```

#### Saturation

Saturation means increasingly large amounts of current are needed to increase
the flux density. The Scaled Quadratic saturation model is currently implemented.
``` math
\begin{aligned}
  k_{sat} =
  \begin{cases}
    S_B(\psi''-S_A)^2 &\text{if } \psi''>S_A \\
    0 &\text{if } \psi''\leq S_A
  \end{cases}
\end{aligned}
```
