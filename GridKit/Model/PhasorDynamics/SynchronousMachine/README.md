# Synchronous Machine Models

## Conventions

![](../../../../docs/Figures/SM1.JPG)

Figure 1: Synchronous Machine. Figure courtesy of
[PowerWorld](https://www.powerworld.com/files/Synchronous-Machines.pdf/)

For the d–q reference frame:
- The q-axis leads the d-axis
- The rotor angle is measured from the q-axis

## Types

- Classical Generator (See [GenClassical](GenClassical/README.md))
- Round Rotor (See [GENROU](GENROU/README.md))
- Salient Rotor/Pole (See [GENSAL](GENSAL/README.md))
- GENPWS
- GENTPF
- GENTPJ
- GENQEC

## Per-Unit Basis

Terminal impedances use the machine base. With a common voltage base,

```math
Z^\mathrm{sys}=Z^\mathrm{mach}\dfrac{S^\mathrm{sys}}{S^\mathrm{base}}
```

$S^\mathrm{sys}$ and $S^\mathrm{base}$ are the system and machine power
bases in megavolt-amperes.

## Saturation

GENROU and GENSAL use the smooth [quadratic ramp](../../../CommonMath.md#quadratic-ramp):

```math
k_\mathrm{sat}=S_Bq(\psi-S_A)
```

Here $\psi=\psi''$ for GENROU and $\psi=E'_q$ for GENSAL. Each model defines
its saturation fit.
