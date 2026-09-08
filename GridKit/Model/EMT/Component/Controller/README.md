# Controller Models

EMT controller models exchange signals with other components. In the switching
inverter examples, PLL supplies angle to the Park operators and frequency to
the controllers through signal ports.

## Grid Following

![Grid-following inverter wiring](../../../../../docs/Figures/EMT/Controller/diagram_gfl.png)

## GFM Voltage Control

![GFM case wiring with PLL](../../../../../docs/Figures/EMT/Controller/diagram_gfm.png)

## Models

- [InnerCurrentControl](InnerCurrentControl/README.md): converter current control in $dq$ coordinates.
- [OuterPowerControl](OuterPowerControl/README.md): active and reactive power control in $dq$ coordinates.
- [OuterVoltageControl](OuterVoltageControl/README.md): filter-capacitor voltage control in $dq$ coordinates.
- [DC Link](DCLink/README.md): capacitor voltage and current balance.
- [IEEET1](IEEET1/README.md)
- [PWM](PWM/README.md)
- [TGOV1](TGOV1/README.md)
- [SEXS-PTI](SEXS-PTI/README.md): simplified excitation system.
- [GASTPTI](GASTPTI/README.md): gas turbine governor with exhaust-temperature limiting.
- [IEEEST](IEEEST/README.md): power system stabilizer.
