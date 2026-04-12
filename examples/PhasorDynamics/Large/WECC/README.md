# Western Electricity Coordinating Council (WECC)

## One-Line Diagram

<div align="center">
   <img align="center" src="wecc.jpg">

  Figure 1: Oneline of the WECC case of [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository)
</div>

## Case Description

This case is adapted from the [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository) WECC Model.

It is configured with GridKit-compatible dynamic models.

Model | Count
---|---
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md) | 243
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md) | 447
[GENROU](../../../../GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/README.md) | 140
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md) | 103
[SEXS-PTI](../../../../GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/README.md) | 103
[IEEEST](../../../../GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/README.md)  | 10
[LoadZIP](../../../../GridKit/Model/PhasorDynamics/LoadZIP/README.md) | 146
[SignalNode](../../../../GridKit/Model/PhasorDynamics/SignalNode/) | 319

## Data Notes

None.


## Events

The following event types are provided for this case.

- Bus fault


## Outstanding

### Dynamics

Only one exciter model is outstanding:
- REECB1
- GAST_PTI
- REGC_A
- REPCA1

The following examples needs to be constructed with this case.
- Line Outage
- Generator Outage
- Forced Oscillations

### Statics

GridKit models Switched Shunts and Transformers with constant impedance, but this erases information and our ability to statically re-initialze the model if initalizing at a different operating point.

- Switched Shunts
- Transformers (LTC and Tap Ratios $\neq 1$)
