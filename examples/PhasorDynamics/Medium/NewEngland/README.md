# **IEEE New England**

## One-Line Diagram

![](newengland.png)

Figure 1: Oneline of the New England IEEE 39-bus case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Case Description

This case is a modified version of the NE 39-bus case, available at the Texas A&M University [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/). It has been modified to include GridKit-compatible generator models.

Model       | Count
------------|--------
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md)         | 39
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md)     | 46
[GENROU](../../../../GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/README.md)       | 10
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md)        | 10
[IEEET1](../../../../GridKit/Model/PhasorDynamics/Exciter/IEEET1/README.md)  | 10
[IEEEST](../../../../GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/README.md)  | 10
[LoadZIP](../../../../GridKit/Model/PhasorDynamics/Load/LoadZIP/README.md) | 19
[SignalNode](../../../../GridKit/Model/PhasorDynamics/SignalNode/README.md) | 40


## Events

The following event types are provided for this case.

- Bus fault


## Outstanding

### Dynamics

Only one exciter model is outstanding:
- EXST1_GE

The following examples needs to be constructed with this case. Currently the only events that have been implemented are faults. This is a practical case to use for modeling outages, or even forced oscillations.
- Line Outage
- Generator Outage
- Forced Oscillations

### Statics

GridKit models Switched Shunts as [LoadZ](../../../../GridKit/Model/PhasorDynamics/Load/LoadZ/README.md) with constant impedance,but this erases information and our ability to statically re-initialze the model if initalizing at a different operating point.

- Transformers (LTC and Tap Ratios $\neq 1$)
