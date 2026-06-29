# Synthetic Texas Case (ACTIVSg2000)

## Case Description

This case is adapted from the ACTIVSg2000 synthetic Texas grid from the Texas A&M Electric Grid Test Case Repository:
https://electricgrids.engr.tamu.edu/electric-grid-test-cases/

It is configured with GridKit-compatible dynamic models.

Model | Count
---|---
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md) | 2000
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md) | 3206
[LoadZ](../../../../GridKit/Model/PhasorDynamics/Load/LoadZ/README.md) | 1507
[GENROU](../../../../GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/README.md) | 544
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md) | 544
[IEEET1](../../../../GridKit/Model/PhasorDynamics/Exciter/IEEET1/README.md) | 544

## Data Notes

- Generator `6215_1` in source data had `Xl == Xdpp == Xqpp = 0.2011`, which can trigger division-by-zero during GENROU initialization.
- In this case, `Xl` was reduced to `0.1911` to satisfy `Xl < Xdpp`.


## Events

The following event types are provided for this case.

- Bus fault


## Outstanding

### Dynamics

Have not been implemented in GridKit:
- ESAC1A
- ESDC1A
- ESST4B
- EXAC2
- EXPIC1
- SCRX
- GGOV1
- HYGOV
- IEEEG1
- GENSAL (In GridKit, not added to this case yet)
- IEEEST (In GridKit, not added to this case yet)

The following examples needs to be constructed with this case.
- Line Outage
- Generator Outage
- Forced Oscillations

### Statics

GridKit models Switched Shunts and Transformers with constant impedance, but this erases information and our ability to statically re-initialze the model if initalizing at a different operating point.

- Switched Shunts
- Transformers (LTC and Tap Ratios $\neq 1$)
