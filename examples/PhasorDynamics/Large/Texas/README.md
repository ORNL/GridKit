# Synthetic Texas Case (ACTIVSg2000)

## Case Description

This case is adapted from the ACTIVSg2000 synthetic Texas grid from the Texas A&M Electric Grid Test Case Repository:
https://electricgrids.engr.tamu.edu/electric-grid-test-cases/

It is configured with GridKit-compatible dynamic models.

Model | Count
---|---
[Bus](../../../../GridKit/Model/PhasorDynamics/Bus/README.md) | 2000
[Branch](../../../../GridKit/Model/PhasorDynamics/Branch/README.md) | 3206
[Load](../../../../GridKit/Model/PhasorDynamics/Load/README.md) | 1507
[GENROU](../../../../GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/README.md) | 544
[TGOV1](../../../../GridKit/Model/PhasorDynamics/Governor/Tgov1/README.md) | 544
[IEEET1](../../../../GridKit/Model/PhasorDynamics/Exciter/IEEET1/README.md) | 544

## Data Notes

- Generator `6215_1` in source data had `Xl == Xdpp == Xqpp = 0.2011`, which can trigger division-by-zero during GENROU initialization.
- In this case, `Xl` was reduced to `0.1911` to satisfy `Xl < Xdpp`.

## Typical Events

- Bus fault (example configuration targets bus `1027`)
- Line outage studies
- Generator outage studies
- Forced oscillation studies
