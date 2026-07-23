# Texas

## Case Description

This case is adapted from the ACTIVSg2000 synthetic Texas grid from the Texas A&M Electric Grid Test Case Repository:
https://electricgrids.engr.tamu.edu/electric-grid-test-cases/

It is configured with GridKit-compatible dynamic models.

## Model Inventory

Model | Count
---|---
Bus | 2000
Branch | 3206
LoadZ | 1507
GENROU | 544
TGOV1 | 544
IEEET1 | 544

## Data Notes

- Generator `6215_1` in source data had `Xl == Xdpp == Xqpp = 0.2011`, which can trigger division-by-zero during GENROU initialization.
- In this case, `Xl` was reduced to `0.1911` to satisfy `Xl < Xdpp`.

## Outstanding

### Dynamics

The following source dynamic models are not yet supported:
- ESAC1A
- ESDC1A
- ESST4B
- EXAC2
- EXPIC1
- SCRX
- GGOV1
- HYGOV
- IEEEG1
- GENSAL (in GridKit, not added to this case yet)
- IEEEST (in GridKit, not added to this case yet)

### Statics

GridKit represents switched shunts and transformers with constant impedance, so the case cannot be re-initialized at a different operating point.

- Switched Shunts
- Transformers (LTC and Tap Ratios $\neq 1$)
