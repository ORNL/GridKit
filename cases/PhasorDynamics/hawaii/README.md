# Hawaii

## One-Line Diagram

![](hawaii.png)

Figure 1: Oneline of the synthetic Hawaii case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Case Description

This case is a modified version of the synthetic Hawaii case, available at the Texas A&M University [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/). It has been modified to include GridKit-compatible generator models.

## Model Inventory

Model | Count
---|---
Bus | 37
Branch | 89
GENROU | 39
TGOV1 | 39
IEEET1 | 39
IEEEST | 14
LoadZIP | 28
SignalNode | 131

## Outstanding

### Dynamics

The following source dynamic models are not yet supported and are represented with surrogates:
- IEEEST (in GridKit, not yet added to this case)
- REGCA
- IEEEG1
- GGOV1
- ESST4B
- REECA1
- EXST1_PTI
- ESST1A

### Statics

GridKit represents switched shunts and transformers with constant impedance, so the case cannot be re-initialized at a different operating point.

- Transformers
- Switched Shunts
