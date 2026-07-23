# New England

## One-Line Diagram

![](newengland.png)

Figure 1: Oneline of the New England IEEE 39-bus case, courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Case Description

This case is a modified version of the NE 39-bus case, available at the Texas A&M University [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/). It has been modified to include GridKit-compatible generator models.

## Model Inventory

Model | Count
---|---
Bus | 39
Branch | 46
GENROU | 10
TGOV1 | 10
IEEET1 | 10
IEEEST | 10
LoadZIP | 19
SignalNode | 40

## Outstanding

### Dynamics

The following source dynamic models are not yet supported:
- EXST1_GE

### Statics

GridKit represents switched shunts and transformers with constant impedance, so the case cannot be re-initialized at a different operating point.

- Transformers (LTC and Tap Ratios $\neq 1$)
