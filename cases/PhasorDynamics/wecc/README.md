# WECC

## One-Line Diagram

![](wecc.jpg)

Figure 1: Oneline of the WECC case of [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository)

## Case Description

This case is adapted from the [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository) WECC Model.

It is configured with GridKit-compatible dynamic models.

## Model Inventory

Model | Count
---|---
Bus | 243
Branch | 447
GENROU | 140
TGOV1 | 103
SEXS-PTI | 103
IEEEST | 10
LoadZIP | 146
SignalNode | 319

## Outstanding

### Dynamics

The following source dynamic models are not yet supported:
- REECB1
- GAST_PTI
- REGCA
- REPCA1

### Statics

GridKit represents switched shunts and transformers with constant impedance, so the case cannot be re-initialized at a different operating point.

- Switched Shunts
- Transformers (LTC and Tap Ratios $\neq 1$)
