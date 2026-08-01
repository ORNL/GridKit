# Western Electricity Coordinating Council (WECC)

## One-Line Diagram

![](wecc.jpg)

Figure 1: Oneline of the WECC case of [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository)

## Case Description

This case is adapted from the [National Laboratory of the Rockies](https://www.nlr.gov/grid/test-case-repository) WECC Model.

It is configured with GridKit-compatible dynamic models.

## Data Notes

None.

## Outstanding

### Dynamics

Only one exciter model is outstanding:
- REECB1
- GAST_PTI
- REGCA
- REPCA1

The following examples needs to be constructed with this case.
- Line Outage
- Generator Outage
- Forced Oscillations

### Statics

GridKit models Switched Shunts and Transformers with constant impedance, but this erases information and our ability to statically re-initialze the model if initalizing at a different operating point.

- Switched Shunts
- Transformers (LTC and Tap Ratios $\neq 1$)
