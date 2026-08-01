# Synthetic Illinois (ACTIVSg200)

## One-Line Diagram

![](oneline.png)

*Figure 1: One-line diagram of ACTIVSg200. Source: Texas A&M University.[^activsg200]*

## Case Description

ACTIVSg200 is a 200-bus synthetic system geographically based on Illinois. It
was built from public information and statistical characteristics of real power
systems; it does not represent the actual Illinois grid. Its dynamic models are
fully represented in GridKit.

## Data Notes

The source `.pwb` file did not include a machine model for
`Bus 197 GIBSON CITY 1 2`. This case adds a `GENROU` model because steady-state
initialization fails if the generator is represented by a negative-impedance
load.

## Outstanding

- Transformers
- Switched shunts

[^activsg200]: Texas A&M University,
    “[Illinois 200-Bus System: ACTIVSg200](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg200/),”
    *Electric Grid Test Case Repository*.
