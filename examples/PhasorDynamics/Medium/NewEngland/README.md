# IEEE New England

## One-Line Diagram

![](oneline.png)

Figure 1: One-line diagram of the IEEE 39-bus New England case. Source:
[PowerWorld](https://www.powerworld.com/WebHelp/).

## Case Description

This case is a modified version of the Texas A&M University New England 39-bus
system[^newengland39] and includes GridKit-compatible generator models.

## Outstanding

- [EXST1_GE](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20EXST1_GE.htm)
- Switched shunts represented as [LoadZ](../../../../docs/GridKit/Model/PhasorDynamics/Load/LoadZ/README.md)
- LTC transformers and non-unit tap ratios

[^newengland39]: Texas A&M University,
    “[New England IEEE 39-Bus System](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/new-england-ieee-39-bus-system/),”
    *Electric Grid Test Case Repository*.
