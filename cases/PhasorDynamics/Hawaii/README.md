# Synthetic Hawaii

## One-Line Diagram

![Hawaii one-line diagram](hawaii.png)

Figure 1: One-line diagram of the synthetic Hawaii case.[^1]

## Development

Unimplemented models:

- `ESST1A`
- [ESST4B](../../../GridKit/Model/PhasorDynamics/Exciter/ESST4B/README.md)
- `EXST1_PTI`
- [GGOV1](../../../GridKit/Model/PhasorDynamics/Governor/GGOV1/README.md)
- [IEEEG1](../../../GridKit/Model/PhasorDynamics/Governor/IEEEG1/README.md)
- `REECA1`

This case was validated against PowerWorld. The comparison results are provided [here](../../../examples/PhasorDynamics/Validation/Hawaii/README.md).

`Hawaii.m` is the MATPOWER export `Hawaii40_20231026.m` of the Hawaii40 case[^2], with the limits and costs for [optimal dispatch](../../../examples/PhasorDynamics/OptimalDispatch/README.md).

[^1]: Texas A&M University, [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/).
[^2]: Texas A&M University, [Hawaii40](https://electricgrids.engr.tamu.edu/hawaii40/).
