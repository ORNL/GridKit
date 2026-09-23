# IEEE 39-Bus New England

## One-Line Diagram

![IEEE39 one-line diagram](IEEE39.png)

Figure 1: One-line diagram of the IEEE 39-Bus case.[^1]

## Development

Unimplemented models:

- `EXST1_GE`

This case was validated against PowerWorld. The comparison results are provided [here](../../../examples/PhasorDynamics/Validation/IEEE39/README.md).

`IEEE39.m` is MATPOWER `case39`[^2], with the limits and costs for [optimal dispatch](../../../examples/PhasorDynamics/OptimalDispatch/README.md).

[^1]: Texas A&M University, [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/).
[^2]: MATPOWER, [`case39.m`](https://github.com/MATPOWER/matpower/blob/master/data/case39.m).
