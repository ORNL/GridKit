# Toy Cases

## Development

Two- and three-bus cases used by the PhasorDynamics examples and integration tests:

- [TwoBusBasic](./TwoBusBasic.case.json): GENROU machine.
- [TwoBusGensal](./TwoBusGensal.case.json): GENSAL machine with TGOV1 and IEEET1 controls.
- [TwoBusIeeet1](./TwoBusIeeet1.case.json): GENROU machine with TGOV1 and IEEET1 controls.
- [TwoBusTgov1](./TwoBusTgov1.case.json): GENROU machine with a TGOV1 governor.
- [ThreeBusBasic](./ThreeBusBasic.case.json): GENROU machines and a constant-impedance load.
- [ThreeBusClassical](./ThreeBusClassical.case.json): Classical machines and a constant-impedance load.
- [ThreeBusConstantSource](./ThreeBusConstantSource.case.json): GENROU machine on a `SignalVoltageOut` bus fed by a constant signal source.
- [ThreeBusRelay](./ThreeBusRelay.case.json): ThreeBusBasic with line `BR_1_2` split at a mid-line bus into BranchBreakers sections protected by a definite-time and an inverse-time relay.
- [ThreeBusZipLoad](./ThreeBusZipLoad.case.json): Classical machines and a ZIP load.
