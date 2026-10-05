# Toy Cases

## Development

Small cases used by the PhasorDynamics examples and integration tests:

- [TwoBusBasic](./TwoBusBasic.case.json): GENROU machine.
- [TwoBusGensal](./TwoBusGensal.case.json): GENSAL machine with TGOV1 and IEEET1 controls.
- [TwoBusIeeet1](./TwoBusIeeet1.case.json): GENROU machine with TGOV1 and IEEET1 controls.
- [TwoBusTgov1](./TwoBusTgov1.case.json): GENROU machine with a TGOV1 governor.
- [ThreeBusBasic](./ThreeBusBasic.case.json): GENROU machines and a constant-impedance load.
- [ThreeBusClassical](./ThreeBusClassical.case.json): Classical machines and a constant-impedance load.
- [ThreeBusConstantSource](./ThreeBusConstantSource.case.json): GENROU machine with a constant signal source and a bus-to-signal adapter.
- [ThreeBusPartitioned](./ThreeBusPartitioned.case.json): Two GENROU machines, with [partition A](./ThreeBusPartitionedA.case.json) and [partition B](./ThreeBusPartitionedB.case.json) for integration testing.
- [ThreeBusZipLoad](./ThreeBusZipLoad.case.json): Classical machines and a ZIP load.
- [TenGenPartitioned](./TenGenPartitioned.case.json): Ten finite GENROU machines and two loads, with [partition A](./TenGenPartitionedA.case.json) and [partition B](./TenGenPartitionedB.case.json) joined at bus 6.
