# TwoArea Hierarchical Test

`TwoAreaHierarchical.case.json` has the same system parameters, buses and devices as the [TwoArea case](../../../../cases/PhasorDynamics/TwoArea/README.md) ([`TwoArea.case.json`](../../../../cases/PhasorDynamics/TwoArea/TwoArea.case.json)), split into the [subsystems](../../../../GridKit/Model/PhasorDynamics/INPUT_FORMAT.md#subsystems) below. The test `PhasorDynamicsIntegrationTwoAreaHierarchical` compares the two files.

System | Buses | Devices
-- | -- | --
Top level | 8 | None
`Area1` | 1, 2, 5, 6, 7 | Generators on buses 1 and 2 with their exciters and governors, the load on bus 7, and the branches between these buses and to bus 8
`Area2` | 3, 4, 9, 10 | Generators on buses 3 and 4 with their exciters and governors, the load on bus 8, and the branches between these buses and to bus 8

Bus 8 connects the two areas, so it stays at the top level, where both areas can use it. The load on bus 8 is in `Area2` because a system with subsystems cannot have devices of its own.

The hierarchical file is a separate copy of the data, so a change to the case must be made in both files.
