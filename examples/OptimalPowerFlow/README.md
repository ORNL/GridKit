# OptimalPowerFlow

| Examples | Purpose |
| --- | --- |
| [EconomicDispatch](EconomicDispatch/README.md) | AC optimal power flow of the cases in `cases/OptimalPowerFlow/`. |

## Running a Study

Build GridKit with Enzyme and Ipopt. For an installation under
`/path/to/install`:

```bash
cd /path/to/install/share/gridkit/examples/OptimalPowerFlow/EconomicDispatch
/path/to/install/bin/EconomicDispatch --case-file IEEE39.case.json
```

`--solver-tolerance` sets the Ipopt tolerance, 1e-8 by default. With
`--reference-cost` in \$/h, the study fails when the cost differs from it by
more than `--tolerance`, relative, 1e-6 by default.
