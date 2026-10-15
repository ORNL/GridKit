# OptimalPowerFlow

The economic dispatch examples
use the reusable MATPOWER cases in `cases/OptimalPowerFlow/`.

| Examples | Purpose |
| --- | --- |
| [EconomicDispatch](EconomicDispatch/README.md) | AC optimal power flow studies. |

## Running a Study

Configure and build GridKit with Enzyme and Ipopt enabled. CMake copies each
case file into the build tree and packages it in the install tree.

For an installation under `/path/to/install`, run:

```bash
cd /path/to/install/share/gridkit/examples/OptimalPowerFlow/EconomicDispatch
/path/to/install/bin/EconomicDispatch IEEE39.m
```

A second argument, a reference cost in \$/h, makes the study fail when the
optimal cost differs from it.
