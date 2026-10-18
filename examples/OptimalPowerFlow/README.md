# OptimalPowerFlow

The economic dispatch examples
use the reusable cases in `cases/OptimalPowerFlow/`.

| Examples | Purpose |
| --- | --- |
| [EconomicDispatch](EconomicDispatch/README.md) | AC optimal power flow studies. |

## Running a Study

Configure and build GridKit with Enzyme and Ipopt enabled. CMake copies each
case file into the build tree and packages it in the install tree.

For an installation under `/path/to/install`, run:

```bash
cd /path/to/install/share/gridkit/examples/OptimalPowerFlow/EconomicDispatch
/path/to/install/bin/EconomicDispatch --case-file IEEE39.case.json
```

`--solver-tolerance` sets the Ipopt convergence tolerance and defaults to 1e-8.

With `--reference-cost` in \$/h, the study fails when the optimal cost differs
from it by more than `--tolerance`, relative, which defaults to 1e-6.
