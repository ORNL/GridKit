# PhasorDynamics

The validation and simulation examples
use the reusable models in `cases/PhasorDynamics/`.

| Examples | Purpose |
| --- | --- |
| [DynamicSimulation](DynamicSimulation/README.md) | Individual simulation studies. |
| [ContingencyAnalysis](ContingencyAnalysis/README.md) | Bus-fault contingency studies. |
| [Validation](Validation/README.md) | System studies compared against PowerWorld reference results. |
| [OptimalDispatch](OptimalDispatch/README.md) | Optimal operating points that start PhasorDynamics in steady state. |

## Running a Study

Configure and build GridKit with SUNDIALS and Enzyme enabled. CMake copies each
study's solver file and the required case and reference files into the build
tree and packages them in the install tree.

For an installation under `/path/to/install`, run:

```bash
cd /path/to/install/share/gridkit/examples/PhasorDynamics/Validation/IEEE39
/path/to/install/bin/DynamicSimulation IEEE39.solver.json
```

The solver file names its case using a local basename, such as
`"system_model_file": "IEEE39.case.json"`. Run the build-tree or installed
study, where CMake has supplied that case file. Monitor output is written
relative to the working directory.

Each validation folder contains the solver file, reference data, and comparison
figures. Its README links to the reusable case description.
