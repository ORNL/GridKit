# PhasorDynamics

The validation and simulation examples
use the reusable models in `cases/PhasorDynamics/`.

| Examples | Purpose |
| --- | --- |
| [Simulation](Simulation/README.md) | Examples using JSON inputs, including a short ACTIVSg10k run without disturbances. |
| `Validation/` | System studies compared against PowerWorld reference results. |

## Running a Study

Configure and build GridKit with SUNDIALS and Enzyme enabled. CMake copies each
study's solver file and the required case and reference files into the build
tree and packages them in the install tree.

From the installed `Validation/IEEE39/` directory, run:

```bash
DynamicSimulation IEEE39.solver.json
```

The solver file names its case using a local basename, such as
`"system_model_file": "IEEE39.case.json"`. Run the build-tree or installed
study, where CMake has supplied that case file. Monitor output is written
relative to the working directory.

Each validation folder contains the solver file, reference data, and comparison
figures. Its README links to the reusable case description.
