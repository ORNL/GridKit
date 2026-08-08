# Architecture

GridKit separates applications, case data, model equations, and numerical
solvers. Cases and models can be reused by different applications.

| Layer | Responsibility | Documentation |
| --- | --- | --- |
| Applications | Configure and run simulations and analyses. | [Apps](../applications/index.md) |
| Cases | Systems and model instances. | [Cases](../cases/index.md) |
| Models | Equations, parameters, ports, initialization, and outputs. | [Models](../models/index.md) |
| Solvers | Dynamic integration and steady-state solution. | [API](../api/index.md) |

The solver-facing boundary is the
[`GridKit::Model::Evaluator`](../api/reference/classGridKit_1_1Model_1_1Evaluator.rst)
interface. Each model family retains its own equations and data contract.
