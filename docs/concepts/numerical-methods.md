# Numerical Methods

GridKit provides dynamic integrators and steady-state nonlinear solvers.

| Task | Implementations and reference |
| --- | --- |
| Differential-algebraic equations (DAEs) | [`IDA`](../api/reference/classAnalysisManager_1_1Sundials_1_1Ida.rst) |
| Native dynamic integration | [`Rosenbrock`](../api/reference/classAnalysisManager_1_1NativeDynamicSolver_1_1Rosenbrock.rst) |
| Step-size policies | [`AdaptiveStep`](../api/reference/classAnalysisManager_1_1NativeDynamicSolver_1_1AdaptiveStep.rst), [`FixedStep`](../api/reference/classAnalysisManager_1_1NativeDynamicSolver_1_1FixedStep.rst) |
| Steady-state nonlinear systems | [`Kinsol`](../api/reference/classAnalysisManager_1_1Sundials_1_1Kinsol.rst) |

See the [API reference](../api.md) for types and parameters and the
[examples](../examples/README.md) for solver configurations.
