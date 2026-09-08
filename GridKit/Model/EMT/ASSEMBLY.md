# Assembly

EMT assembly supports two equation groups. Internal equations fill a model's
DAE rows. Embedded operators contribute to their consumer's equations through
external rows. Each Jacobian block differentiates one equation group.
Electrical devices expose terminal-current signals; only KCL sums those
currents into the bus equations.

## Notation

Symbol | Description
------ | -----------
$\mathbf{y}$ | Internal variables of the model
$\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}$ | Internal variable derivatives
$\mathbf{y}^\mathrm{ext}$ | External variables, internal to a connected model
$\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}$ | External variable derivatives
$\mathbf{f}$ | Internal equation
$\mathbf{f}^\mathrm{ext}$ | External equation
$\alpha$ | Time derivative coefficient supplied by the integrator

## System Form

A model fills its internal rows with its internal equation, and each connected
model $e$ adds its external equation to those same rows. A variable is internal
to exactly one model, and each scalar internal variable corresponds to one
internal row, so the assembled row set of a model is

```math
0 = \mathbf{f}(t,\mathbf{y},\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t},\mathbf{y}^\mathrm{ext}, \dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t})
  + \sum_e \mathbf{f}^\mathrm{ext}_e
```

The [Bus](Component/Bus/README.md) contains `KCL` and Norton sources.
`KCL` owns both current-balance residuals and their Jacobian contributions.
During wiring, `Bus::addCurrent(phase, signal, sign)` registers each current,
with positive sign meaning injection into the bus. KCL reads every registered
signal and appends its gradient with the same sign. Norton terminals register
incident currents positively and shunt currents negatively; other devices
register their existing branch-current outputs. Registration closes at allocation.
Devices never write directly into bus residuals.

## DAE validation

At the executable root, `tagDifferentiable()` derives differential columns
from the assembled `F_yp`, including contributions to variables owned by
other components. EMT derivative coefficients must be independent of time
and state within a fixed topology; zero parameter coefficients and exact
cancellation are removed after accumulation. Classification is then
distributed to every child and embedded operator.

For differential indices `d` and algebraic indices `a`, initialization must
solve for `(yp[d], y[a])` with `y[d]` fixed. The system checks
`M = [F_yp(:, d), F_y(:, a)]` using a sparse maximum matching followed by
equilibrated KLU factorization. Structural slots in `F_y` are retained even
when their current coefficient is zero. A reciprocal pivot-condition estimate
at or below `size() * epsilon` is treated as numerical singularity at working
precision. Diagnostics distinguish structural deficiency from numerical
failure at the supplied state and identify component paths and local indices.

This checks local solvability in the existing coordinates. It neither proves
regularity at all future states nor performs constraint reduction. IDA refreshes
the classification and validation before each restart.

## Operating-point initialization

Components declare the signals they read during initialization, their
initializable outputs, and any outputs whose operating point they require.
The system resolves dependencies through signal gradients and initializes
producers before consumers. Machine mechanical-power and field-voltage
requirements precede their governor and exciter initializers.

`InitialState` reconciles these requirements with prescribed outputs and
declared constants before any state is changed. Each producer then initializes
its own variables using the resolved output values. Unknown state paths,
unknown output names, conflicting requirements, and cyclic initialization
dependencies are errors. Computed outputs own no state; their prescribed values
are checked after the stateful producers initialize.

## Accepted-step history

History belongs to the model or operator that uses it. IDA's
`initializeSimulation(t0)` starts a fresh study and calls `resetHistory()` after
loading its initial state. `restartSimulation(t0)` retains history for an event
at the current time. Restoring saved initial conditions and initializing again
starts another fresh study without reconfiguring the solver. Configuration also
resets history before the initial structural Jacobian evaluation.

`acceptStep(t)` records the consistent initial state and each
accepted internal solver step; trial evaluations and interpolated monitor
samples are never committed. A discontinuity preserves prior history and
commits a new right limit at the event time.

`maximumStepSize()` bounds forward steps, for example by the shortest transport
delay. Containers and components propagate history notifications and take the
smallest child or operator bound. `Delay` supports both this bounded
method-of-steps realization and implicit overlap within a larger adaptive
step; `Propagation` composes rational matrices with that history realization.

`nextDiscontinuityTime(t)` reports the next known delayed arrival. IDA stops
there and calls `beginDiscontinuity(t)` before calculating consistent
right-limit conditions. `Delay` is the only source of such discontinuities:
it implements the `HistoryDiscontinuity` interface, whose private constructor
admits no other model, and the component query is final and only aggregates
that interface through owned operators and container children. Every other
model is smooth by design. Automatic history restarts are supported by
forward DAE simulation; quadrature and adjoint checkpoint replay
do not support these restarts.

## Model Interface

A model implements two member functions. Both read the same inputs and differ
only in where the output lands.

Member function | Output | Placement
--------------- | ------ | ---------
`evaluateInternalResidual(y, yp, y_ext, yp_ext, f)` | $\mathbf{f}$ | Internal rows of the model
`evaluateExternalResidual(y, yp, y_ext, yp_ext, f_ext)` | $\mathbf{f}^\mathrm{ext}$ | Accumulated into internal rows of a connected model

The internal equation fills its output, while the external equation accumulates
its contribution into `f_ext`. Both must be inlinable and reach state only
through their arguments so Enzyme can differentiate them in place.

## Local Jacobian

`assembleJacobian(y_scale, yp_scale)` uses the same local derivatives to form
`y_scale * F_y + yp_scale * F_yp`. The evaluator's `evaluateJacobian()` supplies
`(1, alpha)` for IDA; `(1, 0)` and `(0, 1)` obtain the individual partials.
Changing the active blocks invalidates the cached sparse layout. Contributions
from children, embedded operators, and computed-signal gradients use the same
coefficients.

Rows are grouped by model, the internal rows first, and columns likewise, the
internal variables first. Value and derivative partials share the same block
structure, so $\alpha$ multiplies the second matrix once rather than each of its
terms.

```math
\mathbf{J} =
\begin{bmatrix}
\dfrac{\partial\mathbf{f}}{\partial\mathbf{y}}
&
\dfrac{\partial\mathbf{f}}{\partial\mathbf{y}^\mathrm{ext}}
\\[3ex]
\dfrac{\partial\mathbf{f}^\mathrm{ext}}{\partial\mathbf{y}}
&
\dfrac{\partial\mathbf{f}^\mathrm{ext}}{\partial\mathbf{y}^\mathrm{ext}}
\end{bmatrix}
+\ \alpha
\begin{bmatrix}
\dfrac{\partial\mathbf{f}}{\partial\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}}
&
\dfrac{\partial\mathbf{f}}{\partial\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}}
\\[3ex]
\dfrac{\partial\mathbf{f}^\mathrm{ext}}{\partial\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}}
&
\dfrac{\partial\mathbf{f}^\mathrm{ext}}{\partial\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}}
\end{bmatrix}
```

Region | Equation | Variable
------ | -------- | --------
Upper left | `Equation::Internal` | `Variable::Y`, `Variable::Yp`
Upper right | `Equation::Internal` | `Variable::YExt`, `Variable::YpExt`
Lower left | `Equation::External` | `Variable::Y`, `Variable::Yp`
Lower right | `Equation::External` | `Variable::YExt`, `Variable::YpExt`

Every entry is a contribution of the one model being described, and the
assembled Jacobian is the sum of these contributions over all models.

## Evaluator Design

The evaluator lives in
[SparseJacobian.hpp](../../AutomaticDifferentiation/Enzyme/SparseJacobian.hpp),
and the wrapper is named `ResidualWrapper` to leave the PhasorDynamics `ModelWrapper` intact.

Two independent choices select a block, which equation is differentiated and
which variable it is differentiated against. Both are enums, so one template
covers every block and each of the eight is one pair of enum values.

Variable | `Equation::Internal` | `Equation::External`
-------- | -------------------- | --------------------
`Variable::Y` | $\partial\mathbf{f}/\partial\mathbf{y}$ | $\partial\mathbf{f}^\mathrm{ext}/\partial\mathbf{y}$
`Variable::Yp` | $\partial\mathbf{f}/\partial\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}$ | $\partial\mathbf{f}^\mathrm{ext}/\partial\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}$
`Variable::YExt` | $\partial\mathbf{f}/\partial\mathbf{y}^\mathrm{ext}$ | $\partial\mathbf{f}^\mathrm{ext}/\partial\mathbf{y}^\mathrm{ext}$
`Variable::YpExt` | $\partial\mathbf{f}/\partial\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}$ | $\partial\mathbf{f}^\mathrm{ext}/\partial\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}$

The derivative variables carry $\alpha$ as a value scaling, so the two matrices
above assemble into one Jacobian.

```cpp
enum class Equation
{
  Internal,
  External
};

enum class Variable
{
  Y,
  Yp,
  YExt,
  YpExt
};
```

The wrapper selects the member function from the equation group. Differentiating
against $\dfrac{\mathrm{d}\mathbf{y}^\mathrm{ext}}{\mathrm{d}t}$ requires it to be an input, so both member
functions take it alongside the other three.

```cpp
template <typename ModelT, Equation equation>
struct ResidualWrapper
{
};

template <typename ModelT>
struct ResidualWrapper<ModelT, Equation::Internal>
{
  using ScalarT = typename ModelT::ScalarT;

  static void eval(ModelT*        model,
                   const ScalarT* y,
                   const ScalarT* yp,
                   const ScalarT* y_ext,
                   const ScalarT* yp_ext,
                   ScalarT*       f)
  {
    model->evaluateInternalResidual(y, yp, y_ext, yp_ext, f);
  }
};
```

The `Equation::External` specialization mirrors it and calls
`evaluateExternalResidual`.

The evaluator body is the existing one with two compile-time choices, which
input carries the seed and what the stored value is scaled by. Everything else
is shared.

```cpp
template <typename ModelT, Equation equation, Variable variable>
struct SparseJacobian
{
  using ScalarT = typename ModelT::ScalarT;
  using IdxT    = typename ModelT::IdxT;
  using RealT   = typename ModelT::RealT;

  static void eval(ModelT*        model,
                   const size_t   n_res,
                   const size_t   n_var,
                   const IdxT*    res_indices,
                   const IdxT*    var_indices,
                   const ScalarT* y,
                   const ScalarT* yp,
                   const ScalarT* y_ext,
                   const ScalarT* yp_ext,
                   IdxT*          rows,
                   IdxT*          cols,
                   RealT*         vals,
                   IdxT&          nnz,
                   const RealT    scaling = 1.0)
  {
    if (n_res == 0 || n_var == 0)
    {
      return;
    }

    std::vector<ScalarT> elementary_v(n_var);
    for (size_t var_i = 0; var_i < n_var; ++var_i)
    {
      // Sparse storage. @see LowerSparseStorage.hpp
      ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                   (void*) ident_store<ScalarT, IdxT>,
                                                   var_i);
      ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                     (void*) sparse_store<ScalarT, IdxT>,
                                                     var_i,
                                                     scaling,
                                                     res_indices,
                                                     var_indices,
                                                     rows,
                                                     cols,
                                                     vals,
                                                     &nnz);

      // Elementary vector for Jacobian-vector product
      std::ranges::fill(elementary_v, 0.0);
      elementary_v[var_i] = 1.0;

      void* residual = (void*) ResidualWrapper<ModelT, equation>::eval;
      auto  seed     = elementary_v.data();

      if constexpr (variable == Variable::Y)
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const,
                               model,
                               enzyme_dup,
                               y,
                               output,
                               enzyme_const,
                               yp,
                               enzyme_const,
                               y_ext,
                               enzyme_const,
                               yp_ext,
                               enzyme_dupnoneed,
                               seed,
                               d_output);
      }
      else if constexpr (variable == Variable::Yp)
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const,
                               model,
                               enzyme_const,
                               y,
                               enzyme_dup,
                               yp,
                               output,
                               enzyme_const,
                               y_ext,
                               enzyme_const,
                               yp_ext,
                               enzyme_dupnoneed,
                               seed,
                               d_output);
      }
      else if constexpr (variable == Variable::YExt)
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const,
                               model,
                               enzyme_const,
                               y,
                               enzyme_const,
                               yp,
                               enzyme_dup,
                               y_ext,
                               output,
                               enzyme_const,
                               yp_ext,
                               enzyme_dupnoneed,
                               seed,
                               d_output);
      }
      else
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const,
                               model,
                               enzyme_const,
                               y,
                               enzyme_const,
                               yp,
                               enzyme_const,
                               y_ext,
                               enzyme_dup,
                               yp_ext,
                               output,
                               enzyme_dupnoneed,
                               seed,
                               d_output);
      }
    }
  }
};
```

A model asks for the blocks it has, passing `y_scale` for value variables and
`yp_scale` for derivative variables. Disabled blocks are omitted.

```cpp
using GridKit::Enzyme::Sparse::Equation;
using GridKit::Enzyme::Sparse::SparseJacobian;
using GridKit::Enzyme::Sparse::Variable;

// Lower left of the value matrix, the external equation against internal variables
if (y_scale != 0)
  SparseJacobian<ModelT, Equation::External, Variable::Y>::eval(
      this, n_ext, n_var, ext_indices, var_indices, y, yp, y_ext, yp_ext, rows, cols, vals, nnz, y_scale);

// Same region of the derivative matrix
if (yp_scale != 0)
  SparseJacobian<ModelT, Equation::External, Variable::Yp>::eval(
      this, n_ext, n_var, ext_indices, var_indices, y, yp, y_ext, yp_ext, rows, cols, vals, nnz, yp_scale);
```
