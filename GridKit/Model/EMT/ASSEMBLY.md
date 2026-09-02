# Assembly

An EMT model is defined by two equation groups. The internal equation fills the
model's internal DAE rows. The external equation contributes to the internal
rows of a connected model by accumulation. Each block of the assembled
Jacobian is one derivative of one group.

## Notation

Symbol | Description
------ | -----------
$\mathbf{y}$ | Internal variables of the model
$\dot{\mathbf{y}}$ | Internal variable derivatives
$\mathbf{y}^\text{ext}$ | External variables, internal to a connected model
$\dot{\mathbf{y}}^\text{ext}$ | External variable derivatives
$\mathbf{f}$ | Internal equation
$\mathbf{f}^\text{ext}$ | External equation
$\alpha$ | Time derivative coefficient supplied by the integrator

## System Form

A model fills its internal rows with its internal equation, and each connected
model $e$ adds its external equation to those same rows. A variable is internal
to exactly one model, and each scalar internal variable corresponds to one
internal row, so the assembled row set of a model is

```math
0 = \mathbf{f}(t,\mathbf{y},\dot{\mathbf{y}},\mathbf{y}^\text{ext}, \dot{\mathbf{y}}^\text{ext})
  + \sum_e \mathbf{f}^\text{ext}_e
```

The [Bus](Bus/README.md) is the concrete case. Its internal rows are
current-balance rows into which connected devices accumulate their external
equations. Model documentation marks such a contribution with an arrow,
$\mathbf{f} \leftarrow \dots$.

## Model Interface

A model implements two member functions. Both read the same inputs and differ
only in where the output lands.

Member function | Output | Placement
--------------- | ------ | ---------
`evaluateInternalResidual(y, yp, y_ext, yp_ext, f)` | $\mathbf{f}$ | Internal rows of the model
`evaluateExternalResidual(y, yp, y_ext, yp_ext, f_ext)` | $\mathbf{f}^\text{ext}$ | Accumulated into internal rows of a connected model

The internal equation fills its output, while the external equation accumulates
its contribution into `f_ext`. Both must be inlinable and reach state only
through their arguments so Enzyme can differentiate them in place.

## Local Jacobian

Rows are grouped by model, the internal rows first, and columns likewise, the
internal variables first. Value and derivative partials share the same block
structure, so $\alpha$ multiplies the second matrix once rather than each of its
terms.

```math
\mathbf{J} =
\begin{bmatrix}
\dfrac{\partial\mathbf{f}}{\partial\mathbf{y}}
&
\dfrac{\partial\mathbf{f}}{\partial\mathbf{y}^\text{ext}}
\\[3ex]
\dfrac{\partial\mathbf{f}^\text{ext}}{\partial\mathbf{y}}
&
\dfrac{\partial\mathbf{f}^\text{ext}}{\partial\mathbf{y}^\text{ext}}
\end{bmatrix}
+\ \alpha
\begin{bmatrix}
\dfrac{\partial\mathbf{f}}{\partial\dot{\mathbf{y}}}
&
\dfrac{\partial\mathbf{f}}{\partial\dot{\mathbf{y}}^\text{ext}}
\\[3ex]
\dfrac{\partial\mathbf{f}^\text{ext}}{\partial\dot{\mathbf{y}}}
&
\dfrac{\partial\mathbf{f}^\text{ext}}{\partial\dot{\mathbf{y}}^\text{ext}}
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
`Variable::Y` | $\partial\mathbf{f}/\partial\mathbf{y}$ | $\partial\mathbf{f}^\text{ext}/\partial\mathbf{y}$
`Variable::Yp` | $\partial\mathbf{f}/\partial\dot{\mathbf{y}}$ | $\partial\mathbf{f}^\text{ext}/\partial\dot{\mathbf{y}}$
`Variable::YExt` | $\partial\mathbf{f}/\partial\mathbf{y}^\text{ext}$ | $\partial\mathbf{f}^\text{ext}/\partial\mathbf{y}^\text{ext}$
`Variable::YpExt` | $\partial\mathbf{f}/\partial\dot{\mathbf{y}}^\text{ext}$ | $\partial\mathbf{f}^\text{ext}/\partial\dot{\mathbf{y}}^\text{ext}$

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
against $\dot{\mathbf{y}}^\text{ext}$ requires it to be an input, so both member
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
                               enzyme_const, model,
                               enzyme_dup, y, output,
                               enzyme_const, yp,
                               enzyme_const, y_ext,
                               enzyme_const, yp_ext,
                               enzyme_dupnoneed, seed, d_output);
      }
      else if constexpr (variable == Variable::Yp)
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const, model,
                               enzyme_const, y,
                               enzyme_dup, yp, output,
                               enzyme_const, y_ext,
                               enzyme_const, yp_ext,
                               enzyme_dupnoneed, seed, d_output);
      }
      else if constexpr (variable == Variable::YExt)
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const, model,
                               enzyme_const, y,
                               enzyme_const, yp,
                               enzyme_dup, y_ext, output,
                               enzyme_const, yp_ext,
                               enzyme_dupnoneed, seed, d_output);
      }
      else
      {
        __enzyme_fwddiff<void>(residual,
                               enzyme_const, model,
                               enzyme_const, y,
                               enzyme_const, yp,
                               enzyme_const, y_ext,
                               enzyme_dup, yp_ext, output,
                               enzyme_dupnoneed, seed, d_output);
      }
    }
  }
};
```

A model asks for the blocks it has. The scaling argument is $\alpha$ for the
derivative variables and defaults to one for the others.

```cpp
using GridKit::Enzyme::Sparse::Equation;
using GridKit::Enzyme::Sparse::SparseJacobian;
using GridKit::Enzyme::Sparse::Variable;

// Lower left of the value matrix, the external equation against internal variables
SparseJacobian<LoadZT, Equation::External, Variable::Y>::eval(
    this, n_ext, n_var, ext_indices, var_indices,
    y, yp, y_ext, yp_ext, rows, cols, vals, nnz);

// Same region of the derivative matrix
SparseJacobian<LoadZT, Equation::External, Variable::Yp>::eval(
    this, n_ext, n_var, ext_indices, var_indices,
    y, yp, y_ext, yp_ext, rows, cols, vals, nnz, alpha);
```
