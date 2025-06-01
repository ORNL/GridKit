#pragma once

#include <vector>

extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

template <typename T, typename... Tys>
extern T __enzyme_fwddiff(void*, Tys...) noexcept;

template <typename T, typename... Tys>
extern T __enzyme_todense(Tys...) noexcept;

/// Sparse storage for Enzyme
template <typename ScalarT>
struct Triple
{
  size_t  row;
  size_t  col;
  ScalarT val;
  Triple(Triple&&) = default;

  Triple(size_t row, size_t col, ScalarT val)
    : row(row),
      col(col),
      val(val)
  {
  }
};

__attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(int64_t row, int64_t col, float val, std::vector<Triple<float>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

__attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(int64_t row, int64_t col, double val, std::vector<Triple<double>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

template <typename ScalarT>
__attribute__((always_inline)) static void sparse_store(ScalarT val, int64_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
{
  if (val == 0.0)
    return;
  idx /= sizeof(ScalarT);
  if constexpr (sizeof(ScalarT) == 4)
    inner_storeflt(i, idx, val, triplets);
  else
    inner_storedbl(i, idx, val, triplets);
}

template <typename ScalarT>
__attribute__((always_inline)) static ScalarT sparse_load(int64_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
{
  return 0.0;
}

template <typename ScalarT>
__attribute__((always_inline)) static void ident_store(ScalarT, int64_t idx, size_t i)
{
  assert(0 && "should never load");
}

template <typename ScalarT>
__attribute__((always_inline)) static ScalarT ident_load(int64_t idx, size_t i)
{
  idx /= sizeof(ScalarT);
  return (ScalarT) (idx == i);
}

/// Wrapper around model residual
template <typename T, typename ScalarT>
void residual_wrapper(T* obj, ScalarT* y, ScalarT* f)
{
  obj->evaluateResidualLocally(y, f);
}

/// Function that computes the Jacobian via automatic differentiation
template <typename T, class ScalarT, typename IdxT>
__attribute__((noinline)) void EnzymeSparseModelJacobian(T* model, IdxT n, ScalarT* input, GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
{
  std::vector<Triple<ScalarT>> triplets;
  for (IdxT i = 0; i < n; i++)
  {
    ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, i);
    ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, i, &triplets);

    __enzyme_fwddiff<void>((void*) residual_wrapper<T, ScalarT>,
                           enzyme_const,
                           model,
                           enzyme_dup,
                           input,
                           output,
                           enzyme_dupnoneed,
                           (ScalarT*) 0x1,
                           d_output);
  }

  std::vector<IdxT>    ctemp{};
  std::vector<IdxT>    rtemp{};
  std::vector<ScalarT> valtemp{};
  for (auto& tup : triplets)
  {
    rtemp.push_back(tup.row);
    ctemp.push_back(tup.col);
    valtemp.push_back(tup.val);
  }
  jac.setValues(rtemp, ctemp, valtemp);
}

