#include <assert.h>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <vector>

#include <LinearAlgebra/SparseMatrix/COO_Matrix.hpp>

/**
 * @brief Standalone example that computes the sparse Jacobian of a vector-valued function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Modify the sparse storage provided to Enzyme to directly operate on std::vector and GridKit::LinearAlgebra::COO_Matrix
 * TODO: Convert this into a unit test.
 */

using SparseMatrix = GridKit::LinearAlgebra::COO_Matrix<double, size_t>;
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

/// Vector-valued function to differentiate
template <class ScalarT, typename IdxT>
__attribute__((always_inline)) static void f(size_t n, ScalarT* input, ScalarT* output)
{
  for (IdxT i = 0; i < n; i++)
  {
    output[i] = input[i] * input[i];
  }
}

/// Reference Jacobian
template <class ScalarT, typename IdxT>
void jac_f_ref(std::vector<ScalarT> x, std::vector<ScalarT> y, SparseMatrix& jac)
{
  std::vector<IdxT>    ctemp{};
  std::vector<IdxT>    rtemp{};
  std::vector<ScalarT> valtemp{};
  for (IdxT idy = 0; idy < y.size(); ++idy)
  {
    for (IdxT idx = 0; idx < x.size(); ++idx)
    {
      if (idx == idy)
      {
        rtemp.push_back(idx);
        ctemp.push_back(idy);
        valtemp.push_back(2.0 * x[idx]);
      }
    }
  }
  jac.setValues(rtemp, ctemp, valtemp);
}

/// Function that computes the Jacobian via automatic differentiation
template <class ScalarT, typename IdxT>
__attribute__((noinline)) void jac_f(IdxT n, ScalarT* input, SparseMatrix& jac)
{
  std::vector<Triple<ScalarT>> triplets;
  for (IdxT i = 0; i < n; i++)
  {
    ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, i);
    ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, i, &triplets);

    __enzyme_fwddiff<void>((void*) f<ScalarT, IdxT>,
                           enzyme_const,
                           n,
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

/// Compare two sparse matrices
template <class ScalarT, typename IdxT>
void check(SparseMatrix matrix_1, SparseMatrix matrix_2, int& fail)
{
  std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> entries_1 = matrix_1.getEntries();
  const auto [rcord_1, ccord_1, vals_1]                                               = entries_1;
  std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> entries_2 = matrix_2.getEntries();
  const auto [rcord_2, ccord_2, vals_2]                                               = entries_2;
  for (IdxT ind = 0; ind < vals_1.size(); ++ind)
  {
    if (rcord_1[ind] != rcord_2[ind])
      fail++;
    if (ccord_1[ind] != ccord_2[ind])
      fail++;
    if (std::abs(vals_1[ind] - vals_2[ind]) > std::numeric_limits<ScalarT>::epsilon())
      fail++;
  }
}

int main()
{
  /// Vector and matrix declarations
  size_t              n = 5;
  std::vector<double> x(n);
  std::vector<double> sq(n);
  SparseMatrix        dsq     = SparseMatrix(n, n);
  SparseMatrix        dsq_ref = SparseMatrix(n, n);

  /// Input initialization
  double val = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    x[i]  = val;
    val  += 1.0;
  }

  /// Function evaluation
  f<double, size_t>(x.size(), x.data(), sq.data());

  /// Reference Jacobian
  jac_f_ref<double, size_t>(x, sq, dsq_ref);

  /// Enzyme Jacobian
  jac_f<double, size_t>(n, x.data(), dsq);

  /// Check
  int fail = 0;
  check<double, size_t>(dsq, dsq_ref, fail);
  bool verbose = true;
  if (verbose)
  {
    dsq.printMatrix("Autodiff Jacobian");
    dsq_ref.printMatrix("Reference Jacobian");
  }

  return fail;
}
