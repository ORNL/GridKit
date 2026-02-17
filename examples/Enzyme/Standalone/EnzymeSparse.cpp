#include <assert.h>
#include <stdio.h>

#include <cmath>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>
#include <GridKit/Testing/Testing.hpp>

/**
 * @brief Standalone example that computes the sparse Jacobian of a vector-valued function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Modify the sparse storage provided to Enzyme to directly operate on std::vector and COO_Matrix
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

[[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(size_t row, size_t col, float val, std::vector<Triple<float>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

[[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(size_t row, size_t col, double val, std::vector<Triple<double>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

template <typename ScalarT>
__attribute__((always_inline)) static void sparse_store(ScalarT val, size_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
{
  if (val == 0.0)
    return;
  idx /= sizeof(ScalarT);
  if constexpr (sizeof(ScalarT) == 4)
    inner_storeflt(idx, i, val, triplets);
  else
    inner_storedbl(idx, i, val, triplets);
}

template <typename ScalarT>
__attribute__((always_inline)) static ScalarT sparse_load(size_t, size_t, std::vector<Triple<ScalarT>>&)
{
  return 0.0;
}

template <typename ScalarT>
__attribute__((always_inline)) static void ident_store(ScalarT, size_t, size_t)
{
  assert(0 && "should never store");
}

template <typename ScalarT>
__attribute__((always_inline)) static ScalarT ident_load(size_t idx, size_t i)
{
  idx /= sizeof(ScalarT);
  return (ScalarT) (idx == i);
}

/// Vector-valued function to differentiate
template <typename ScalarT>
__attribute__((always_inline)) static void f(size_t N, ScalarT* input, ScalarT* output)
{
  for (size_t idx = 0; idx < N; ++idx)
  {
    output[idx] = 0.0;
    for (size_t idy = 0; idy <= idx; idy++)
    {
      output[idx] += input[idy] * input[idy];
    }
  }
}

/// Reference Jacobian
template <typename ScalarT>
void jac_f_ref(std::vector<ScalarT> x, std::vector<ScalarT> y, SparseMatrix& jac)
{
  std::vector<size_t>  ctemp{};
  std::vector<size_t>  rtemp{};
  std::vector<ScalarT> valtemp{};
  for (size_t idy = 0; idy < y.size(); ++idy)
  {
    for (size_t idx = 0; idx < x.size(); ++idx)
    {
      if (idy <= idx)
      {
        rtemp.push_back(idx);
        ctemp.push_back(idy);
        valtemp.push_back(2.0 * x[idy]);
      }
    }
  }
  jac.setValues(rtemp, ctemp, valtemp);
}

/// Function that computes the Jacobian via automatic differentiation
template <typename ScalarT>
__attribute__((noinline)) void jac_f(size_t N, ScalarT* input, SparseMatrix& jac)
{
  std::vector<Triple<ScalarT>> triplets;
  for (size_t i = 0; i < N; i++)
  {
    ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, i);
    ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, i, &triplets);

    __enzyme_fwddiff<void>((void*) f<ScalarT>,
                           enzyme_const,
                           N,
                           enzyme_dup,
                           input,
                           output,
                           enzyme_dupnoneed,
                           (ScalarT*) 0x1,
                           d_output);
  }

  std::vector<size_t>  ctemp{};
  std::vector<size_t>  rtemp{};
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
void check(SparseMatrix matrix_1, SparseMatrix matrix_2, int& fail)
{
  std::tuple<std::vector<size_t>&, std::vector<size_t>&, std::vector<double>&> entries_1 = matrix_1.getEntries();
  const auto [rcord_1, ccord_1, vals_1]                                                  = entries_1;
  std::tuple<std::vector<size_t>&, std::vector<size_t>&, std::vector<double>&> entries_2 = matrix_2.getEntries();
  const auto [rcord_2, ccord_2, vals_2]                                                  = entries_2;
  for (size_t ind = 0; ind < vals_1.size(); ++ind)
  {
    if (rcord_1[ind] != rcord_2[ind])
      fail++;
    if (ccord_1[ind] != ccord_2[ind])
      fail++;
    if (!GridKit::Testing::isEqual(vals_1[ind], vals_2[ind]))
      fail++;
  }
}

int main()
{
  /// Vector and matrix declarations
  size_t              N = 5;
  std::vector<double> x(N);
  std::vector<double> sq(N);
  SparseMatrix        dsq     = SparseMatrix(N, N);
  SparseMatrix        dsq_ref = SparseMatrix(N, N);

  /// Input initialization
  double val = 0.0;
  for (size_t i = 0; i < N; ++i)
  {
    x[i]  = val;
    val  += 1.0;
  }

  /// Function evaluation
  f(x.size(), x.data(), sq.data());

  /// Reference Jacobian
  jac_f_ref(x, sq, dsq_ref);

  /// Enzyme Jacobian
  jac_f<double>(N, x.data(), dsq);

  /// Check
  int fail = 0;
  check(dsq, dsq_ref, fail);
  bool verbose = true;
  if (verbose)
  {
    dsq.printMatrix("Autodiff Jacobian");
    dsq_ref.printMatrix("Reference Jacobian");
  }

  return fail;
}
