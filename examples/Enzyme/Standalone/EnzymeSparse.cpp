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
template <typename T>
struct Triple
{
  size_t row;
  size_t col;
  T      val;
  Triple(Triple&&) = default;

  Triple(size_t row, size_t col, T val)
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

template <typename T>
__attribute__((always_inline)) static void sparse_store(T val, int64_t idx, size_t i, std::vector<Triple<T>>& triplets)
{
  if (val == 0.0)
    return;
  idx /= sizeof(T);
  if constexpr (sizeof(T) == 4)
    inner_storeflt(i, idx, val, triplets);
  else
    inner_storedbl(i, idx, val, triplets);
}

template <typename T>
__attribute__((always_inline)) static T sparse_load(int64_t idx, size_t i, std::vector<Triple<T>>& triplets)
{
  return 0.0;
}

template <typename T>
__attribute__((always_inline)) static void ident_store(T, int64_t idx, size_t i)
{
  assert(0 && "should never load");
}

template <typename T>
__attribute__((always_inline)) static T ident_load(int64_t idx, size_t i)
{
  idx /= sizeof(T);
  return (T) (idx == i);
}

/// Vector-valued function to differentiate
template <class ScalarT, typename IdxT>
__attribute__((always_inline)) static void f(size_t N, ScalarT* input, ScalarT* output)
{
  for (IdxT i = 0; i < N; i++)
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
__attribute__((noinline)) void jac_f(IdxT N, ScalarT* input, SparseMatrix& jac)
{
  std::vector<Triple<ScalarT>> triplets;
  for (IdxT i = 0; i < N; i++)
  {
    ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, i);
    ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, i, &triplets);

    __enzyme_fwddiff<void>((void*) f<ScalarT, IdxT>,
                           enzyme_const,
                           N,
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
  f<double, size_t>(x.size(), x.data(), sq.data());

  /// Reference Jacobian
  jac_f_ref<double, size_t>(x, sq, dsq_ref);

  /// Enzyme Jacobian
  jac_f<double, size_t>(N, x.data(), dsq);

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
