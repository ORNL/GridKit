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
template <typename T>
__attribute__((always_inline)) static void f(size_t N, T* input, T* output)
{
  for (size_t i = 0; i < N; i++)
  {
    output[i] = input[i] * input[i];
  }
}

/// Reference Jacobian
template <typename T>
void jac_f_ref(std::vector<T> x, std::vector<T> y, SparseMatrix& jac)
{
  std::vector<size_t> ctemp{};
  std::vector<size_t> rtemp{};
  std::vector<T>      valtemp{};
  for (int idy = 0; idy < y.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
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
template <typename T>
__attribute__((noinline)) void jac_f(size_t N, T* input, SparseMatrix& jac)
{
  std::vector<Triple<T>> triplets;
  for (size_t i = 0; i < N; i++)
  {
    T* output   = __enzyme_todense<T*>((void*) ident_load<T>, (void*) ident_store<T>, i);
    T* d_output = __enzyme_todense<T*>((void*) sparse_load<T>, (void*) sparse_store<T>, i, &triplets);

    __enzyme_fwddiff<void>((void*) f<T>,
                           enzyme_const,
                           N,
                           enzyme_dup,
                           input,
                           output,
                           enzyme_dupnoneed,
                           (T*) 0x1,
                           d_output);
  }

  std::vector<size_t> ctemp{};
  std::vector<size_t> rtemp{};
  std::vector<T>      valtemp{};
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
  for (int ind = 0; ind < vals_1.size(); ++ind)
  {
    if (rcord_1[ind] != rcord_2[ind])
      fail++;
    if (ccord_1[ind] != ccord_2[ind])
      fail++;
    if (std::abs(vals_1[ind] - vals_2[ind]) > std::numeric_limits<double>::epsilon())
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
  for (int i = 0; i < N; ++i)
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
