#include <assert.h>
#include <stdio.h>

#include <cmath>
#include <vector>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/LinearAlgebra/MemoryUtils.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>
#include <GridKit/Testing/Testing.hpp>

/**
 * @brief Standalone example that computes the sparse Jacobian of a vector-valued function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

using SparseMatrix = GridKit::LinearAlgebra::CooMatrix<double, size_t>;
using namespace GridKit::Enzyme::Sparse;

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
SparseMatrix* jac_f_ref(std::vector<ScalarT> x, std::vector<ScalarT> y)
{
  assert(x.size() == y.size());
  size_t   N           = x.size();
  size_t   nnz         = N * (N + 1) / 2; // lower triangular matrix, given the function above
  size_t*  rows        = new size_t[nnz];
  size_t*  cols        = new size_t[nnz];
  ScalarT* vals        = new ScalarT[nnz];
  size_t   current_nnz = 0;
  for (size_t idy = 0; idy < y.size(); ++idy)
  {
    for (size_t idx = 0; idx < x.size(); ++idx)
    {
      if (idy <= idx && current_nnz < nnz)
      {
        rows[current_nnz] = idx;
        cols[current_nnz] = idy;
        vals[current_nnz] = 2.0 * x[idy];
        current_nnz++;
      }
    }
  }
  // Hijacking constructor sets rows, cols and vals to nullptr
  return new SparseMatrix(N, N, current_nnz, &rows, &cols, &vals);
}

/// Function that computes the Jacobian via automatic differentiation
template <typename ScalarT>
__attribute__((noinline)) SparseMatrix* jac_f(size_t N, ScalarT* input)
{
  size_t* index_maps = new size_t[N];
  for (size_t i = 0; i < N; i++)
  {
    index_maps[i] = i;
  }

  size_t*  rows_buffer = new size_t[N * N];
  size_t*  cols_buffer = new size_t[N * N];
  ScalarT* vals_buffer = new ScalarT[N * N];
  size_t   current_nnz = 0;
  for (size_t i = 0; i < N; i++)
  {
    ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, size_t>,
                                                 (void*) ident_store<ScalarT, size_t>,
                                                 i);
    ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, size_t>,
                                                   (void*) sparse_store<ScalarT, size_t>,
                                                   i,
                                                   1.0,
                                                   index_maps,
                                                   index_maps,
                                                   rows_buffer,
                                                   cols_buffer,
                                                   vals_buffer,
                                                   &current_nnz);

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

  size_t*  rows = new size_t[current_nnz];
  size_t*  cols = new size_t[current_nnz];
  ScalarT* vals = new ScalarT[current_nnz];
  for (size_t ind = 0; ind < current_nnz; ++ind)
  {
    rows[ind] = rows_buffer[ind];
    cols[ind] = cols_buffer[ind];
    vals[ind] = vals_buffer[ind];
  }
  delete[] index_maps;
  delete[] rows_buffer;
  delete[] cols_buffer;
  delete[] vals_buffer;
  return new SparseMatrix(N, N, current_nnz, &rows, &cols, &vals);
}

/// Compare two sparse matrices
void check(SparseMatrix* matrix_1, SparseMatrix* matrix_2, int& fail)
{
  size_t  nnz_1  = matrix_1->getNnz();
  size_t* rows_1 = matrix_1->getRowData(GridKit::LinearAlgebra::memory::HOST);
  size_t* cols_1 = matrix_1->getColData(GridKit::LinearAlgebra::memory::HOST);
  double* vals_1 = matrix_1->getValues(GridKit::LinearAlgebra::memory::HOST);

  size_t  nnz_2  = matrix_2->getNnz();
  size_t* rows_2 = matrix_2->getRowData(GridKit::LinearAlgebra::memory::HOST);
  size_t* cols_2 = matrix_2->getColData(GridKit::LinearAlgebra::memory::HOST);
  double* vals_2 = matrix_2->getValues(GridKit::LinearAlgebra::memory::HOST);

  if (nnz_1 != nnz_2)
    fail++;

  for (size_t ind = 0; ind < nnz_1; ++ind)
  {
    std::cout << rows_1 << ", " << cols_1 << ", " << vals_1 << "\n";
    std::cout << rows_2 << ", " << cols_2 << ", " << vals_2 << "\n";
    if (rows_1[ind] != rows_2[ind])
      fail++;
    if (cols_1[ind] != cols_2[ind])
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
  SparseMatrix* dsq_ref = jac_f_ref(x, sq);

  /// Enzyme Jacobian
  SparseMatrix* dsq = jac_f<double>(N, x.data());

  /// Check
  int fail = 0;
  check(dsq, dsq_ref, fail);
  bool verbose = true;
  if (verbose)
  {
    std::cout << "Autodiff Jacobian\n";
    dsq->print();
    std::cout << "Reference Jacobian\n";
    dsq_ref->print();
  }

  delete dsq_ref;
  delete dsq;

  return fail;
}
