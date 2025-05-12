#include <iostream>
#include <limits>
#include <vector>

#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>

/**
 * @brief Standalone example that computes the Jacobian of a vector-valued function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

using DenseMatrix = GridKit::LinearAlgebra::DenseMatrix<double, size_t>;
int enzyme_dupnoneed;
int enzyme_dup;
int enzyme_const;
template <class ScalarT, typename IdxT>
void __enzyme_fwddiff(void*, int, std::vector<ScalarT>, std::vector<ScalarT>, int, std::vector<ScalarT>, std::vector<ScalarT>*);

template <class ScalarT>
inline ScalarT square_scalar(ScalarT x)
{
  return x * x;
}

template <class ScalarT>
inline ScalarT dsquare_ref_scalar(ScalarT x)
{
  return 2.0 * x;
}

/// Vector-valued function to differentiate
template <class ScalarT, typename IdxT>
void square(std::vector<ScalarT> x, std::vector<ScalarT>& y)
{
  for (IdxT idx = 0; idx < x.size(); ++idx)
  {
    y[idx] = square_scalar<ScalarT>(x[idx]);
  }
}

/// Reference Jacobian
template <class ScalarT, typename IdxT>
void dsquare_ref(std::vector<ScalarT> x, std::vector<ScalarT> y, DenseMatrix& dy)
{
  for (IdxT idy = 0; idy < y.size(); ++idy)
  {
    for (IdxT idx = 0; idx < x.size(); ++idx)
    {
      if (idx == idy)
        dy.setValue(idx, idy, dsquare_ref_scalar<ScalarT>(x[idx]));
    }
  }
}

/// Function that computes the Jacobian via automatic differentiation
template <class ScalarT, typename IdxT>
void dsquare(std::vector<ScalarT> x, std::vector<ScalarT> y, DenseMatrix& dy)
{
  std::vector<ScalarT> v(x.size());
  std::vector<ScalarT> d_y(y.size());
  for (IdxT idy = 0; idy < y.size(); ++idy)
  {
    /// Elementary vector for Jacobian-vector product
    for (IdxT idx = 0; idx < x.size(); ++idx)
    {
      v[idx] = 0.0;
    }
    v[idy] = 1.0;

    /// Autodiff
    __enzyme_fwddiff<ScalarT, IdxT>((void*) square<ScalarT, IdxT>, 
                                    enzyme_dup, x, v, 
                                    enzyme_dupnoneed, y, &d_y);

    /// Store result
    for (IdxT idx = 0; idx < x.size(); ++idx)
    {
      dy.setValue(idx, idy, d_y[idx]);
    }
  }
}

int main()
{
  /// Vector and matrix declarations
  constexpr size_t       N = 10;
  std::vector<double> x(N);
  std::vector<double> sq(N);
  DenseMatrix         dsq     = DenseMatrix(N, N);
  DenseMatrix         dsq_ref = DenseMatrix(N, N);

  /// Random input values
  srand(time(NULL));
  for (size_t idx = 0; idx < x.size(); ++idx)
  {
    x[idx] = rand();
  }

  /// Function evaluation
  square<double, size_t>(x, sq);

  /// Reference Jacobian
  dsquare_ref<double, size_t>(x, sq, dsq_ref);

  /// Enzyme Jacobian
  dsquare<double, size_t>(x, sq, dsq);

  /// Check
  int  fail    = 0;
  bool verbose = true;
  for (size_t idy = 0; idy < sq.size(); ++idy)
  {
    for (size_t idx = 0; idx < x.size(); ++idx)
    {
      if (std::abs(dsq.getValue(idx, idy) - dsq_ref.getValue(idx, idy)) > std::numeric_limits<double>::epsilon())
      {
        fail++;
        if (verbose)
        {
          std::cout << "Result incorrect at line = " << idy << ", column = " << idx << "\n";
          std::cout << "x = " << x[idx] << ", x^2 = " << sq[idx] << ", d(x^2)/dx = " << dsq.getValue(idx, idy) << "\n";
        }
      }
    }
  }
  if (verbose)
  {
    dsq.printMatrix("Autodiff Jacobian");
    dsq_ref.printMatrix("Reference Jacobian");
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
