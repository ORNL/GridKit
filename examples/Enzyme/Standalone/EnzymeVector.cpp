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
int  enzyme_dupnoneed;
int  enzyme_dup;
int  enzyme_const;
void __enzyme_fwddiff(void*, int, std::vector<double>, std::vector<double>, int, std::vector<double>, std::vector<double>*);

inline double square_scalar(double x)
{
  return x * x;
}

inline double dsquare_ref_scalar(double x)
{
  return 2.0 * x;
}

// Vector-valued function to differentiate
void square(std::vector<double> x, std::vector<double>& y)
{
  for (int idx = 0; idx < x.size(); ++idx)
  {
    y[idx] = 0.0;
    for (int idy = 0; idy <= idx; idy++)
    {
      y[idx] += square_scalar(x[idy]);
    }
  }
}

// Reference Jacobian
void dsquare_ref(std::vector<double> x, std::vector<double> y, DenseMatrix& dy)
{
  for (int idy = 0; idy < y.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
    {
      if (idy <= idx)
        dy.setValue(idx, idy, dsquare_ref_scalar(x[idy]));
    }
  }
}

// Function that computes the Jacobian via automatic differentiation
void dsquare(std::vector<double> x, std::vector<double> y, DenseMatrix& dy)
{
  std::vector<double> v(x.size());
  std::vector<double> d_y(y.size());
  for (int idy = 0; idy < y.size(); ++idy)
  {
    // Elementary vector for Jacobian-vector product
    for (int idx = 0; idx < x.size(); ++idx)
    {
      v[idx] = 0.0;
    }
    v[idy] = 1.0;

    // Autodiff
    __enzyme_fwddiff((void*) square, enzyme_dup, x, v, enzyme_dupnoneed, y, &d_y);

    // Store result
    for (int idx = 0; idx < x.size(); ++idx)
    {
      dy.setValue(idx, idy, d_y[idx]);
    }
  }
}

int main()
{
  // Vector and matrix declarations
  constexpr int       N = 10;
  std::vector<double> x(N);
  std::vector<double> sq(N);
  DenseMatrix         dsq     = DenseMatrix(N, N);
  DenseMatrix         dsq_ref = DenseMatrix(N, N);

  // Random input values
  srand(time(NULL));
  for (int idx = 0; idx < x.size(); ++idx)
  {
    x[idx] = rand();
  }

  // Function evaluation
  square(x, sq);

  // Reference Jacobian
  dsquare_ref(x, sq, dsq_ref);

  // Enzyme Jacobian
  dsquare(x, sq, dsq);

  // Check
  int  fail    = 0;
  bool verbose = true;
  for (int idy = 0; idy < sq.size(); ++idy)
  {
    for (int idx = 0; idx < x.size(); ++idx)
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
