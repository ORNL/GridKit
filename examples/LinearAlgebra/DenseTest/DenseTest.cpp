
#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>

int main()
{
  int fail = 0;

  size_t                                              m = 4;
  size_t                                              n = 4;
  GridKit::LinearAlgebra::DenseMatrix<double, size_t> A =
      GridKit::LinearAlgebra::DenseMatrix<double, size_t>(m, n);

  double val = 0.0;
  for (size_t j = 0; j < n; ++j)
  {
    for (size_t i = 0; i < m; ++i)
    {
      A.setValue(i, j, val);
      val += 1.0;
    }
  }
  A.printMatrix("Dense matrix test output");
  (A.getValuesCOO())->printMatrix();

  return fail;
}
