#include <iostream>

#include <GridKit/Utilities/Testing.hpp>

#include "VectorModel.hpp"

/**
 * @brief Example that computes the Jacobian of a vector-valued residual
 * (implemented as the member function of a class and operating directly on class members)
 * by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

inline double dsquare_ref_scalar(double x)
{
  return 2.0 * x;
}

// Reference Jacobian
DenseMatrix dsquare_ref(std::vector<double> x, std::vector<double> y)
{
  DenseMatrix jac(x.size(), y.size());
  for (size_t idy = 0; idy < y.size(); ++idy)
  {
    for (size_t idx = 0; idx < x.size(); ++idx)
    {
      if (idy <= idx)
        jac.setValue(idx, idy, dsquare_ref_scalar(x[idy]));
    }
  }
  return jac;
}

int main()
{
  // Size and variable declarations
  constexpr size_t    n = 10;
  std::vector<double> var(n);

  // Random input values
  srand(static_cast<unsigned int>(time(NULL)));
  for (size_t idx = 0; idx < var.size(); ++idx)
  {
    var[idx] = rand();
  }

  // Model
  VectorModel* vector_model = new VectorModel(n);
  vector_model->setVariable(var);
  vector_model->evalResidual();
  vector_model->evalJacobian();
  std::vector<double> var_temp = vector_model->getVariable();
  std::vector<double> res      = vector_model->getResidual();
  DenseMatrix         jac      = vector_model->getJacobian();

  // Reference Jacobian
  DenseMatrix jac_ref = dsquare_ref(var, res);

  // Check
  int  fail    = 0;
  bool verbose = true;
  for (size_t idy = 0; idy < res.size(); ++idy)
  {
    for (size_t idx = 0; idx < var.size(); ++idx)
    {
      if (!GridKit::Testing::isEqual(jac.getValue(idx, idy), jac_ref.getValue(idx, idy)))
      {
        fail++;
        if (verbose)
        {
          std::cout << "Result incorrect at line = " << idy << ", column = " << idx << "\n";
          std::cout << "x = " << var_temp[idx] << ", x^2 = " << res[idx] << ", d(x^2)/dx = " << jac.getValue(idx, idy) << "\n";
        }
      }
    }
  }
  if (verbose)
  {
    jac.printMatrix("Autodiff Jacobian");
    jac_ref.printMatrix("Reference Jacobian");
  }
  std::cout << "Status: " << fail << "\n";

  // Cleanup
  delete vector_model;

  return fail;
}
