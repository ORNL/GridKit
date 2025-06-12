#pragma once

#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>

using DenseMatrix = GridKit::LinearAlgebra::DenseMatrix<double, size_t>;

/**
 * @brief Class providing methods to evaluate a vector-valued residual and its Jacobian.
 * This is used to test automatic differentiation.
 */
class VectorModel
{
private:
  std::vector<double> x_, f_;
  DenseMatrix         df_dx_;
  inline double       square_scalar(double);
  void                square(std::vector<double>&, std::vector<double>&);

public:
  VectorModel(size_t);
  void                 setVariable(std::vector<double>);
  void                 evalResidual();
  void                 evalJacobian();
  std::vector<double>& getVariable();
  std::vector<double>& getResidual();
  DenseMatrix&         getJacobian();
  ~VectorModel();
};
