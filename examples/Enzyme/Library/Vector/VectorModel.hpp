#pragma once

#include <LinearAlgebra/DenseMatrix/DenseMatrix.hpp>

/**
 * @brief Class providing methods to evaluate a vector-valued residual and its Jacobian.
 * This is used to test automatic differentiation.
 */
template <class ScalarT, typename IdxT>
class VectorModel
{
private:
  using DenseMatrix = GridKit::LinearAlgebra::DenseMatrix<ScalarT, IdxT>;

  std::vector<ScalarT> x_, f_;
  DenseMatrix          df_dx_;
  inline ScalarT       square_scalar(ScalarT);
  void                 square(std::vector<ScalarT>&, std::vector<ScalarT>&);

public:
  VectorModel(IdxT);
  void                  setVariable(std::vector<ScalarT>);
  void                  evaluateResidual();
  void                  evaluateJacobian();
  std::vector<ScalarT>& getVariable();
  std::vector<ScalarT>& getResidual();
  DenseMatrix&          getJacobian();
  ~VectorModel();
};
