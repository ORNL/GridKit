
#include "Branch.hpp"

#include <cmath>
#include <iostream>

#include <Model/PowerFlow/Bus/BaseBus.hpp>
#include <Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{

  /*!
   * @brief Constructor for a pi-model branch
   *
   * Arguments passed to ModelEvaluatorImpl:
   * - Number of equations = 0
   * - Number of independent variables = 0
   * - Number of quadratures = 0
   * - Number of optimization parameters = 0
   */

  template <class ScalarT, typename IdxT>
  Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2)
    : R_(0.0),
      X_(0.01),
      G_(0.0),
      B_(0.0),
      fbusID_(0),
      tbusID_(0),
      bus1_(bus1),
      bus2_(bus2)
  {
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  Branch<ScalarT, IdxT>::Branch(real_type R, real_type X, real_type G, real_type B, bus_type* bus1, bus_type* bus2)
    : R_(R),
      X_(X),
      G_(G),
      B_(B),
      fbusID_(0),
      tbusID_(0),
      bus1_(bus1),
      bus2_(bus2)
  {
  }

  template <class ScalarT, typename IdxT>
  Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, BranchData& data)
    : R_(data.r),
      X_(data.x),
      G_(0.0),
      B_(data.b),
      fbusID_(data.fbus),
      tbusID_(data.tbus),
      bus1_(bus1),
      bus2_(bus2)
  {
    size_ = 0;
  }

  template <class ScalarT, typename IdxT>
  Branch<ScalarT, IdxT>::~Branch()
  {
    // std::cout << "Destroy Branch..." << std::endl;
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::allocate()
  {
    // std::cout << "Allocate Branch..." << std::endl;
    return 0;
  }

  /**
   * Initialization of the branch model
   *
   */
  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /**
   * \brief Identify differential variables.
   */
  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * \brief Residual contribution of the branch is pushed to the
   * two terminal buses.
   *
   * @todo Add and verify conductance to ground (B and G)
   */
  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluating branch residual ...\n";
    real_type b      = -X_ / (R_ * R_ + X_ * X_);
    real_type g      = R_ / (R_ * R_ + X_ * X_);
    ScalarT   dtheta = theta1() - theta2();

    P1() -= (g + 0.5 * G_) * V1() * V1() + V1() * V2() * (-g * cos(dtheta) - b * sin(dtheta));
    Q1() -= (-b - 0.5 * B_) * V1() * V1() + V1() * V2() * (-g * sin(dtheta) + b * cos(dtheta));
    P2() -= (g + 0.5 * G_) * V2() * V2() + V1() * V2() * (-g * cos(dtheta) + b * sin(dtheta));
    Q2() -= (-b - 0.5 * B_) * V2() * V2() + V1() * V2() * (g * sin(dtheta) + b * cos(dtheta));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::evaluateJacobian()
  {
    std::cout << "Evaluate Jacobian for Branch..." << std::endl;
    std::cout << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::evaluateIntegrand()
  {
    // std::cout << "Evaluate Integrand for Branch..." << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::initializeAdjoint()
  {
    // std::cout << "Initialize adjoint for Branch..." << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    // std::cout << "Evaluate adjoint residual for Branch..." << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Branch<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    // std::cout << "Evaluate adjoint Integrand for Branch..." << std::endl;
    return 0;
  }

  // Available template instantiations
  template class Branch<double, long int>;
  template class Branch<double, size_t>;

} // namespace GridKit
