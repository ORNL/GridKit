
#include "Branch.hpp"

#include <cmath>
#include <iostream>

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

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

  template <typename scalar_type, typename index_type>
  Branch<scalar_type, index_type>::Branch(BusT* bus1, BusT* bus2)
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

  template <typename scalar_type, typename index_type>
  Branch<scalar_type, index_type>::Branch(RealT R, RealT X, RealT G, RealT B, BusT* bus1, BusT* bus2)
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

  template <typename scalar_type, typename index_type>
  Branch<scalar_type, index_type>::Branch(BusT* bus1, BusT* bus2, BranchData& data)
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

  template <typename scalar_type, typename index_type>
  Branch<scalar_type, index_type>::~Branch()
  {
    // std::cout << "Destroy Branch..." << std::endl;
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::allocate()
  {
    // std::cout << "Allocate Branch..." << std::endl;
    return 0;
  }

  /**
   * Initialization of the branch model
   *
   */
  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /**
   * \brief Identify differential variables.
   */
  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Compute the absolute tolerance for each variable in the model
   *
   * @param rel_tol The relative tolerance which can be used to pick the
   *        absolute tolerance.
   * @return int 0 if successful, non-zero otherwise.
   *
   * This represents a "noise" level close to zero for which pure relative
   * error cannot be used.
   */
  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::setAbsoluteTolerance(RealT)
  {
    return 0;
  }

  /**
   * \brief Residual contribution of the branch is pushed to the
   * two terminal buses.
   *
   * @todo Add and verify conductance to ground (B and G)
   */
  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::evaluateResidual()
  {
    // std::cout << "Evaluating branch residual ...\n";
    RealT   b      = -X_ / (R_ * R_ + X_ * X_);
    RealT   g      = R_ / (R_ * R_ + X_ * X_);
    ScalarT dtheta = theta1() - theta2();

    P1() -= (g + 0.5 * G_) * V1() * V1() + V1() * V2() * (-g * cos(dtheta) - b * sin(dtheta));
    Q1() -= (-b - 0.5 * B_) * V1() * V1() + V1() * V2() * (-g * sin(dtheta) + b * cos(dtheta));
    P2() -= (g + 0.5 * G_) * V2() * V2() + V1() * V2() * (-g * cos(dtheta) + b * sin(dtheta));
    Q2() -= (-b - 0.5 * B_) * V2() * V2() + V1() * V2() * (g * sin(dtheta) + b * cos(dtheta));

    if (bus1_->size() > 0)
    {
      bus1_->getResidual().setDataUpdated();
    }
    if (bus2_->size() > 0)
    {
      bus2_->getResidual().setDataUpdated();
    }

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::evaluateJacobian()
  {
    std::cout << "Evaluate Jacobian for Branch..." << std::endl;
    std::cout << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::evaluateIntegrand()
  {
    // std::cout << "Evaluate Integrand for Branch..." << std::endl;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::initializeAdjoint()
  {
    // std::cout << "Initialize adjoint for Branch..." << std::endl;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::evaluateAdjointResidual()
  {
    // std::cout << "Evaluate adjoint residual for Branch..." << std::endl;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Branch<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    // std::cout << "Evaluate adjoint Integrand for Branch..." << std::endl;
    return 0;
  }

  // Available template instantiations
  template class Branch<double, long int>;
  template class Branch<double, size_t>;

} // namespace GridKit
