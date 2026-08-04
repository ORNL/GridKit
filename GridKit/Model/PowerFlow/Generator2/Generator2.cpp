
#include "Generator2.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>

#include <GridKit/Model/PowerFlow/Bus/BusSlack.hpp>

namespace GridKit
{

  /*!
   * @brief Constructor for a simple generator model
   *
   * Arguments passed to ModelEvaluatorImpl:
   * - Number of equations = 2
   * - Number of quadratures = 1
   * - Number of optimization parameters = 1
   */
  template <class ScalarT, typename IdxT>
  Generator2<ScalarT, IdxT>::Generator2(bus_type* bus)
    : ModelEvaluatorImpl<ScalarT, IdxT>(2, 1, 1),
      H_(5.0),
      D_(0.005),
      Pm_(0.7),
      Xdp_(0.5),
      Eqp_(0.93),
      omega_s_(1.0),
      omega_b_(2.0 * 60.0 * std::numbers::pi_v<RealT>),
      omega_up_(omega_s_ + 0.0002),
      omega_lo_(omega_s_ - 0.0002),
      theta_s_(1.0),
      c_(10000.0),
      beta_(2),
      bus_(bus)
  {
  }

  template <class ScalarT, typename IdxT>
  Generator2<ScalarT, IdxT>::~Generator2()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::allocate()
  {
    tag_.resize(static_cast<size_t>(size_));
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::tagDifferentiable()
  {
    tag_[0] = true;
    tag_[1] = true;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::initialize()
  {
    auto* y        = y_.getData();
    auto* yp       = yp_.getData();
    auto* param    = param_.getData();
    auto* param_up = param_up_.getData();
    auto* param_lo = param_lo_.getData();

    // Set optimization parameter value and bounds
    param[0]    = Pm_;
    param_up[0] = 1.5;
    param_lo[0] = 0.5;

    y[0]  = asin((Pm_ * Xdp_) / (Eqp_ * V())) + theta(); // <~ asin(Pm/Pmax)
    y[1]  = omega_s_;
    yp[0] = 0.0;
    yp[1] = 0.0;

    y_.setDataUpdated();
    yp_.setDataUpdated();
    param_.setDataUpdated();
    param_up_.setDataUpdated();
    param_lo_.setDataUpdated();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::evaluateResidual()
  {
    const auto* y     = y_.getData();
    const auto* yp    = yp_.getData();
    auto*       f     = f_.getData();
    const auto* param = param_.getData();

    f[0] = -yp[0] + omega_b_ * (y[1] - omega_s_);
    f[1] = -yp[1] + omega_s_ / (2.0 * H_) * (param[0] - Eqp_ / Xdp_ * V() * std::sin(y[0] - theta()) - D_ * (y[1] - omega_s_));
    f_.setDataUpdated();
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::evaluateJacobian()
  {
    std::cout << "Evaluate Jacobian for Gen2..." << std::endl;
    std::cout << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::evaluateIntegrand()
  {
    const auto* y = y_.getData();
    auto*       g = g_.getData();

    g[0] = frequencyPenalty(y[1]);
    g_.setDataUpdated();
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::initializeAdjoint()
  {
    const auto* y   = y_.getData();
    auto*       yB  = yB_.getData();
    auto*       ypB = ypB_.getData();

    yB[0]  = 0.0;
    yB[1]  = 0.0;
    ypB[0] = 0.0;
    ypB[1] = frequencyPenaltyDer(y[1]);

    yB_.setDataUpdated();
    ypB_.setDataUpdated();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    const auto* y   = y_.getData();
    const auto* yB  = yB_.getData();
    const auto* ypB = ypB_.getData();
    auto*       fB  = fB_.getData();

    fB[0] = -ypB[0] + omega_s_ / (2.0 * H_) * Eqp_ / Xdp_ * V() * std::cos(y[0] - theta()) * yB[1];
    fB[1] = -ypB[1] + omega_s_ / (2.0 * H_) * D_ * yB[1] - omega_b_ * yB[0] + frequencyPenaltyDer(y[1]);
    fB_.setDataUpdated();
    return 0;
  }

  // template <class ScalarT, typename IdxT>
  // int Generator2<ScalarT, IdxT>::evaluateAdjointJacobian()
  // {
  //     std::cout << "Evaluate adjoint Jacobian for Gen2..." << std::endl;
  //     std::cout << "Adjoint Jacobian evaluation not implemented!" << std::endl;
  //     return 0;
  // }

  template <class ScalarT, typename IdxT>
  int Generator2<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    // std::cout << "Evaluate adjoint Integrand for Gen2..." << std::endl;
    const auto* yB = yB_.getData();
    auto*       gB = gB_.getData();

    gB[0] = -omega_s_ / (2.0 * H_) * yB[1];
    gB_.setDataUpdated();
    return 0;
  }

  //
  // Private functions
  //

  /**
   * Frequency penalty is used as the objective function for the generator model.
   */
  template <class ScalarT, typename IdxT>
  ScalarT Generator2<ScalarT, IdxT>::frequencyPenalty(ScalarT omega)
  {
    return c_ * std::pow(std::max(0.0, std::max(omega - omega_up_, omega_lo_ - omega)), beta_);
  }

  /**
   * Derivative of frequency penalty cannot be written in terms of min/max functions.
   * Need to expand conditional statements instead.
   */
  template <class ScalarT, typename IdxT>
  ScalarT Generator2<ScalarT, IdxT>::frequencyPenaltyDer(ScalarT omega)
  {
    if (omega > omega_up_)
    {
      return beta_ * c_ * std::pow(omega - omega_up_, beta_ - 1.0);
    }
    else if (omega < omega_lo_)
    {
      return beta_ * c_ * std::pow(omega - omega_lo_, beta_ - 1.0);
    }
    else
    {
      return 0.0;
    }
  }

  // Available template instantiations
  template class Generator2<double, long int>;
  template class Generator2<double, size_t>;

} // namespace GridKit
