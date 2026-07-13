
#define _USE_MATH_DEFINES
#include "Generator4Param.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>

namespace GridKit
{

  /*!
   * @brief Constructor for a simple generator model
   *
   * Arguments passed to ModelEvaluatorImpl:
   * - Number of equations = 4 differential + 2 algebraic = 6
   * - Number of quadratures = 1
   * - Number of optimization parameters = 1
   */
  template <class ScalarT, typename IdxT>
  Generator4Param<ScalarT, IdxT>::Generator4Param(bus_type* bus, ScalarT P0, ScalarT Q0)
    : ModelEvaluatorImpl<ScalarT, IdxT>(6, 1, 1),
      H_(5.0),
      D_(0.04),
      Xq_(0.85),
      Xd_(1.05),
      Xqp_(0.35),
      Xdp_(0.35),
      Rs_(0.01),
      Tq0p_(1.0), // [s]
      Td0p_(8.0), // [s]
      Ef_(1.45),
      Pm_(1.0),
      omega_s_(1.0),
      omega_b_(2.0 * 60.0 * M_PI),
      P0_(P0),
      Q0_(Q0),
      bus_(bus)
  {
  }

  template <class ScalarT, typename IdxT>
  Generator4Param<ScalarT, IdxT>::~Generator4Param()
  {
  }

  /*!
   * @brief This function will be used to allocate sparse Jacobian matrices.
   *
   */
  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::allocate()
  {
    // std::cout << "Allocate Generator4Param..." << std::endl;
    tag_.resize(static_cast<size_t>(size_));
    return 0;
  }

  /**
   * @brief Initialization of the generator model
   *
   * Initialization equations are derived from example 9.2 in Power System
   * Modeling and Scripting, Federico Milano, Chapter 9, p. 225:
   * \f{eqnarray*}{
   * &~& \omega_0 = 0, \\
   * &~& \delta_0 = \tan^{-1} \left(\frac{X_q P_0 - R_s Q_0}{V_0^2 + R_s P_0 + X_q Q_0} \right) + \theta_0, \\
   * &~& \phi_0   = \delta_0 - \theta_0 + \tan^{-1} \left( \frac{Q_0}{P_0} \right), \\
   * &~& I_{d0}   = \frac{\sqrt{P_0^2 + Q_0^2}}{V_0} \sin(\phi_0), \\
   * &~& I_{q0}   = \frac{\sqrt{P_0^2 + Q_0^2}}{V_0} \cos(\phi_0), \\
   * &~& E_{d0}'  = V_0 \sin(\delta_0 - \theta_0) + R_s I_{d0} - X_q' I_{q0}, \\
   * &~& E_{q0}'  = V_0 \cos(\delta_0 - \theta_0) + R_s I_{q0} + X_d' I_{d0}
   * \f}
   *
   * The input from exciter and governor is set to the steady state value:
   * \f{eqnarray*}{
   * &~& E_{f0} = E_{q0}' + (X_d - X_d') I_{d0}, \\
   * &~& P_{m0} = E_{d0}' I_{d0} + E_{q0}' I_{q0} + ( X_q' - X_d') I_{d0} I_{q0}
   * \f}
   *
   */
  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::initialize()
  {
    // std::cout << "Initialize Generator4Param..." << std::endl;

    // Compute initial guess for the generator voltage phase
    const ScalarT delta = atan((Xq_ * P0_ - Rs_ * Q0_) / (V() * V() + Rs_ * P0_ + Xq_ * Q0_)) + theta();

    // Compute initial guess for the generator current phase
    const ScalarT phi = theta() - delta - atan(Q0_ / P0_);

    // Compute initial gueses for generator currents and potentials in d-q frame
    const ScalarT Id = std::sqrt(P0_ * P0_ + Q0_ * Q0_) / V() * sin(phi);
    const ScalarT Iq = std::sqrt(P0_ * P0_ + Q0_ * Q0_) / V() * cos(phi);
    const ScalarT Ed = V() * sin(theta() - delta) + Rs_ * Id + Xqp_ * Iq;
    const ScalarT Eq = V() * cos(theta() - delta) + Rs_ * Iq - Xdp_ * Id;

    auto* y  = y_.getData();
    auto* yp = yp_.getData();

    y[0]  = delta;
    y[1]  = omega_s_;
    y[2]  = Ed;
    y[3]  = Eq;
    y[4]  = Id;
    y[5]  = Iq;
    yp[0] = 0.0;
    yp[1] = 0.0;
    yp[2] = 0.0;
    yp[3] = 0.0;
    yp[4] = 0.0;
    yp[5] = 0.0;

    // Set control parameter values here.
    Ef_ = Eq - (Xd_ - Xdp_) * Id;                      // <~ set to steady state value
    Pm_ = Ed * Id + Eq * Iq + (Xdp_ - Xqp_) * Id * Iq; // <~ set to steady state value

    // Initialize optimization parameters
    auto* param    = param_.getData();
    auto* param_up = param_up_.getData();
    auto* param_lo = param_lo_.getData();

    param[0]    = H_;
    param_up[0] = 10.0;
    param_lo[0] = 2.0;

    return 0;
  }

  /**
   * \brief Identify differential variables.
   */
  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::tagDifferentiable()
  {
    tag_[0] = true;
    tag_[1] = true;
    tag_[2] = true;
    tag_[3] = true;

    for (IdxT i = 4; i < size_; ++i)
    {
      tag_[static_cast<size_t>(i)] = false;
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    std::fill(abs_tol_.getData(), abs_tol_.getData() + abs_tol_.getSize(), rel_tol);
    return 0;
  }

  /**
   * @brief Computes residual vector for the generator model.
   *
   * Residual equations are given per model in Power System Modeling and
   * Scripting, Federico Milano, Chapter 15, p. 334:
   * \f{eqnarray*}{
   * f_0: &~& \dot{\delta} -\omega_b (\omega - \omega_s), \\
   * f_1: &~& 2H/\omega_s \dot{\omega} - L_m(P_m) + E_q' I_q + E_d' I_d + (X_q' - X_d')I_d I_q  + D (\omega - \omega_s), \\
   * f_2: &~& T_{q0}' \dot{E}_d' + E_d' - (X_q - X_q')I_q, \\
   * f_3: &~& T_{d0}' \dot{E}_q' + E_q' + (X_d - X_d')I_d - E_f, \\
   * f_4: &~& R_s I_d - X_q' I_q + V \sin(\delta - \theta) - E_d', \\
   * f_5: &~& R_s I_q + X_d' I_d + V \cos(\delta - \theta) - E_q',
   * \f}
   * where \f$ \Omega_b \f$ is the synchronous frequency in [rad/s], and
   * overdot denotes time derivative.
   *
   * Generator injection active and reactive power are
   * \f{eqnarray*}{
   * P_g &=& E_d' I_d + E_q' I_q + (X_q' - X_d') I_d I_q - R_s (I_d^2 + I_q^2), \\
   * Q_q &=& E_q' I_d - E_d' I_q - X_q' I_q^2 - X_d' I_d^2, \\
   * \f}
   * respectively.
   *
   * State variables are:
   * \f$ y_0 = \omega \f$, \f$ y_1 = \delta \f$, \f$ y_2 = E_d' \f$, \f$ y_3 = E_q' \f$,
   * \f$ y_4 = I_d \f$, \f$ y_5 = I_q \f$.
   *
   */
  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::evaluateResidual()
  {
    // std::cout << "Evaluate residual for Generator4Param..." << std::endl;
    auto* f = f_.getData();

    f[0] = dotDelta() - omega_b_ * (omega() - omega_s_);
    f[1] = (2.0 * H()) / omega_s_ * dotOmega() - Pm() + Eqp() * Iq() + Edp() * Id() + (-Xdp_ + Xqp_) * Id() * Iq() + D_ * (omega() - omega_s_);
    f[2] = Tq0p_ * dotEdp() + Edp() - (Xq_ - Xqp_) * Iq();
    f[3] = Td0p_ * dotEqp() + Eqp() + (Xd_ - Xdp_) * Id() - Ef();
    f[4] = Rs_ * Id() - Xqp_ * Iq() + V() * sin(delta() - theta()) - Edp();
    f[5] = Xdp_ * Id() + Rs_ * Iq() + V() * cos(delta() - theta()) - Eqp();

    // Compute active and reactive load provided by the infinite bus.
    P() += Pg();
    Q() += Qg();

    // std::cout << "Residual: t = " << time_ << std::endl;

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::evaluateJacobian()
  {
    std::cerr << "Evaluate Jacobian for Generator4Param..." << std::endl;
    std::cerr << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::evaluateIntegrand()
  {
    // std::cout << "Evaluate Integrand for Generator4Param..." << std::endl;
    auto* g = g_.getData();

    g[0] = trajectoryPenalty(time_);
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::initializeAdjoint()
  {
    // std::cout << "Initialize adjoint for Generator4Param..." << std::endl;
    auto* yB  = yB_.getData();
    auto* ypB = ypB_.getData();

    for (IdxT i = 0; i < size_; ++i)
    {
      yB[static_cast<size_t>(i)]  = 0.0;
      ypB[static_cast<size_t>(i)] = 0.0;
    }
    ypB[2] = -trajectoryPenaltyDerEdp(time_) / Tq0p_;
    ypB[3] = -trajectoryPenaltyDerEqp(time_) / Td0p_;

    return 0;
  }

  /**
   * @brief Computes adjoint residual vector for the generator model.
   *
   * Adjoint residual equations are given as:
   * \f{eqnarray*}{
   * f_{B0}: &~& \dot{y}_{B0} - y_{B4} V \cos(\delta - \theta) + y_{B5} V \sin(\delta - \theta), \\
   * f_{B1}: &~& 2H/\omega_s \dot{y}_{B1} + y_{B0} \omega_b - y_{B1} D + y_{B9} (1 - T_2/T_1) - y_{B10} K T_2/T_1 + g_{\omega}(\omega), \\
   * f_{B2}: &~& T_{q0}' \dot{y}_{B2} - y_{B1} I_d - y_{B2} + y_{B4} + y_{B6} I_d - y_{B7} I_q, \\
   * f_{B3}: &~& T_{d0}' \dot{y}_{B3} - y_{B1} I_q - y_{B3} + y_{B5} + y_{B6} I_q + y_{B7} I_d, \\
   * f_{B4}: &~& -y_{B1} (E_d' + (-X_d'+X_q') I_q) - y_{B3} (X_d - X_d') - y_{B4} R_s - y_{B5} X_d' + y_{B6} (E_d' + (X_q' - X_d') I_q - 2 R_s I_d) + y_{B7} (E_q' - 2 X_d' I_d), \\
   * f_{B5}: &~& -y_{B1} (E_q' + (-X_d'+X_q') I_d) + y_{B2} (X_q - X_q') + y_{B4} X_q' - y_{B5} R_s + y_{B6} (E_q' + (X_q' - X_d') I_d - 2 R_s I_q) - y_{B7} (E_d' + 2 X_q' I_q). \\
   * \f}
   *
   */
  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    // std::cout << "Evaluate adjoint residual for Generator4Param..." << std::endl;
    ScalarT sinPhi = sin(delta() - theta());
    ScalarT cosPhi = cos(delta() - theta());

    auto* yB  = yB_.getData();
    auto* ypB = ypB_.getData();
    auto* fB  = fB_.getData();

    // Generator adjoint
    fB[0] = ypB[0] - yB[4] * V() * cosPhi + yB[5] * V() * sinPhi;
    fB[1] = 2.0 * H() / omega_s_ * ypB[1] + yB[0] * omega_b_ - yB[1] * D_; //+ frequencyPenaltyDer(omega());
    fB[2] = Tq0p_ * ypB[2] - yB[1] * Id() - yB[2] + yB[4] + trajectoryPenaltyDerEdp(time_);
    fB[3] = Td0p_ * ypB[3] - yB[1] * Iq() - yB[3] + yB[5] + trajectoryPenaltyDerEqp(time_);
    fB[4] = -yB[1] * (Edp() + (Xqp_ - Xdp_) * Iq()) - yB[3] * (Xd_ - Xdp_) - yB[4] * Rs_ - yB[5] * Xdp_;
    fB[5] = -yB[1] * (Eqp() + (Xqp_ - Xdp_) * Id()) + yB[2] * (Xq_ - Xqp_) + yB[4] * Xqp_ - yB[5] * Rs_;

    return 0;
  }

  // template <class ScalarT, typename IdxT>
  // int Generator4Param<ScalarT, IdxT>::evaluateAdjointJacobian()
  // {
  //     std::cout << "Evaluate adjoint Jacobian for Generator4Param..." << std::endl;
  //     std::cout << "Adjoint Jacobian evaluation not implemented!" << std::endl;
  //     return 0;
  // }

  template <class ScalarT, typename IdxT>
  int Generator4Param<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    // std::cout << "Evaluate adjoint Integrand for Generator4Param..." << std::endl;
    auto* yB = yB_.getData();
    auto* gB = gB_.getData();

    gB[0] = -2.0 * yB[1] * dotOmega() / omega_s_;

    return 0;
  }

  //
  // Private functions
  //

  /**
   * Generator active power Pg.
   *
   * \f[ P_g = E_q' I_q + E_d' I_d + (X_q' - X_d') I_q I_d - R_a (I_d^2 + I_q^2) \f]
   *
   */
  template <class ScalarT, typename IdxT>
  ScalarT Generator4Param<ScalarT, IdxT>::Pg()
  {
    auto* y = y_.getData();
    return y[5] * V() * cos(theta() - y[0]) + y[4] * V() * sin(theta() - y[0]);
  }

  /**
   * Generator reactive power Qg.
   *
   * \f[ Q_g = E_q' I_d - E_d' I_q - X_d' I_d^2 - X_q' I_q^2 \f]
   */
  template <class ScalarT, typename IdxT>
  ScalarT Generator4Param<ScalarT, IdxT>::Qg()
  {
    auto* y = y_.getData();
    return y[5] * V() * sin(theta() - y[0]) - y[4] * V() * cos(theta() - y[0]);
  }

  /**
   * @brief Difference between computed system state and look-up table value.
   *
   * @todo Look-up table should probably live outside the generator model.
   */
  template <class ScalarT, typename IdxT>
  ScalarT Generator4Param<ScalarT, IdxT>::trajectoryPenalty(ScalarT t) const
  {
    size_t N  = table_.size();
    double ti = table_[0][0];
    double tf = table_[N - 1][0];
    double dt = (tf - ti) / static_cast<double>(N - 1);
    int    n  = static_cast<int>(std::trunc(t / tf * static_cast<double>(N - 1)));

    double Edp_est = 0.0;
    double Eqp_est = 0.0;

    if (t >= ti && t < tf)
    {
      // Interpolate from look-up table
      Edp_est = (table_[static_cast<size_t>(n + 1)][3] - table_[static_cast<size_t>(n)][3]) / (table_[static_cast<size_t>(n + 1)][0] - table_[static_cast<size_t>(n)][0]) * (t - table_[static_cast<size_t>(n)][0]) + table_[static_cast<size_t>(n)][3];
      Eqp_est = (table_[static_cast<size_t>(n + 1)][4] - table_[static_cast<size_t>(n)][4]) / (table_[static_cast<size_t>(n + 1)][0] - table_[static_cast<size_t>(n)][0]) * (t - table_[static_cast<size_t>(n)][0]) + table_[static_cast<size_t>(n)][4];
    }
    else
    {
      if (tf <= t && t < tf + dt)
      {
        // Extrapolate from look-up table
        Edp_est = (table_[static_cast<size_t>(n)][3] - table_[static_cast<size_t>(n - 1)][3]) / (table_[static_cast<size_t>(n)][0] - table_[static_cast<size_t>(n - 1)][0]) * (t - table_[static_cast<size_t>(n - 1)][0]) + table_[static_cast<size_t>(n - 1)][3];
        Eqp_est = (table_[static_cast<size_t>(n)][4] - table_[static_cast<size_t>(n - 1)][4]) / (table_[static_cast<size_t>(n)][0] - table_[static_cast<size_t>(n - 1)][0]) * (t - table_[static_cast<size_t>(n - 1)][0]) + table_[static_cast<size_t>(n - 1)][4];
      }
      else
      {
        // Too far away to extrapolate
        std::cerr << "Trajectory penalty: Out of time bounds at time " << t << "\n";
        return -1.0;
      }
    }
    double d = (Edp() - Edp_est);
    double q = (Eqp() - Eqp_est);
    return (d * d + q * q);
  }

  template <class ScalarT, typename IdxT>
  ScalarT Generator4Param<ScalarT, IdxT>::trajectoryPenaltyDerEdp(ScalarT t) const
  {
    size_t N       = table_.size();
    double ti      = table_[0][0];
    double tf      = table_[N - 1][0];
    double dt      = (tf - ti) / static_cast<double>(N - 1);
    int    n       = static_cast<int>(std::trunc(t / tf * static_cast<double>(N - 1)));
    double Edp_est = 0.0;

    if (t >= ti && t < tf)
    {
      Edp_est = (table_[static_cast<size_t>(n + 1)][3] - table_[static_cast<size_t>(n)][3]) / (table_[static_cast<size_t>(n + 1)][0] - table_[static_cast<size_t>(n)][0]) * (t - table_[static_cast<size_t>(n)][0]) + table_[static_cast<size_t>(n)][3];
    }
    else
    {
      if (tf <= t && t < tf + dt)
      {
        Edp_est = (table_[static_cast<size_t>(n)][3] - table_[static_cast<size_t>(n - 1)][3]) / (table_[static_cast<size_t>(n)][0] - table_[static_cast<size_t>(n - 1)][0]) * (t - table_[static_cast<size_t>(n - 1)][0]) + table_[static_cast<size_t>(n - 1)][3];
      }
      else
      {
        std::cerr << "Trajectory penalty: Out of time bounds at time " << t << "\n";
        return -1.0;
      }
    }
    double d = (Edp() - Edp_est);

    return 2.0 * d;
  }

  template <class ScalarT, typename IdxT>
  ScalarT Generator4Param<ScalarT, IdxT>::trajectoryPenaltyDerEqp(ScalarT t) const
  {
    size_t N       = table_.size();
    double ti      = table_[0][0];
    double tf      = table_[N - 1][0];
    double dt      = (tf - ti) / static_cast<double>(N - 1);
    int    n       = static_cast<int>(std::trunc(t / tf * static_cast<double>(N - 1)));
    double Eqp_est = 0.0;

    if (t >= ti && t < tf)
    {
      Eqp_est = (table_[static_cast<size_t>(n + 1)][4] - table_[static_cast<size_t>(n)][4]) / (table_[static_cast<size_t>(n + 1)][0] - table_[static_cast<size_t>(n)][0]) * (t - table_[static_cast<size_t>(n)][0]) + table_[static_cast<size_t>(n)][4];
    }
    else
    {
      if (tf <= t && t < tf + dt)
      {
        Eqp_est = (table_[static_cast<size_t>(n)][4] - table_[static_cast<size_t>(n - 1)][4]) / (table_[static_cast<size_t>(n)][0] - table_[static_cast<size_t>(n - 1)][0]) * (t - table_[static_cast<size_t>(n - 1)][0]) + table_[static_cast<size_t>(n - 1)][4];
      }
      else
      {
        std::cerr << "Trajectory penalty: Out of time bounds at time " << t << "\n";
        return -1.0;
      }
    }
    double q = (Eqp() - Eqp_est);

    return 2.0 * q;
  }

  template class Generator4Param<double, long int>;
  template class Generator4Param<double, size_t>;

} // namespace GridKit
