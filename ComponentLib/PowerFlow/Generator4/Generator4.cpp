/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit™. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */


#include <iostream>
#include <cmath>
#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>
#include "Generator4.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a simple generator model
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 4 differential + 2 algebraic = 6
 * - Number of quadratures = 1
 * - Number of optimization parameters = 2
 */
template <class ScalarT, typename IdxT>
Generator4<ScalarT, IdxT>::Generator4(bus_type* bus, ScalarT P0, ScalarT Q0)
  : ModelEvaluatorImpl<ScalarT, IdxT>(6, 1, 2),
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
    omega_b_(2.0*60.0*M_PI),
    omega_up_(omega_s_ + 0.0001),
    omega_lo_(omega_s_ - 0.0001),
    c_(10000.0),
    beta_(2),
    P0_(P0),
    Q0_(Q0),
    bus_(bus)
{
}

template <class ScalarT, typename IdxT>
Generator4<ScalarT, IdxT>::~Generator4()
{
}

/*!
 * @brief This function will be used to allocate sparse Jacobian matrices.
 *
 */
template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::allocate()
{
    //std::cout << "Allocate Generator4..." << std::endl;
    tag_.resize(size_);

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
int Generator4<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize Generator4..." << std::endl;

    // Compute initial guess for the generator voltage phase
    const ScalarT delta = atan((Xq_*P0_ - Rs_*Q0_) / (V()*V() + Rs_*P0_ + Xq_*Q0_)) + theta();

    // Compute initial guess for the generator current phase
    const ScalarT phi   = theta() - delta - atan(Q0_/P0_);

    // Compute initial gueses for generator currents and potentials in d-q frame
    const ScalarT Id = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * sin(phi);
    const ScalarT Iq = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * cos(phi);
    const ScalarT Ed = V()*sin(theta() - delta) + Rs_*Id + Xqp_*Iq;
    const ScalarT Eq = V()*cos(theta() - delta) + Rs_*Iq - Xdp_*Id;

    y_[0] =  delta;
    y_[1] =  omega_s_;
    y_[2] =  Ed;
    y_[3] =  Eq;
    y_[4] =  Id;
    y_[5] =  Iq;
    yp_[0] = 0.0;
    yp_[1] = 0.0;
    yp_[2] = 0.0;
    yp_[3] = 0.0;
    yp_[4] = 0.0;
    yp_[5] = 0.0;

    // Set control parameter values here.
    Ef_ = Eq - (Xd_ - Xdp_)*Id;                // <~ set to steady state value
    Pm_ = Ed*Id + Eq*Iq + (Xdp_ - Xqp_)*Id*Iq; // <~ set to steady state value

    // Initialize optimization parameters
    param_[0] = Pm_;
    param_up_[0] = 1.5;
    param_lo_[0] = 0.0;

    param_[1] = Ef_;
    param_up_[1] = 1.7;
    param_lo_[1] = 0.0;

    return 0;
}

/**
 * \brief Identify differential variables.
 */
template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::tagDifferentiable()
{
    tag_[0] = true;
    tag_[1] = true;
    tag_[2] = true;
    tag_[3] = true;

    for (IdxT i=4; i < size_; ++i)
    {
        tag_[i] = false;
    }

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
int Generator4<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluate residual for Generator4..." << std::endl;
    f_[0] = dotDelta() -  omega_b_* (omega() - omega_s_);
    f_[1] = (2.0*H_)/omega_s_*dotOmega() - Pm() + Eqp()*Iq() + Edp()*Id() + (- Xdp_ + Xqp_)*Id()*Iq() + D_*(omega() - omega_s_);
    f_[2] = Tq0p_*dotEdp() + Edp() - (Xq_ - Xqp_)*Iq();
    f_[3] = Td0p_*dotEqp() + Eqp() + (Xd_ - Xdp_)*Id() - Ef();
    f_[4] =  Rs_*Id() - Xqp_*Iq() + V()*sin(delta() - theta()) - Edp();
    f_[5] =  Xdp_*Id() + Rs_*Iq() + V()*cos(delta() - theta()) - Eqp();

    // Compute active and reactive load provided by the infinite bus.
    P() += Pg();
    Q() += Qg();

    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateJacobian()
{
    std::cerr << "Evaluate Jacobian for Generator4..." << std::endl;
    std::cerr << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateIntegrand()
{
    // std::cout << "Evaluate Integrand for Generator4..." << std::endl;
    g_[0] = frequencyPenalty(y_[1]);
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::initializeAdjoint()
{
    //std::cout << "Initialize adjoint for Generator4..." << std::endl;
    for (IdxT i=0; i<size_; ++i)
    {
        yB_[i] = 0.0;
        ypB_[i] = 0.0;
    }
    ypB_[1] = frequencyPenaltyDer(y_[1]);

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
int Generator4<ScalarT, IdxT>::evaluateAdjointResidual()
{
    // std::cout << "Evaluate adjoint residual for Generator4..." << std::endl;
    ScalarT sinPhi = sin(delta() - theta());
    ScalarT cosPhi = cos(delta() - theta());

    // Generator adjoint
    fB_[0] = ypB_[0] - yB_[4]*V()*cosPhi + yB_[5]*V()*sinPhi;
    fB_[1] = 2.0*H_/omega_s_*ypB_[1] + yB_[0]*omega_b_ - yB_[1]*D_ + frequencyPenaltyDer(omega());
    fB_[2] = Tq0p_*ypB_[2] - yB_[1]*Id() - yB_[2] + yB_[4];
    fB_[3] = Td0p_*ypB_[3] - yB_[1]*Iq() - yB_[3] + yB_[5];
    fB_[4] = -yB_[1]*(Edp() + (Xqp_ - Xdp_)*Iq()) - yB_[3]*(Xd_ - Xdp_) - yB_[4]*Rs_ - yB_[5]*Xdp_;
    fB_[5] = -yB_[1]*(Eqp() + (Xqp_ - Xdp_)*Id()) + yB_[2]*(Xq_ - Xqp_) + yB_[4]*Xqp_ - yB_[5]*Rs_;

    return 0;
}

// template <class ScalarT, typename IdxT>
// int Generator4<ScalarT, IdxT>::evaluateAdjointJacobian()
// {
//     std::cout << "Evaluate adjoint Jacobian for Generator4..." << std::endl;
//     std::cout << "Adjoint Jacobian evaluation not implemented!" << std::endl;
//     return 0;
// }

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    // std::cout << "Evaluate adjoint Integrand for Generator4..." << std::endl;
    gB_[0] = yB_[1];
    gB_[1] = yB_[3];

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
ScalarT Generator4<ScalarT, IdxT>::Pg()
{
    return y_[5]*V()*cos(theta() - y_[0]) + y_[4]*V()*sin(theta() - y_[0]);
}

/**
 * Generator reactive power Qg.
 *
 * \f[ Q_g = E_q' I_d - E_d' I_q - X_d' I_d^2 - X_q' I_q^2 \f]
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::Qg()
{
    return y_[5]*V()*sin(theta() - y_[0]) - y_[4]*V()*cos(theta() - y_[0]);
}

/**
 * Frequency penalty is used as the objective function for the generator model.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::frequencyPenalty(ScalarT omega)
{
    return c_ * pow(std::max(0.0, std::max(omega - omega_up_, omega_lo_ - omega)), beta_);
}

/**
 * Derivative of frequency penalty cannot be written in terms of min/max functions.
 * Need to expand conditional statements instead.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::frequencyPenaltyDer(ScalarT omega)
{
    if (omega > omega_up_)
    {
        return beta_ * c_ * pow(omega - omega_up_, beta_ - 1.0);
    }
    else if (omega < omega_lo_)
    {
        return beta_ * c_ * pow(omega - omega_lo_, beta_ - 1.0);
    }
    else
    {
        return 0.0;
    }
}



// Available template instantiations
template class Generator4<double, long int>;
template class Generator4<double, size_t>;


} // namespace ModelLib

