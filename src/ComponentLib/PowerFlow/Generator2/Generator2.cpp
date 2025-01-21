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
#include <ComponentLib/PowerFlow/Bus/BusSlack.hpp>
#include "Generator2.hpp"

namespace ModelLib {

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
    omega_b_(2.0*60.0*M_PI),
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
    tag_.resize(size_);
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
int Generator2<ScalarT, IdxT>::initialize()
{
    // Set optimization parameter value and bounds
    param_[0] = Pm_;
    param_up_[0] = 1.5;
    param_lo_[0] = 0.5;

    y_[0] = asin((Pm_*Xdp_)/(Eqp_*V())) + theta(); // <~ asin(Pm/Pmax)
    y_[1] = omega_s_;
    yp_[0] = 0.0;
    yp_[1] = 0.0;

    return 0;
}

template <class ScalarT, typename IdxT>
int Generator2<ScalarT, IdxT>::evaluateResidual()
{
    f_[0] = -yp_[0] + omega_b_*(y_[1]-omega_s_);
    f_[1] = -yp_[1] + omega_s_/(2.0*H_)*( param_[0] - Eqp_/Xdp_*V()*sin(y_[0] - theta()) - D_*(y_[1]-omega_s_) );
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
    g_[0] = frequencyPenalty(y_[1]);
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator2<ScalarT, IdxT>::initializeAdjoint()
{
    yB_[0] = 0.0;
    yB_[1] = 0.0;
    ypB_[0] = 0.0;
    ypB_[1] = frequencyPenaltyDer(y_[1]);

    return 0;
}

template <class ScalarT, typename IdxT>
int Generator2<ScalarT, IdxT>::evaluateAdjointResidual()
{
    fB_[0]  = -ypB_[0] + omega_s_/(2.0*H_)*Eqp_/Xdp_*V()*cos(y_[0] - theta()) * yB_[1];
    fB_[1]  = -ypB_[1] + omega_s_/(2.0*H_)*D_ * yB_[1] - omega_b_*yB_[0] + frequencyPenaltyDer(y_[1]);
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
    gB_[0] = -omega_s_/(2.0*H_) * yB_[1];
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
    return c_ * pow(std::max(0.0, std::max(omega - omega_up_, omega_lo_ - omega)), beta_);
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
template class Generator2<double, long int>;
template class Generator2<double, size_t>;


} // namespace ModelLib