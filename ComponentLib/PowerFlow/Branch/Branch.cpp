/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov> and Duan Nan <duan4@llnl.gov>.
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
#include <PowerSystemData.hpp>

#include "Branch.hpp"

namespace ModelLib {

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
    //std::cout << "Destroy Branch..." << std::endl;
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int Branch<ScalarT, IdxT>::allocate()
{
    //std::cout << "Allocate Branch..." << std::endl;
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
    real_type b = -X_/(R_*R_ + X_*X_);
    real_type g =  R_/(R_*R_ + X_*X_);
    ScalarT dtheta = theta1() - theta2();

    P1() -= ( g + 0.5*G_)*V1()*V1() + V1()*V2()*(-g*cos(dtheta) - b*sin(dtheta));
    Q1() -= (-b - 0.5*B_)*V1()*V1() + V1()*V2()*(-g*sin(dtheta) + b*cos(dtheta));
    P2() -= ( g + 0.5*G_)*V2()*V2() + V1()*V2()*(-g*cos(dtheta) + b*sin(dtheta));
    Q2() -= (-b - 0.5*B_)*V2()*V2() + V1()*V2()*( g*sin(dtheta) + b*cos(dtheta));

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
    //std::cout << "Initialize adjoint for Branch..." << std::endl;
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

} //namespace ModelLib
