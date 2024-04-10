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
#include "BusPV.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a PV bus
 *
 * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
 * - Number of equations = 1 (size_)
 * - Number of variables = 1 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusPV<ScalarT, IdxT>::BusPV()
  : BaseBus<ScalarT, IdxT>(0), V_(0.0), theta0_(0.0)
{
    //std::cout << "Create BusPV..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
}

/*!
 * @brief Constructor for a PV bus
 *
 * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
 * - Number of equations = 1 (size_)
 * - Number of variables = 1 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusPV<ScalarT, IdxT>::BusPV(ScalarT V, ScalarT theta0)
  : BaseBus<ScalarT, IdxT>(0), V_(V), theta0_(theta0)
{
    //std::cout << "Create BusPV..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
}

template <class ScalarT, typename IdxT>
BusPV<ScalarT, IdxT>::BusPV(BusData& data)
  : BaseBus<ScalarT, IdxT>(data.bus_i), V_(data.Vm), theta0_(data.Va)
{
    //std::cout << "Create BusPV ..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 1;
}

template <class ScalarT, typename IdxT>
BusPV<ScalarT, IdxT>::~BusPV()
{
    //std::cout << "Destroy Gen2..." << std::endl;
}

/*!
 * @brief allocate method resizes local solution and residual vectors.
 */
template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::allocate()
{
    //std::cout << "Allocate PV bus ..." << std::endl;
    f_.resize(size_);
    y_.resize(size_);
    yp_.resize(size_);
    tag_.resize(size_);

    fB_.resize(size_);
    yB_.resize(size_);
    ypB_.resize(size_);

    return 0;
}


template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::tagDifferentiable()
{
    tag_[0] = false;
    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize BusPV..." << std::endl;
    theta() = theta0_;
    yp_[0] = 0.0;

    return 0;
}

/*!
 * @brief PV bus does not compute residuals, so here we just reset residual values.
 *
 * @warning This implementation assumes bus residuals are always evaluated
 * _before_ component model residuals.
 *
 */
template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a PV bus ...\n";
    P() = 0.0; // <-- Residual P
    Q() = 0.0; // <-- Output Qg, the reactive power generator needs to supply

    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::initializeAdjoint()
{
    // std::cout << "Initialize BusPV..." << std::endl;
    yB_[0] = 0.0;
    ypB_[0] = 0.0;

    return 0;
}

template <class ScalarT, typename IdxT>
int BusPV<ScalarT, IdxT>::evaluateAdjointResidual()
{
    fB_[0] = 0.0;

    return 0;
}

// Available template instantiations
template class BusPV<double, long int>;
template class BusPV<double, size_t>;


} // namespace ModelLib

