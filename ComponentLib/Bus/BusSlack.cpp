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
#include "BusSlack.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a slack bus
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusSlack<ScalarT, IdxT>::BusSlack()
    : BaseBus<ScalarT, IdxT>(0), V_(0.0), theta_(0.0), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
{
    //std::cout << "Create BusSlack..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 0;
}

/*!
 * @brief BusSlack constructor.
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusSlack<ScalarT, IdxT>::BusSlack(ScalarT V, ScalarT theta)
    : BaseBus<ScalarT, IdxT>(0), V_(V), theta_(theta), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
{
    //std::cout << "Create BusSlack..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;
    P() = 0.0;
    Q() = 0.0;
    size_ = 0;
}

template <class ScalarT, typename IdxT>
BusSlack<ScalarT, IdxT>::BusSlack(BusData& data)
    : BaseBus<ScalarT, IdxT>(data.bus_i), V_(data.Vm), theta_(data.Va)
{
    //std::cout << "Create BusSlack..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;
    P() = 0.0;
    Q() = 0.0;
    size_ = 0;
}

template <class ScalarT, typename IdxT>
BusSlack<ScalarT, IdxT>::~BusSlack()
{
}

template <class ScalarT, typename IdxT>
int BusSlack<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a slack bus ...\n";
    P() = 0.0;
    Q() = 0.0;
    return 0;
}

template <class ScalarT, typename IdxT>
int BusSlack<ScalarT, IdxT>::evaluateAdjointResidual()
{
    PB() = 0.0;
    QB() = 0.0;
    return 0;
}


// Available template instantiations
template class BusSlack<double, long int>;
template class BusSlack<double, size_t>;


} // namespace ModelLib

