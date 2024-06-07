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
#include <vector>
#include "Load.hpp"
#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>

namespace ModelLib {

/*!
 * @brief Constructor for a constant load model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Load<ScalarT, IdxT>::Load(bus_type* bus, ScalarT P, ScalarT Q)
  : P_(P),
    Q_(Q),
    busID_(0),
    bus_(bus)
{
    //std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
}

template <class ScalarT, typename IdxT>
Load<ScalarT, IdxT>::Load(bus_type* bus, LoadData& data)
  : P_(data.Pd),
    Q_(data.Qd),
    busID_(data.bus_i),
    bus_(bus)
{
    //std::cout << "Create a load model with " << size_ << " variables ...\n";
    size_ = 0;
}

template <class ScalarT, typename IdxT>
Load<ScalarT, IdxT>::~Load()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::allocate()
{
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating load residual ...\n";
    bus_->P() -= P_;
    bus_->Q() -= Q_;
    return 0;
}

template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Load<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}




// Available template instantiations
template class Load<double, long int>;
template class Load<double, size_t>;


} //namespace ModelLib

