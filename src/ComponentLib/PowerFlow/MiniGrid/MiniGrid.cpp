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
#include "MiniGrid.hpp"
#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>

namespace ModelLib {

/*!
 * @brief Constructor for a constant load model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
MiniGrid<ScalarT, IdxT>::MiniGrid()
  : ModelEvaluatorImpl<ScalarT, IdxT>(3, 0, 0),
    Pl2_(  2.5),
    Ql2_( -0.8),
    Pg3_(  2.0),
    V1_ (  1.0),
    th1_(  0.0),
    V3_ (  1.1),
    B12_( 10.0),
    B13_( 15.0),
    B22_(-22.0),
    B23_( 12.0)
{
    //std::cout << "Create a load model with " << size_ << " variables ...\n";
    rtol_ = 1e-5;
    atol_ = 1e-5;
}

template <class ScalarT, typename IdxT>
MiniGrid<ScalarT, IdxT>::~MiniGrid()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int MiniGrid<ScalarT, IdxT>::allocate()
{
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int MiniGrid<ScalarT, IdxT>::initialize()
{
    th2() = 0.0; // th2
    V2()  = 1.0; // V2
    th3() = 0.0; // th3
    return 0;
}


/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int MiniGrid<ScalarT, IdxT>::evaluateResidual()
{
    f_[0] = -Pl2_ - V2()*(V1_*B12_*sin(th2()-th1_) + V3_*B23_*sin(th2() - th3()));
    f_[1] = -Ql2_ + V2()*(V1_*B12_*cos(th2()-th1_) + B22_*V2() + V3_*B23_*cos(th2() - th3()));
    f_[2] =  Pg3_ - V3_ *(V1_*B13_*sin(th3()-th1_) + V2()*B23_*sin(th3() - th2()));

    return 0;
}

template <class ScalarT, typename IdxT>
int MiniGrid<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}

// Available template instantiations
template class MiniGrid<double, long int>;
template class MiniGrid<double, size_t>;


} //namespace ModelLib

