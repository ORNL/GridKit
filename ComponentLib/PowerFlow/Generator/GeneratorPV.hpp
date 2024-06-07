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

#pragma once

#include <vector>
#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include "GeneratorBase.hpp"

namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}


namespace ModelLib
{
     /*!
     * @brief Implementation of a PV generator.
     *
     */
    template  <class ScalarT, typename IdxT>
    class GeneratorPV : public GeneratorBase<ScalarT, IdxT>
    {
        using GeneratorBase<ScalarT, IdxT>::size_;
        using GeneratorBase<ScalarT, IdxT>::nnz_;
        using GeneratorBase<ScalarT, IdxT>::time_;
        using GeneratorBase<ScalarT, IdxT>::alpha_;
        using GeneratorBase<ScalarT, IdxT>::y_;
        using GeneratorBase<ScalarT, IdxT>::yp_;
        using GeneratorBase<ScalarT, IdxT>::tag_;
        using GeneratorBase<ScalarT, IdxT>::f_;
        using GeneratorBase<ScalarT, IdxT>::g_;
        using GeneratorBase<ScalarT, IdxT>::yB_;
        using GeneratorBase<ScalarT, IdxT>::ypB_;
        using GeneratorBase<ScalarT, IdxT>::fB_;
        using GeneratorBase<ScalarT, IdxT>::gB_;
        using GeneratorBase<ScalarT, IdxT>::param_;

        using bus_type = BaseBus<ScalarT, IdxT>;
        using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
        using GenData = GridKit::PowerSystemData::GenData<real_type, IdxT>;

    public:
        GeneratorPV(bus_type* bus, GenData& data);
        virtual ~GeneratorPV();

        int allocate();
        int initialize();
        int tagDifferentiable();
        int evaluateResidual();
        int evaluateJacobian();
        int evaluateIntegrand();

        int initializeAdjoint();
        int evaluateAdjointResidual();
        //int evaluateAdjointJacobian();
        int evaluateAdjointIntegrand();
        
        virtual ScalarT& P()
        {
            return P_;
        }

        virtual const ScalarT& P() const
        {
            return P_;
        }

        /// @brief  Reactive power excess on PV bus
        /// @return reference to negative PV generator reactive power
        virtual ScalarT& Q()
        {
            return (bus_->Q());
        }

        /// @brief  Reactive power excess on PV bus
        /// @return const reference to negative PV generator reactive power
        virtual const ScalarT& Q() const
        {
            return (bus_->Q());
        }

    private:
        ScalarT P_;
        bus_type* bus_;
    };
}

