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

#ifndef _BUS_PV_HPP_
#define _BUS_PV_HPP_

#include <cassert>
#include "BaseBus.hpp"
#include <PowerSystemData.hpp>

namespace ModelLib
{
    /*!
     * @brief Implementation of a PV bus.
     *
     * Voltage _V_ and phase _theta_ are variables in PV bus model.
     * Active and reactive power, _P_ and _Q_, are residual components.
     *
     *
     */
    template  <class ScalarT, typename IdxT>
    class BusPV : public BaseBus<ScalarT, IdxT>
    {
        using BaseBus<ScalarT, IdxT>::size_;
        using BaseBus<ScalarT, IdxT>::y_;
        using BaseBus<ScalarT, IdxT>::yp_;
        using BaseBus<ScalarT, IdxT>::yB_;
        using BaseBus<ScalarT, IdxT>::ypB_;
        using BaseBus<ScalarT, IdxT>::f_;
        using BaseBus<ScalarT, IdxT>::fB_;
        using BaseBus<ScalarT, IdxT>::tag_;

    public:
        using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
        using BusData = GridKit::PowerSystemData::BusData<real_type, IdxT>;

        BusPV();
        BusPV(ScalarT V, ScalarT theta0);
        BusPV(BusData& data);
        virtual ~BusPV();

        virtual int allocate();
        virtual int tagDifferentiable();
        virtual int initialize();
        virtual int evaluateResidual();
        virtual int initializeAdjoint();
        virtual int evaluateAdjointResidual();

        virtual ScalarT& V()
        {
            return V_;
        }

        virtual const ScalarT& V() const
        {
            return V_;
        }

        virtual ScalarT& theta()
        {
            return y_[0];
        }

        virtual const ScalarT& theta() const
        {
            return y_[0];
        }

        virtual ScalarT& P()
        {
            return f_[0];
        }

        virtual const ScalarT& P() const
        {
            return f_[0];
        }


        virtual ScalarT& Q()
        {
            return Q_;
        }

        virtual const ScalarT& Q() const
        {
            return Q_;
        }

        virtual ScalarT& lambdaP()
        {
            assert(false);
            return thetaB_;
        }

        virtual const ScalarT& lambdaP() const
        {
            assert(false);
            return thetaB_;
        }

        virtual ScalarT& lambdaQ()
        {
            assert(false);
            return VB_;
        }

        virtual const ScalarT& lambdaQ() const
        {
            assert(false);
            return VB_;
        }

        virtual ScalarT& PB()
        {
            assert(false);
            return PB_;
        }

        virtual const ScalarT& PB() const
        {
            assert(false);
            return PB_;
        }

        virtual ScalarT& QB()
        {
            assert(false);
            return QB_;
        }

        virtual const ScalarT& QB() const
        {
            assert(false);
            return QB_;
        }

        virtual const int BusType() const
        {
            return BaseBus<ScalarT, IdxT>::BusType::PV;
        }

    private:
        ScalarT V_;      ///< Bus voltage magnitude
        ScalarT theta0_; ///< Default initial value for phase
        ScalarT Q_;      ///< Reactive power that generator needs to provide

        ScalarT VB_;
        ScalarT thetaB_;
        ScalarT PB_;
        ScalarT QB_;
    };

} // namespace ModelLib


#endif // _BUS_PV_HPP_
