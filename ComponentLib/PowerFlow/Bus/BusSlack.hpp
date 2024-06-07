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

#ifndef _BUS_SLACK_HPP_
#define _BUS_SLACK_HPP_

#include "BaseBus.hpp"
#include <PowerSystemData.hpp>

namespace ModelLib
{
    /*!
     * @brief Implementation of a slack bus.
     *
     * Slack bus sets voltage _V_ and phase _theta_ as constants.
     * Active and reactive power, _P_ and _Q_, are component model outputs,
     * but are computed outside the BusSlack class.
     *
     *
     */
    template  <class ScalarT, typename IdxT>
    class BusSlack : public BaseBus<ScalarT, IdxT>
    {
        using BaseBus<ScalarT, IdxT>::size_;
        using BaseBus<ScalarT, IdxT>::y_;
        using BaseBus<ScalarT, IdxT>::yp_;
        using BaseBus<ScalarT, IdxT>::f_;
        using BaseBus<ScalarT, IdxT>::g_;
        using BaseBus<ScalarT, IdxT>::atol_;
        using BaseBus<ScalarT, IdxT>::rtol_;

    public:
        using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
        using BusData = GridKit::PowerSystemData::BusData<real_type, IdxT>;

        BusSlack();
        BusSlack(ScalarT V, ScalarT theta);
        BusSlack(BusData& data);
        virtual ~BusSlack();
        virtual int evaluateResidual();
        virtual int evaluateAdjointResidual();

        /// @todo Should slack bus allow changing voltage?
        virtual ScalarT& V()
        {
            return V_;
        }

        virtual const ScalarT& V() const
        {
            return V_;
        }

        /// @todo Should slack bus allow changing phase?
        virtual ScalarT& theta()
        {
            return theta_;
        }

        virtual const ScalarT& theta() const
        {
            return theta_;
        }

        virtual ScalarT& P()
        {
            return P_;
        }

        virtual const ScalarT& P() const
        {
            return P_;
        }

        virtual ScalarT& Q()
        {
            return Q_;
        }

        virtual const ScalarT& Q() const
        {
            return Q_;
        }

        /// @todo Should slack bus allow changing voltage?
        virtual ScalarT& lambdaP()
        {
            return thetaB_;
        }

        virtual const ScalarT& lambdaP() const
        {
            return thetaB_;
        }

        /// @todo Should slack bus allow changing phase?
        virtual ScalarT& lambdaQ()
        {
            return VB_;
        }

        virtual const ScalarT& lambdaQ() const
        {
            return VB_;
        }

        virtual ScalarT& PB()
        {
            return PB_;
        }

        virtual const ScalarT& PB() const
        {
            return PB_;
        }

        virtual ScalarT& QB()
        {
            return QB_;
        }

        virtual const ScalarT& QB() const
        {
            return QB_;
        }

        virtual const int BusType() const
        {
            return BaseBus<ScalarT, IdxT>::BusType::Slack;
        }

    private:
        ScalarT V_;
        ScalarT theta_;
        ScalarT P_;
        ScalarT Q_;

        ScalarT VB_;
        ScalarT thetaB_;
        ScalarT PB_;
        ScalarT QB_;

    }; // class BusSlack

} // namespace ModelLib


#endif // _BUS_SLACK_HPP_
