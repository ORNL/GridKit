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

#ifndef _BASE_BUS_HPP_
#define _BASE_BUS_HPP_

#include <ModelEvaluatorImpl.hpp>

namespace ModelLib
{
    /*!
     * @brief Base class for all power flow buses.
     *
     * Derived bus types:
     *   0 - swing bus (V and theta are constants)
     *   1 - PV bus    (P and V are constants)
     *   2 - PQ bus    (P and Q are constants)
     *
     * @todo Consider static instead of dynamic polymorphism for
     * bus types. Create Bus class that takes template parameter
     * BusType.
     */
    template  <class ScalarT, typename IdxT>
    class BaseBus : public ModelEvaluatorImpl<ScalarT, IdxT>
    {
    protected:
        using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::rtol_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::atol_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

    public:
        typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;

        enum BusType{PQ=1, PV, Slack, Isolated};

        BaseBus(IdxT id) : busID_(id) {}
        virtual ~BaseBus(){}

        // Set defaults for ModelEvaluator methods
        virtual int allocate() { return 0;}
        virtual int initialize() { return 0;}
        virtual int tagDifferentiable() { return 0;}
        virtual int evaluateResidual() { return 0;}
        virtual int evaluateJacobian() { return 0;}
        virtual int evaluateIntegrand() { return 0;}

        virtual int initializeAdjoint() { return 0;}
        virtual int evaluateAdjointResidual() { return 0;}
        //virtual int evaluateAdjointJacobian() { return 0;}
        virtual int evaluateAdjointIntegrand() { return 0;}
        virtual void updateTime(real_type, real_type) {} // <- throw exception here

        // Pure virtual methods specific to Bus types
        virtual ScalarT& V() = 0;
        virtual const ScalarT& V() const = 0;
        virtual ScalarT& theta() = 0;
        virtual const ScalarT& theta() const = 0;
        virtual ScalarT& P() = 0;
        virtual const ScalarT& P() const = 0;
        virtual ScalarT& Q() = 0;
        virtual const ScalarT& Q() const = 0;

        virtual ScalarT& lambdaP() = 0;
        virtual const ScalarT& lambdaP() const = 0;
        virtual ScalarT& lambdaQ() = 0;
        virtual const ScalarT& lambdaQ() const = 0;
        virtual ScalarT& PB() = 0;
        virtual const ScalarT& PB() const = 0;
        virtual ScalarT& QB() = 0;
        virtual const ScalarT& QB() const = 0;

        virtual const int BusType() const = 0;

        virtual const IdxT BusID() const
        {
            return busID_;
        }

    protected:
        const IdxT busID_;
    }; // class BaseBus

} // namespace ModelLib


#endif // _BASE_BUS_HPP_
