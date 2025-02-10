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

#include <Model/PhasorDynamics/BusBase.hpp>


// Forward declaration of BusData structure
namespace GridKit
{
namespace PowerSystemData
{
    template <typename RealT, typename IdxT>
    struct BusData;
}
}

namespace GridKit
{
namespace PhasorDynamics
{
    /*!
     * @brief Implementation of a PQ bus.
     *
     * Voltage _V_ and phase _theta_ are variables in PQ bus model.
     * Active and reactive power, _P_ and _Q_, are residual components.
     *
     *
     */
    template  <class ScalarT, typename IdxT>
    class Bus : public BusBase<ScalarT, IdxT>
    {
        using BusBase<ScalarT, IdxT>::size_;
        using BusBase<ScalarT, IdxT>::y_;
        using BusBase<ScalarT, IdxT>::yp_;
        using BusBase<ScalarT, IdxT>::yB_;
        using BusBase<ScalarT, IdxT>::ypB_;
        using BusBase<ScalarT, IdxT>::f_;
        using BusBase<ScalarT, IdxT>::fB_;
        using BusBase<ScalarT, IdxT>::tag_;

    public:
        using real_type = typename BusBase<ScalarT, IdxT>::real_type;
        using BusData   = GridKit::PowerSystemData::BusData<real_type, IdxT>;

        Bus();
        Bus(ScalarT Vr, ScalarT Vi);
        Bus(BusData& data);
        virtual ~Bus();

        virtual int allocate();
        virtual int tagDifferentiable();
        virtual int initialize();
        virtual int evaluateResidual();
        virtual int initializeAdjoint();
        virtual int evaluateAdjointResidual();

        virtual ScalarT& Vr()
        {
            return y_[0];
        }

        virtual const ScalarT& Vr() const
        {
            return y_[0];
        }

        virtual ScalarT& Vi()
        {
            return y_[1];
        }

        virtual const ScalarT& Vi() const
        {
            return y_[1];
        }

        virtual ScalarT& Ir()
        {
            return f_[0];
        }

        virtual const ScalarT& Ir() const
        {
            return f_[0];
        }

        virtual ScalarT& Ii()
        {
            return f_[1];
        }

        virtual const ScalarT& Ii() const
        {
            return f_[1];
        }

        virtual ScalarT& VrB()
        {
            return yB_[0];
        }

        virtual const ScalarT& VrB() const
        {
            return yB_[0];
        }

        virtual ScalarT& ViB()
        {
            return yB_[1];
        }

        virtual const ScalarT& ViB() const
        {
            return yB_[1];
        }

        virtual ScalarT& IrB()
        {
            return fB_[0];
        }

        virtual const ScalarT& IrB() const
        {
            return fB_[0];
        }

        virtual ScalarT& IiB()
        {
            return fB_[1];
        }

        virtual const ScalarT& IiB() const
        {
            return fB_[1];
        }

        // virtual const int BusType() const
        // {
        //     return BaseBus<ScalarT, IdxT>::BusType::PQ;
        // }

        // virtual const IdxT BusID() const
        // {
        //     return busID_;
        // }

    private:
        // Default initial values for voltage and phase on PQ bus
        // const IdxT busID_{static_cast<IdxT>(-1)};
        ScalarT Vr0_{0.0};
        ScalarT Vi0_{0.0};

    };

} // PhasorDynamics
} // namespace GridKit
