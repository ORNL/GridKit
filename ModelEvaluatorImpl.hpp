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

#ifndef _MODEL_EVALUATOR_IMPL_HPP_
#define _MODEL_EVALUATOR_IMPL_HPP_

#include <vector>
#include <ModelEvaluator.hpp>

namespace ModelLib
{

    /*!
     * @brief Model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class ModelEvaluatorImpl : public ModelEvaluator<ScalarT, IdxT>
    {
    public:
        typedef typename ModelEvaluator<ScalarT, IdxT>::real_type real_type;

        ModelEvaluatorImpl()
          : size_(0),
            size_quad_(0),
            size_opt_(0)
        {}

        ModelEvaluatorImpl(IdxT size, IdxT size_quad, IdxT size_opt)
          : size_(size),
            size_quad_(size_quad),
            size_opt_(size_opt),
            y_(size_),
            yp_(size_),
            f_(size_),
            g_(size_quad_),
            yB_(size_),
            ypB_(size_),
            fB_(size_),
            gB_(size_opt_),
            J_(COO_Matrix<ScalarT, IdxT>()),
            param_(size_opt_),
            param_up_(size_opt_),
            param_lo_(size_opt_)
        {
        }

        virtual IdxT size()
        {
            return size_;
        }

        virtual IdxT nnz()
        {
            return nnz_;
        }

        virtual bool hasJacobian()
        {
            return false;
        }

        virtual IdxT size_quad()
        {
            return size_quad_;
        }

        virtual IdxT size_opt()
        {
            return size_opt_;
        }

        // virtual void updateTime(real_type t, real_type a)
        // {
        //     time_ = t;
        //     alpha_ = a;
        //     std::cout << "updateTime: t = " << time_ << "\n";
        // }

        virtual void setTolerances(real_type& rtol, real_type& atol) const
        {
            rtol = rtol_;
            atol = atol_;
        }

        virtual void setMaxSteps(IdxT& msa) const
        {
            msa = max_steps_;
        }

        std::vector<ScalarT>& y()
        {
            return y_;
        }

        const std::vector<ScalarT>& y() const
        {
            return y_;
        }

        std::vector<ScalarT>& yp()
        {
            return yp_;
        }

        const std::vector<ScalarT>& yp() const
        {
            return yp_;
        }

        std::vector<bool>& tag()
        {
            return tag_;
        }

        const std::vector<bool>& tag() const
        {
            return tag_;
        }

        std::vector<ScalarT>& yB()
        {
            return yB_;
        }

        const std::vector<ScalarT>& yB() const
        {
            return yB_;
        }

        std::vector<ScalarT>& ypB()
        {
            return ypB_;
        }

        const std::vector<ScalarT>& ypB() const
        {
            return ypB_;
        }

        std::vector<ScalarT>& param()
        {
            return param_;
        }

        const std::vector<ScalarT>& param() const
        {
            return param_;
        }

        std::vector<ScalarT>& param_up()
        {
            return param_up_;
        }

        const std::vector<ScalarT>& param_up() const
        {
            return param_up_;
        }

        std::vector<ScalarT>& param_lo()
        {
            return param_lo_;
        }

        const std::vector<ScalarT>& param_lo() const
        {
            return param_lo_;
        }

        std::vector<ScalarT>& getResidual()
        {
            return f_;
        }

        const std::vector<ScalarT>& getResidual() const
        {
            return f_;
        }

        COO_Matrix<ScalarT, IdxT>& getJacobian()
        {
            return J_;
        }

        const COO_Matrix<ScalarT, IdxT>& getJacobian() const
        {
            return J_;
        }

        std::vector<ScalarT>& getIntegrand()
        {
            return g_;
        }

        const std::vector<ScalarT>& getIntegrand() const
        {
            return g_;
        }

        std::vector<ScalarT>& getAdjointResidual()
        {
            return fB_;
        }

        const std::vector<ScalarT>& getAdjointResidual() const
        {
            return fB_;
        }

        std::vector<ScalarT>& getAdjointIntegrand()
        {
            return gB_;
        }

        const std::vector<ScalarT>& getAdjointIntegrand() const
        {
            return gB_;
        }

        //@todo Fix ID naming
        IdxT getIDcomponent()
        {
            return idc_;
        }



    protected:
        IdxT size_;
        IdxT nnz_;
        IdxT size_quad_;
        IdxT size_opt_;

        std::vector<ScalarT> y_;
        std::vector<ScalarT> yp_;
        std::vector<bool> tag_;
        std::vector<ScalarT> f_;
        std::vector<ScalarT> g_;

        std::vector<ScalarT> yB_;
        std::vector<ScalarT> ypB_;
        std::vector<ScalarT> fB_;
        std::vector<ScalarT> gB_;

        COO_Matrix<ScalarT, IdxT> J_;

        std::vector<ScalarT> param_;
        std::vector<ScalarT> param_up_;
        std::vector<ScalarT> param_lo_;

        real_type time_;
        real_type alpha_;

        real_type rtol_;
        real_type atol_;

        IdxT max_steps_;

        IdxT idc_;

    };


} // namespace ModelLib

#endif // _MODEL_EVALUATOR_IMPL_HPP_