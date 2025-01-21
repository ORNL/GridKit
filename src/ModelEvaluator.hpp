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

#ifndef _MODEL_EVALUATOR_HPP_
#define _MODEL_EVALUATOR_HPP_

#include <vector>
#include <ScalarTraits.hpp>
#include <SparseMatrix/COO_Matrix.hpp>

namespace ModelLib
{

    /*!
     * @brief Abstract class describing a model.
     *
     */
    template <class ScalarT, typename IdxT>
    class ModelEvaluator
    {
    public:
        typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

        ModelEvaluator(){}
        virtual ~ModelEvaluator(){}

        virtual int allocate() = 0;
        virtual int initialize() = 0;
        virtual int tagDifferentiable() = 0;
        virtual int evaluateResidual() = 0;
        virtual int evaluateJacobian() = 0;
        virtual int evaluateIntegrand() = 0;

        virtual int initializeAdjoint() = 0;
        virtual int evaluateAdjointResidual() = 0;
        // virtual int evaluateAdjointJacobian() = 0;
        virtual int evaluateAdjointIntegrand() = 0;

        virtual IdxT size() = 0;
        virtual IdxT nnz() = 0;

        /**
         * @brief Is the Jacobian defined. Used in IDA to determine wether DQ is used or not
         * 
         * @return true 
         * @return false 
         */
        virtual bool hasJacobian() = 0;
        
        virtual IdxT size_quad() = 0;
        virtual IdxT size_opt() = 0;
        virtual void updateTime(real_type t, real_type a) = 0;
        virtual void setTolerances(real_type& rtol, real_type& atol) const = 0;
        virtual void setMaxSteps(IdxT& msa) const = 0;

        virtual std::vector<ScalarT>& y() = 0;
        virtual const std::vector<ScalarT>& y() const = 0;

        virtual std::vector<ScalarT>& yp() = 0;
        virtual const std::vector<ScalarT>& yp() const = 0;

        virtual std::vector<bool>& tag() = 0;
        virtual const std::vector<bool>& tag() const = 0;

        virtual std::vector<ScalarT>& yB() = 0;
        virtual const std::vector<ScalarT>& yB() const = 0;

        virtual std::vector<ScalarT>& ypB() = 0;
        virtual const std::vector<ScalarT>& ypB() const = 0;

        virtual std::vector<ScalarT>& param() = 0;
        virtual const std::vector<ScalarT>& param() const = 0;

        virtual std::vector<ScalarT>& param_up() = 0;
        virtual const std::vector<ScalarT>& param_up() const = 0;

        virtual std::vector<ScalarT>& param_lo() = 0;
        virtual const std::vector<ScalarT>& param_lo() const = 0;

        virtual std::vector<ScalarT>& getResidual() = 0;
        virtual const std::vector<ScalarT>& getResidual() const = 0;


        virtual COO_Matrix<ScalarT, IdxT>& getJacobian() = 0;
        virtual const COO_Matrix<ScalarT, IdxT>& getJacobian() const = 0;

        virtual std::vector<ScalarT>& getIntegrand() = 0;
        virtual const std::vector<ScalarT>& getIntegrand() const = 0;

        virtual std::vector<ScalarT>& getAdjointResidual() = 0;
        virtual const std::vector<ScalarT>& getAdjointResidual() const = 0;

        virtual std::vector<ScalarT>& getAdjointIntegrand() = 0;
        virtual const std::vector<ScalarT>& getAdjointIntegrand() const = 0;

    };


} // namespace ModelLib

#endif // _MODEL_EVALUATOR_HPP_