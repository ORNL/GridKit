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



#ifndef _IPOPT_DYNAMIC_CONSTRAINT_HPP_
#define _IPOPT_DYNAMIC_CONSTRAINT_HPP_

#include <IpTNLP.hpp>
#include "OptimizationSolver.hpp"

namespace AnalysisManager {

    namespace IpoptInterface {

        /**
         * Implementation of Ipopt's pure virtual TNLP class.
         *
         * TNLP defines Ipopt's interface to the model. This is in fact
         * the model evaluator interface to Ipopt. In this case however,
         * the model evaluator calls dynamic solver to compute the objective
         * and the gradient.
         *
         * \note This clas is based on Cosmin's reformulation of the dynamic
         * constrained optimization problem. For now it is hard-wired to
         * 1-parameter optimization problems.
         *
         */
        template <class ScalarT, typename IdxT>
        class DynamicConstraint : public Ipopt::TNLP, public OptimizationSolver<ScalarT, IdxT>
        {
            using OptimizationSolver<ScalarT, IdxT>::integrator_;
            using OptimizationSolver<ScalarT, IdxT>::model_;

            typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

            typedef Ipopt::Index Index;
            typedef Ipopt::Number Number;
            typedef Ipopt::SolverReturn SolverReturn;
            typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;
            typedef Ipopt::IpoptData IpoptData;

        public:
            DynamicConstraint(Sundials::Ida<ScalarT, IdxT>* integrator);
            virtual ~DynamicConstraint();

            /// Returns sizes of the model components
            virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                      Index& nnz_h_lag, IndexStyleEnum& index_style);

            /// Returns problem bounds
            virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                         Index m, Number* g_l, Number* g_u);

            /// Initialize optimization
            virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                            bool init_z, Number* z_L, Number* z_U,
                                            Index m, bool init_lambda,
                                            Number* lambda);

            /// Evaluate objective
            virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

            /// Evaluate objective gradient
            virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

            /// Evaluate constraint residuals (not used here)
            virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

            /// Evaluate Jacobian (not used here)
            virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                                    Index m, Index nele_jac, Index* iRow, Index *jCol,
                                    Number* values);

            /// Evaluate Hessian (have Ipopt estimate Hessian)
            virtual bool eval_h(Index n, const Number* x, bool new_x,
                                Number obj_factor, Index m, const Number* lambda,
                                bool new_lambda, Index nele_hess, Index* iRow,
                                Index* jCol, Number* values);

            /// Postprocessing of the results (not used here)
            virtual void finalize_solution(SolverReturn status,
                                           Index n, const Number* x, const Number* z_L, const Number* z_U,
                                           Index m, const Number* g, const Number* lambda,
                                           Number obj_value,
                                           const IpoptData* ip_data,
                                           IpoptCalculatedQuantities* ip_cq);
        private:
            real_type t_init_;
            real_type t_final_;
            int nout_;
        };

    } // namespace IpoptInterface
} // namespace AnalysisManager

#endif // _IPOPT_DYNAMIC_CONSTRAINT_HPP_
