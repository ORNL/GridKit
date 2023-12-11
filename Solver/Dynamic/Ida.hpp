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


#ifndef _IDA_HPP_
#define _IDA_HPP_

#include <iostream>
#include <exception>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>       /* access to KLU linear solver          */
#include <sunlinsol/sunlinsol_dense.h>     /* access to dense linear solver        */

#include "ModelEvaluator.hpp"
#include "DynamicSolver.hpp"

namespace AnalysisManager
{
    namespace Sundials
    {
        template <class ScalarT, typename IdxT>
        class Ida : public DynamicSolver<ScalarT, IdxT>
        {
            using DynamicSolver<ScalarT, IdxT>::model_;

            typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

        public:
            Ida(ModelLib::ModelEvaluator<ScalarT, IdxT>* model);
            ~Ida();

            int configureSimulation();
            int configureLinearSolver();
            int getDefaultInitialCondition();
            int setIntegrationTime(real_type t_init, real_type t_final, int nout);
            int initializeSimulation(real_type t0, bool findConsistent=false);
            int runSimulation(real_type tf, int nout=1);
            int deleteSimulation();

            int configureQuadrature();
            int initializeQuadrature();
            int runSimulationQuadrature(real_type tf, int nout=1);
            int deleteQuadrature();

            int configureAdjoint();
            int configureLinearSolverBackward();
            int initializeAdjoint(IdxT steps = 100);
            int initializeBackwardSimulation(real_type tf);
            int runForwardSimulation(real_type tf, int nout=1);
            int runBackwardSimulation(real_type t0);
            int deleteAdjoint();


            int saveInitialCondition()
            {
                N_VScale(1.0, yy_, yy0_);
                N_VScale(1.0, yp_, yp0_);
                return 0;
            }

            int getSavedInitialCondition()
            {
                N_VScale(1.0, yy0_, yy_);
                N_VScale(1.0, yp0_, yp_);
                return 0;
            }

            real_type getInitialTime()
            {
                return t_init_;
            }

            real_type getFinalTime()
            {
                return t_final_;
            }

            int getNumberOutputTimes()
            {
                return nout_;
            }

            const real_type* getIntegral() const
            {
                return NV_DATA_S(q_);
            }

            real_type* getIntegral()
            {
                return NV_DATA_S(q_);
            }

            const real_type* getAdjointIntegral() const
            {
                return NV_DATA_S(qB_);
            }

            real_type* getAdjointIntegral()
            {
                return NV_DATA_S(qB_);
            }

            void printOutput(realtype t);
            void printSpecial(realtype t, N_Vector x);
            void printFinalStats();

        private:
            static int Residual(realtype t,
                                N_Vector yy, N_Vector yp,
                                N_Vector rr, void *user_data);
            
            static int Jac(realtype t, realtype cj,
               N_Vector yy, N_Vector yp, N_Vector resvec,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

            static int Integrand(realtype t,
                                 N_Vector yy,   N_Vector yp,
                                 N_Vector rhsQ, void *user_data);

            static int adjointResidual(realtype t,
                                       N_Vector yy,  N_Vector yp,
                                       N_Vector yyB, N_Vector ypB,
                                       N_Vector rrB, void *user_data);

            static int adjointIntegrand(realtype t,
                                        N_Vector yy,  N_Vector yp,
                                        N_Vector yyB, N_Vector ypB,
                                        N_Vector rhsQB, void *user_data);

        private:
            void* solver_;
            SUNContext context_;
            SUNMatrix JacobianMat_;
            SUNMatrix JacobianMatB_;
            SUNLinearSolver linearSolver_;
            SUNLinearSolver linearSolverB_;

            real_type t_init_;
            real_type t_final_;
            int nout_; ///< Number of integration outputs

            N_Vector yy_;  ///< Solution vector
            N_Vector yp_;  ///< Solution derivatives vector
            N_Vector tag_; ///< Tags differential variables
            N_Vector q_;   ///< Integrand vector

            N_Vector yy0_; ///< Storage for initial values
            N_Vector yp0_; ///< Storage for initial derivatives

            N_Vector yyB_; ///< Adjoint solution vector
            N_Vector ypB_; ///< Adjoint solution derivatives vector
            N_Vector qB_;  ///< Backward integrand vector

            int backwardID_;

        private:
            //static void copyMat(ModelEvaluator::Mat& J, SlsMat Jida);
            static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
            static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
            static void copyVec(const std::vector<bool>& x, N_Vector y);

            //int check_flag(void *flagvalue, const char *funcname, int opt);
            inline void checkAllocation(void* v, const char* functionName);
            inline void checkOutput(int retval, const char* functionName);

        };

        /// Simple exception to use within Ida class.
        class SundialsException : public std::exception
        {
            virtual const char* what() const throw()
            {
                return "Method in Ida class failed!\n";
            }
        };


    } // namespace Sundials


} // namespace AnalysisManager


#endif // _IDA_HPP_
