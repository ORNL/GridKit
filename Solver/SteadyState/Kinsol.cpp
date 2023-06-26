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

/**
 * @file Kinsol.cpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 * 
 * Contains definition of interface to KINSOL nonlinear solver from
 * SUNDIALS library.
 * 
 */

#include <iostream>
#include <iomanip>

#include <kinsol/kinsol.h>             // access to KINSOL func., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector      
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix      
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver

#include "ModelEvaluator.hpp"
#include "Kinsol.hpp"


namespace AnalysisManager
{

namespace Sundials
{

    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::Kinsol(ModelLib::ModelEvaluator<ScalarT, IdxT>* model) 
        : SteadyStateSolver<ScalarT, IdxT>(model)
    {
        int retval = 0;

        // Create the SUNDIALS context that all SUNDIALS objects require
        retval = SUNContext_Create(NULL, &context_);
        checkOutput(retval, "SUNContext");

        solver_ = KINCreate(context_);
        tag_ = NULL;
    }

    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::~Kinsol()
    {
        SUNContext_Free(&context_);
        KINFree(&solver_);

        N_VDestroy_Serial(this->yy_);
        N_VDestroy_Serial(this->yy0_);
        N_VDestroy_Serial(this->scale_);

        SUNMatDestroy(this->JacobianMat_);
        SUNLinSolFree_Dense(this->linearSolver_);
        
        solver_ = 0;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::configureSimulation()
    {
        int retval = 0;

        // Allocate solution vectors
        yy_ = N_VNew_Serial(model_->size(), context_);
        checkAllocation((void*) yy_, "N_VNew_Serial");

        // Allocate scaling vector
        scale_ = N_VClone(yy_);
        checkAllocation((void*) scale_, "N_VClone");

        // Create vectors to store restart initial condition
        yy0_ = N_VClone(yy_);
        checkAllocation((void*) yy0_, "N_VClone");

        // Allocate and initialize KIN workspace
        retval = KINInit(solver_, this->Residual, yy_);
        checkOutput(retval, "KINInit");

        // Set pointer to model data
        retval = KINSetUserData(solver_, model_);
        checkOutput(retval, "KINSetUserData");

        // Set output verbosity level
        retval = KINSetPrintLevel(solver_, 0);
        checkOutput(retval, "KINSetPrintLevel");

        // Set tolerances
        realtype fnormtol;  ///< Residual tolerance
        realtype scsteptol; ///< Scaled step tolerance

        model_->setTolerances(fnormtol, scsteptol); ///< \todo Function name should be "getTolerances"!
        retval = KINSetFuncNormTol(solver_, fnormtol);
        checkOutput(retval, "KINSetFuncNormTol");

        retval = KINSetScaledStepTol(solver_, scsteptol);
        checkOutput(retval, "KINSetScaledStepTol");

        // Set up linear solver
        JacobianMat_ = SUNDenseMatrix(model_->size(), model_->size(), context_);
        checkAllocation((void*) JacobianMat_, "SUNDenseMatrix");

        linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
        checkAllocation((void*) linearSolver_, "SUNLinSol_Dense");

        retval = KINSetLinearSolver(solver_, linearSolver_, JacobianMat_);
        checkOutput(retval, "KINSetLinearSolver");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::configureLinearSolver()
    {
        int retval = 0;

        // Set up linear solver
        JacobianMat_ = SUNDenseMatrix(model_->size(), model_->size(), context_);
        checkAllocation((void*) JacobianMat_, "SUNDenseMatrix");

        linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
        checkAllocation((void*) linearSolver_, "SUNLinSol_Dense");

        retval = KINSetLinearSolver(solver_, linearSolver_, JacobianMat_);
        checkOutput(retval, "KINSetLinearSolver");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::getDefaultInitialCondition()
    {
        model_->initialize();

        copyVec(model_->y(), yy_);

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::runSimulation()
    {
        int retval = 0;
        N_VConst(1.0, scale_);
        retval = KINSol(solver_, yy_, KIN_LINESEARCH, scale_, scale_);
        checkOutput(retval, "KINSol");
        //printOutput(tout); 
        //std::cout << "\n";
        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::deleteSimulation()
    {
        SUNLinSolFree(linearSolver_);
        KINFree(&solver_);
        N_VDestroy(yy_);
        N_VDestroy(scale_);
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::Residual(N_Vector yy, N_Vector rr, void *user_data)
    {
        ModelLib::ModelEvaluator<ScalarT, IdxT>* model = 
            static_cast<ModelLib::ModelEvaluator<ScalarT, IdxT>*>(user_data);

        copyVec(yy, model->y());

        model->evaluateResidual();
        const std::vector<ScalarT>& f = model->getResidual();
        copyVec(f, rr);

        return 0;
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const N_Vector x, std::vector< ScalarT >& y)
    {
        const ScalarT* xdata = NV_DATA_S(x);
        for(unsigned int i = 0; i < y.size(); ++i)
        {
            y[i] = xdata[i];
        }
    }


    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const std::vector< ScalarT >& x, N_Vector y)
    {
        ScalarT* ydata = NV_DATA_S(y);
        for(unsigned int i = 0; i < x.size(); ++i)
        {
            ydata[i] = x[i];
        }
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const std::vector< bool >& x, N_Vector y)
    {
        ScalarT* ydata = NV_DATA_S(y);
        for(unsigned int i = 0; i < x.size(); ++i)
        {
            if (x[i])
                ydata[i] = 1.0;
            else
                ydata[i] = 0.0;
        }
    }


    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printOutput()
    {
        realtype *yval  = N_VGetArrayPointer_Serial(yy_);

        std::cout << std::setprecision(5) << std::setw(7);
        for (IdxT i = 0; i < model_->size(); ++i)
        {
            std::cout << yval[i] << " ";
        }
        std::cout << "\n";
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printSpecial(realtype t, N_Vector y)
    {
        realtype *yval = N_VGetArrayPointer_Serial(y);
        IdxT N = N_VGetLength_Serial(y);
        std::cout << "{";
        std::cout << std::setprecision(5) << std::setw(7) << t;
        for (IdxT i = 0; i < N; ++i)
        {
            std::cout << ", " << yval[i];
        }
        std::cout << "},\n";
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printFinalStats()
    {
        int retval = 0;
        void* mem = solver_;
        long int nni;
        long int nfe;
        long int nje;
        long int nlfe;

        // retval = KINGetNumSteps(mem, &nst);
        // checkOutput(retval, "KINGetNumSteps");
        retval = KINGetNumNonlinSolvIters(mem, &nni);
        checkOutput(retval, "KINGetNumNonlinSolvIters");
        retval = KINGetNumFuncEvals(mem, &nfe);
        checkOutput(retval, "KINGetNumFuncEvals");
        retval = KINGetNumJacEvals(mem, &nje);
        checkOutput(retval, "KINGetNumJacEvals");
        retval = KINGetNumLinFuncEvals(mem, &nlfe);
        checkOutput(retval, "KINGetNumLinFuncEvals");

        // std::cout << "\nFinal Run Statistics: \n\n";
        std::cout << "Number of nonlinear iterations     = " <<  nni  << "\n";
        std::cout << "Number of function evaluations     = " <<  nfe  << "\n";
        std::cout << "Number of Jacobian evaluations     = " <<  nje  << "\n";
        std::cout << "Number of linear function evals.   = " <<  nlfe << "\n";
    }


    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::checkAllocation(void* v, const char* functionName)
    {
        if (v == NULL)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n";
            throw SundialsException();
        }
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::checkOutput(int retval, const char* functionName)
    {
        if (retval < 0)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n";
            throw SundialsException();
        }
    }

    // Compiler will prevent building modules with data type incompatible with realtype
    template class Kinsol<realtype, long int>;
    template class Kinsol<realtype, int>;
    template class Kinsol<realtype, size_t>;

} // namespace Sundials
} // namespace AnalysisManager
