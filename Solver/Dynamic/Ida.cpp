/*
 * 
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 * 
 * This file is part of GridKit. For details, see github.com/LLNL/GridKit 
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


#include <iostream>
#include <iomanip>

#include "ModelEvaluator.hpp"
#include "Ida.hpp"
#include <klu.h>
#include <idas/idas.h>
#include <idas/idas_klu.h>
#include <idas/idas_dense.h>


namespace AnalysisManager
{

namespace Sundials 
{

    template <typename ScalarT, typename T, typename I>
    Ida<ScalarT, T, I>::Ida(ModelLib::ModelEvaluator<ScalarT, T, I>* model) : DynamicSolver<ScalarT, T, I>(model)
    {
        solver_ = IDACreate();
        tag_ = NULL;
    }
    
    template <typename ScalarT, typename T, typename I>
    Ida<ScalarT, T, I>::~Ida()
    {
    }
            
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::configureSimulation()
    {
        int retval = 0;
        
        // Allocate solution vectors
        yy_ = N_VNew_Serial(model_->size());
        checkVectorAllocation(yy_, "N_VNew_Serial");
        yp_ = N_VClone(yy_);
        checkVectorAllocation(yp_, "N_VClone");
        
        // Create vectors to store restart initial condition
        yy0_ = N_VClone(yy_);
        checkVectorAllocation(yy0_, "N_VClone");
        yp0_ = N_VClone(yy_);
        checkVectorAllocation(yp0_, "N_VClone");
        
        // Dummy initial time; will be overridden. 
        const realtype t0 = RCONST(0.0);
  
        // Allocate and initialize IDA workspace
        retval = IDAInit(solver_, this->Residual, t0, yy_, yp_);
        checkOutput(retval, "IDAInit");
        
        // Set pointer to model data
        retval = IDASetUserData(solver_, model_);
        checkOutput(retval, "IDASetUserData");
        
        // Set tolerances        
        realtype rtol;
        realtype atol;

        model_->setTolerances(rtol, atol); ///< \todo Function name should be "getTolerances"!
        retval = IDASStolerances(solver_, rtol, atol);
        checkOutput(retval, "IDASStolerances");
                
        // Tag differential variables
        std::vector<bool>& tag = model_->tag();
        if (tag.size() == model_->size())
        {
            tag_ = N_VClone(yy_);
            checkVectorAllocation(tag_, "N_VClone");
            model_->tagDifferentiable();
            copyVec(tag, tag_);
            
            retval = IDASetId(solver_, tag_);
            checkOutput(retval, "IDASetId");
            retval = IDASetSuppressAlg(solver_, TRUE);
            checkOutput(retval, "IDASetSuppressAlg");
        }
        
        // Set up linear solver
        retval = IDADense(solver_, model_->size());
        checkOutput(retval, "IDADense");
        
        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::getDefaultInitialCondition()
    {
        model_->initialize();
        
        copyVec(model_->y(), yy_);
        copyVec(model_->yp(), yp_);
        
        return 0;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::setIntegrationTime(T t_init, T t_final, int nout)
    {
        t_init_  = t_init;
        t_final_ = t_final;
        nout_    = nout;
        return 0;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::initializeSimulation(T t0, bool findConsistent)
    {
        int retval = 0;
        
        // Need to reinitialize IDA to set to get correct initial conditions
        retval = IDAReInit(solver_, t0, yy_, yp_);
        checkOutput(retval, "IDAReInit");
        
        // Find a consistent set of initial conditions for DAE
        if (findConsistent)
        {
            int initType = IDA_Y_INIT;
            
            if (tag_)
                initType = IDA_YA_YDP_INIT;
            
            retval = IDACalcIC(solver_, initType, 0.1);
            checkOutput(retval, "IDACalcIC");
        }

        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::runSimulation(T tf, int nout)
    {
        int retval = 0;
        int iout = 0;
        T tret;
        T dt = tf/nout;
        T tout = dt;
        
        /* In loop, call IDASolve, print results, and test for error.
         *     Break out of loop when NOUT preset output times have been reached. */
        while(nout > iout) 
        {
            retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
            checkOutput(retval, "IDASolve");
            
            if (retval == IDA_SUCCESS) 
            {
                ++iout;
                tout += dt;
            }
        }
        
        return retval;        
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::deleteSimulation()
    {
        IDAFree(&solver_);
        N_VDestroy(yy_);
        N_VDestroy(yp_);
        return 0;
    }
    
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::configureQuadrature()
    {
        int retval = 0;
        
        // Create and initialize quadratures 
        q_ = N_VNew_Serial(model_->size_quad());
        checkVectorAllocation(q_, "N_VNew_Serial");
        
        // Set integrand function and allocate quadrature workspace
        retval = IDAQuadInit(solver_, this->Integrand, q_);
        checkOutput(retval, "IDAQuadInit");
        
        // Set tolerances and error control for quadratures
        T rtol, atol;
        model_->setTolerances(rtol, atol);
        
        // Set tolerances for quadrature stricter than for integration
        retval = IDAQuadSStolerances(solver_, rtol*0.1, atol*0.1);
        checkOutput(retval, "IDAQuadSStolerances");
        
        // Include quadrature in eror checking
        retval = IDASetQuadErrCon(solver_, TRUE);
        checkOutput(retval, "IDASetQuadErrCon");
        
        return retval;
    }
    
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::initializeQuadrature()
    {
        int retval = 0;
        
        // Set all quadratures to zero
        N_VConst(RCONST(0.0), q_);
    
        // Initialize quadratures
        retval = IDAQuadReInit(solver_, q_);
        checkOutput(retval, "IDAQuadInit");
        
        return retval;
    }
    
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::runSimulationQuadrature(T tf, int nout)
    {
        int retval = 0;
        T tret;
        
        //std::cout << "Forward integration for initial value problem ... \n";

        T dt = tf/nout;
        T tout = dt;
        for(int i = 0; i < nout; ++i) 
        {
            retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
            checkOutput(retval, "IDASolve");
            
            if (retval == IDA_SUCCESS) 
            {
                tout += dt;
            }
            
            retval = IDAGetQuad(solver_, &tret, q_);
            checkOutput(retval, "IDASolve");
        }
    
        return retval;
    }
    
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::deleteQuadrature()
    {
        IDAQuadFree(solver_);
        N_VDestroy(q_);
  
        return 0;
    }
    
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::configureAdjoint()
    {
        // Allocate adjoint vector, derivatives and quadrature
        yyB_ = N_VNew_Serial(model_->size());
        checkVectorAllocation(yyB_, "N_VNew_Serial");
        
        ypB_ = N_VClone(yyB_);
        checkVectorAllocation(ypB_, "N_VClone");
        
        qB_ = N_VNew_Serial(model_->size_quad());
        checkVectorAllocation(qB_, "N_VNew_Serial");
        
        return 0;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::initializeAdjoint(I steps)
    {
        int retval = 0;
        
        // Create adjoint workspace
        retval = IDAAdjInit(solver_, steps, IDA_HERMITE);
        checkOutput(retval, "IDAAdjInit");
        
        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::initializeBackwardSimulation(T tf)
    {
        int retval = 0;
        realtype rtol;
        realtype atol;
        
        model_->initializeAdjoint();
        
        copyVec(model_->yB(),  yyB_);
        copyVec(model_->ypB(), ypB_);
        N_VConst(0.0, qB_);
        
        retval = IDACreateB(solver_, &backwardID_);
        checkOutput(retval, "IDACreateB");
        
        // IDAInitB must be called after forward simulation run. 
        retval = IDAInitB(solver_, backwardID_, this->adjointResidual, tf, yyB_, ypB_);
        checkOutput(retval, "IDAInitB");
        
        model_->setTolerances(rtol, atol);
        retval = IDASStolerancesB(solver_, backwardID_, rtol, atol);
        checkOutput(retval, "IDASStolerancesB");
        
        retval = IDASetUserDataB(solver_, backwardID_, model_);
        checkOutput(retval, "IDASetUserDataB");
        
        retval = IDASetMaxNumStepsB(solver_, backwardID_, 1000);
        
        retval = IDADenseB(solver_, backwardID_, model_->size());
        checkOutput(retval, "IDADenseB");
        
        // Also reinitialize quadratures. 
        retval = IDAQuadInitB(solver_, backwardID_, this->adjointIntegrand, qB_);
        checkOutput(retval, "IDAQuadInitB");
        
        retval = IDAQuadSStolerancesB(solver_, backwardID_, rtol*0.1, atol*0.1);
        checkOutput(retval, "IDAQuadSStolerancesB");
        
        // Include quadratures in error control
        retval = IDASetQuadErrConB(solver_, backwardID_, TRUE);
        checkOutput(retval, "IDASetQuadErrConB");
        
        
        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::runForwardSimulation(T tf, int nout)
    {
        int retval = 0;
        int ncheck;
        T time;
        
        //std::cout << "Forward integration for adjoint analysis ... \n";
        
        T dt = tf/nout;
        T tout = dt;
        for(int i = 0; i < nout; ++i) 
        {
            retval = IDASolveF(solver_, tout, &time, yy_, yp_, IDA_NORMAL, &ncheck);
            checkOutput(retval, "IDASolveF");
            
            if (retval == IDA_SUCCESS) 
            {
                tout += dt;
            }
            
            retval = IDAGetQuad(solver_, &time, q_);
            checkOutput(retval, "IDASolve");
        }
        
        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::runBackwardSimulation(T t_init)
    {
        int retval = 0;
        long int nstB;
        T time;
        
        //std::cout << "Backward integration for adjoint analysis ... ";
        
        retval = IDASolveB(solver_, t_init, IDA_NORMAL);
        checkOutput(retval, "IDASolveB");
        
        IDAGetNumSteps(IDAGetAdjIDABmem(solver_, backwardID_), &nstB);
        //std::cout << "done ( nst = " << nstB << " )\n";
        
        retval = IDAGetB(solver_, backwardID_, &time, yyB_, ypB_);
        checkOutput(retval, "IDAGetB");
        
        retval = IDAGetQuadB(solver_, backwardID_, &time, qB_);
        checkOutput(retval, "IDAGetB");
        
        return retval;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::deleteAdjoint()
    {
        IDAAdjFree(solver_);
        return 0;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::Residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    {
        ModelLib::ModelEvaluator<ScalarT, T, I>* model = static_cast<ModelLib::ModelEvaluator<ScalarT, T, I>*>(user_data);
        
        model->updateTime(tres, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        
        model->evaluateResidual();
        const std::vector<ScalarT>& f = model->getResidual();
        copyVec(f, rr);

        return 0;
    }
    

    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::Integrand(realtype tt, N_Vector yy, N_Vector yp, N_Vector rhsQ, void *user_data)
    {
        ModelLib::ModelEvaluator<ScalarT, T, I>* model = static_cast<ModelLib::ModelEvaluator<ScalarT, T, I>*>(user_data);
        
        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        
        model->evaluateIntegrand();
        const std::vector<ScalarT>& g = model->getIntegrand();
        copyVec(g, rhsQ);

        return 0;
    }
    
    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::adjointResidual(realtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void *user_data)
    {
        ModelLib::ModelEvaluator<ScalarT, T, I>* model = static_cast<ModelLib::ModelEvaluator<ScalarT, T, I>*>(user_data);
        
        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        copyVec(yyB, model->yB());
        copyVec(ypB, model->ypB());
        
        model->evaluateAdjointResidual();
        const std::vector<ScalarT>& fB = model->getAdjointResidual();
        copyVec(fB, rrB);
        
        return 0;
    }
    

    template <typename ScalarT, typename T, typename I>
    int Ida<ScalarT, T, I>::adjointIntegrand(realtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rhsQB, void *user_data)
    {
        ModelLib::ModelEvaluator<ScalarT, T, I>* model = static_cast<ModelLib::ModelEvaluator<ScalarT, T, I>*>(user_data);
        
        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        copyVec(yyB, model->yB());
        copyVec(ypB, model->ypB());
        
        model->evaluateAdjointIntegrand();
        const std::vector<ScalarT>& gB = model->getAdjointIntegrand();
        copyVec(gB, rhsQB);

        return 0;
    }
    

    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::copyVec(const N_Vector x, std::vector< ScalarT >& y)
    {
        const ScalarT* xdata = NV_DATA_S(x);
        for(unsigned int i = 0; i < y.size(); ++i)
        {
            y[i] = xdata[i];
        }
    }
    
    
    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::copyVec(const std::vector< ScalarT >& x, N_Vector y)
    {
        ScalarT* ydata = NV_DATA_S(y);
        for(unsigned int i = 0; i < x.size(); ++i)
        {
            ydata[i] = x[i];
        }
    }
    
    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::copyVec(const std::vector< bool >& x, N_Vector y)
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
    
    
    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::printOutput(realtype t)
    {
        realtype *yval = N_VGetArrayPointer_Serial(yy_);
        
        std::cout << std::setprecision(5) << std::setw(7) << t << " ";
        for (I i = 0; i < model_->size(); ++i)
        {
            std::cout << yval[i] << " ";
        }
        std::cout << "\n";
    }
    
    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::printFinalStats()
    {
        int retval = 0;
        void* mem = solver_;
        long int nst;
        long int nre;
        long int nje;
        long int nni;
        long int netf;
        long int ncfn;
        
        retval = IDAGetNumSteps(mem, &nst);
        checkOutput(retval, "IDAGetNumSteps");
        retval = IDAGetNumResEvals(mem, &nre);
        checkOutput(retval, "IDAGetNumResEvals");
        retval = IDASlsGetNumJacEvals(mem, &nje);
        checkOutput(retval, "IDASlsGetNumJacEvals");
        retval = IDAGetNumNonlinSolvIters(mem, &nni);
        checkOutput(retval, "IDAGetNumNonlinSolvIters");
        retval = IDAGetNumErrTestFails(mem, &netf);
        checkOutput(retval, "IDAGetNumErrTestFails");
        retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
        checkOutput(retval, "IDAGetNumNonlinSolvConvFails");
        
        std::cout << "\nFinal Run Statistics: \n\n";
        std::cout << "Number of steps                    = " <<  nst  << "\n";
        std::cout << "Number of residual evaluations     = " <<  nre  << "\n";
        std::cout << "Number of Jacobian evaluations     = " <<  nje  << "\n";
        std::cout << "Number of nonlinear iterations     = " <<  nni  << "\n";
        std::cout << "Number of error test failures      = " <<  netf << "\n";
        std::cout << "Number of nonlinear conv. failures = " <<  ncfn << "\n";
    }
    
    
    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::checkVectorAllocation(N_Vector v, const char* functionName)
    {
        if (v == NULL)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n"; 
            throw SundialsException();
        }
    }

    template <typename ScalarT, typename T, typename I>
    void Ida<ScalarT, T, I>::checkOutput(int retval, const char* functionName)
    {
        if (retval < 0)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n"; 
            throw SundialsException();
        }
    }

    // Compiler will prevent building modules with data type incompatible with realtype
    template class Ida<realtype, realtype, long int>;
    template class Ida<realtype, realtype, size_t>;
    
} // namespace Sundials
} // namespace AnalysisManager


