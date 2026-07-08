
#include "Ida.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

#include <idas/idas.h>
#include <idas/idas_ls.h>

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{

  namespace Sundials
  {

    template <class ScalarT, typename IdxT>
    Ida<ScalarT, IdxT>::Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : DynamicSolver<ScalarT, IdxT>(model)
    {
      int retval = 0;

      // Create the SUNDIALS context that all SUNDIALS objects require
      retval = SUNContext_Create(SUN_COMM_NULL, &context_);
      checkOutput(retval, "SUNContext");
      solver_ = IDACreate(context_);
    }

    /**
     * @brief Destroy the Ida< Scalar T,  Idx T>:: Ida object
     *
     * @note if sysmodel is freed before this will fail. May want something agnostic to this
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    Ida<ScalarT, IdxT>::~Ida()
    {
      deleteQuadrature();
      deleteAdjoint();
      deleteSimulation();
      deleteBackwardSimulation();
      SUNContext_Free(&context_);
    }

    /**
     * @brief Configure the simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureSimulation()
    {
      int retval = 0;

      // Allocate solution vectors
      yy_ = N_VNew_Serial(static_cast<sunindextype>(model_->size()), context_);
      checkAllocation((void*) yy_, "N_VNew_Serial");
      yp_ = N_VClone(yy_);
      checkAllocation((void*) yp_, "N_VClone");

      // get intial conditions
      this->getDefaultInitialCondition();

      // Create vectors to store restart initial condition
      yy0_ = N_VClone(yy_);
      checkAllocation((void*) yy0_, "N_VClone");
      yp0_ = N_VClone(yp_);
      checkAllocation((void*) yp0_, "N_VClone");

      // Dummy initial time; will be overridden.
      const RealT t0 = 0.0;

      // Allocate and initialize IDA workspace
      retval = IDAInit(solver_, this->Residual, t0, yy_, yp_);
      checkOutput(retval, "IDAInit");

      // Set pointer to model data
      retval = IDASetUserData(solver_, model_);
      checkOutput(retval, "IDASetUserData");

      // Tag differential variables
      tag_ = N_VClone(yy_);
      checkAllocation((void*) tag_, "N_VClone");
      model_->tagDifferentiable();
      copyVec(model_->tag(), tag_);

      retval = IDASetId(solver_, tag_);
      checkOutput(retval, "IDASetId");

      setIDAOptions(solver_, time_step_, rel_tol_, abs_tol_override_, max_steps_, suppress_alg_);

      // Set up linear solver
      return this->configureLinearSolver();
    }

    /**
     * @brief Configure the linear solver
     *
     * @note This currently uses pre-processor directives to set dense or sparse
     * linear solvers
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolver()
    {
      int retval = 0;

/// Todo - Implement a cleaner way to handle the Sparse versus Dense solvers
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
      if (model_->hasJacobian())
      {
        this->configureLinearSolverSparse();
      }
      else
      {
        this->configureLinearSolverDense();
      }
#else
      /// Todo - Improve error handling capabilities and hasJacobian_ ownership
      if (model_->hasJacobian())
      {
        Log::warning() << "SUNDIALS is not configured with KLU, but the model has a (sparse) Jacobian. "
                       << "Falling back to dense Jacobian.\n";
      }
      this->configureLinearSolverDense();
#endif

      return retval;
    }

#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
    /**
     * @brief Configure a sparse linear solver
     *
     * @note This method is only available if SUNDIALS is configured with KLU
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolverSparse()
    {
      int retval = 0;

      sunindextype n   = static_cast<sunindextype>(model_->size());
      sunindextype nnz = static_cast<sunindextype>(model_->getCsrJacobian()->getNnz());

      JacobianMat_ = SUNSparseMatrix(n,
                                     n,
                                     nnz,
                                     CSR_MAT,
                                     context_);
      checkAllocation((void*) JacobianMat_, "SUNSparseMatrix");

      linearSolver_ = SUNLinSol_KLU(yy_, JacobianMat_, context_);
      checkAllocation((void*) linearSolver_, "SUNLinSol_KLU");

      retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
      checkOutput(retval, "IDASetLinearSolver");

      retval = IDASetJacFn(solver_, this->Jac);
      checkOutput(retval, "IDASetJacFn");

      return retval;
    }
#endif

    /**
     * @brief Configure a dense linear solver
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolverDense()
    {
      int retval = 0;

      JacobianMat_ = SUNDenseMatrix(static_cast<sunindextype>(model_->size()),
                                    static_cast<sunindextype>(model_->size()),
                                    context_);
      checkAllocation((void*) JacobianMat_, "SUNDenseMatrix");

      linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
      checkAllocation((void*) linearSolver_, "SUNLinSol_Dense");

      retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
      checkOutput(retval, "IDASetLinearSolver");

      return retval;
    }

    /**
     * @brief Get default initial condition
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::getDefaultInitialCondition()
    {
      model_->initialize();

      copyVec(model_->y(), yy_);
      copyVec(model_->yp(), yp_);

      return 0;
    }

    /**
     * @brief Set integration time
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::setIntegrationTime(RealT t_init, RealT t_final, int nout)
    {
      t_init_  = t_init;
      t_final_ = t_final;
      nout_    = nout;
      return 0;
    }

    /**
     * @brief Initialize the simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeSimulation(RealT t0, bool findConsistent)
    {
      int retval = 0;

      t_init_ = t0;

      // Need to reinitialize IDA to set to get correct initial conditions
      retval = IDAReInit(solver_, t0, yy_, yp_);
      checkOutput(retval, "IDAReInit");

      // Find a consistent set of initial conditions for DAE
      if (findConsistent)
      {
        int initType = IDA_Y_INIT;

        if (tag_)
          initType = IDA_YA_YDP_INIT;

        retval = IDACalcIC(solver_, initType, t0 + 0.1);
        checkOutput(retval, "IDACalcIC");

        retval = IDAGetConsistentIC(solver_, yy_, yp_);
        checkOutput(retval, "IDAGetConsistentIC");

        copyVec(yy_, model_->y());
        copyVec(yp_, model_->yp());
      }

      return retval;
    }

    /**
     * @brief Run the IDA solver on the given model and produce a solution at
     * the given final time.
     *
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Matrix and vector index data type
     * @param tf The final simulation time.
     * @param nout The number of integration segmentstimes.
     * @param step_callback An optional callback which, if provided, will be
     * called after each time the IDA solver has been invoked with the value
     * of `t` that IDA has calculated the last step at. The provided model will
     * be updated with the latest values of `y` and `yp` before the callback is
     * invoked.
     * @return int zero if successful, error code otherwise.
     *
     * @note The actual time of the final IDA solution should be somewhat
     * close to `tf`, however due to rounding error the precise final time may
     * be before or after `tf`.
     *
     * @todo Consider adding initial time as the function argument, as well.
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulation(RealT tf, int nout, const std::optional<std::function<void(RealT)>> step_callback)
    {
      int   retval = 0;
      int   iout   = 0;
      RealT tret;
      RealT dt   = (tf - t_init_) / static_cast<RealT>(nout);
      RealT tout = t_init_ + dt;

      // In loop, call IDASolve, print results, and test for error.
      //  Break out of loop when NOUT preset output times have been reached.
      // printOutput(0.0);
      while (nout > iout)
      {
        retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
        checkOutput(retval, "IDASolve");

        if (step_callback.has_value() || model_->monitoring())
        {
          // The callback may try to observe upated values in the model, so we
          // should update them here (At this point, the model's values are one
          // internal integrator step out of date)
          copyVec(yy_, model_->y());
          copyVec(yp_, model_->yp());
          model_->updateTime(tret, 0.0);

          if (model_->monitoring())
          {
            model_->printMonitoredVariables();
          }
          if (step_callback.has_value())
          {
            (*step_callback)(tret);
          }
        }

        if (retval == IDA_SUCCESS)
        {
          ++iout;
          tout += dt;
        }
      }

      // Final copy out. No guarantee last residual evaluation is final step.
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());
      model_->updateTime(tf, 0.0);
      // if (model_->monitoring())
      // {
      //   model_->printMonitoredVariables();
      // }

      // std::cout << "\n";
      return retval;
    }

    /**
     * @brief Delete the simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteSimulation()
    {
      N_VDestroy(yy_);
      N_VDestroy(yp_);
      N_VDestroy(tag_);
      N_VDestroy(yy0_);
      N_VDestroy(yp0_);
      SUNLinSolFree(linearSolver_);
      SUNMatDestroy(JacobianMat_);
      IDAFree(&solver_);
      return 0;
    }

    /**
     * @brief Configure quadrature
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureQuadrature()
    {
      int retval = 0;

      // Create and initialize quadratures
      q_ = N_VNew_Serial(static_cast<sunindextype>(model_->sizeQuadrature()), context_);
      checkAllocation((void*) q_, "N_VNew_Serial");

      // Set integrand function and allocate quadrature workspace
      retval = IDAQuadInit(solver_, this->Integrand, q_);
      checkOutput(retval, "IDAQuadInit");

      // Set tolerances and error control for quadratures
      setQuadratureTolerance(solver_, quadrature_rel_tol_, quadrature_abs_tol_override_);

      // Include quadrature in eror checking
      retval = IDASetQuadErrCon(solver_, SUNTRUE);
      checkOutput(retval, "IDASetQuadErrCon");

      return retval;
    }

    /**
     * @brief Initialize quadrature
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeQuadrature()
    {
      int retval = 0;

      // Set all quadratures to zero
      N_VConst(0.0, q_);

      // Initialize quadratures
      retval = IDAQuadReInit(solver_, q_);
      checkOutput(retval, "IDAQuadInit");

      return retval;
    }

    /**
     * @brief Run simulation with quadrature
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulationQuadrature(RealT tf, int nout)
    {
      int   retval = 0;
      RealT tret;

      // std::cout << "Forward integration for initial value problem ... \n";

      RealT dt   = tf / static_cast<RealT>(nout);
      RealT tout = dt;
      // printOutput(0.0);
      // printSpecial(0.0, yy_);
      for (int i = 0; i < nout; ++i)
      {
        retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
        checkOutput(retval, "IDASolve");
        // printSpecial(tout, yy_);
        // printOutput(tout);

        if (retval == IDA_SUCCESS)
        {
          tout += dt;
        }

        retval = IDAGetQuad(solver_, &tret, q_);
        checkOutput(retval, "IDAGetQuad");
      }

      // Final copy out. No gaurentee last residual evaluation is final step.
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());
      model_->updateTime(tf, 0.0);

      return retval;
    }

    /**
     * @brief Delete quadrature
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteQuadrature()
    {
      IDAQuadFree(solver_);
      N_VDestroy(q_);

      return 0;
    }

    /**
     * @brief Configure adjoint
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureAdjoint()
    {
      // Allocate adjoint vector, derivatives and quadrature
      yyB_ = N_VNew_Serial(static_cast<sunindextype>(model_->size()), context_);
      checkAllocation((void*) yyB_, "N_VNew_Serial");

      ypB_ = N_VClone(yyB_);
      checkAllocation((void*) ypB_, "N_VClone");

      qB_ = N_VNew_Serial(static_cast<sunindextype>(model_->sizeParams()), context_);
      checkAllocation((void*) qB_, "N_VNew_Serial");

      return 0;
    }

    /**
     * @brief Initialize adjoint
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeAdjoint(IdxT steps)
    {
      int retval = 0;

      // Create adjoint workspace
      retval = IDAAdjInit(solver_, static_cast<long>(steps), IDA_HERMITE);
      checkOutput(retval, "IDAAdjInit");

      return retval;
    }

    /**
     * @brief Initialize backward simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeBackwardSimulation(RealT tf)
    {
      int retval = 0;

      model_->initializeAdjoint();

      copyVec(model_->yB(), yyB_);
      copyVec(model_->ypB(), ypB_);
      N_VConst(0.0, qB_);

      retval = IDACreateB(solver_, &backwardID_);
      checkOutput(retval, "IDACreateB");

      // IDAInitB must be called after forward simulation run.
      retval = IDAInitB(solver_, backwardID_, this->adjointResidual, tf, yyB_, ypB_);
      checkOutput(retval, "IDAInitB");

      setIDAOptions(IDAGetAdjIDABmem(solver_, backwardID_),
                    backward_time_step_,
                    backward_rel_tol_,
                    backward_abs_tol_override_,
                    backward_max_steps_,
                    backward_suppress_alg_);

      retval = IDASetUserDataB(solver_, backwardID_, model_);
      checkOutput(retval, "IDASetUserDataB");

      // Allocate Jacobian matrix, if not already
      if (JacobianMatB_ == nullptr)
      {
        JacobianMatB_ = SUNDenseMatrix(static_cast<sunindextype>(model_->size()),
                                       static_cast<sunindextype>(model_->size()),
                                       context_);
        checkAllocation((void*) JacobianMatB_, "SUNDenseMatrix");
      }

      // Allocate linear solver, if not already
      if (linearSolverB_ == nullptr)
      {
        linearSolverB_ = SUNLinSol_Dense(yyB_, JacobianMatB_, context_);
        checkAllocation((void*) linearSolverB_, "SUNLinSol_Dense");
      }

      // Setup linear solver (only dense supported at this time)
      retval = IDASetLinearSolverB(solver_, backwardID_, linearSolverB_, JacobianMatB_);
      checkOutput(retval, "IDASetLinearSolverB");

      // Also reinitialize quadratures.
      retval = IDAQuadInitB(solver_, backwardID_, this->adjointIntegrand, qB_);
      checkOutput(retval, "IDAQuadInitB");

      setQuadratureTolerance(IDAGetAdjIDABmem(solver_, backwardID_), backward_quadrature_rel_tol_, backward_quadrature_abs_tol_override_);

      // Include quadratures in error control
      retval = IDASetQuadErrConB(solver_, backwardID_, SUNTRUE);
      checkOutput(retval, "IDASetQuadErrConB");

      return retval;
    }

    /**
     * @brief Configure linear solver for backward simulation
     *
     * @note This only supports dense linear solvers at the moment
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolverBackward()
    {
      int retval = 0;

      // Create Jacobian matrix
      JacobianMatB_ = SUNDenseMatrix(static_cast<sunindextype>(model_->size()),
                                     static_cast<sunindextype>(model_->size()),
                                     context_);
      checkAllocation((void*) JacobianMatB_, "SUNDenseMatrix");

      // Create linear solver
      linearSolverB_ = SUNLinSol_Dense(yyB_, JacobianMatB_, context_);
      checkAllocation((void*) linearSolverB_, "SUNLinSol_Dense");

      // Attach linear solver to IDA
      retval = IDASetLinearSolverB(solver_, backwardID_, linearSolverB_, JacobianMatB_);
      checkOutput(retval, "IDASetLinearSolverB");

      return retval;
    }

    /**
     * @brief Run forward simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runForwardSimulation(RealT tf, int nout)
    {
      int   retval = 0;
      int   ncheck;
      RealT time;

      // std::cout << "Forward integration for adjoint analysis ... \n";

      RealT dt   = tf / static_cast<RealT>(nout);
      RealT tout = dt;
      for (int i = 0; i < nout; ++i)
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

      // Final copy out. No gaurentee last residual evaluation is final step.
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());
      model_->updateTime(tf, 0.0);

      return retval;
    }

    /**
     * @brief Run backward simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runBackwardSimulation(RealT t_init)
    {
      int      retval = 0;
      long int nstB;
      RealT    time;

      // std::cout << "Backward integration for adjoint analysis ... ";

      retval = IDASolveB(solver_, t_init, IDA_NORMAL);
      checkOutput(retval, "IDASolveB");

      IDAGetNumSteps(IDAGetAdjIDABmem(solver_, backwardID_), &nstB);
      // std::cout << "done ( nst = " << nstB << " )\n";

      retval = IDAGetB(solver_, backwardID_, &time, yyB_, ypB_);
      checkOutput(retval, "IDAGetB");

      // Copy back into model
      copyVec(yyB_, model_->yB());
      copyVec(ypB_, model_->ypB());

      retval = IDAGetQuadB(solver_, backwardID_, &time, qB_);
      checkOutput(retval, "IDAGetQuadB");

      return retval;
    }

    /**
     * @brief Delete adjoint
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteAdjoint()
    {
      IDAAdjFree(solver_);
      return 0;
    }

    /**
     * @brief Delete backward simulation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteBackwardSimulation()
    {
      N_VDestroy(yyB_);
      N_VDestroy(ypB_);
      N_VDestroy(qB_);
      SUNLinSolFree(linearSolverB_);
      SUNMatDestroy(JacobianMatB_);

      return 0;
    }

    /**
     * @brief Residual evaluation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Residual(RealT tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());
      copyVec(yp, model->yp());
      model->updateTime(tres, 0.0);

      model->evaluateResidual();
      copyVec(model->getResidual(), rr);

      return 0;
    }

    /**
     * @brief Jacobian evaluation
     *
     * @note The model Jacobian is stored in CSR format.
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Jac(RealT t, RealT cj, N_Vector yy, N_Vector yp, N_Vector, SUNMatrix J, void* user_data, N_Vector, N_Vector, N_Vector)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());
      copyVec(yp, model->yp());
      model->updateTime(t, cj);

      model->evaluateJacobian();

      using CsrMatrixT = GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>;
      CsrMatrixT* Jac  = model->getCsrJacobian();

      SUNMatZero(J);

      sunindextype* sun_row_ptrs = SUNSparseMatrix_IndexPointers(J);
      sunindextype* sun_cols     = SUNSparseMatrix_IndexValues(J);
      RealT*        sun_vals     = SUNSparseMatrix_Data(J);

      IdxT n   = Jac->getNumRows();
      IdxT nnz = Jac->getNnz();

      // Get reference to the jacobian entries
      IdxT*  row_ptrs = Jac->getRowData();
      IdxT*  cols     = Jac->getColData();
      RealT* vals     = Jac->getValues();

      // Copy data from model jac to sundials
      std::copy(row_ptrs, row_ptrs + n + 1, sun_row_ptrs);
      std::copy(cols, cols + nnz, sun_cols);
      std::copy(vals, vals + nnz, sun_vals);

      return 0;
    }

    /**
     * @brief Integrand evaluation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Integrand(RealT tt, N_Vector yy, N_Vector yp, N_Vector rhsQ, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());
      copyVec(yp, model->yp());
      model->updateTime(tt, 0.0);

      model->evaluateIntegrand();
      copyVec(model->getIntegrand(), rhsQ);

      return 0;
    }

    /**
     * @brief Adjoint residual evaluation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointResidual(RealT tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());
      copyVec(yp, model->yp());
      copyVec(yyB, model->yB());
      copyVec(ypB, model->ypB());
      model->updateTime(tt, 0.0);

      model->evaluateAdjointResidual();
      copyVec(model->getAdjointResidual(), rrB);

      return 0;
    }

    /**
     * @brief Adjoint integrand evaluation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointIntegrand(RealT tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rhsQB, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());
      copyVec(yp, model->yp());
      copyVec(yyB, model->yB());
      copyVec(ypB, model->ypB());
      model->updateTime(tt, 0.0);

      model->evaluateAdjointIntegrand();
      copyVec(model->getAdjointIntegrand(), rhsQB);

      return 0;
    }

    /**
     * @brief Copy SUNDIALS N_Vector to Vector
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const N_Vector x, VectorT& y)
    {
      const auto xsize = static_cast<size_t>(N_VGetLength(x));
      if (xsize != y.size())
      {
        std::cerr << "\nN_Vector size (" << xsize << ") does not match vector size ("
                  << y.size() << ").\n\n";
        throw SundialsException();
      }

      const ScalarT* xdata = N_VGetArrayPointer(x);
      std::copy_n(xdata, y.size(), y.getData(GridKit::memory::HOST));
    }

    /**
     * @brief Copy Vector to SUNDIALS N_Vector
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const VectorT& x, N_Vector y)
    {
      const auto ysize = static_cast<size_t>(N_VGetLength(y));
      if (x.size() != ysize)
      {
        std::cerr << "\nvector size (" << x.size() << ") does not match N_Vector size ("
                  << ysize << ").\n\n";
        throw SundialsException();
      }

      ScalarT* ydata = N_VGetArrayPointer(y);
      std::copy_n(x.getData(GridKit::memory::HOST), x.size(), ydata);
    }

    /**
     * @brief Print output
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printOutput(RealT t) const
    {
      RealT* yval  = N_VGetArrayPointer(yy_);
      RealT* ypval = N_VGetArrayPointer(yp_);

      std::cout << std::setprecision(5) << std::setw(7) << t << " ";
      for (IdxT i = 0; i < model_->size(); ++i)
      {
        std::cout << yval[i] << " ";
      }
      for (IdxT i = 0; i < model_->size(); ++i)
      {
        std::cout << ypval[i] << " ";
      }
      std::cout << "\n";
    }

    /**
     * @brief Special print
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printSpecial(RealT t, N_Vector y) const
    {
      RealT* yval = N_VGetArrayPointer(y);
      IdxT   N    = static_cast<IdxT>(N_VGetLength(y));
      std::cout << "{";
      std::cout << std::setprecision(5) << std::setw(7) << t;
      for (IdxT i = 0; i < N; ++i)
      {
        std::cout << ", " << yval[i];
      }
      std::cout << "},\n";
    }

    /**
     * @brief Print final stats
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printFinalStats() const
    {
      int retval = IDAPrintAllStats(solver_, stdout, SUN_OUTPUTFORMAT_TABLE);
      checkOutput(retval, "IDAPrintAllStats");
    }

    /**
     * @brief Accumulate another stats object into this one, allowing for stats to be kept
     *        across multiple simulations with IDA
     */
    IdaStats& IdaStats::operator+=(const IdaStats& other)
    {
      num_steps_                       += other.num_steps_;
      num_residual_evals_              += other.num_residual_evals_;
      num_linear_decompositions_       += other.num_linear_decompositions_;
      num_error_test_fails_            += other.num_error_test_fails_;
      num_nonlinear_iters_             += other.num_nonlinear_iters_;
      num_nonlinear_convergence_fails_ += other.num_nonlinear_convergence_fails_;

      return *this;
    }

    /**
     * @brief Generate a string containing all of the stats in a formatted report.
     *        All columns are aligned. To change the width of a column,
     *        modify `label_width` or `stat_width`.
     */
    std::string IdaStats::report() const
    {
      int               label_width = 30;
      int               stat_width  = 12;
      std::stringstream out;

      out << std::setw(label_width) << "Steps" << " : " << std::setw(stat_width) << num_residual_evals_ << '\n'
          << std::setw(label_width) << "Residual evals" << " : " << std::setw(stat_width) << num_linear_decompositions_ << '\n'
          << std::setw(label_width) << "Linear decompositions" << " : " << std::setw(stat_width) << num_linear_decompositions_ << '\n'
          << std::setw(label_width) << "Error test failures" << " : " << std::setw(stat_width) << num_error_test_fails_ << '\n'
          << std::setw(label_width) << "Nonlinear iterations" << " : " << std::setw(stat_width) << num_nonlinear_iters_ << '\n'
          << std::setw(label_width) << "Nonlinear convergence failures" << " : " << std::setw(stat_width) << num_nonlinear_convergence_fails_;

      return out.str();
    }

    /**
     * @brief Construct and return an `IdaStats` object containing the statistics of the current IDA workspace.
     *        Several statistics returned by IDA are ignored because they are about the current state of IDA,
     *        rather than about the simulation at large.
     */
    template <class ScalarT, typename IdxT>
    IdaStats Ida<ScalarT, IdxT>::getStats() const
    {
      IdaStats stats;

      // Dummies for ignoring stats
      int         dummy;
      sunrealtype dummy2;

      int retval = IDAGetIntegratorStats(solver_,
                                         &stats.num_steps_,
                                         &stats.num_residual_evals_,
                                         &stats.num_linear_decompositions_,
                                         &stats.num_error_test_fails_,
                                         &dummy,
                                         &dummy,
                                         &dummy2,
                                         &dummy2,
                                         &dummy2,
                                         &dummy2);
      checkOutput(retval, "IDAGetIntegratorStats");

      retval = IDAGetNonlinSolvStats(solver_, &stats.num_nonlinear_iters_, &stats.num_nonlinear_convergence_fails_);
      checkOutput(retval, "IDAGetNonlinSolvStats");

      return stats;
    }

    /**
     * @brief Check SUNDIALS allocation
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkAllocation(void* v, const char* functionName)
    {
      if (v == NULL)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n";
        throw SundialsException();
      }
    }

    /**
     * @brief Check SUNDIALS output
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkOutput(int retval, const char* functionName)
    {
      if (retval < 0)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n";
        throw SundialsException();
      }
    }

    /**
     * @brief Set fixed step size and tolerances for the nonlinear solver
     *
     * @param time_step The fixed step size to use or 0 for adaptive
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setFixedStep(ScalarT time_step)
    {
      time_step_ = time_step;
    }

    /**
     * @brief Set fixed step size and tolerances for the nonlinear solver for
     *        the backward simulation
     *
     * @param time_step The fixed step size to use or 0 for adaptive
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setBackwardFixedStep(ScalarT time_step)
    {
      backward_time_step_ = time_step;
    }

    /**
     * @brief Set the relative tolerance and optionally override the absolute
     *        tolerance
     *
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setTolerance(ScalarT rel_tol,
                                          ScalarT abs_tol_override)
    {
      rel_tol_          = rel_tol;
      abs_tol_override_ = abs_tol_override;
    }

    /**
     * @brief Set the relative tolerance and optionally override the absolute
     *        tolerance for the backward simulation
     *
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setBackwardTolerance(ScalarT rel_tol,
                                                  ScalarT abs_tol_override)
    {
      backward_rel_tol_          = rel_tol;
      backward_abs_tol_override_ = abs_tol_override;
    }

    /**
     * @brief Set the quadrature relative tolerance and optionally override the
     *        absolute tolerance
     *
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setQuadratureTolerance(ScalarT rel_tol,
                                                    ScalarT abs_tol_override)
    {
      quadrature_rel_tol_          = rel_tol;
      quadrature_abs_tol_override_ = abs_tol_override;
    }

    /**
     * @brief Set the quadrature relative tolerance and optionally override the
     *        absolute tolerance for the backward simulation
     *
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setBackwardQuadratureTolerance(ScalarT rel_tol,
                                                            ScalarT abs_tol_override)
    {
      backward_quadrature_rel_tol_          = rel_tol;
      backward_quadrature_abs_tol_override_ = abs_tol_override;
    }

    /**
     * @brief Set whether IDA suppresses local error tests on algebraic variables
     *
     * @param suppress If true, algebraic variables are excluded from IDA's
     *        local error test
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setSuppressAlgebraicErrors(bool suppress)
    {
      suppress_alg_ = suppress;
    }

    /**
     * @brief Set whether IDA suppresses local error tests on algebraic
     *        variables for the backward simulation
     *
     * @param suppress If true, algebraic variables are excluded from IDA's
     *        local error test
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setBackwardSuppressAlgebraicErrors(bool suppress)
    {
      backward_suppress_alg_ = suppress;
    }

    /**
     * @brief Set the maximum number of steps
     *
     * @param max_steps The maximum number of steps
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setMaxSteps(IdxT max_steps)
    {
      max_steps_ = max_steps;
    }

    /**
     * @brief Set the maximum number of steps for the backward simulation
     *
     * @param max_steps The maximum number of steps
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setBackwardMaxSteps(IdxT max_steps)
    {
      backward_max_steps_ = max_steps;
    }

    /**
     * @brief A helper function to set common IDA options
     *
     * @param mem The IDA memory (either forward or backward)
     * @param time_step The fixed step size to use
     * @param rel_tol The relative tolerance to use for the nonlinear solver
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance for the nonlinear solver rather than the
     *        model's default absolute tolerance
     * @param max_steps The maximum number of steps
     * @param suppress_alg If true, algebraic variables are excluded from IDA's
     *        local error test
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setIDAOptions(void*   mem,
                                           ScalarT time_step,
                                           ScalarT rel_tol,
                                           ScalarT abs_tol_override,
                                           IdxT    max_steps,
                                           bool    suppress_alg)
    {
      int retval = 0;
      retval     = IDASetMinStep(mem, time_step);
      checkOutput(retval, "IDASetMinStep");
      retval = IDASetMaxStep(mem, time_step);
      checkOutput(retval, "IDASetMaxStep");
      retval = IDASetMaxNumSteps(mem, static_cast<long int>(max_steps));
      checkOutput(retval, "IDASetMaxNumSteps");
      retval = IDASetSuppressAlg(mem, suppress_alg ? SUNTRUE : SUNFALSE);
      checkOutput(retval, "IDASetSuppressAlg");

      if (time_step == 0)
      {
        setTolerance(mem, rel_tol, abs_tol_override);
      }
      else
      {
        /* Since the starting procedure is first order, the maximum global order
         * of convergence is two */
        retval = IDASetMaxOrd(mem, 2);
        checkOutput(retval, "IDASetMaxOrd");

        /* Enable more nonlinear iterations because a failed nonlinear solve
         * causes a failed integration with fixed steps */
        static constexpr int FIXED_STEP_NONLIN_ITRS = 16;
        retval                                      = IDASetMaxNonlinIters(mem, FIXED_STEP_NONLIN_ITRS);
        checkOutput(retval, "IDASetMaxNonlinIters");

        // Set a large tolerance so the error test will never fail
        static constexpr RealT FIXED_STEP_TOL_FAC = 1e100;
        setTolerance(mem,
                     FIXED_STEP_TOL_FAC * rel_tol,
                     FIXED_STEP_TOL_FAC * abs_tol_override,
                     FIXED_STEP_TOL_FAC);

        /* We want the nonlinear solver tolerance to be ~rel_tol, but the with
         * the large tolerances set above, we need to choose this tolerance to
         * "undo" the fac scaling. */
        retval = IDASetNonlinConvCoef(mem, 1 / FIXED_STEP_TOL_FAC);
        checkOutput(retval, "IDASetNonlinConvCoef");
      }
    }

    /**
     * @brief A helper function to set the relative tolerance and optionally
     *        override the absolute tolerance
     *
     * @param mem The IDA memory (either forward or backward)
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @param abs_tol_fac A factor to apply to the absolute tolerance if not
     *        overridden
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setTolerance(void*   mem,
                                          ScalarT rel_tol,
                                          ScalarT abs_tol_override,
                                          ScalarT abs_tol_fac)
    {
      int retval = 0;

      if (abs_tol_override > 0)
      {
        retval = IDASStolerances(mem, rel_tol, abs_tol_override);
        checkOutput(retval, "IDASStolerances");
        return;
      }

      N_Vector abs_tol_vec = N_VClone(yy_);
      checkAllocation((void*) abs_tol_vec, "N_VClone");
      model_->setAbsoluteTolerance(rel_tol);
      copyVec(model_->absoluteTolerance(), abs_tol_vec);
      N_VScale(abs_tol_fac, abs_tol_vec, abs_tol_vec);

      retval = IDASVtolerances(mem, rel_tol, abs_tol_vec);
      checkOutput(retval, "IDASVtolerances");

      N_VDestroy(abs_tol_vec);
    }

    /**
     * @brief A helper function to set the quadrature relative tolerance and
     *        optionally override the absolute tolerance
     *
     * @param mem The IDA memory (either forward or backward)
     * @param rel_tol The relative tolerance to use
     * @param abs_tol_override If positive, this value will be used as the
     *        absolute tolerance rather than the model's default absolute
     *        tolerance
     * @param abs_tol_fac A factor to apply to the absolute tolerance if not
     *        overridden
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     */
    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::setQuadratureTolerance(void*   mem,
                                                    ScalarT rel_tol,
                                                    ScalarT abs_tol_override)
    {
      int retval = 0;

      if (abs_tol_override > 0)
      {
        retval = IDAQuadSStolerances(mem, rel_tol, abs_tol_override);
        checkOutput(retval, "IDAQuadSStolerances");
        return;
      }

      N_Vector abs_tol_vec = N_VClone(yy_);
      checkAllocation((void*) abs_tol_vec, "N_VClone");
      model_->setAbsoluteTolerance(rel_tol);
      copyVec(model_->absoluteTolerance(), abs_tol_vec);

      retval = IDAQuadSVtolerances(mem, rel_tol, abs_tol_vec);
      checkOutput(retval, "IDAQuadSVtolerances");

      N_VDestroy(abs_tol_vec);
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class Ida<sunrealtype, long int>;
    template class Ida<sunrealtype, int>;
    template class Ida<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
