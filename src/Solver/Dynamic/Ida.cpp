
#include "Ida.hpp"

#include <iomanip>
#include <iostream>

#include <idas/idas.h>
#include <idas/idas_ls.h>

#include "Model/Evaluator.hpp"

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
      const real_type t0 = 0.0;

      // Allocate and initialize IDA workspace
      retval = IDAInit(solver_, this->Residual, t0, yy_, yp_);
      checkOutput(retval, "IDAInit");

      // Set pointer to model data
      retval = IDASetUserData(solver_, model_);
      checkOutput(retval, "IDASetUserData");

      // Set tolerances
      real_type rel_tol;
      real_type abs_tol;

      model_->setTolerances(rel_tol, abs_tol); ///< \todo Function name should be "getTolerances"!
      retval = IDASStolerances(solver_, rel_tol, abs_tol);
      checkOutput(retval, "IDASStolerances");

      IdxT msa;
      model_->setMaxSteps(msa);

      /// \todo Need to set max number of steps based on user input!
      retval = IDASetMaxNumSteps(solver_, static_cast<long>(msa));
      checkOutput(retval, "IDASetMaxNumSteps");

      // Tag differential variables
      std::vector<bool>& tag = model_->tag();
      if (static_cast<IdxT>(tag.size()) == model_->size())
      {
        tag_ = N_VClone(yy_);
        checkAllocation((void*) tag_, "N_VClone");
        model_->tagDifferentiable();
        copyVec(tag, tag_);

        retval = IDASetId(solver_, tag_);
        checkOutput(retval, "IDASetId");
        retval = IDASetSuppressAlg(solver_, SUNTRUE);
        checkOutput(retval, "IDASetSuppressAlg");
      }

      // Set up linear solver
      return this->configureLinearSolver();
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolver()
    {
      int retval = 0;
      if (model_->hasJacobian())
      {
        JacobianMat_ = SUNSparseMatrix(static_cast<sunindextype>(model_->size()),
                                       static_cast<sunindextype>(model_->size()),
                                       static_cast<sunindextype>(model_->size() * model_->size()),
                                       CSR_MAT,
                                       context_);
        checkAllocation((void*) JacobianMat_, "SUNSparseMatrix");

        linearSolver_ = SUNLinSol_KLU(yy_, JacobianMat_, context_);
        checkAllocation((void*) linearSolver_, "SUNLinSol_KLU");

        retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
        checkOutput(retval, "IDASetLinearSolver");

        retval = IDASetJacFn(solver_, this->Jac);
        checkOutput(retval, "IDASetJacFn");
      }
      else
      {
        JacobianMat_ = SUNDenseMatrix(static_cast<sunindextype>(model_->size()),
                                      static_cast<sunindextype>(model_->size()),
                                      context_);
        checkAllocation((void*) JacobianMat_, "SUNDenseMatrix");

        linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
        checkAllocation((void*) linearSolver_, "SUNLinSol_Dense");

        retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
        checkOutput(retval, "IDASetLinearSolver");
      }

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::getDefaultInitialCondition()
    {
      model_->initialize();

      copyVec(model_->y(), yy_);
      copyVec(model_->yp(), yp_);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::setIntegrationTime(real_type t_init, real_type t_final, int nout)
    {
      t_init_  = t_init;
      t_final_ = t_final;
      nout_    = nout;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeSimulation(real_type t0, bool findConsistent)
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

        retval = IDACalcIC(solver_, initType, 0.1);
        checkOutput(retval, "IDACalcIC");

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
    int Ida<ScalarT, IdxT>::runSimulation(real_type tf, int nout, const std::optional<std::function<void(real_type)>> step_callback)
    {
      int       retval = 0;
      int       iout   = 0;
      real_type tret;
      real_type dt   = (tf - t_init_) / static_cast<real_type>(nout);
      real_type tout = t_init_ + dt;

      // In loop, call IDASolve, print results, and test for error.
      //  Break out of loop when NOUT preset output times have been reached.
      // printOutput(0.0);
      while (nout > iout)
      {
        retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
        checkOutput(retval, "IDASolve");

        if (step_callback.has_value())
        {
          // The callback may try to observe upated values in the model, so we
          // should update them here (At this point, the model's values are one
          // internal integrator step out of date)
          model_->updateTime(tret, 0.0);
          copyVec(yy_, model_->y());
          copyVec(yp_, model_->yp());

          (*step_callback)(tret);
        }

        if (retval == IDA_SUCCESS)
        {
          ++iout;
          tout += dt;
        }
      }

      // Final copy out. No guarantee last residual evaluation is final step.
      model_->updateTime(tf, 0.0);
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());

      // std::cout << "\n";
      return retval;
    }

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
      real_type rel_tol, abs_tol;
      model_->setTolerances(rel_tol, abs_tol);

      // Set tolerances for quadrature stricter than for integration
      retval = IDAQuadSStolerances(solver_, rel_tol * 0.1, abs_tol * 0.1);
      checkOutput(retval, "IDAQuadSStolerances");

      // Include quadrature in eror checking
      retval = IDASetQuadErrCon(solver_, SUNTRUE);
      checkOutput(retval, "IDASetQuadErrCon");

      return retval;
    }

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulationQuadrature(real_type tf, int nout)
    {
      int       retval = 0;
      real_type tret;

      // std::cout << "Forward integration for initial value problem ... \n";

      real_type dt   = tf / static_cast<real_type>(nout);
      real_type tout = dt;
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
      model_->updateTime(tf, 0.0);
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteQuadrature()
    {
      IDAQuadFree(solver_);
      N_VDestroy(q_);

      return 0;
    }

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeAdjoint(IdxT steps)
    {
      int retval = 0;

      // Create adjoint workspace
      retval = IDAAdjInit(solver_, static_cast<long>(steps), IDA_HERMITE);
      checkOutput(retval, "IDAAdjInit");

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeBackwardSimulation(real_type tf)
    {
      int       retval = 0;
      real_type rel_tol;
      real_type abs_tol;

      model_->initializeAdjoint();

      copyVec(model_->yB(), yyB_);
      copyVec(model_->ypB(), ypB_);
      N_VConst(0.0, qB_);

      retval = IDACreateB(solver_, &backwardID_);
      checkOutput(retval, "IDACreateB");

      // IDAInitB must be called after forward simulation run.
      retval = IDAInitB(solver_, backwardID_, this->adjointResidual, tf, yyB_, ypB_);
      checkOutput(retval, "IDAInitB");

      model_->setTolerances(rel_tol, abs_tol);
      retval = IDASStolerancesB(solver_, backwardID_, rel_tol, abs_tol);
      checkOutput(retval, "IDASStolerancesB");

      retval = IDASetUserDataB(solver_, backwardID_, model_);
      checkOutput(retval, "IDASetUserDataB");

      /// \todo Need to set max number of steps based on user input!
      retval = IDASetMaxNumStepsB(solver_, backwardID_, 2000);
      checkOutput(retval, "IDASetMaxNumSteps");

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

      // retval = IDAQuadSStolerancesB(solver_, backwardID_, rel_tol*1.1, abs_tol*1.1);
      retval = IDAQuadSStolerancesB(solver_, backwardID_, rel_tol * 0.1, abs_tol * 0.1);
      checkOutput(retval, "IDAQuadSStolerancesB");

      // Include quadratures in error control
      retval = IDASetQuadErrConB(solver_, backwardID_, SUNTRUE);
      checkOutput(retval, "IDASetQuadErrConB");

      return retval;
    }

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runForwardSimulation(real_type tf, int nout)
    {
      int       retval = 0;
      int       ncheck;
      real_type time;

      // std::cout << "Forward integration for adjoint analysis ... \n";

      real_type dt   = tf / static_cast<real_type>(nout);
      real_type tout = dt;
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
      model_->updateTime(tf, 0.0);
      copyVec(yy_, model_->y());
      copyVec(yp_, model_->yp());

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runBackwardSimulation(real_type t_init)
    {
      int       retval = 0;
      long int  nstB;
      real_type time;

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteAdjoint()
    {
      IDAAdjFree(solver_);
      return 0;
    }

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Residual(real_type tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      model->updateTime(tres, 0.0);
      copyVec(yy, model->y());
      copyVec(yp, model->yp());

      model->evaluateResidual();
      const std::vector<ScalarT>& f = model->getResidual();
      copyVec(f, rr);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Jac(real_type t, real_type cj, N_Vector yy, N_Vector yp, N_Vector, SUNMatrix J, void* user_data, N_Vector, N_Vector, N_Vector)
    {

      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      model->updateTime(t, cj);
      copyVec(yy, model->y());
      copyVec(yp, model->yp());

      model->evaluateJacobian();
      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> Jac = model->getJacobian();

      // Get reference to the jacobian entries
      std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> tpm = Jac.getEntries();
      const auto [r, c, val]                                                        = tpm;

      // get the CSR row pointers from COO matrix
      std::vector<IdxT> csrrowdata = Jac.getCSRRowData();

      SUNMatZero(J);

      // Set row pointers
      sunindextype* rowptrs = SUNSparseMatrix_IndexPointers(J);
      std::copy(csrrowdata.cbegin(), csrrowdata.cend(), rowptrs);

      sunindextype* colvals = SUNSparseMatrix_IndexValues(J);
      real_type*    data    = SUNSparseMatrix_Data(J);
      // Copy data from model jac to sundials
      std::copy(c.cbegin(), c.cend(), colvals);
      std::copy(val.cbegin(), val.cend(), data);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Integrand(real_type tt, N_Vector yy, N_Vector yp, N_Vector rhsQ, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      model->updateTime(tt, 0.0);
      copyVec(yy, model->y());
      copyVec(yp, model->yp());

      model->evaluateIntegrand();
      const std::vector<ScalarT>& g = model->getIntegrand();
      copyVec(g, rhsQ);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointResidual(real_type tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

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

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointIntegrand(real_type tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rhsQB, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

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

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const N_Vector x, std::vector<ScalarT>& y)
    {
      const ScalarT* xdata = N_VGetArrayPointer(x);
      std::copy_n(xdata, y.size(), y.begin());
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const std::vector<ScalarT>& x, N_Vector y)
    {
      ScalarT* ydata = N_VGetArrayPointer(y);
      std::copy(x.cbegin(), x.cend(), ydata);
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const std::vector<bool>& x, N_Vector y)
    {
      ScalarT* ydata = N_VGetArrayPointer(y);
      std::copy(x.cbegin(), x.cend(), ydata);
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printOutput(real_type t)
    {
      real_type* yval  = N_VGetArrayPointer(yy_);
      real_type* ypval = N_VGetArrayPointer(yp_);

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

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printSpecial(real_type t, N_Vector y)
    {
      real_type* yval = N_VGetArrayPointer(y);
      IdxT       N    = static_cast<IdxT>(N_VGetLength(y));
      std::cout << "{";
      std::cout << std::setprecision(5) << std::setw(7) << t;
      for (IdxT i = 0; i < N; ++i)
      {
        std::cout << ", " << yval[i];
      }
      std::cout << "},\n";
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printFinalStats()
    {
      int retval = IDAPrintAllStats(solver_, stdout, SUN_OUTPUTFORMAT_TABLE);
      checkOutput(retval, "IDAPrintAllStats");
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkAllocation(void* v, const char* functionName)
    {
      if (v == NULL)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n";
        throw SundialsException();
      }
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkOutput(int retval, const char* functionName)
    {
      if (retval < 0)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n";
        throw SundialsException();
      }
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class Ida<sunrealtype, long int>;
    template class Ida<sunrealtype, int>;
    template class Ida<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
