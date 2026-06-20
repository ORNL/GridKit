#include "Arkode.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

#include <sunlinsol/sunlinsol_klu.h>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    Arkode<ScalarT, IdxT>::Arkode(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : DynamicSolver<ScalarT, IdxT>(model)
    {
      int retval = SUNContext_Create(SUN_COMM_NULL, &context_);
      checkOutput(retval, "SUNContext");
    }

    template <class ScalarT, typename IdxT>
    Arkode<ScalarT, IdxT>::~Arkode()
    {
      deleteSimulation();
      SUNContext_Free(&context_);
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::configureSimulation()
    {
      yy_ = N_VNew_Serial(static_cast<sunindextype>(model_->size()), context_);
      checkAllocation(static_cast<void*>(yy_), "N_VNew_Serial");

      getDefaultInitialCondition();

      yy0_ = N_VClone(yy_);
      checkAllocation(static_cast<void*>(yy0_), "N_VClone");
      N_VScale(1.0, yy_, yy0_);

      solver_ = createStepper(RealT{0}, yy_);
      checkAllocation(solver_, "createStepper");

      int retval = ARKodeSetUserData(solver_, this);
      checkOutput(retval, "ARKodeSetUserData");

      setArkodeOptions(solver_, time_step_, rel_tol_, abs_tol_override_, max_steps_);
      return configureLinearSolver();
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::configureLinearSolver()
    {
      if (!model_->hasJacobian() || model_->getCsrJacobian() == nullptr)
      {
        Log::error() << "Arkode requires an analytic sparse Jacobian.\n";
        throw SundialsException();
      }

      if (requiresIdentityMass() && !massMatrixIsIdentity(t_init_))
      {
        Log::error() << "ErkStep requires an identity mass matrix.\n";
        throw SundialsException();
      }

      int retval = 0;
      if (usesImplicitSolver())
        retval = setupJacobianSolver();
      if (usesMassMatrix())
        retval = setupMassSolver();
      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::setupJacobianSolver()
    {
      const sunindextype n   = static_cast<sunindextype>(model_->size());
      const sunindextype nnz = static_cast<sunindextype>(model_->getCsrJacobian()->getNnz());

      JacobianMat_ = SUNSparseMatrix(n, n, nnz, CSR_MAT, context_);
      checkAllocation(static_cast<void*>(JacobianMat_), "SUNSparseMatrix");

      linearSolver_ = SUNLinSol_KLU(yy_, JacobianMat_, context_);
      checkAllocation(static_cast<void*>(linearSolver_), "SUNLinSol_KLU");

      int retval = ARKodeSetLinearSolver(solver_, linearSolver_, JacobianMat_);
      checkOutput(retval, "ARKodeSetLinearSolver");

      retval = ARKodeSetJacFn(solver_, this->Jac);
      checkOutput(retval, "ARKodeSetJacFn");

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::setupMassSolver()
    {
      const sunindextype n   = static_cast<sunindextype>(model_->size());
      const sunindextype nnz = static_cast<sunindextype>(model_->getCsrJacobian()->getNnz());

      MassMat_ = SUNSparseMatrix(n, n, nnz, CSR_MAT, context_);
      checkAllocation(static_cast<void*>(MassMat_), "SUNSparseMatrix");

      massLinearSolver_ = SUNLinSol_KLU(yy_, MassMat_, context_);
      checkAllocation(static_cast<void*>(massLinearSolver_), "SUNLinSol_KLU");

      int retval = ARKodeSetMassLinearSolver(solver_, massLinearSolver_, MassMat_, SUNFALSE);
      checkOutput(retval, "ARKodeSetMassLinearSolver");

      retval = ARKodeSetMassFn(solver_, this->Mass);
      checkOutput(retval, "ARKodeSetMassFn");

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::getDefaultInitialCondition()
    {
      model_->initialize();
      copyVec(model_->y(), yy_);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::setIntegrationTime(RealT t_init, RealT t_final, int nout)
    {
      t_init_  = t_init;
      t_final_ = t_final;
      nout_    = nout;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::initializeSimulation(RealT t0)
    {
      t_init_    = t0;
      int retval = ARKodeReset(solver_, t0, yy_);
      checkOutput(retval, "ARKodeReset");
      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::runSimulation(RealT tf, int nout, std::optional<std::function<void(RealT)>> step_callback)
    {
      int   retval = 0;
      int   iout   = 0;
      RealT tret   = t_init_;
      RealT dt     = (tf - t_init_) / static_cast<RealT>(nout);
      RealT tout   = t_init_ + dt;

      while (nout > iout)
      {
        retval = ARKodeEvolve(solver_, tout, yy_, &tret, ARK_NORMAL);
        checkOutput(retval, "ARKodeEvolve");

        if (step_callback.has_value() || model_->monitoring())
        {
          copyVec(yy_, model_->y());
          model_->updateTime(tret, 0.0);

          if (model_->monitoring())
            model_->printMonitoredVariables();
          if (step_callback.has_value())
            (*step_callback)(tret);
        }

        if (retval >= 0)
        {
          ++iout;
          tout += dt;
        }
      }

      copyVec(yy_, model_->y());
      model_->updateTime(tf, 0.0);
      t_final_ = tf;
      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::deleteSimulation()
    {
      if (yy_)
      {
        N_VDestroy(yy_);
        yy_ = nullptr;
      }
      if (yy0_)
      {
        N_VDestroy(yy0_);
        yy0_ = nullptr;
      }
      if (linearSolver_)
      {
        SUNLinSolFree(linearSolver_);
        linearSolver_ = nullptr;
      }
      if (massLinearSolver_)
      {
        SUNLinSolFree(massLinearSolver_);
        massLinearSolver_ = nullptr;
      }
      if (JacobianMat_)
      {
        SUNMatDestroy(JacobianMat_);
        JacobianMat_ = nullptr;
      }
      if (MassMat_)
      {
        SUNMatDestroy(MassMat_);
        MassMat_ = nullptr;
      }
      if (solver_)
        ARKodeFree(&solver_);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::Rhs(RealT t, N_Vector yy, N_Vector ydot, void* user_data)
    {
      auto* solver = static_cast<Arkode<ScalarT, IdxT>*>(user_data);
      auto* model  = solver->model_;

      copyVec(yy, model->y());
      std::vector<ScalarT>& yp = model->yp();
      std::fill(yp.begin(), yp.end(), ScalarT{0});

      model->updateTime(t, 0.0);
      model->evaluateResidual();
      copyVec(model->getResidual(), ydot);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::Jac(RealT t, N_Vector yy, N_Vector, SUNMatrix J, void* user_data, N_Vector, N_Vector, N_Vector)
    {
      auto* solver = static_cast<Arkode<ScalarT, IdxT>*>(user_data);
      auto* model  = solver->model_;

      copyVec(yy, model->y());
      std::vector<ScalarT>& yp = model->yp();
      std::fill(yp.begin(), yp.end(), ScalarT{0});

      model->updateTime(t, 0.0);
      model->evaluateJacobian();
      copyCsrToSun(model->getCsrJacobian(), J);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::Mass(RealT t, SUNMatrix M, void* user_data, N_Vector, N_Vector, N_Vector)
    {
      auto* solver = static_cast<Arkode<ScalarT, IdxT>*>(user_data);
      return solver->recoverMassMatrix(t, M);
    }

    template <class ScalarT, typename IdxT>
    int Arkode<ScalarT, IdxT>::recoverMassMatrix(RealT t, SUNMatrix M)
    {
      model_->updateTime(t, 1.0);
      model_->evaluateJacobian();
      CsrMatrixT*        J   = model_->getCsrJacobian();
      const auto         nnz = static_cast<size_t>(J->getNnz());
      std::vector<RealT> j1(J->getValues(), J->getValues() + nnz);

      model_->updateTime(t, 0.0);
      model_->evaluateJacobian();
      J = model_->getCsrJacobian();

      SUNMatZero(M);
      sunindextype* sun_row_ptrs = SUNSparseMatrix_IndexPointers(M);
      sunindextype* sun_cols     = SUNSparseMatrix_IndexValues(M);
      sunrealtype*  sun_vals     = SUNSparseMatrix_Data(M);

      const auto   n       = static_cast<size_t>(J->getNumRows());
      const IdxT*  row_ptr = J->getRowData();
      const IdxT*  cols    = J->getColData();
      const RealT* j0      = J->getValues();

      for (size_t i = 0; i < n + 1; ++i)
        sun_row_ptrs[i] = static_cast<sunindextype>(row_ptr[i]);
      for (size_t k = 0; k < nnz; ++k)
      {
        sun_cols[k] = static_cast<sunindextype>(cols[k]);
        sun_vals[k] = static_cast<sunrealtype>(j0[k] - j1[k]);
      }

      return 0;
    }

    template <class ScalarT, typename IdxT>
    bool Arkode<ScalarT, IdxT>::massMatrixIsIdentity(RealT t)
    {
      model_->updateTime(t, 1.0);
      model_->evaluateJacobian();
      CsrMatrixT*        J   = model_->getCsrJacobian();
      const auto         nnz = static_cast<size_t>(J->getNnz());
      std::vector<RealT> j1(J->getValues(), J->getValues() + nnz);

      model_->updateTime(t, 0.0);
      model_->evaluateJacobian();
      J = model_->getCsrJacobian();

      const auto   n       = static_cast<size_t>(J->getNumRows());
      const IdxT*  row_ptr = J->getRowData();
      const IdxT*  cols    = J->getColData();
      const RealT* j0      = J->getValues();

      const RealT       tol = static_cast<RealT>(100) * std::numeric_limits<RealT>::epsilon();
      std::vector<bool> diagonal_seen(n, false);

      for (size_t row = 0; row < n; ++row)
      {
        for (IdxT k = row_ptr[row]; k < row_ptr[row + 1]; ++k)
        {
          const auto  col   = static_cast<size_t>(cols[k]);
          const RealT value = j0[k] - j1[static_cast<size_t>(k)];
          const RealT ref   = (col == row) ? RealT{1} : RealT{0};
          const RealT scale = RealT{1} + std::abs(ref);
          if (std::abs(value - ref) > tol * scale)
            return false;
          if (col == row)
            diagonal_seen[row] = true;
        }
      }

      return std::all_of(diagonal_seen.begin(), diagonal_seen.end(), [](bool seen)
                         { return seen; });
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setArkodeOptions(void* mem, RealT time_step, RealT rel_tol, RealT abs_tol_override, IdxT max_steps)
    {
      int retval = ARKodeSetMaxNumSteps(mem, static_cast<long int>(max_steps));
      checkOutput(retval, "ARKodeSetMaxNumSteps");

      if (order_ > 0)
      {
        retval = ARKodeSetOrder(mem, order_);
        checkOutput(retval, "ARKodeSetOrder");
      }

      setTolerance(mem, rel_tol, abs_tol_override);

      if (time_step != RealT{0})
      {
        retval = ARKodeSetFixedStep(mem, time_step);
        checkOutput(retval, "ARKodeSetFixedStep");
      }
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setTolerance(void* mem, RealT rel_tol, RealT abs_tol_override)
    {
      if (abs_tol_override > RealT{0})
      {
        int retval = ARKodeSStolerances(mem, rel_tol, abs_tol_override);
        checkOutput(retval, "ARKodeSStolerances");
        return;
      }

      N_Vector abs_tol_vec = N_VClone(yy_);
      checkAllocation(static_cast<void*>(abs_tol_vec), "N_VClone");
      model_->setAbsoluteTolerance(rel_tol);
      copyVec(model_->absoluteTolerance(), abs_tol_vec);

      int retval = ARKodeSVtolerances(mem, rel_tol, abs_tol_vec);
      checkOutput(retval, "ARKodeSVtolerances");

      N_VDestroy(abs_tol_vec);
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setFixedStep(ScalarT time_step)
    {
      time_step_ = static_cast<RealT>(time_step);
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setOrder(int order)
    {
      order_ = order;
      if (solver_ && order_ > 0)
      {
        int retval = ARKodeSetOrder(solver_, order_);
        checkOutput(retval, "ARKodeSetOrder");
      }
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setTolerance(ScalarT rel_tol, ScalarT abs_tol_override)
    {
      rel_tol_          = static_cast<RealT>(rel_tol);
      abs_tol_override_ = static_cast<RealT>(abs_tol_override);
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::setMaxSteps(IdxT max_steps)
    {
      max_steps_ = max_steps;
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::copyVec(const N_Vector x, std::vector<ScalarT>& y)
    {
      const auto xsize = static_cast<size_t>(N_VGetLength(x));
      if (xsize != y.size())
      {
        std::cerr << "\nN_Vector size (" << xsize << ") does not match std::vector size ("
                  << y.size() << ").\n\n";
        throw SundialsException();
      }

      const sunrealtype* xdata = N_VGetArrayPointer(x);
      for (size_t i = 0; i < y.size(); ++i)
        y[i] = static_cast<ScalarT>(xdata[i]);
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::copyVec(const std::vector<ScalarT>& x, N_Vector y)
    {
      const auto ysize = static_cast<size_t>(N_VGetLength(y));
      if (x.size() != ysize)
      {
        std::cerr << "\nstd::vector size (" << x.size() << ") does not match N_Vector size ("
                  << ysize << ").\n\n";
        throw SundialsException();
      }

      sunrealtype* ydata = N_VGetArrayPointer(y);
      for (size_t i = 0; i < x.size(); ++i)
        ydata[i] = static_cast<sunrealtype>(x[i]);
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::copyCsrToSun(CsrMatrixT* src, SUNMatrix dst)
    {
      SUNMatZero(dst);
      sunindextype* sun_row_ptrs = SUNSparseMatrix_IndexPointers(dst);
      sunindextype* sun_cols     = SUNSparseMatrix_IndexValues(dst);
      sunrealtype*  sun_vals     = SUNSparseMatrix_Data(dst);

      const auto   n       = static_cast<size_t>(src->getNumRows());
      const auto   nnz     = static_cast<size_t>(src->getNnz());
      const IdxT*  row_ptr = src->getRowData();
      const IdxT*  cols    = src->getColData();
      const RealT* vals    = src->getValues();

      for (size_t i = 0; i < n + 1; ++i)
        sun_row_ptrs[i] = static_cast<sunindextype>(row_ptr[i]);
      for (size_t k = 0; k < nnz; ++k)
      {
        sun_cols[k] = static_cast<sunindextype>(cols[k]);
        sun_vals[k] = static_cast<sunrealtype>(vals[k]);
      }
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::checkAllocation(void* v, const char* functionName)
    {
      if (v == nullptr)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n";
        throw SundialsException();
      }
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::checkOutput(int retval, const char* functionName)
    {
      if (retval < 0)
      {
        std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n";
        throw SundialsException();
      }
    }

    template <class ScalarT, typename IdxT>
    void Arkode<ScalarT, IdxT>::printFinalStats() const
    {
      int retval = ARKodePrintAllStats(solver_, stdout, SUN_OUTPUTFORMAT_TABLE);
      checkOutput(retval, "ARKodePrintAllStats");
    }

    template class Arkode<sunrealtype, long int>;
    template class Arkode<sunrealtype, int>;
    template class Arkode<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
