/**
 * @file Kinsol.cpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains definition of interface to KINSOL nonlinear solver from
 * SUNDIALS library.
 */

#include <algorithm>
#include <format>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

#include <kinsol/kinsol.h>             // access to KINSOL func., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix

#include <GridKit/Model/Evaluator.hpp>

#include "Kinsol.hpp"

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::Kinsol(
      GridKit::Model::Evaluator<ScalarT, IdxT>* model
    )
      : Kinsol(model->size())
    {
      model_ = model;
      residual_ =
        [model](std::span<const ScalarT> values, std::span<ScalarT> residual)
      {
        std::copy(values.begin(), values.end(), model->y().getData());
        model->y().setDataUpdated();
        const int status = model->evaluateResidual();
        std::copy_n(
          model->getResidual().getData(), residual.size(), residual.begin()
        );
        return status;
      };
    }

    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::Kinsol(IdxT size)
      : SteadyStateSolver<ScalarT, IdxT>(nullptr), size_(size)
    {
      int retval = 0;

      // Create the SUNDIALS context that all SUNDIALS objects require
      retval = SUNContext_Create(SUN_COMM_NULL, &context_);
      checkOutput(retval, "SUNContext");

      solver_ = KINCreate(context_);
    }

    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::~Kinsol()
    {
      deleteSimulation();
      SUNContext_Free(&context_);
      solver_ = nullptr;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::configureSimulation()
    {
      if (model_)
      {
        size_ = model_->size();
      }
      if (size_ <= 0)
      {
        throw std::invalid_argument("KINSOL requires a positive problem size");
      }

      int retval = 0;

      // Allocate solution vectors
      yy_ = N_VNew_Serial(static_cast<sunindextype>(size_), context_);
      checkAllocation(context_, (void*) yy_, "N_VNew_Serial");

      // Allocate scaling vector
      scale_ = N_VClone(yy_);
      checkAllocation(context_, (void*) scale_, "N_VClone");

      // Create vectors to store restart initial condition
      yy0_ = N_VClone(yy_);
      checkAllocation(context_, (void*) yy0_, "N_VClone");

      // Allocate and initialize KIN workspace
      retval = KINInit(solver_, this->Residual, yy_);
      checkOutput(retval, "KINInit");

      // Set the callback context
      retval = KINSetUserData(solver_, this);
      checkOutput(retval, "KINSetUserData");

      // Set up linear solver
      return this->configureLinearSolver();
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::configureLinearSolver()
    {
      int retval = 0;

      // Set up linear solver
      JacobianMat_ = SUNDenseMatrix(
        static_cast<sunindextype>(size_),
        static_cast<sunindextype>(size_),
        context_
      );
      checkAllocation(context_, (void*) JacobianMat_, "SUNDenseMatrix");

      linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
      checkAllocation(context_, (void*) linearSolver_, "SUNLinSol_Dense");

      retval = KINSetLinearSolver(solver_, linearSolver_, JacobianMat_);
      checkOutput(retval, "KINSetLinearSolver");

      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::getDefaultInitialCondition()
    {
      if (!model_)
      {
        throw std::logic_error(
          "KINSOL requires an explicit initial guess when used with a callback"
        );
      }

      model_->initialize();

      copyVec(model_->y(), yy_);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::runSimulation()
    {
      N_VConst(1.0, scale_);
      return solve();
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::runSimulation(
      std::span<ScalarT> values, ResidualFunction residual, bool line_search
    )
    {
      if (values.size() != static_cast<std::size_t>(size_) || !residual)
      {
        throw std::invalid_argument(
          "KINSOL requires a matching initial guess and a residual callback"
        );
      }

      // ensure the model callback is restored on scope exit
      struct Restore
      {
        ResidualFunction& active;
        ResidualFunction  previous;

        ~Restore()
        {
          active.swap(previous);
        }
      } restore{residual_, std::move(residual_)};

      residual_ = std::move(residual);
      std::copy(values.begin(), values.end(), N_VGetArrayPointer(yy_));
      N_VConst(1.0, scale_);
      const int status = solve(line_search);
      std::copy_n(N_VGetArrayPointer(yy_), values.size(), values.begin());
      return status;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::solve(bool line_search)
    {
      if (!residual_)
      {
        throw std::logic_error("KINSOL requires a residual callback");
      }

      callback_error_  = nullptr;
      const int retval = KINSol(
        solver_, yy_, line_search ? KIN_LINESEARCH : KIN_NONE, scale_, scale_
      );
      if (callback_error_)
      {
        std::rethrow_exception(callback_error_);
      }
      checkOutput(retval, "KINSol");
      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::deleteSimulation()
    {
      KINFree(&solver_);
      N_VDestroy(this->yy_);
      N_VDestroy(this->yy0_);
      N_VDestroy(this->scale_);
      SUNMatDestroy(this->JacobianMat_);
      SUNLinSolFree_Dense(this->linearSolver_);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::Residual(
      N_Vector yy, N_Vector rr, void* user_data
    )
    {
      auto& self = *static_cast<Kinsol*>(user_data);
      try
      {
        return self.residual_(
          {N_VGetArrayPointer(yy), static_cast<std::size_t>(self.size_)},
          {N_VGetArrayPointer(rr), static_cast<std::size_t>(self.size_)}
        );
      }
      catch (...)
      {
        self.callback_error_ = std::current_exception();
        return -1;
      }
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const N_Vector x, VectorT& y)
    {
      const ScalarT* xdata = N_VGetArrayPointer(x);
      std::copy_n(xdata, static_cast<size_t>(y.getSize()), y.getData());
      y.setDataUpdated();
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const VectorT& x, N_Vector y)
    {
      const auto* xdata = x.getData();
      auto*       ydata = N_VGetArrayPointer(y);
      std::copy_n(xdata, static_cast<size_t>(x.getSize()), ydata);
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printOutput() const
    {
      sunrealtype* yval = N_VGetArrayPointer(yy_);

      std::cout << std::setprecision(5) << std::setw(7);
      for (IdxT i = 0; i < size_; ++i)
      {
        std::cout << yval[i] << " ";
      }
      std::cout << "\n";
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printSpecial(sunrealtype t, N_Vector y) const
    {
      sunrealtype* yval = N_VGetArrayPointer_Serial(y);
      IdxT         N    = static_cast<IdxT>(N_VGetLength_Serial(y));
      std::cout << "{";
      std::cout << std::setprecision(5) << std::setw(7) << t;
      for (IdxT i = 0; i < N; ++i)
      {
        std::cout << ", " << yval[i];
      }
      std::cout << "},\n";
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printFinalStats() const
    {
      int retval = KINPrintAllStats(solver_, stdout, SUN_OUTPUTFORMAT_TABLE);
      checkOutput(retval, "KINPrintAllStats");
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::checkAllocation(
      SUNContext c, void* v, const char* functionName
    )
    {
      if (v == NULL)
      {
        // NOTE: this will only work if called immediately after the function
        auto errorCode = SUNContext_GetLastError(c);
        Logger::error() << std::format(
          "Function {} failed: {}\n", functionName, SUNGetErrMsg(errorCode)
        );
        throw SundialsException();
      }
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::checkOutput(
      int retval, const char* functionName
    )
    {
      if (retval < 0)
      {
        Logger::error() << std::format(
          "Function {} failed: {}\n", functionName, SUNGetErrMsg(retval)
        );
        throw SundialsException();
      }
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::setTolerance(ScalarT tol)
    {
      setTolerance(tol, tol);
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::setTolerance(
      ScalarT residual_tol, ScalarT step_tol
    )
    {
      int retval = KINSetFuncNormTol(solver_, residual_tol);
      checkOutput(retval, "KINSetFuncNormTol");

      retval = KINSetScaledStepTol(solver_, step_tol);
      checkOutput(retval, "KINSetScaledStepTol");
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::setFunctionError(RealT error)
    {
      checkOutput(KINSetRelErrFunc(solver_, error), "KINSetRelErrFunc");
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::setJacobianUpdateInterval(long int iterations)
    {
      checkOutput(
        KINSetMaxSetupCalls(solver_, iterations), "KINSetMaxSetupCalls"
      );
    }

    // Compiler will prevent building modules with data type incompatible with
    // sunrealtype
    template class Kinsol<sunrealtype, long int>;
    template class Kinsol<sunrealtype, int>;
    template class Kinsol<sunrealtype, size_t>;
  } // namespace Sundials
} // namespace AnalysisManager
