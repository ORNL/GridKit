
/**
 * @file Kinsol.cpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains definition of interface to KINSOL nonlinear solver from
 * SUNDIALS library.
 *
 */

#include "Kinsol.hpp"

#include <iomanip>
#include <iostream>

#include <kinsol/kinsol.h>             // access to KINSOL func., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{

  namespace Sundials
  {

    template <class ScalarT, typename IdxT>
    Kinsol<ScalarT, IdxT>::Kinsol(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : SteadyStateSolver<ScalarT, IdxT>(model)
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
      int retval = 0;

      // Allocate solution vectors
      yy_ = N_VNew_Serial(static_cast<sunindextype>(model_->size()), context_);
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

      // Set up linear solver
      return this->configureLinearSolver();
    }

    template <class ScalarT, typename IdxT>
    int Kinsol<ScalarT, IdxT>::configureLinearSolver()
    {
      int retval = 0;

      // Set up linear solver
      JacobianMat_ = SUNDenseMatrix(static_cast<sunindextype>(model_->size()),
                                    static_cast<sunindextype>(model_->size()),
                                    context_);
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
      // printOutput(tout);
      // std::cout << "\n";
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
    int Kinsol<ScalarT, IdxT>::Residual(N_Vector yy, N_Vector rr, void* user_data)
    {
      GridKit::Model::Evaluator<ScalarT, IdxT>* model =
          static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

      copyVec(yy, model->y());

      model->evaluateResidual();
      copyVec(model->getResidual(), rr);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const N_Vector x, VectorT& y)
    {
      const ScalarT* xdata = N_VGetArrayPointer(x);
      std::copy_n(xdata, static_cast<size_t>(y.getSize()), y.getData());
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::copyVec(const VectorT& x, N_Vector y)
    {
      ScalarT* ydata = N_VGetArrayPointer(y);
      std::copy_n(x.getData(), static_cast<size_t>(x.getSize()), ydata);
    }

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::printOutput() const
    {
      sunrealtype* yval = N_VGetArrayPointer(yy_);

      std::cout << std::setprecision(5) << std::setw(7);
      for (IdxT i = 0; i < model_->size(); ++i)
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

    template <class ScalarT, typename IdxT>
    void Kinsol<ScalarT, IdxT>::setTolerance(ScalarT tol)
    {
      int retval = KINSetFuncNormTol(solver_, tol);
      checkOutput(retval, "KINSetFuncNormTol");

      retval = KINSetScaledStepTol(solver_, tol);
      checkOutput(retval, "KINSetScaledStepTol");
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class Kinsol<sunrealtype, long int>;
    template class Kinsol<sunrealtype, int>;
    template class Kinsol<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
