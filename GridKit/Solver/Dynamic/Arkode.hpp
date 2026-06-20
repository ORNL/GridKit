#pragma once

#include <functional>
#include <optional>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/SundialsException.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <class, typename>
    class CsrMatrix;
  }
} // namespace GridKit

namespace AnalysisManager
{
  namespace Sundials
  {
    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class Arkode : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

    protected:
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using CsrMatrixT = GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>;

    public:
      Arkode(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      virtual ~Arkode();

      int configureSimulation();
      int configureLinearSolver();
      int getDefaultInitialCondition();
      int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      int initializeSimulation(RealT t0);
      int runSimulation(RealT tf, int nout = 1, std::optional<std::function<void(RealT)>> step_callback = {});
      int deleteSimulation();

      RealT getInitialTime()
      {
        return t_init_;
      }

      RealT getFinalTime()
      {
        return t_final_;
      }

      int getNumberOutputTimes()
      {
        return nout_;
      }

      void printFinalStats() const;

      void setFixedStep(ScalarT time_step);
      void setOrder(int order);
      using DynamicSolver<ScalarT, IdxT>::setTolerance;
      void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) override;
      void setMaxSteps(IdxT max_steps) override;

    protected:
      virtual void* createStepper(RealT t0, N_Vector y0) = 0;

      virtual bool usesImplicitSolver() const
      {
        return true;
      }

      virtual bool usesMassMatrix() const
      {
        return true;
      }

      virtual bool requiresIdentityMass() const
      {
        return false;
      }

      int  setupJacobianSolver();
      int  setupMassSolver();
      int  recoverMassMatrix(RealT t, SUNMatrix M);
      bool massMatrixIsIdentity(RealT t);

      static int Rhs(RealT t, N_Vector yy, N_Vector ydot, void* user_data);
      static int Jac(RealT t, N_Vector yy, N_Vector fy, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
      static int Mass(RealT t, SUNMatrix M, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

      void setArkodeOptions(void* mem, RealT time_step, RealT rel_tol, RealT abs_tol_override, IdxT max_steps);
      void setTolerance(void* mem, RealT rel_tol, RealT abs_tol_override);

      static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
      static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
      static void copyCsrToSun(CsrMatrixT* src, SUNMatrix dst);

      static void checkAllocation(void* v, const char* functionName);
      static void checkOutput(int retval, const char* functionName);

    protected:
      static constexpr RealT DEFAULT_REL_TOL = 1e-5;

      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNMatrix       MassMat_{};
      SUNLinearSolver linearSolver_{};
      SUNLinearSolver massLinearSolver_{};

      RealT t_init_{};
      RealT t_final_{};
      int   nout_{};

      N_Vector yy_{};
      N_Vector yy0_{};

      RealT time_step_{};
      RealT rel_tol_{DEFAULT_REL_TOL};
      RealT abs_tol_override_{};
      IdxT  max_steps_{};
      int   order_{};
    };

  } // namespace Sundials
} // namespace AnalysisManager
