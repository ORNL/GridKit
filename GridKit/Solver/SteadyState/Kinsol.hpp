
/**
 * @file Kinsol.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains declaration of interface to KINSOL nonlinear solver from
 * SUNDIALS library.
 */

#pragma once

#include <exception>
#include <functional>
#include <iostream>
#include <span>

#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/SteadyState/SteadyStateSolver.hpp>
#include <GridKit/Solver/SundialsException.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class Kinsol : public SteadyStateSolver<ScalarT, IdxT>
    {
      using SteadyStateSolver<ScalarT, IdxT>::model_;
      using Logger = ::GridKit::Utilities::Logger;
      using RealT  = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT =
        typename GridKit::Model::Evaluator<ScalarT, IdxT>::VectorT;

      static_assert(
        std::is_same_v<RealT, sunrealtype>,
        "RealT must be the same type as sunrealtype"
      );

    public:
      using ResidualFunction =
        std::function<int(std::span<const ScalarT>, std::span<ScalarT>)>;

      Kinsol(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      explicit Kinsol(IdxT size);

      ~Kinsol();

      int configureSimulation();
      int configureLinearSolver();
      int getDefaultInitialCondition();
      // int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      // int initializeSimulation();
      int runSimulation();

      /// Use the given values as the initial guess and return the accepted
      /// solution in the same vector.
      ///
      /// @param line_search whether or not to use KINSOL's backtracking line
      ///                    search
      int runSimulation(
        std::span<ScalarT> values,
        ResidualFunction   residual,
        bool               line_search = true
      );

      int deleteSimulation();

      void setTolerance(ScalarT tol) override;
      void setTolerance(ScalarT residual_tol, ScalarT step_tol);

      /// Estimated relative error of residual evaluations.
      ///
      /// This is used to choose finite-difference perturbations.
      void setFunctionError(RealT error);

      /// Maximum number of Newton iterations between Jacobian evaluations.
      void setJacobianUpdateInterval(long int iterations);

      // int configureQuadrature();
      // int initializeQuadrature();
      // int runSimulationQuadrature(RealT tf, int nout=1);
      // int deleteQuadrature();

      // int configureAdjoint();
      // int configureLinearSolverBackward();
      // int initializeAdjoint(IdxT steps = 100);
      // int initializeBackwardSimulation(RealT tf);
      // int runForwardSimulation(RealT tf, int nout=1);
      // int runBackwardSimulation(RealT t0);
      // int deleteAdjoint();

      int saveInitialCondition()
      {
        N_VScale(1.0, yy_, yy0_);
        return 0;
      }

      int getSavedInitialCondition()
      {
        N_VScale(1.0, yy0_, yy_);
        return 0;
      }

      // RealT getInitialTime()
      // {
      //     return t_init_;
      // }

      // RealT getFinalTime()
      // {
      //     return t_final_;
      // }

      // int getNumberOutputTimes()
      // {
      //     return nout_;
      // }

      // const RealT* getIntegral() const
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // RealT* getIntegral()
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // const RealT* getAdjointIntegral() const
      // {
      //     return N_VGetArrayPointer(qB_);
      // }

      // RealT* getAdjointIntegral()
      // {
      //     return N_VGetArrayPointer(qB_);
      // }

      void printOutput() const;
      void printSpecial(sunrealtype t, N_Vector x) const;
      void printFinalStats() const;

    private:
      static int Residual(N_Vector yy, N_Vector rr, void* user_data);
      int        solve(bool line_search = true);

      // static int Integrand(sunrealtype t,
      //                      N_Vector yy,   N_Vector yp,
      //                      N_Vector rhsQ, void *user_data);

      // static int adjointResidual(sunrealtype t,
      //                            N_Vector yy,  N_Vector yp,
      //                            N_Vector yyB, N_Vector ypB,
      //                            N_Vector rrB, void *user_data);

      // static int adjointIntegrand(sunrealtype t,
      //                             N_Vector yy,  N_Vector yp,
      //                             N_Vector yyB, N_Vector ypB,
      //                             N_Vector rhsQB, void *user_data);

    private:
      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNLinearSolver linearSolver_{};

      IdxT               size_;
      ResidualFunction   residual_;
      std::exception_ptr callback_error_;

      N_Vector yy_{};    ///< Solution vector
      N_Vector scale_{}; ///< Scaling vector
      N_Vector q_{};     ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      static void copyVec(const N_Vector x, VectorT& y);
      static void copyVec(const VectorT& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      static inline void checkAllocation(
        SUNContext ctx, void* v, const char* functionName
      );
      static inline void checkOutput(int retval, const char* functionName);
    };

  } // namespace Sundials

} // namespace AnalysisManager
