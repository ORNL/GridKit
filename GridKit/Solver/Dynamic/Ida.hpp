
#pragma once

#include <exception>
#include <functional>
#include <iostream>
#include <optional>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense linear solver        */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */

#include <GridKit/Definitions.hpp>

#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
#include <sunlinsol/sunlinsol_klu.h> /* access to KLU linear solver          */
#endif

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    using Log = ::GridKit::Utilities::Logger;

    struct IdaStats
    {
      long int num_steps_                       = 0;
      long int num_residual_evals_              = 0;
      long int num_linear_decompositions_       = 0;
      long int num_error_test_fails_            = 0;
      long int num_nonlinear_iters_             = 0;
      long int num_nonlinear_convergence_fails_ = 0;

      IdaStats&   operator+=(const IdaStats& other);
      std::string report() const;
    };

    template <class ScalarT, typename IdxT>
    class Ida : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT    = typename EvaluatorT::VectorT;

    public:
      Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Ida();

      int configureSimulation();
      int configureLinearSolver();
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
      int configureLinearSolverSparse();
#endif
      int configureLinearSolverDense();
      int getDefaultInitialCondition();
      int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      int initializeSimulation(RealT t0, bool findConsistent = true);

      int runSimulation(RealT tf, RealT dt_monitor = 0, std::optional<std::function<void(RealT)>> step_callback = {});
      int deleteSimulation();

      int configureQuadrature();
      int initializeQuadrature();
      int runSimulationQuadrature(RealT tf, RealT dt_monitor = 0);
      int deleteQuadrature();

      int configureAdjoint();
      int configureLinearSolverBackward();
      int initializeAdjoint(IdxT steps = 100);
      int initializeBackwardSimulation(RealT tf);
      int runForwardSimulation(RealT tf, RealT dt_monitor = 0);
      int runBackwardSimulation(RealT t0);
      int deleteAdjoint();
      int deleteBackwardSimulation();

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

      RealT getInitialTime()
      {
        return t_init_;
      }

      const RealT* getIntegral() const
      {
        return N_VGetArrayPointer(q_);
      }

      RealT* getIntegral()
      {
        return N_VGetArrayPointer(q_);
      }

      const RealT* getAdjointIntegral() const
      {
        return N_VGetArrayPointer(qB_);
      }

      RealT* getAdjointIntegral()
      {
        return N_VGetArrayPointer(qB_);
      }

      void printOutput(RealT t) const;
      void printSpecial(RealT t, N_Vector x) const;
      void printFinalStats() const;

      void setFixedStep(ScalarT time_step);
      void setBackwardFixedStep(ScalarT time_step);
      using DynamicSolver<ScalarT, IdxT>::setTolerance;
      void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) override;
      void setBackwardTolerance(ScalarT rel_tol, ScalarT abs_tol_override = 0);
      void setQuadratureTolerance(ScalarT rel_tol,
                                  ScalarT abs_tol_override = 0);
      void setBackwardQuadratureTolerance(ScalarT rel_tol,
                                          ScalarT abs_tol_override = 0);
      void setSuppressAlgebraicErrors(bool suppress);
      void setBackwardSuppressAlgebraicErrors(bool suppress);
      void setMaxSteps(IdxT maxSteps) override;
      void setBackwardMaxSteps(IdxT maxSteps);

      IdaStats getStats() const;

    private:
      static int Residual(RealT    t,
                          N_Vector yy,
                          N_Vector yp,
                          N_Vector rr,
                          void*    user_data);

      static int Jac(RealT     t,
                     RealT     cj,
                     N_Vector  yy,
                     N_Vector  yp,
                     N_Vector  resvec,
                     SUNMatrix J,
                     void*     user_data,
                     N_Vector  tmp1,
                     N_Vector  tmp2,
                     N_Vector  tmp3);

      static int Integrand(RealT    t,
                           N_Vector yy,
                           N_Vector yp,
                           N_Vector rhsQ,
                           void*    user_data);

      static int adjointResidual(RealT    t,
                                 N_Vector yy,
                                 N_Vector yp,
                                 N_Vector yyB,
                                 N_Vector ypB,
                                 N_Vector rrB,
                                 void*    user_data);

      static int adjointIntegrand(RealT    t,
                                  N_Vector yy,
                                  N_Vector yp,
                                  N_Vector yyB,
                                  N_Vector ypB,
                                  N_Vector rhsQB,
                                  void*    user_data);

      int   getMonitorStepCount(RealT tf, RealT dt_monitor) const;
      RealT getMonitorTime(RealT tf, RealT dt_monitor, int step, int nsteps) const;
      void  updateModelState(RealT t);

    private:
      static constexpr ScalarT DEFAULT_REL_TOL = 1e-5;

      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNMatrix       JacobianMatB_{};
      SUNLinearSolver linearSolver_{};
      SUNLinearSolver linearSolverB_{};

      RealT t_init_{};

      N_Vector yy_{};  ///< Solution vector
      N_Vector yp_{};  ///< Solution derivatives vector
      N_Vector tag_{}; ///< Tags differential variables
      N_Vector q_{};   ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values
      N_Vector yp0_{}; ///< Storage for initial derivatives

      N_Vector yyB_{}; ///< Adjoint solution vector
      N_Vector ypB_{}; ///< Adjoint solution derivatives vector
      N_Vector qB_{};  ///< Backward integrand vector

      int backwardID_{};

      RealT time_step_{};
      RealT rel_tol_{DEFAULT_REL_TOL};
      RealT abs_tol_override_{};
      IdxT  max_steps_{};
      bool  suppress_alg_{false};

      RealT backward_time_step_{};
      RealT backward_rel_tol_{DEFAULT_REL_TOL};
      RealT backward_abs_tol_override_{};
      IdxT  backward_max_steps_{};
      bool  backward_suppress_alg_{false};

      RealT quadrature_rel_tol_{0.1 * DEFAULT_REL_TOL};
      RealT quadrature_abs_tol_override_{};

      RealT backward_quadrature_rel_tol_{0.1 * DEFAULT_REL_TOL};
      RealT backward_quadrature_abs_tol_override_{};

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      static void copyVec(const N_Vector x, VectorT& y);
      static void copyVec(const VectorT& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      static void checkAllocation(void* v, const char* functionName);
      static void checkOutput(int retval, const char* functionName);

      void setIDAOptions(void*   mem,
                         ScalarT time_step,
                         ScalarT rel_tol,
                         ScalarT abs_tol_override,
                         IdxT    max_steps,
                         bool    suppress_alg);
      void setTolerance(void*   mem,
                        ScalarT rel_tol,
                        ScalarT abs_tol_override,
                        ScalarT abs_tol_fac = 1);
      void setQuadratureTolerance(void*   mem,
                                  ScalarT rel_tol,
                                  ScalarT abs_tol_override);
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
