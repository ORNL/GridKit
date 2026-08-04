
#pragma once

#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

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
      std::string sundials_version_;
      int         sundials_logging_level_          = 0;
      long int    num_steps_                       = 0;
      long int    num_residual_evals_              = 0;
      long int    num_linear_solver_setups_        = 0;
      long int    num_error_test_fails_            = 0;
      long int    num_backtrack_operations_        = 0;
      long int    num_nonlinear_iters_             = 0;
      long int    num_nonlinear_convergence_fails_ = 0;
      long int    num_nonlinear_step_fails_        = 0;
      long int    num_jacobian_evals_              = 0;
      long int    last_jacobian_step_              = 0;
      long int    num_linear_iters_                = 0;
      long int    num_linear_convergence_fails_    = 0;
      long int    num_linear_residual_evals_       = 0;
      long int    num_preconditioner_evals_        = 0;
      long int    num_preconditioner_solves_       = 0;
      long int    num_jtimes_setup_evals_          = 0;
      long int    num_jtimes_evals_                = 0;
      long int    last_linear_flag_                = 0;
      std::string last_linear_flag_name_;
      int         last_order_          = 0;
      int         current_order_       = 0;
      sunrealtype actual_initial_step_ = 0.0;
      sunrealtype last_step_           = 0.0;
      sunrealtype current_step_        = 0.0;
      sunrealtype current_time_        = 0.0;
      sunrealtype current_cj_          = 0.0;
      sunrealtype jacobian_time_       = 0.0;
      sunrealtype jacobian_cj_         = 0.0;

      IdaStats&   operator+=(const IdaStats& other);
      std::string report() const;
    };

    enum class IdaConsistentICType
    {
      Y,
      YA_YDP
    };

    enum class IdaLogLevel
    {
      Error   = 1,
      Warning = 2
    };

    struct IdaLogOptions
    {
      std::filesystem::path file;
      IdaLogLevel           level{IdaLogLevel::Warning};
    };

    template <class ScalarT, typename IdxT>
    class Ida : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT    = typename EvaluatorT::VectorT;

    public:
      using OutputCallback       = std::function<void(RealT)>;
      using InternalStepCallback = std::function<void(const IdaStats&)>;

      Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model,
          std::optional<IdaLogOptions>              log_options = {});
      ~Ida();

      int configureSimulation();
      int configureLinearSolver();
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
      int configureLinearSolverSparse();
#endif
      int configureLinearSolverDense();
      int getDefaultInitialCondition();
      int setMaxOrder(int max_order);
      int setMaxStep(RealT hmax);
      int initializeSimulation(RealT                t0,
                               bool                 findConsistent  = true,
                               std::optional<RealT> consistent_tout = {});

      int runSimulation(RealT                         tf,
                        std::optional<RealT>          dt_monitor    = RealT{0},
                        std::optional<OutputCallback> step_callback = {});
      int runSimulationWithStepHistory(RealT                         tf,
                                       std::optional<RealT>          dt_monitor,
                                       const InternalStepCallback&   internal_step_callback,
                                       std::optional<OutputCallback> step_callback = {});
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
      void setConsistentICType(IdaConsistentICType consistent_ic_type);
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
      int   getIDAConsistentICType() const;
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

      N_Vector yy_step_{}; ///< Storage for the solver step output during output interpolation
      N_Vector yp_step_{}; ///< Storage for the solver step derivatives during output interpolation

      N_Vector yyB_{}; ///< Adjoint solution vector
      N_Vector ypB_{}; ///< Adjoint solution derivatives vector
      N_Vector qB_{};  ///< Backward integrand vector

      int       backwardID_{};
      SUNLogger logger_{};

      RealT               time_step_{};
      RealT               rel_tol_{DEFAULT_REL_TOL};
      RealT               abs_tol_override_{};
      IdxT                max_steps_{};
      bool                suppress_alg_{false};
      IdaConsistentICType consistent_ic_type_{IdaConsistentICType::YA_YDP};

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
      void        configureLogger(const IdaLogOptions& options);
      static void copyVec(const N_Vector x, VectorT& y);
      static void copyVec(const VectorT& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);
      void        interpolateSolution(RealT t);
      void        publishOutput(RealT t, const std::optional<OutputCallback>& step_callback);

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
