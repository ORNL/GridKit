
#pragma once

#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <vector>

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
      long int num_jacobian_evals_              = 0;
      long int num_linear_decompositions_       = 0;
      long int num_error_test_fails_            = 0;
      long int num_nonlinear_iters_             = 0;
      long int num_nonlinear_convergence_fails_ = 0;

      IdaStats&   operator+=(const IdaStats& other);
      std::string report() const;
    };

    enum class IdaConsistentICType
    {
      Y,
      YA_YDP
    };

    /**
     * @brief Integrator state recorded after one IDA internal step
     *
     * Counters are cumulative within an IDA segment. Differencing adjacent
     * records in the same segment gives the work performed by one step.
     */
    struct IdaStepRecord
    {
      int         segment                         = 0;
      sunrealtype t                               = 0.0;
      sunrealtype h                               = 0.0;
      sunrealtype h_next                          = 0.0;
      int         order                           = 0;
      int         order_next                      = 0;
      long int    num_steps                       = 0;
      long int    num_residual_evals              = 0;
      long int    num_jacobian_evals              = 0;
      long int    num_error_test_fails            = 0;
      long int    num_nonlinear_iters             = 0;
      long int    num_nonlinear_convergence_fails = 0;
    };

    enum class KluOrdering
    {
      AMD     = 0,
      COLAMD  = 1,
      NATURAL = 2
    };

    /**
     * @brief Optional forward-simulation controls passed directly to IDA
     */
    template <class RealT>
    struct IdaOptions
    {
      RealT rel_tol{1e-5};
      RealT abs_tol{};

      std::optional<RealT> fixed_step;
      std::optional<RealT> init_step;
      std::optional<RealT> min_step;
      std::optional<RealT> max_step;

      std::optional<int>      max_order;
      std::optional<long int> max_num_steps;
      std::optional<int>      max_err_test_fails;
      std::optional<bool>     suppress_alg;

      std::optional<int>   max_nonlin_iters;
      std::optional<int>   max_conv_fails;
      std::optional<RealT> nonlin_conv_coef;

      std::optional<int>   max_num_steps_ic;
      std::optional<int>   max_num_jacs_ic;
      std::optional<int>   max_num_iters_ic;
      std::optional<int>   max_backs_ic;
      std::optional<bool>  line_search_off_ic;
      std::optional<RealT> nonlin_conv_coef_ic;
      std::optional<RealT> step_tolerance_ic;

      std::optional<bool>        linear_solution_scaling;
      std::optional<RealT>       delta_cj_lsetup;
      std::optional<KluOrdering> klu_ordering;
    };

    template <class ScalarT, typename IdxT>
    class Ida : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT    = typename EvaluatorT::VectorT;

    public:
      using Options = IdaOptions<RealT>;

      Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Ida();

      int configureSimulation();
      int configureLinearSolver();
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
      int configureLinearSolverSparse();
#endif
      int configureLinearSolverDense();
      int getDefaultInitialCondition();
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
      void printPerformanceStats() const;

      void setFixedStep(ScalarT time_step);
      void setBackwardFixedStep(ScalarT time_step);
      using DynamicSolver<ScalarT, IdxT>::setTolerance;
      void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) override;
      void setBackwardTolerance(ScalarT rel_tol, ScalarT abs_tol_override = 0);
      void setQuadratureTolerance(ScalarT rel_tol,
                                  ScalarT abs_tol_override = 0);
      void setBackwardQuadratureTolerance(ScalarT rel_tol,
                                          ScalarT abs_tol_override = 0);
      void setOptions(const Options& options);
      void setSuppressAlgebraicErrors(bool suppress);
      void setBackwardSuppressAlgebraicErrors(bool suppress);
      void setConsistentICType(IdaConsistentICType consistent_ic_type);
      void setKluOrdering(KluOrdering ordering);
      void setMaxSteps(IdxT maxSteps) override;
      void setBackwardMaxSteps(IdxT maxSteps);

      IdaStats getStats() const;

      void enableStepTrace(bool enable = true)
      {
        trace_enabled_ = enable;
      }

      void clearStepTrace()
      {
        step_trace_.clear();
      }

      void setTraceSegment(int segment)
      {
        trace_segment_ = segment;
      }

      const std::vector<IdaStepRecord>& getStepTrace() const
      {
        return step_trace_;
      }

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
      int   runSimulationTraced(
            RealT                                            tf,
            RealT                                            dt_monitor,
            const std::optional<std::function<void(RealT)>>& step_callback);
      void recordStep(RealT t);

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

      Options             options_{};
      IdaConsistentICType consistent_ic_type_{IdaConsistentICType::YA_YDP};

      RealT    backward_time_step_{};
      RealT    backward_rel_tol_{DEFAULT_REL_TOL};
      RealT    backward_abs_tol_override_{};
      long int backward_max_steps_{};
      bool     backward_suppress_alg_{false};

      RealT quadrature_rel_tol_{0.1 * DEFAULT_REL_TOL};
      RealT quadrature_abs_tol_override_{};

      RealT backward_quadrature_rel_tol_{0.1 * DEFAULT_REL_TOL};
      RealT backward_quadrature_abs_tol_override_{};

      std::vector<IdaStepRecord> step_trace_{};
      bool                       trace_enabled_{false};
      int                        trace_segment_{0};

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      static void copyVec(const N_Vector x, VectorT& y);
      static void copyVec(const VectorT& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      static void checkAllocation(void* v, const char* functionName);
      static void checkOutput(int retval, const char* functionName);

      void        setIDAOptions(void*    mem,
                                ScalarT  time_step,
                                ScalarT  rel_tol,
                                ScalarT  abs_tol_override,
                                long int max_steps,
                                bool     suppress_alg);
      void        applyIDAOptions(void* mem, const Options& options);
      void        applyLinearSolverOptions(void* mem, const Options& options);
      static void validateOptions(const Options& options);
      void        setTolerance(void*   mem,
                               ScalarT rel_tol,
                               ScalarT abs_tol_override,
                               ScalarT abs_tol_fac = 1);
      void        setQuadratureTolerance(void*   mem,
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
