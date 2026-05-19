
#pragma once

#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <optional>
#include <string>

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
      long int    num_jacobian_eval_steps_         = 0;
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

      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    public:
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
      int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      int initializeSimulation(RealT                t0,
                               bool                 findConsistent  = false,
                               std::optional<RealT> consistent_tout = {});

      int runSimulation(RealT tf, int nout = 1, std::optional<std::function<void(RealT)>> step_callback = {});
      int deleteSimulation();

      int configureQuadrature();
      int initializeQuadrature();
      int runSimulationQuadrature(RealT tf, int nout = 1);
      int deleteQuadrature();

      int configureAdjoint();
      int configureLinearSolverBackward();
      int initializeAdjoint(IdxT steps = 100);
      int initializeBackwardSimulation(RealT tf);
      int runForwardSimulation(RealT tf, int nout = 1);
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

      RealT getFinalTime()
      {
        return t_final_;
      }

      int getNumberOutputTimes()
      {
        return nout_;
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

      void printOutput(RealT t);
      void printSpecial(RealT t, N_Vector x);
      void printFinalStats();

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

    private:
      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNMatrix       JacobianMatB_{};
      SUNLinearSolver linearSolver_{};
      SUNLinearSolver linearSolverB_{};

      RealT t_init_{};
      RealT t_final_{};
      int   nout_{}; ///< Number of integration outputs

      N_Vector yy_{};  ///< Solution vector
      N_Vector yp_{};  ///< Solution derivatives vector
      N_Vector tag_{}; ///< Tags differential variables
      N_Vector q_{};   ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values
      N_Vector yp0_{}; ///< Storage for initial derivatives

      N_Vector yyB_{}; ///< Adjoint solution vector
      N_Vector ypB_{}; ///< Adjoint solution derivatives vector
      N_Vector qB_{};  ///< Backward integrand vector

      int       backwardID_{};
      SUNLogger logger_{};

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      void        configureLogger(const IdaLogOptions& options);
      static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
      static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      static void checkAllocation(void* v, const char* functionName);
      static void checkOutput(int retval, const char* functionName);
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
