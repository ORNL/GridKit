#pragma once

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>

#include <arkode/arkode_pdaestep.h>

namespace AnalysisManager
{
  namespace Sundials
  {
    // TODO: remove once model supports partitioning
    using Mask = std::vector<bool>;

    struct PartitionedStats
    {
      long int num_steps_                        = 0;
      long int num_component_residual_evals_     = 0;
      long int num_component_linear_setups_      = 0;
      long int num_component_error_test_fails_   = 0;
      long int num_component_nonlinear_iters_    = 0;
      long int num_component_nonlinear_failures_ = 0;
      long int num_coupling_linear_setups_       = 0;
      long int num_coupling_nonlinear_iters_     = 0;
      long int num_coupling_nonlinear_failures_  = 0;

      std::string report() const;
    };

    // TODO: unify with SundialsException
    class PartitionedException : public std::runtime_error
    {
    public:
      explicit PartitionedException(const std::string& message)
        : std::runtime_error(message)
      {
      }
    };

    template <class ScalarT, typename IdxT>
    class Partitioned : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT    = typename EvaluatorT::VectorT;

    public:
      Partitioned(EvaluatorT*       model,
                  std::vector<Mask> component_masks,
                  Mask              coupling_mask);
      ~Partitioned() override;

      int configureSimulation();
      int initializeSimulation(RealT t0);
      int runSimulation(
          RealT                                     tf,
          RealT                                     dt_monitor    = 0,
          std::optional<std::function<void(RealT)>> step_callback = {});
      int deleteSimulation();

      void setFixedStep(ScalarT time_step);
      using DynamicSolver<ScalarT, IdxT>::setTolerance;
      void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) override;
      void setMaxSteps(IdxT max_steps) override;

      PartitionedStats getStats() const;

    private:
      static int ComponentResidual(int      partition,
                                   RealT    t,
                                   N_Vector y,
                                   N_Vector w,
                                   N_Vector yp,
                                   N_Vector residual,
                                   void*    user_data);
      static int ComponentJacTimes(int      partition,
                                   RealT    t,
                                   N_Vector y,
                                   N_Vector w,
                                   N_Vector yp,
                                   N_Vector residual,
                                   N_Vector v,
                                   N_Vector Jv,
                                   RealT    cj,
                                   void*    user_data,
                                   N_Vector tmp1,
                                   N_Vector tmp2);
      static int AlgebraicResidual(RealT    t,
                                   N_Vector y,
                                   N_Vector w,
                                   N_Vector residual,
                                   void*    user_data);
      static int AlgebraicJacTimes(RealT    t,
                                   N_Vector y,
                                   N_Vector w,
                                   N_Vector v,
                                   N_Vector Jv,
                                   void*    user_data,
                                   N_Vector tmp);

      int componentResidual(int      partition,
                            RealT    t,
                            N_Vector y,
                            N_Vector w,
                            N_Vector yp,
                            N_Vector residual);
      int componentJacTimes(int      partition,
                            RealT    t,
                            RealT    cj,
                            N_Vector y,
                            N_Vector w,
                            N_Vector yp,
                            N_Vector v,
                            N_Vector Jv);
      int algebraicResidual(RealT    t,
                            N_Vector y,
                            N_Vector w,
                            N_Vector residual);
      int algebraicJacTimes(RealT    t,
                            N_Vector y,
                            N_Vector w,
                            N_Vector v,
                            N_Vector Jv);

      void validateAndBuildIndices();
      void allocateStateVectors();
      void createSolver(RealT t0);
      void configurePartitionSolvers();
      void configureTolerances();
      void updateModelState(RealT t);
      void scatterFullState(N_Vector y);
      void scatterAlgebraicState(N_Vector y, N_Vector w);
      void scatterComponentState(int      partition,
                                 N_Vector y,
                                 N_Vector w,
                                 N_Vector yp);
      void multiplyMaskedJacobian(const std::vector<IdxT>& indices,
                                  N_Vector                 v,
                                  N_Vector                 Jv);

      void        gather(const VectorT&           source,
                         const std::vector<IdxT>& indices,
                         N_Vector                 destination) const;
      void        scatter(N_Vector                 source,
                          const std::vector<IdxT>& indices,
                          VectorT&                 destination) const;
      static void checkOutput(int retval, const char* function_name);
      static void checkAllocation(const void* pointer, const char* function_name);
      int         getMonitorStepCount(RealT tf, RealT dt_monitor) const;
      RealT       getMonitorTime(RealT tf,
                                 RealT dt_monitor,
                                 int   step,
                                 int   steps) const;

      std::vector<Mask>              component_masks_;
      Mask                           coupling_mask_;
      std::vector<std::vector<IdxT>> differential_indices_;
      std::vector<std::vector<IdxT>> algebraic_indices_;
      std::vector<std::vector<IdxT>> component_indices_;
      std::vector<IdxT>              coupling_indices_;
      std::vector<IdxT>              all_algebraic_indices_;

      SUNContext                   context_{};
      void*                        solver_{};
      std::vector<N_Vector>        x_;
      std::vector<N_Vector>        z_;
      std::vector<N_Vector>        xp_;
      std::vector<N_Vector>        zp_;
      N_Vector                     w_{};
      N_Vector                     wp_{};
      N_Vector                     y_{};
      N_Vector                     yp_{};
      SUNLinearSolver              algebraic_linear_solver_{};
      std::vector<SUNLinearSolver> component_linear_solvers_;

      RealT t_init_{};
      RealT time_step_{};
      RealT rel_tol_{1.0e-5};
      RealT abs_tol_override_{};
      IdxT  max_steps_{};
    };
  } // namespace Sundials
} // namespace AnalysisManager
