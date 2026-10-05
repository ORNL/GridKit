#pragma once

#include <functional>
#include <memory>
#include <vector>

#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_stepper.h>

#include <GridKit/Model/Coupling.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    /**
     * @brief Lie–Trotter splitting of coupled models (ARKODE SplittingStep).
     *
     * The state is a ManyVector with one block per partition. Each partition's
     * integrator advances its own block, reading its coupling inputs from the
     * other blocks when its stage starts. Outputs are each partition's accepted
     * state with the inputs it was solved with.
     *
     * With a coupling tolerance, a step is accepted only if no coupling input
     * changed by more than the tolerance during it, and a SUNDIALS PI
     * controller sizes the next step from that change.
     */
    template <class ScalarT, typename IdxT>
    class SplittingStep
    {
    public:
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using SolverT    = DynamicSolver<ScalarT, IdxT>;
      using CouplingT  = GridKit::Model::Coupling<ScalarT, IdxT>;

      SplittingStep();
      ~SplittingStep();
      SplittingStep(const SplittingStep&)            = delete;
      SplittingStep& operator=(const SplittingStep&) = delete;

      /// The solver and its model must outlive this object.
      void addPartition(SolverT& solver, std::vector<CouplingT> couplings);

      /// Step size; with a coupling tolerance, the first step after each initialization.
      void setFixedStep(RealT step);
      /// Largest change of a coupling input within a step, relative to 1 + |input| (0: fixed steps).
      void setCouplingTolerance(RealT tolerance);
      /// Consistent-coupling tolerances, as for Ida (abs_tol <= 0: each variable's own).
      void setTolerance(RealT rel_tol, RealT abs_tol);
      void setOutput(std::function<void(RealT)> output);

      void configureSimulation();
      void initializeSimulation(RealT t0);
      void runSimulation(RealT tf, RealT dt_monitor = 0);

      long numSteps() const;
      long numRejectedSteps() const;

    private:
      struct Input
      {
        sunindextype block;
        IdxT         index;
        ScalarT*     value;
        RealT        abs_tol; ///< The source variable's own absolute tolerance
      };

      struct Partition
      {
        SolverT*               solver{};
        SUNStepper             stepper{}; ///< The solver over its own state
        SUNStepper             block{};   ///< The solver over block `index`, given to ARKODE
        sunindextype           index{};
        std::vector<CouplingT> couplings;
        std::vector<Input>     inputs;
      };

      static Partition& content(SUNStepper block);
      static void       couple(const Partition& partition, N_Vector y);
      static SUNErrCode resetBlock(SUNStepper block, sunrealtype t, N_Vector y);
      static int        evolveBlock(SUNStepper block, sunrealtype tout, N_Vector y, sunrealtype* tret);
      static SUNErrCode setBlockStopTime(SUNStepper block, sunrealtype tstop);
      static SUNErrCode setBlockStepDirection(SUNStepper block, sunrealtype direction);

      RealT couplingMismatch(N_Vector y, RealT rel_tol, RealT abs_tol) const;
      void  advance(RealT tout);

      SUNContext                              context_{};
      std::vector<std::unique_ptr<Partition>> partitions_;
      std::vector<N_Vector>                   blocks_; ///< Owned; the ManyVector does not own them
      N_Vector                                y_{};
      N_Vector                                y_saved_{}; ///< Step start, restored on rejection
      void*                                   arkode_{};
      SUNAdaptController                      controller_{};
      std::function<void(RealT)>              output_;
      RealT                                   first_step_{};
      RealT                                   step_{};
      RealT                                   coupling_tol_{};
      RealT                                   rel_tol_{};
      RealT                                   abs_tol_{};
      RealT                                   t_{};
      long                                    steps_{};
      long                                    rejected_steps_{};
    };
  } // namespace Sundials
} // namespace AnalysisManager
