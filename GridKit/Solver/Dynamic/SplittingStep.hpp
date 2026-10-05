#pragma once

#include <exception>
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
     * @brief Multicolor Gauss-Seidel splitting of coupled partitions (ARKODE SplittingStep).
     *
     * Partitions that share no coupling form a color. ARKODE advances the
     * colors in turn, and a color advances its partitions concurrently. Inputs
     * ramp toward partitions already advanced and extrapolate the others. A
     * SUNDIALS PI controller bounds the endpoint prediction mismatch. Output
     * follows acceptance of the complete step.
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

      /// The solver, its model, and its private context must outlive this object.
      void addPartition(SolverT& solver, std::vector<CouplingT> couplings, SUNContext context);

      /// Partitions of one color advance on up to this many threads (OpenMP builds).
      void setNumThreads(int count);

      /// Step size; with a coupling tolerance, the first step after each initialization.
      void setFixedStep(RealT step);
      /// Endpoint prediction mismatch relative to 1 + |input| (0: fixed steps).
      void setCouplingTolerance(RealT tolerance);
      /// Consistent-coupling tolerances, as for Ida (abs_tol <= 0: each variable's own).
      void setTolerance(RealT rel_tol, RealT abs_tol);
      void setOutput(std::function<void(RealT)> output);

      void configureSimulation();
      void initializeSimulation(RealT t0);
      void runSimulation(RealT tf, RealT dt_monitor = 0);

      long  numSteps() const;
      long  numRejectedSteps() const;
      /// Wall time spent advancing a partition, in seconds.
      RealT partitionTime(std::size_t partition) const;

    private:
      using InputT = GridKit::Model::Input<ScalarT>;

      struct Partition;

      /// A coupling resolved to its source partition
      struct Link
      {
        const Partition* source;
        IdxT             index;
        InputT*          input;
        RealT            abs_tol; ///< The source variable's own absolute tolerance
        RealT            start{}; ///< Source value at the step start
        RealT            rate{};  ///< Source rate over the last step
      };

      struct Partition
      {
        SolverT*               solver{};
        SUNContext             context{}; ///< Borrowed from the regional solver
        SUNStepper             stepper{}; ///< The solver over its own state
        sunindextype           index{};
        RealT                  time{};    ///< Time of the state in its block
        RealT                  seconds{}; ///< Wall time spent advancing it
        std::vector<CouplingT> couplings;
        std::vector<Link>      links;
      };

      /// Partitions that share no coupling; one ARKODE partition
      struct Color
      {
        std::vector<Partition*> partitions;
        SUNStepper              stepper{};
        int                     threads{1};
        std::exception_ptr      failure; ///< Rethrown by advance(), since ARKODE is C
      };

      static Color& content(SUNStepper stepper);
      static RealT  value(const Link& link, N_Vector y);
      static void   couple(const Partition& partition, N_Vector y, RealT t);
      static int    evolvePartition(Partition& partition, sunrealtype tout, N_Vector y, sunrealtype* tret);
      template <class Function>
      static std::exception_ptr concurrently(const Color& color, Function&& function);
      static SUNErrCode         resetColor(SUNStepper stepper, sunrealtype t, N_Vector y);
      static int                evolveColor(SUNStepper stepper, sunrealtype tout, N_Vector y, sunrealtype* tret);
      static SUNErrCode         setColorStopTime(SUNStepper stepper, sunrealtype tstop);
      static SUNErrCode         setColorStepDirection(SUNStepper stepper, sunrealtype direction);

      void  colorPartitions();
      RealT couplingMismatch(N_Vector y, RealT t, RealT rel_tol, RealT abs_tol) const;
      void  advance(RealT tout);

      int                                     num_threads_{1};
      SUNContext                              context_{};
      std::vector<std::unique_ptr<Partition>> partitions_;
      std::vector<std::unique_ptr<Color>>     colors_; ///< In ARKODE's order
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
