#pragma once

#include <atomic>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <string>
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
     * @brief Multicolor Gauss-Seidel splitting of coupled partitions.
     *
     * Partitions that share no coupling form a color, and a coupled pair
     * advances in color order. A partition starts once the partitions it
     * follows have finished, so independent partitions advance concurrently.
     * Inputs ramp toward partitions already advanced and extrapolate the
     * others. A SUNDIALS PI controller bounds the endpoint prediction mismatch,
     * and KINSOL accelerates the consistent-coupling sweeps. Output comes from
     * each partition's interpolant once the complete step is accepted.
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

      /// Partitions advance on up to this many threads (OpenMP builds).
      void setNumThreads(int count);

      /// Step size; with a coupling tolerance, the first step after each initialization.
      void setFixedStep(RealT step);
      /// Endpoint prediction mismatch relative to 1 + |input| (0: fixed steps).
      void setCouplingTolerance(RealT tolerance);
      /// Consistent-coupling tolerances, as for Ida (abs_tol <= 0: each variable's own).
      void setTolerance(RealT rel_tol, RealT abs_tol);
      /// Called with each output time while every model holds its state at that time.
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
        SolverT*                solver{};
        SUNContext              context{}; ///< Borrowed from the regional solver
        SUNStepper              stepper{}; ///< The solver over its own state
        sunindextype            index{};
        RealT                   time{};    ///< Time of the state in its block
        RealT                   seconds{}; ///< Wall time spent advancing it
        std::vector<CouplingT>  couplings;
        std::vector<Link>       links;
        std::vector<Partition*> successors;     ///< Coupled partitions of later colors
        int                     predecessors{}; ///< Coupled partitions of earlier colors
        std::atomic<int>        waiting{};      ///< Predecessors not yet finished
        std::vector<RealT>      y_start;        ///< State at the step start, restored on rejection
        std::vector<RealT>      outputs;        ///< States at the output times inside the step
      };

      static RealT& value(const Link& link, N_Vector y);
      static void   couple(const Partition& partition, N_Vector y, RealT t);
      static int    sweep(N_Vector u, N_Vector g, void* user_data);

      void colorPartitions();
      template <class Function>
      void        traverse(Function&& function);
      void        evolve(Partition& partition, RealT end);
      void        emitOutputs();
      void        consistentCoupling();
      RealT       couplingMismatch(N_Vector y, RealT t, RealT rel_tol, RealT abs_tol, const Link** worst = nullptr) const;
      std::string describe(const Link& link) const;
      void        advance(RealT tout);

      int                                     num_threads_{1};
      SUNContext                              context_{};
      std::vector<std::unique_ptr<Partition>> partitions_;
      std::vector<Partition*>                 roots_;  ///< Partitions of the first color
      std::vector<Link*>                      links_;  ///< Every link, in the order of u_
      std::vector<N_Vector>                   blocks_; ///< Owned; the ManyVector does not own them
      N_Vector                                y_{};
      N_Vector                                u_{}; ///< Coupling sources, the consistency unknowns
      N_Vector                                u_scale_{};
      N_Vector                                ones_{};
      void*                                   kinsol_{};
      std::exception_ptr                      failure_; ///< Raised inside KINSOL, rethrown after it
      const Link*                             worst_{}; ///< Largest mismatch of the last step
      SUNAdaptController                      controller_{};
      std::function<void(RealT)>              output_;
      std::deque<RealT>                       output_times_; ///< Pending, ascending
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
