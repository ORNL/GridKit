#include "SplittingStep.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include <kinsol/kinsol.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#include <GridKit/Definitions.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    namespace
    {
      using Log = GridKit::Utilities::Logger;

      // ARKODE defaults (arkode_adapt_impl.h, arkode_impl.h)
      constexpr double SAFETY         = 0.9;  ///< Fraction of the controller's step taken
      constexpr double GROWTH         = 20.0; ///< Largest step growth
      constexpr double ETAMIN         = 0.1;  ///< Largest step reduction
      constexpr double ETAMXF         = 0.3;  ///< Largest step after a rejection, relative to it
      constexpr int    MAX_REJECTIONS = 7;    ///< MAXNEF: rejections allowed per step

      /// Enough Gauss-Seidel sweeps for a contraction of 0.8 to reduce a mismatch of 1e9 to one.
      constexpr long   MAX_SWEEPS        = 100;
      /// Previous sweeps that Anderson acceleration combines.
      constexpr long   ANDERSON_DEPTH    = 5;
      /// Shortest step, relative to the first, before the coupling counts as collapsed.
      constexpr double MIN_STEP_FRACTION = 1.0e-4;

      void check(int flag, const char* call)
      {
        if (flag < 0)
        {
          Log::error() << "SplittingStep: " << call << " failed with flag " << flag << '\n';
          throw std::runtime_error(std::string("SplittingStep: ") + call + " failed");
        }
      }

      template <typename T>
      T allocated(T object, const char* call)
      {
        if (object == nullptr)
        {
          throw std::runtime_error(std::string("SplittingStep: ") + call + " failed");
        }
        return object;
      }
    } // namespace

    template <class ScalarT, typename IdxT>
    SplittingStep<ScalarT, IdxT>::SplittingStep()
    {
      check(SUNContext_Create(SUN_COMM_NULL, &context_), "SUNContext_Create");
    }

    template <class ScalarT, typename IdxT>
    SplittingStep<ScalarT, IdxT>::~SplittingStep()
    {
      KINFree(&kinsol_);
      SUNAdaptController_Destroy(controller_);
      N_VDestroy(ones_);
      N_VDestroy(u_scale_);
      N_VDestroy(u_);
      N_VDestroy(y_);
      for (auto block : blocks_)
      {
        N_VDestroy(block);
      }
      for (auto& partition : partitions_)
      {
        // SUNStepper_Destroy dereferences the stepper unchecked.
        if (partition->stepper)
        {
          SUNStepper_Destroy(&partition->stepper);
        }
      }
      SUNContext_Free(&context_);
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::addPartition(SolverT& solver, std::vector<CouplingT> couplings, SUNContext context)
    {
      if (y_)
      {
        throw std::logic_error("SplittingStep: partitions are added before configureSimulation");
      }
      // Owned before the stepper exists, so the destructor releases it if the call below throws.
      auto& partition     = *partitions_.emplace_back(std::make_unique<Partition>());
      partition.solver    = &solver;
      partition.context   = context;
      partition.index     = static_cast<sunindextype>(partitions_.size() - 1);
      partition.couplings = std::move(couplings);
      partition.stepper   = solver.createSUNStepper();
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::setNumThreads(int count)
    {
      if (y_)
        throw std::logic_error("SplittingStep: threads are selected before configuration");
      if (count < 1)
        throw std::invalid_argument("SplittingStep: threads must be positive");
#ifndef GRIDKIT_ENABLE_OPENMP
      if (count != 1)
        throw std::invalid_argument("SplittingStep: multiple threads require OpenMP");
#endif
      num_threads_ = count;
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::setFixedStep(RealT step)
    {
      first_step_ = step;
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::setCouplingTolerance(RealT tolerance)
    {
      coupling_tol_ = tolerance;
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::setTolerance(RealT rel_tol, RealT abs_tol)
    {
      rel_tol_ = rel_tol;
      abs_tol_ = abs_tol;
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::setOutput(std::function<void(RealT)> output)
    {
      output_ = std::move(output);
    }

    template <class ScalarT, typename IdxT>
    long SplittingStep<ScalarT, IdxT>::numSteps() const
    {
      return steps_;
    }

    template <class ScalarT, typename IdxT>
    long SplittingStep<ScalarT, IdxT>::numRejectedSteps() const
    {
      return rejected_steps_;
    }

    template <class ScalarT, typename IdxT>
    typename SplittingStep<ScalarT, IdxT>::RealT SplittingStep<ScalarT, IdxT>::partitionTime(std::size_t partition) const
    {
      return partitions_.at(partition)->seconds;
    }

    template <class ScalarT, typename IdxT>
    typename SplittingStep<ScalarT, IdxT>::RealT& SplittingStep<ScalarT, IdxT>::value(const Link& link, N_Vector y)
    {
      return N_VGetSubvectorArrayPointer_ManyVector(y, link.source->index)[link.index];
    }

    /// Inputs over a stage from t: a source already past t is interpolated to
    /// its new value; one still at t is extrapolated at its last rate.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::couple(const Partition& partition, N_Vector y, RealT t)
    {
      for (const auto& link : partition.links)
      {
        const RealT current = value(link, y);
        link.input->value   = current;
        link.input->rate    = link.rate;
        link.input->start   = t;
        if (link.source->time > t)
        {
          link.input->value = link.start;
          link.input->rate  = (current - link.start) / (link.source->time - t);
        }
      }
    }

    /**
     * @brief Runs `function` on every partition, each once its predecessors
     * have finished, so independent partitions run concurrently.
     *
     * After a failure the remaining partitions are skipped; the first failure
     * is rethrown once all have finished.
     */
    template <class ScalarT, typename IdxT>
    template <class Function>
    void SplittingStep<ScalarT, IdxT>::traverse(Function&& function)
    {
      for (auto& partition : partitions_)
      {
        partition->waiting = partition->predecessors;
      }
      std::vector<std::exception_ptr> failures(partitions_.size());
      std::atomic<bool>               failed{false};
      std::function<void(Partition&)> run;
      run = [&](Partition& partition)
      {
        if (!failed)
        {
          Log::ScopedOutput output;
          try
          {
            function(partition);
          }
          catch (...)
          {
            failures[static_cast<std::size_t>(partition.index)] = std::current_exception();
            failed                                              = true;
          }
        }
        for (auto* successor : partition.successors)
        {
          if (--successor->waiting == 0)
          {
#ifdef GRIDKIT_ENABLE_OPENMP
#pragma omp task
#endif
            run(*successor);
          }
        }
      };
#ifdef GRIDKIT_ENABLE_OPENMP
#pragma omp parallel num_threads(num_threads_) if (num_threads_ > 1)
#pragma omp single
#endif
      for (auto* root : roots_)
      {
#ifdef GRIDKIT_ENABLE_OPENMP
#pragma omp task
#endif
        run(*root);
      }
      for (const auto& failure : failures)
      {
        if (failure)
        {
          std::rethrow_exception(failure);
        }
      }
    }

    /// Advance one partition from t_ to `end`, keeping its state at each output time inside the step.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::evolve(Partition& partition, RealT end)
    {
      const auto  start = std::chrono::steady_clock::now();
      N_Vector    block = blocks_[static_cast<std::size_t>(partition.index)];
      const auto* state = N_VGetArrayPointer(block);
      const auto  size  = static_cast<std::size_t>(N_VGetLength(block));

      partition.y_start.assign(state, state + size);
      partition.outputs.clear();

      couple(partition, y_, t_);
      check(SUNStepper_Reset(partition.stepper, t_, block), "SUNStepper_Reset");
      check(SUNStepper_SetStopTime(partition.stepper, end), "SUNStepper_SetStopTime");
      RealT reached = t_;
      for (const RealT t : output_times_)
      {
        if (!(t < end))
        {
          break;
        }
        check(SUNStepper_Evolve(partition.stepper, t, block, &reached), "SUNStepper_Evolve");
        partition.outputs.insert(partition.outputs.end(), state, state + size);
      }
      check(SUNStepper_Evolve(partition.stepper, end, block, &reached), "SUNStepper_Evolve");
      partition.time     = reached;
      partition.seconds += std::chrono::duration<RealT>(std::chrono::steady_clock::now() - start).count();
    }

    /**
     * @brief G(u) for KINSOL: one Gauss-Seidel sweep of consistent states,
     * starting from the coupling sources u and returning the sources it
     * produces. Failures are kept for initializeSimulation, since KINSOL is C.
     */
    template <class ScalarT, typename IdxT>
    int SplittingStep<ScalarT, IdxT>::sweep(N_Vector u, N_Vector g, void* user_data)
    {
      auto& self = *static_cast<SplittingStep*>(user_data);
      try
      {
        const auto* sources = N_VGetArrayPointer(u);
        for (std::size_t i = 0; i < self.links_.size(); ++i)
        {
          value(*self.links_[i], self.y_) = sources[i];
        }
        self.traverse([&self](Partition& partition)
                      {
          couple(partition, self.y_, self.t_);
          auto& model = *partition.solver->getModel();
          auto* block = N_VGetSubvectorArrayPointer_ManyVector(self.y_, partition.index);
          std::copy_n(block, model.size(), model.y().getData());
          model.y().setDataUpdated();
          check(partition.solver->computeConsistentState(self.t_, self.t_ + self.first_step_), "computeConsistentState");
          std::copy_n(model.y().getData(), model.size(), block); });
        auto* produced = N_VGetArrayPointer(g);
        for (std::size_t i = 0; i < self.links_.size(); ++i)
        {
          produced[i] = value(*self.links_[i], self.y_);
        }
        return 0;
      }
      catch (...)
      {
        self.failure_ = std::current_exception();
        return -1;
      }
    }

    /**
     * @brief Colors partitions so that none shares a color with one it reads or
     * that reads it, and orders coupled partitions by color.
     *
     * Greedy first fit, most-coupled partitions first (ties in registration
     * order). Neither the coloring nor the order depends on the thread count.
     */
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::colorPartitions()
    {
      const auto                         count = partitions_.size();
      std::vector<std::set<std::size_t>> neighbors(count);
      for (std::size_t p = 0; p < count; ++p)
      {
        for (const auto& link : partitions_[p]->links)
        {
          const auto source = static_cast<std::size_t>(link.source->index);
          neighbors[p].insert(source);
          neighbors[source].insert(p);
        }
      }

      std::vector<std::size_t> order(count);
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::stable_sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q)
                       { return neighbors[p].size() > neighbors[q].size(); });

      // First fit: the smallest color none of the neighbors colored so far holds.
      std::vector<std::size_t> color(count, count);
      std::size_t              colors = 0;
      for (const auto p : order)
      {
        std::set<std::size_t> taken;
        for (const auto q : neighbors[p])
        {
          taken.insert(color[q]);
        }
        color[p] = 0;
        while (taken.count(color[p]))
        {
          ++color[p];
        }
        colors = std::max(colors, color[p] + 1);
      }
      if (colors < 2)
      {
        throw std::invalid_argument("SplittingStep: partitions must be coupled");
      }

      for (std::size_t p = 0; p < count; ++p)
      {
        for (const auto q : neighbors[p])
        {
          if (color[q] > color[p])
          {
            partitions_[p]->successors.push_back(partitions_[q].get());
            ++partitions_[q]->predecessors;
          }
        }
        if (color[p] == 0)
        {
          roots_.push_back(partitions_[p].get());
        }
      }
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::configureSimulation()
    {
      if (y_)
      {
        throw std::logic_error("SplittingStep: already configured");
      }
      if (!std::isfinite(first_step_) || !(first_step_ > 0.0))
      {
        throw std::invalid_argument("SplittingStep: the step must be positive");
      }

      if (partitions_.size() < 2)
        throw std::invalid_argument("SplittingStep: at least two partitions are required");

      std::set<SUNContext>                          contexts;
      std::map<const EvaluatorT*, const Partition*> partition_of;
      for (const auto& partition : partitions_)
      {
        auto& model = *partition->solver->getModel();
        if (!partition_of.emplace(&model, partition.get()).second)
        {
          throw std::invalid_argument("SplittingStep: a model is added twice");
        }
        if (!partition->context || !contexts.insert(partition->context).second)
          throw std::invalid_argument("SplittingStep: each partition needs its own solver context");
        auto block = allocated(N_VNew_Serial(static_cast<sunindextype>(model.size()), partition->context), "N_VNew_Serial");
        blocks_.push_back(block);
        std::copy_n(model.y().getData(), model.size(), N_VGetArrayPointer(block));
      }

      std::set<InputT*> inputs;
      for (const auto& partition : partitions_)
      {
        for (const auto& coupling : partition->couplings)
        {
          const auto source = partition_of.find(coupling.source);
          if (source == partition_of.end()
              || std::cmp_less(coupling.index, 0)
              || std::cmp_greater_equal(coupling.index, N_VGetLength(blocks_[static_cast<std::size_t>(source->second->index)]))
              || coupling.input == nullptr || !inputs.insert(coupling.input).second)
          {
            throw std::invalid_argument("SplittingStep: a coupling needs a registered source, an index in it, and its own input");
          }
          RealT abs_tol = abs_tol_;
          if (abs_tol <= 0.0)
          {
            abs_tol = coupling.source->absoluteTolerance().getData()[coupling.index];
          }
          partition->links.push_back({source->second, coupling.index, coupling.input, abs_tol});
        }
      }
      for (const auto& partition : partitions_)
      {
        for (auto& link : partition->links)
        {
          links_.push_back(&link);
        }
      }
      colorPartitions();

      y_ = allocated(N_VNew_ManyVector(static_cast<sunindextype>(blocks_.size()), blocks_.data(), context_), "N_VNew_ManyVector");
      if (coupling_tol_ > 0.0)
      {
        controller_ = allocated(SUNAdaptController_PI(context_), "SUNAdaptController_PI");
      }

      u_       = allocated(N_VNew_Serial(static_cast<sunindextype>(links_.size()), context_), "N_VNew_Serial");
      u_scale_ = allocated(N_VClone(u_), "N_VClone");
      ones_    = allocated(N_VClone(u_), "N_VClone");
      N_VConst(1.0, ones_);
      kinsol_ = allocated(KINCreate(context_), "KINCreate");
      check(KINSetMAA(kinsol_, ANDERSON_DEPTH), "KINSetMAA");
      check(KINInit(kinsol_, sweep, u_), "KINInit");
      check(KINSetUserData(kinsol_, this), "KINSetUserData");
      check(KINSetNumMaxIters(kinsol_, MAX_SWEEPS), "KINSetNumMaxIters");
      // The scaled max norm of a sweep's change, as couplingMismatch measures it.
      check(KINSetFuncNormTol(kinsol_, 1.0), "KINSetFuncNormTol");
    }

    /// Largest change any input would see if re-read now, in units of
    /// abs_tol + rel_tol |value| (abs_tol <= 0: each input's own).
    template <class ScalarT, typename IdxT>
    typename SplittingStep<ScalarT, IdxT>::RealT
    SplittingStep<ScalarT, IdxT>::couplingMismatch(N_Vector y, RealT t, RealT rel_tol, RealT abs_tol, const Link** worst) const
    {
      RealT mismatch = 0.0;
      for (const auto* link : links_)
      {
        const RealT source = value(*link, y);
        RealT       weight = abs_tol + rel_tol * std::abs(source);
        if (abs_tol <= 0.0)
        {
          weight = link->abs_tol + rel_tol * std::abs(source);
        }
        const RealT change = std::abs(link->input->at(t) - source) / weight;
        if (!std::isfinite(change))
        {
          throw std::runtime_error("SplittingStep: coupling value is not finite");
        }
        if (change > mismatch && worst != nullptr)
        {
          *worst = link;
        }
        mismatch = std::max(mismatch, change);
      }
      return mismatch;
    }

    /// The partition a link feeds, its place among that partition's couplings, and its source.
    template <class ScalarT, typename IdxT>
    std::string SplittingStep<ScalarT, IdxT>::describe(const Link& link) const
    {
      for (const auto& partition : partitions_)
      {
        const auto& links = partition->links;
        if (!links.empty() && &link >= links.data() && &link < links.data() + links.size())
        {
          return "partition " + std::to_string(partition->index) + " coupling " + std::to_string(&link - links.data())
                 + " (source partition " + std::to_string(link.source->index) + ", variable " + std::to_string(link.index) + ")";
        }
      }
      return "an unknown coupling";
    }

    /**
     * @brief Consistent coupling at t_, differential states fixed.
     *
     * KINSOL iterates Gauss-Seidel sweeps, with Anderson acceleration, until
     * every partition was solved with the coupling values it now has. Each
     * step then starts from the source values this leaves.
     */
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::consistentCoupling()
    {
      auto* sources = N_VGetArrayPointer(u_);
      auto* scale   = N_VGetArrayPointer(u_scale_);
      for (std::size_t i = 0; i < links_.size(); ++i)
      {
        sources[i]    = value(*links_[i], y_);
        RealT abs_tol = abs_tol_;
        if (abs_tol <= 0.0)
        {
          abs_tol = links_[i]->abs_tol;
        }
        scale[i] = 1.0 / (abs_tol + rel_tol_ * std::abs(sources[i]));
      }
      const int flag = KINSol(kinsol_, u_, KIN_FP, ones_, u_scale_);
      if (failure_)
      {
        std::rethrow_exception(std::exchange(failure_, nullptr));
      }
      if (flag < 0)
      {
        throw std::runtime_error("SplittingStep: consistency iteration failed at t = " + std::to_string(t_)
                                 + " (KINSOL flag " + std::to_string(flag) + ")");
      }
      for (auto* link : links_)
      {
        link->start = value(*link, y_);
      }
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::initializeSimulation(RealT t0)
    {
      // Inputs may jump here (an event), so no rate carries over.
      t_ = t0;
      for (auto& partition : partitions_)
      {
        partition->time = t0;
        for (auto& link : partition->links)
        {
          link.rate = 0.0;
        }
      }
      consistentCoupling();
      step_ = first_step_;
      if (controller_)
      {
        check(SUNAdaptController_Reset(controller_), "SUNAdaptController_Reset");
      }
      if (output_)
      {
        output_(t0);
      }
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::runSimulation(RealT tf, RealT dt_monitor)
    {
      output_times_.clear();
      if (output_)
      {
        const RealT t0     = t_;
        const int   nsteps = monitorStepCount(t0, tf, dt_monitor);
        for (int i = 1; i <= nsteps; ++i)
        {
          output_times_.push_back(monitorTime(t0, tf, dt_monitor, i, nsteps));
        }
      }
      advance(tf);
    }

    /// Writes the outputs the accepted step reached: those inside it from the
    /// states each partition kept, one at its end from the current state.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::emitOutputs()
    {
      std::size_t kept = 0;
      for (; !output_times_.empty() && output_times_.front() < t_; output_times_.pop_front(), ++kept)
      {
        for (auto& partition : partitions_)
        {
          auto&      model = *partition->solver->getModel();
          const auto size  = static_cast<std::size_t>(model.size());
          std::copy_n(partition->outputs.data() + kept * size, size, model.y().getData());
          model.y().setDataUpdated();
        }
        output_(output_times_.front());
      }
      if (kept > 0)
      {
        for (auto& partition : partitions_)
        {
          auto& model = *partition->solver->getModel();
          std::copy_n(N_VGetArrayPointer(blocks_[static_cast<std::size_t>(partition->index)]), model.size(), model.y().getData());
          model.y().setDataUpdated();
        }
      }
      if (!output_times_.empty() && output_times_.front() == t_)
      {
        output_(t_);
        output_times_.pop_front();
      }
    }

    /// Accept or reject all regional solves together.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::advance(RealT tout)
    {
      int rejections = 0;
      while (t_ < tout)
      {
        if (step_ < MIN_STEP_FRACTION * first_step_)
        {
          throw std::runtime_error("SplittingStep: coupling step collapsed at t = " + std::to_string(t_)
                                   + "; largest mismatch on " + describe(*worst_));
        }
        const RealT t_start = t_;
        const RealT end     = std::min(t_start + step_, tout); // Lands exactly on tout
        const RealT h       = end - t_start;
        if (!std::isfinite(h) || !(h > 0.0))
          throw std::runtime_error("SplittingStep: step cannot advance time");
        traverse([&](Partition& partition)
                 { evolve(partition, end); });

        if (controller_)
        {
          // How far each partition's inputs ended from their sources; O(h^2) once rates are known.
          const RealT error = couplingMismatch(y_, end, coupling_tol_, coupling_tol_, &worst_);
          RealT       next  = h;
          check(SUNAdaptController_EstimateStep(controller_, h, 1, error, &next), "SUNAdaptController_EstimateStep");
          if (error > 1.0)
          {
            if (++rejections == MAX_REJECTIONS)
            {
              throw std::runtime_error("SplittingStep: coupling tolerance not met at t = " + std::to_string(t_start)
                                       + "; largest mismatch on " + describe(*worst_));
            }
            // Back to the accepted state, made consistent: the coupling it was
            // accepted with is off by up to the tolerance, which no shorter
            // retry could remove.
            for (auto& partition : partitions_)
            {
              std::copy(partition->y_start.begin(), partition->y_start.end(), N_VGetArrayPointer(blocks_[static_cast<std::size_t>(partition->index)]));
              partition->time = t_start;
            }
            consistentCoupling();
            step_ = std::clamp(SAFETY * next, ETAMIN * h, ETAMXF * h);
            ++rejected_steps_;
            continue;
          }
          check(SUNAdaptController_UpdateH(controller_, h, error), "SUNAdaptController_UpdateH");
          // Bound growth from the intended step: one shortened to reach tout does not limit the next.
          step_      = std::clamp(SAFETY * next, ETAMIN * step_, GROWTH * step_);
          rejections = 0;
        }
        t_ = end;
        for (auto* link : links_)
        {
          const RealT current = value(*link, y_);
          link->rate          = (current - link->start) / h;
          link->start         = current;
        }
        ++steps_;
        emitOutputs();
      }
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class SplittingStep<sunrealtype, long int>;
    template class SplittingStep<sunrealtype, int>;
    template class SplittingStep<sunrealtype, size_t>;
  } // namespace Sundials
} // namespace AnalysisManager
