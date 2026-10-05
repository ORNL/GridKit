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

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#include <GridKit/Definitions.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include <arkode/arkode.h>
#include <arkode/arkode_splittingstep.h>

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
      constexpr int MAX_SWEEPS = 100;

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
      ARKodeFree(&arkode_);
      SUNAdaptController_Destroy(controller_);
      N_VDestroy(y_saved_);
      N_VDestroy(y_);
      for (auto block : blocks_)
      {
        N_VDestroy(block);
      }
      // SUNStepper_Destroy dereferences the stepper unchecked.
      for (auto& color : colors_)
      {
        if (color->stepper)
        {
          SUNStepper_Destroy(&color->stepper);
        }
      }
      for (auto& partition : partitions_)
      {
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
      // Owned before the steppers exist, so the destructor releases them if a call below throws.
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
    typename SplittingStep<ScalarT, IdxT>::Color& SplittingStep<ScalarT, IdxT>::content(SUNStepper stepper)
    {
      void* content = nullptr;
      SUNStepper_GetContent(stepper, &content);
      return *static_cast<Color*>(content);
    }

    template <class ScalarT, typename IdxT>
    typename SplittingStep<ScalarT, IdxT>::RealT SplittingStep<ScalarT, IdxT>::value(const Link& link, N_Vector y)
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

    template <class ScalarT, typename IdxT>
    int SplittingStep<ScalarT, IdxT>::evolvePartition(Partition& partition, sunrealtype tout, N_Vector y, sunrealtype* tret)
    {
      const auto start   = std::chrono::steady_clock::now();
      const int  flag    = SUNStepper_Evolve(partition.stepper, tout, y, tret);
      partition.seconds += std::chrono::duration<RealT>(std::chrono::steady_clock::now() - start).count();
      if (flag == SUN_SUCCESS)
        partition.time = *tret;
      return flag;
    }

    /// Runs `function` on every partition of a color at once. Failures are
    /// returned, not thrown, since ARKODE calls this through C.
    template <class ScalarT, typename IdxT>
    template <class Function>
    std::exception_ptr SplittingStep<ScalarT, IdxT>::concurrently(const Color& color, Function&& function)
    {
      const auto                      count = static_cast<std::ptrdiff_t>(color.partitions.size());
      std::vector<std::exception_ptr> failures(color.partitions.size());
#ifdef GRIDKIT_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(color.threads) if (color.threads > 1)
#endif
      for (std::ptrdiff_t i = 0; i < count; ++i)
      {
        Log::ScopedOutput output;
        try
        {
          function(*color.partitions[static_cast<std::size_t>(i)]);
        }
        catch (...)
        {
          failures[static_cast<std::size_t>(i)] = std::current_exception();
        }
      }
      for (auto& failure : failures)
      {
        if (failure)
        {
          return failure;
        }
      }
      return nullptr;
    }

    // ARKODE passes the same vector to Reset and Evolve, and a partition only
    // writes its own block, so the other blocks need no saving.
    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::resetColor(SUNStepper stepper, sunrealtype t, N_Vector y)
    {
      for (const auto* partition : content(stepper).partitions)
      {
        couple(*partition, y, t);
        const SUNErrCode err = SUNStepper_Reset(partition->stepper, t, N_VGetSubvector_ManyVector(y, partition->index));
        if (err != SUN_SUCCESS)
        {
          return err;
        }
      }
      return SUN_SUCCESS;
    }

    template <class ScalarT, typename IdxT>
    int SplittingStep<ScalarT, IdxT>::evolveColor(SUNStepper stepper, sunrealtype tout, N_Vector y, sunrealtype* tret)
    {
      auto& color   = content(stepper);
      color.failure = concurrently(color, [&](Partition& partition)
                                   {
        RealT reached = tout;
        check(evolvePartition(partition, tout, N_VGetSubvector_ManyVector(y, partition.index), &reached), "SUNStepper_Evolve"); });

      *tret = tout;
      if (color.failure)
      {
        return SUN_ERR_OP_FAIL;
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::setColorStopTime(SUNStepper stepper, sunrealtype tstop)
    {
      for (const auto* partition : content(stepper).partitions)
      {
        const SUNErrCode err = SUNStepper_SetStopTime(partition->stepper, tstop);
        if (err != SUN_SUCCESS)
        {
          return err;
        }
      }
      return SUN_SUCCESS;
    }

    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::setColorStepDirection(SUNStepper stepper, sunrealtype direction)
    {
      for (const auto* partition : content(stepper).partitions)
      {
        const SUNErrCode err = SUNStepper_SetStepDirection(partition->stepper, direction);
        if (err != SUN_SUCCESS)
        {
          return err;
        }
      }
      return SUN_SUCCESS;
    }

    /**
     * @brief Colors partitions so that none shares a color with one it reads or
     * that reads it.
     *
     * Greedy first fit, most-coupled partitions first (ties in registration
     * order). The coloring does not depend on the thread count.
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

      std::vector<const Color*> color_of(count, nullptr);
      for (const auto p : order)
      {
        auto color = std::find_if(colors_.begin(), colors_.end(), [&](const auto& candidate)
                                  { return std::none_of(neighbors[p].begin(), neighbors[p].end(), [&](std::size_t q)
                                                        { return color_of[q] == candidate.get(); }); });
        if (color == colors_.end())
        {
          color = colors_.insert(colors_.end(), std::make_unique<Color>());
        }
        color_of[p] = color->get();
        (*color)->partitions.push_back(partitions_[p].get());
      }

      for (const auto& color : colors_)
      {
        color->threads = static_cast<int>(std::min(color->partitions.size(), static_cast<std::size_t>(num_threads_)));
        check(SUNStepper_Create(context_, &color->stepper), "SUNStepper_Create");
        check(SUNStepper_SetContent(color->stepper, color.get()), "SUNStepper_SetContent");
        check(SUNStepper_SetResetFn(color->stepper, resetColor), "SUNStepper_SetResetFn");
        check(SUNStepper_SetEvolveFn(color->stepper, evolveColor), "SUNStepper_SetEvolveFn");
        check(SUNStepper_SetStopTimeFn(color->stepper, setColorStopTime), "SUNStepper_SetStopTimeFn");
        check(SUNStepper_SetStepDirectionFn(color->stepper, setColorStepDirection), "SUNStepper_SetStepDirectionFn");
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

      colorPartitions();
      if (colors_.size() < 2)
      {
        throw std::invalid_argument("SplittingStep: partitions must be coupled");
      }
      std::vector<SUNStepper> steppers;
      for (const auto& color : colors_)
      {
        steppers.push_back(color->stepper);
      }

      y_       = allocated(N_VNew_ManyVector(static_cast<sunindextype>(blocks_.size()), blocks_.data(), context_), "N_VNew_ManyVector");
      y_saved_ = allocated(N_VClone(y_), "N_VClone");
      // Lie-Trotter over the colors; higher orders need negative substeps, which DAE partitions reject.
      arkode_  = allocated(SplittingStepCreate(steppers.data(), static_cast<int>(steppers.size()), 0.0, y_, context_), "SplittingStepCreate");
      // Outputs come from accepted step ends, never interpolated states.
      check(ARKodeSetInterpolantType(arkode_, ARK_INTERP_NONE), "ARKodeSetInterpolantType");
      if (coupling_tol_ > 0.0)
      {
        controller_ = allocated(SUNAdaptController_PI(context_), "SUNAdaptController_PI");
      }
    }

    /// Largest change any input would see if re-read now, in units of
    /// abs_tol + rel_tol |value| (abs_tol <= 0: each input's own).
    template <class ScalarT, typename IdxT>
    typename SplittingStep<ScalarT, IdxT>::RealT
    SplittingStep<ScalarT, IdxT>::couplingMismatch(N_Vector y, RealT t, RealT rel_tol, RealT abs_tol) const
    {
      RealT mismatch = 0.0;
      for (const auto& partition : partitions_)
      {
        for (const auto& link : partition->links)
        {
          const RealT source = value(link, y);
          RealT       weight = abs_tol + rel_tol * std::abs(source);
          if (abs_tol <= 0.0)
          {
            weight = link.abs_tol + rel_tol * std::abs(source);
          }
          const RealT change = std::abs(link.input->at(t) - source) / weight;
          if (!std::isfinite(change))
          {
            throw std::runtime_error("SplittingStep: coupling value is not finite");
          }
          mismatch = std::max(mismatch, change);
        }
      }
      return mismatch;
    }

    /**
     * @brief Consistent coupling at t0, differential states fixed.
     *
     * Block Gauss-Seidel until every partition was solved with the coupling
     * values it now has.
     */
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

      RealT mismatch = std::numeric_limits<RealT>::infinity();
      for (int sweep = 0; sweep < MAX_SWEEPS && mismatch > 1.0; ++sweep)
      {
        for (const auto& color : colors_)
        {
          for (const auto* partition : color->partitions)
          {
            couple(*partition, y_, t0);
          }
          const auto failure = concurrently(*color, [&](Partition& partition)
                                            {
            auto& model = *partition.solver->getModel();
            auto* block = N_VGetSubvectorArrayPointer_ManyVector(y_, partition.index);
            std::copy_n(block, model.size(), model.y().getData());
            model.y().setDataUpdated();
            check(partition.solver->computeConsistentState(t0, t0 + first_step_), "computeConsistentState");
            std::copy_n(model.y().getData(), model.size(), block); });
          if (failure)
          {
            std::rethrow_exception(failure);
          }
        }
        mismatch = couplingMismatch(y_, t0, rel_tol_, abs_tol_);
      }
      if (mismatch > 1.0)
      {
        throw std::runtime_error("SplittingStep: consistency iteration failed at t = " + std::to_string(t0)
                                 + " (mismatch " + std::to_string(mismatch) + ")");
      }

      for (auto& partition : partitions_)
      {
        for (auto& link : partition->links)
        {
          link.start = value(link, y_);
        }
      }
      check(ARKodeReset(arkode_, t0, y_), "ARKodeReset");
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
      const RealT t0     = t_;
      const int   nsteps = monitorStepCount(t0, tf, dt_monitor);
      for (int i = 1; i <= nsteps; ++i)
      {
        advance(monitorTime(t0, tf, dt_monitor, i, nsteps));
        if (output_)
        {
          output_(t_);
        }
      }
    }

    /// Accept or reject all regional solves together.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::advance(RealT tout)
    {
      int rejections = 0;
      while (t_ < tout)
      {
        const RealT h       = std::min(step_, tout - t_);
        const RealT t_start = t_;
        if (!std::isfinite(h) || !(t_start + h > t_start))
          throw std::runtime_error("SplittingStep: step cannot advance time");
        N_VScale(1.0, y_, y_saved_);
        RealT trial_time = t_start + h;
        check(ARKodeSetFixedStep(arkode_, h), "ARKodeSetFixedStep");
        check(ARKodeSetStopTime(arkode_, tout), "ARKodeSetStopTime");
        const int flag = ARKodeEvolve(arkode_, tout, y_, &trial_time, ARK_ONE_STEP);
        for (auto& color : colors_)
        {
          if (color->failure)
          {
            std::rethrow_exception(std::exchange(color->failure, nullptr));
          }
        }
        check(flag, "ARKodeEvolve");
        if (!(trial_time > t_start))
        {
          throw std::runtime_error("SplittingStep: no progress at t = " + std::to_string(t_start));
        }

        if (controller_)
        {
          // How far each partition's inputs ended from their sources; O(h^2) once rates are known.
          const RealT error = couplingMismatch(y_, trial_time, coupling_tol_, coupling_tol_);
          RealT       next  = h;
          check(SUNAdaptController_EstimateStep(controller_, h, 1, error, &next), "SUNAdaptController_EstimateStep");
          if (error > 1.0)
          {
            if (++rejections == MAX_REJECTIONS)
            {
              throw std::runtime_error("SplittingStep: coupling tolerance not met at t = " + std::to_string(t_start));
            }
            N_VScale(1.0, y_saved_, y_);
            for (auto& partition : partitions_)
            {
              partition->time = t_start;
            }
            check(ARKodeReset(arkode_, t_, y_), "ARKodeReset");
            step_ = std::clamp(SAFETY * next, ETAMIN * h, ETAMXF * h);
            ++rejected_steps_;
            continue;
          }
          check(SUNAdaptController_UpdateH(controller_, h, error), "SUNAdaptController_UpdateH");
          // Bound growth from the intended step: one shortened to reach tout does not limit the next.
          step_      = std::clamp(SAFETY * next, ETAMIN * step_, GROWTH * step_);
          rejections = 0;
        }
        t_ = trial_time;
        for (auto& partition : partitions_)
        {
          for (auto& link : partition->links)
          {
            const RealT current = value(link, y_);
            link.rate           = (current - link.start) / (t_ - t_start);
            link.start          = current;
          }
        }
        ++steps_;
      }
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class SplittingStep<sunrealtype, long int>;
    template class SplittingStep<sunrealtype, int>;
    template class SplittingStep<sunrealtype, size_t>;
  } // namespace Sundials
} // namespace AnalysisManager
