#include "SplittingStep.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

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
      for (auto& partition : partitions_)
      {
        // SUNStepper_Destroy dereferences the stepper unchecked.
        if (partition->block)
        {
          SUNStepper_Destroy(&partition->block);
        }
        if (partition->stepper)
        {
          SUNStepper_Destroy(&partition->stepper);
        }
      }
      SUNContext_Free(&context_);
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::addPartition(SolverT& solver, std::vector<CouplingT> couplings)
    {
      if (arkode_)
      {
        throw std::logic_error("SplittingStep: partitions are added before configureSimulation");
      }
      // Owned before the steppers exist, so the destructor releases them if a call below throws.
      auto& partition     = *partitions_.emplace_back(std::make_unique<Partition>());
      partition.solver    = &solver;
      partition.index     = static_cast<sunindextype>(partitions_.size() - 1);
      partition.couplings = std::move(couplings);
      partition.stepper   = solver.createSUNStepper();
      check(SUNStepper_Create(context_, &partition.block), "SUNStepper_Create");
      check(SUNStepper_SetContent(partition.block, &partition), "SUNStepper_SetContent");
      check(SUNStepper_SetResetFn(partition.block, resetBlock), "SUNStepper_SetResetFn");
      check(SUNStepper_SetEvolveFn(partition.block, evolveBlock), "SUNStepper_SetEvolveFn");
      check(SUNStepper_SetStopTimeFn(partition.block, setBlockStopTime), "SUNStepper_SetStopTimeFn");
      check(SUNStepper_SetStepDirectionFn(partition.block, setBlockStepDirection), "SUNStepper_SetStepDirectionFn");
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
    typename SplittingStep<ScalarT, IdxT>::Partition& SplittingStep<ScalarT, IdxT>::content(SUNStepper block)
    {
      void* content = nullptr;
      SUNStepper_GetContent(block, &content);
      return *static_cast<Partition*>(content);
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

    // ARKODE passes the same vector to Reset and Evolve, and a partition only
    // writes its own block, so the other blocks need no saving.
    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::resetBlock(SUNStepper block, sunrealtype t, N_Vector y)
    {
      const auto& partition = content(block);
      couple(partition, y, t);
      return SUNStepper_Reset(partition.stepper, t, N_VGetSubvector_ManyVector(y, partition.index));
    }

    template <class ScalarT, typename IdxT>
    int SplittingStep<ScalarT, IdxT>::evolveBlock(SUNStepper block, sunrealtype tout, N_Vector y, sunrealtype* tret)
    {
      auto&     partition = content(block);
      const int flag      = SUNStepper_Evolve(partition.stepper, tout, N_VGetSubvector_ManyVector(y, partition.index), tret);
      partition.time      = *tret;
      return flag;
    }

    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::setBlockStopTime(SUNStepper block, sunrealtype tstop)
    {
      return SUNStepper_SetStopTime(content(block).stepper, tstop);
    }

    template <class ScalarT, typename IdxT>
    SUNErrCode SplittingStep<ScalarT, IdxT>::setBlockStepDirection(SUNStepper block, sunrealtype direction)
    {
      return SUNStepper_SetStepDirection(content(block).stepper, direction);
    }

    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::configureSimulation()
    {
      if (arkode_)
      {
        throw std::logic_error("SplittingStep: already configured");
      }
      if (!(first_step_ > 0.0))
      {
        throw std::invalid_argument("SplittingStep: the step must be positive");
      }

      std::map<const EvaluatorT*, const Partition*> partition_of;
      std::vector<SUNStepper>                       steppers;
      for (const auto& partition : partitions_)
      {
        auto& model = *partition->solver->getModel();
        if (!partition_of.emplace(&model, partition.get()).second)
        {
          throw std::invalid_argument("SplittingStep: a model is added twice");
        }
        auto block = allocated(N_VNew_Serial(static_cast<sunindextype>(model.size()), context_), "N_VNew_Serial");
        blocks_.push_back(block);
        std::copy_n(model.y().getData(), model.size(), N_VGetArrayPointer(block));
        steppers.push_back(partition->block);
      }

      std::set<InputT*> inputs;
      for (const auto& partition : partitions_)
      {
        for (const auto& coupling : partition->couplings)
        {
          const auto source = partition_of.find(coupling.source);
          if (source == partition_of.end()
              || static_cast<sunindextype>(coupling.index) >= N_VGetLength(blocks_[static_cast<std::size_t>(source->second->index)])
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

      y_       = allocated(N_VNew_ManyVector(static_cast<sunindextype>(blocks_.size()), blocks_.data(), context_), "N_VNew_ManyVector");
      y_saved_ = allocated(N_VClone(y_), "N_VClone");
      // Lie–Trotter is SplittingStep's default; higher orders need negative substeps, which DAE partitions reject.
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
    SplittingStep<ScalarT, IdxT>::couplingMismatch(N_Vector y, RealT rel_tol, RealT abs_tol) const
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
          const RealT change = std::abs(link.input->at(t_) - source) / weight;
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
        for (const auto& partition : partitions_)
        {
          couple(*partition, y_, t0);
          auto& model = *partition->solver->getModel();
          auto* block = N_VGetSubvectorArrayPointer_ManyVector(y_, partition->index);
          std::copy_n(block, model.size(), model.y().getData());
          model.y().setDataUpdated();
          check(partition->solver->computeConsistentState(t0, t0 + first_step_), "computeConsistentState");
          std::copy_n(model.y().getData(), model.size(), block);
        }
        mismatch = couplingMismatch(y_, rel_tol_, abs_tol_);
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

    /// Step to exactly tout, one splitting step per ARKODE call.
    template <class ScalarT, typename IdxT>
    void SplittingStep<ScalarT, IdxT>::advance(RealT tout)
    {
      int rejections = 0;
      while (t_ < tout)
      {
        const RealT h       = std::min(step_, tout - t_);
        const RealT t_start = t_;
        N_VScale(1.0, y_, y_saved_);
        check(ARKodeSetFixedStep(arkode_, h), "ARKodeSetFixedStep");
        check(ARKodeSetStopTime(arkode_, tout), "ARKodeSetStopTime");
        check(ARKodeEvolve(arkode_, tout, y_, &t_, ARK_ONE_STEP), "ARKodeEvolve");
        if (!(t_ > t_start))
        {
          throw std::runtime_error("SplittingStep: no progress at t = " + std::to_string(t_start));
        }

        if (controller_)
        {
          // How far each partition's inputs ended from their sources; O(h^2) once rates are known.
          const RealT error = couplingMismatch(y_, coupling_tol_, coupling_tol_);
          RealT       next  = h;
          check(SUNAdaptController_EstimateStep(controller_, h, 1, error, &next), "SUNAdaptController_EstimateStep");
          if (error > 1.0)
          {
            if (++rejections == MAX_REJECTIONS)
            {
              throw std::runtime_error("SplittingStep: coupling tolerance not met at t = " + std::to_string(t_start));
            }
            N_VScale(1.0, y_saved_, y_);
            t_ = t_start;
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
