#pragma once

#include <functional>

#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "AnalysisUtilities.hpp"

namespace GridKit::EMT
{
  /// Validated study actions, applied together at each event time.
  template <typename ScalarT, typename IdxT>
  class EventSchedule
  {
    using System = SystemModel<ScalarT, IdxT>;
    using Ida    = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;
    using RealT  = typename System::RealT;
    using Stats  = AnalysisManager::Sundials::IdaStats;

  public:
    EventSchedule(System& system, const StudyData& study)
      : system_(system), study_(study)
    {
      validateEventTimes(study_.events, study_.tmax);
      for (const auto& event : study_.events)
      {
        if (const auto* command = std::get_if<SwitchEvent>(&event.action))
          system_.getSwitch(command->element_id);
        else
        {
          const auto& step = std::get<SignalStep>(event.action);
          if (!std::isfinite(step.value) || !system_.signal(step.signal_id).constant())
            throw std::invalid_argument("signal_step requires a declared constant and finite value: " + step.signal_id);
        }
      }
    }

    /// Apply time-zero actions after model initialization, before IDA sees it.
    void configure(Ida& ida)
    {
      if (!study_.events.empty() && study_.events.front().time == 0.0)
        applyGroup();
      ida.configureSimulation();
    }

    Stats run(Ida& ida, std::function<void(RealT)> record = {})
    {
      const auto monitor = [&](RealT time)
      { if (record) record(time); };
      Stats total;
      RealT time = 0.0;
      ida.initializeSimulation(time);
      monitor(time);
      while (next_ < study_.events.size())
      {
        time = static_cast<RealT>(study_.events[next_].time);
        ida.runSimulation(time, study_.dt_monitor, monitor);
        total += ida.getStats();

        if (applyGroup())
        {
          system_.updateTime(time, RealT{1});
          if (system_.evaluateResidual() != 0 || system_.evaluateJacobian() != 0)
            throw std::runtime_error("EMT event model assembly failed");
          ida.configureLinearSolver();
        }
        // A discontinuity changes algebraic values and derivatives, not states.
        ida.setConsistentICType(AnalysisManager::Sundials::IdaConsistentICType::YA_YDP);
        ida.initializeSimulation(time);
        monitor(time);
      }
      if (study_.tmax > time)
        ida.runSimulation(study_.tmax, study_.dt_monitor, monitor);
      total += ida.getStats();
      return total;
    }

  private:
    bool applyGroup()
    {
      const double time             = study_.events[next_].time;
      bool         topology_changed = false;
      do
      {
        const auto& action = study_.events[next_++].action;
        if (const auto* command = std::get_if<SwitchEvent>(&action))
        {
          system_.getSwitch(command->element_id)->setOpen(command->open);
          topology_changed = true;
        }
        else
        {
          const auto& step = std::get<SignalStep>(action);
          system_.signal(step.signal_id).setConstantValue(step.value);
        }
      } while (next_ < study_.events.size() && study_.events[next_].time == time);
      if (topology_changed)
        system_.resetJacobianStructure();
      return topology_changed;
    }

    System&          system_;
    const StudyData& study_;
    size_t           next_{0};
  };
} // namespace GridKit::EMT
