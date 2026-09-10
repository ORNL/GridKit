#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>

#include "AnalysisUtilities.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {

    using SimulationSystem = SystemModel<double, std::size_t>;
    using SimulationSolver =
      AnalysisManager::Sundials::Ida<double, std::size_t>;
    using Signals        = std::map<std::string, double>;
    using SignalBindings = std::map<std::string, std::size_t>;

    inline void configureSolver(SimulationSolver& ida, const StudyData& data)
    {
      ida.setTolerance(data.rel_tol, data.abs_tol);
      ida.setFixedStep(data.dt_fixed);
      ida.setMaxSteps(data.max_steps);
      ida.setMaxOrder(data.max_order);
      ida.setConsistentICType(data.consistent_ic_type);
      ida.configureSimulation();
    }

    inline void validateStudy(const StudyData& data)
    {
      if (!std::isfinite(data.tmax) || data.tmax <= 0.0
          || !std::isfinite(data.dt_monitor) || data.dt_monitor < 0.0)
      {
        throw std::invalid_argument(
          "Study requires a positive tmax and a nonnegative monitor interval"
        );
      }

      double previous = 0.0;
      for (const auto& event : data.events)
      {
        if (!std::isfinite(event.time) || event.time < previous
            || event.time > data.tmax)
        {
          throw std::invalid_argument(
            "Events must be ordered and within the simulation interval"
          );
        }

        previous = event.time;

        if ((event.type == SystemEvent::Type::FAULT_ON
             || event.type == SystemEvent::Type::FAULT_OFF)
            && event.bus == INVALID_INDEX<std::size_t>)
        {
          throw std::invalid_argument("Fault events require a bus number");
        }
      }
    }

    inline void applyEvent(SimulationSystem& sys, const SystemEvent& event)
    {
      int ret = 0;
      switch (event.type)
      {
      case SystemEvent::Type::FAULT_ON:
        ret = sys.getBus(event.bus)->setFault(true, event.R, event.X);
        break;
      case SystemEvent::Type::FAULT_OFF:
        ret = sys.getBus(event.bus)->clearFault();
        break;
      case SystemEvent::Type::IN_SERVICE:
        sys.getBranch(event.device)->setInService(true);
        break;
      case SystemEvent::Type::OUT_OF_SERVICE:
        sys.getBranch(event.device)->setInService(false);
        break;
      }
      if (ret != 0)
      {
        throw std::runtime_error(
          "Event failed at time " + std::to_string(event.time)
        );
      }
    }

    /// Run a standalone study with event handling and monitoring.
    inline Testing::TestStatus runStandaloneStudy(
      const StudyData& data,
      bool             consistent_initial = true,
      bool             print_results      = true
    )
    {
      validateStudy(data);
      SimulationSystem sys(data.model_data);
      sys.allocate();
      SimulationSolver ida(&sys);
      configureSolver(ida, data);
      ida.initializeSimulation(0.0, consistent_initial);
      for (const auto& event : data.events)
      {
        ida.runSimulation(event.time, data.dt_monitor);
        applyEvent(sys, event);
        ida.initializeSimulation(event.time);
      }
      ida.runSimulation(data.tmax, data.dt_monitor);
      sys.stopMonitor();
      return checkErrors(data, print_results);
    }

    inline bool sameTime(double a, double b)
    {
      return std::abs(a - b) <= 32.0 * std::numeric_limits<double>::epsilon()
                                  * std::max({1.0, std::abs(a), std::abs(b)});
    }

    inline void suffixOutputs(StudyData& data, const std::string& suffix)
    {
      const auto rename = [&](const std::filesystem::path& path)
      {
        return path.parent_path()
               / (path.stem().string() + suffix + path.extension().string());
      };

      if (!data.output_file.empty())
      {
        data.output_file = rename(data.output_file);
      }

      for (auto& sink : data.model_data.monitor_sink)
      {
        sink.file_name = rename(sink.file_name).string();
        if (!data.output_file.empty()
            && sink.format == Model::VariableMonitorFormat::CSV
            && sink.delim == ",")
        {
          data.output_file = sink.file_name;
        }
      }
    }

    /// Simulation partition with externally supplied signal inputs.
    class CoupledSimulation
    {
      struct Input
      {
        double      value{0.0};
        std::size_t index{INVALID_INDEX<std::size_t>};
      };

      const StudyData&                  data_;
      SignalBindings                    outputs_;
      std::map<std::string, Input>      inputs_;
      SimulationSystem                  sys_;
      std::unique_ptr<SimulationSolver> ida_;
      double                            time_{0.0};
      std::size_t                       event_{0};
      std::size_t                       monitor_{1};
      bool                              pending_monitor_{false};

      bool applyEvents()
      {
        bool changed = false;
        while (event_ < data_.events.size()
               && sameTime(data_.events[event_].time, time_))
        {
          applyEvent(sys_, data_.events[event_++]);
          changed = true;
        }
        return changed;
      }

    public:
      CoupledSimulation(
        const StudyData&      data,
        const SignalBindings& inputs,
        const SignalBindings& outputs,
        const Signals*        initial_inputs = nullptr
      )
        : data_(data), outputs_(outputs), sys_(data.model_data)
      {
        validateStudy(data);
        if (data.consistent_ic_type
            != AnalysisManager::Sundials::IdaConsistentICType::YA_YDP)
        {
          throw std::invalid_argument(
            "co-simulation requires consistent_ic_type=ya_ydp to preserve "
            "differential states during interface solves"
          );
        }

        for (const auto& [name, id] : inputs)
        {
          auto* signal = sys_.getSignalNode(id);
          if (signal->assigned() || signal->linked())
          {
            throw std::invalid_argument(
              "co-simulation input has a model-owned source: " + name
            );
          }

          auto& input = inputs_[name];
          for (const auto& bus : data.model_data.bus)
          {
            // TODO: this doesn't seem like the cleanest way to initialize the
            //       voltages
            for (auto component : {BusSignalInputs::vr, BusSignalInputs::vi})
            {
              const auto found = bus.signal_inputs.find(component);
              if (found != bus.signal_inputs.end() && found->second == id)
              {
                input.value =
                  component == BusSignalInputs::vr ? bus.Vr0 : bus.Vi0;
              }
            }
          }
          signal->setAssigned();
          signal->link(&input.value, &input.index);
        }
        if (initial_inputs)
        {
          setInputs(*initial_inputs);
        }
        sys_.allocate();
        for (const auto& [name, id] : outputs_)
        {
          if (!sys_.getSignalNode(id)->linked())
          {
            throw std::invalid_argument(
              "co-simulation output has no source: " + name
            );
          }
        }
        ida_ = std::make_unique<SimulationSolver>(&sys_);
        configureSolver(*ida_, data);
      }

      void setInputs(const Signals& signals)
      {
        if (signals.size() != inputs_.size())
        {
          throw std::invalid_argument(
            "peer signal names do not match the configured inputs"
          );
        }

        for (auto& [name, input] : inputs_)
        {
          const auto found = signals.find(name);

          if (found == signals.end() || !std::isfinite(found->second))
          {
            throw std::invalid_argument(
              "missing or non-finite input signal: " + name
            );
          }

          input.value = found->second;
        }
      }

      Signals outputs()
      {
        // refresh current outputs accumulated during residual evaluation
        sys_.evaluateResidual();

        Signals values;
        for (const auto& [name, id] : outputs_)
        {
          values.emplace(name, sys_.getSignalNode(id)->read());
        }

        return values;
      }

      void start()
      {
        applyEvents();
        ida_->initializeSimulation(0.0, true, false);
        pending_monitor_ = true;
      }

      void refresh()
      {
        ida_->initializeSimulation(time_, true, false);
      }

      void commit()
      {
        if (pending_monitor_ && sys_.monitoring())
        {
          sys_.evaluateResidual();
          sys_.printMonitoredVariables();
        }
        pending_monitor_ = false;
      }

      double time() const
      {
        return time_;
      }

      void advance(double target)
      {
        if (!std::isfinite(target) || target <= time_ || target > data_.tmax)
        {
          throw std::invalid_argument("Invalid co-simulation step time");
        }

        // restart after an external input change
        ida_->initializeSimulation(time_, true, false);
        while (time_ < target && !sameTime(time_, target))
        {
          double next = target;
          if (event_ < data_.events.size())
          {
            next = std::min(next, data_.events[event_].time);
          }

          const double monitor_time =
            data_.dt_monitor > 0.0
              ? static_cast<double>(monitor_) * data_.dt_monitor
              : data_.tmax;
          next = std::min(next, monitor_time);
          const bool print =
            sameTime(next, monitor_time) || sameTime(next, data_.tmax);
          const bool at_target = sameTime(next, target);
          ida_->runSimulationStep(next, 0.0, 1, 1, {}, print && !at_target);
          pending_monitor_ = print && at_target;
          time_            = next;
          if (sameTime(time_, monitor_time))
          {
            monitor_++;
          }
          if (applyEvents())
          {
            ida_->initializeSimulation(time_, true, false);
          }
        }
        time_ = target;
      }

      Testing::TestStatus finish(bool print_results)
      {
        sys_.stopMonitor();
        return checkErrors(data_, print_results);
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
