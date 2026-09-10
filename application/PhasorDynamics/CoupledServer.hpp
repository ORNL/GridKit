#pragma once

#include <GridKit/Definitions.hpp>

#ifdef GRIDKIT_ENABLE_ZMQ

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include "CoSimulationTypes.hpp"
#include "SimulationOptions.hpp"
#include "SimulationRunner.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Implementation of a server coupled through co-simulation.
    class CoupledServer
    {
      const StudyData&                   base_;
      const SimulationOptions&           options_;
      StudyData                          data_;
      std::unique_ptr<CoupledSimulation> simulation_;
      CoSimulationMessage                outgoing_;
      std::uint64_t                      session_{0};
      Testing::TestStatus                status_ = true;

      void beginStudy(const CoSimulationMessage& message)
      {
        if (simulation_ || message.session <= session_ || message.time != 0.0
            || !sameTime(message.final_time, base_.tmax)
            || !sameTime(message.interval, options_.interval))
        {
          throw std::runtime_error(
            "invalid co-simulation initialization or incompatible study timing"
          );
        }

        session_ = message.session;
        data_    = base_;

        if (options_.serve)
        {
          suffixOutputs(data_, "_session_" + std::to_string(session_));
        }

        simulation_ = std::make_unique<CoupledSimulation>(
          data_, options_.inputs, options_.outputs, &message.signals
        );
        simulation_->start();
        outgoing_.signals    = simulation_->outputs();
        outgoing_.final_time = data_.tmax;
        outgoing_.interval   = options_.interval;
        for (const auto& event : data_.events)
        {
          outgoing_.event_times.push_back(event.time);
        }
      }

      bool finishStudy(const CoSimulationMessage& message)
      {
        if (!sameTime(message.time, data_.tmax)
            || !sameTime(simulation_->time(), data_.tmax))
        {
          throw std::runtime_error("co-simulation ended before the final time");
        }

        status_ *= simulation_->finish(false);
        simulation_.reset();
        return options_.serve;
      }

    public:
      CoupledServer(const StudyData& base, const SimulationOptions& options)
        : base_(base), options_(options)
      {
      }

      CoSimulationMessage extract() const
      {
        return outgoing_;
      }

      Testing::TestStatus status() const
      {
        return status_;
      }

      /// Whether no simulations are running.
      bool finished() const noexcept
      {
        return !simulation_;
      }

      bool handleRequest(const CoSimulationMessage& message)
      {
        outgoing_ =
          CoSimulationMessage(message.type, message.session, message.time);
        if (message.type == CoSim::MessageType::INIT)
        {
          beginStudy(message);
          return true;
        }

        if (message.type == CoSim::MessageType::SHUTDOWN)
        {
          if (simulation_ || message.session != session_)
          {
            throw std::runtime_error(
              "shutdown requires a completed co-simulation session"
            );
          }
          return false;
        }

        if (!simulation_ || message.session != session_)
        {
          throw std::runtime_error(
            "co-simulation message has no matching active session"
          );
        }

        switch (message.type)
        {
        case CoSim::MessageType::STEP:
          if (message.time - simulation_->time() > options_.interval
              && !sameTime(
                message.time - simulation_->time(), options_.interval
              ))
          {
            throw std::runtime_error(
              "co-simulation step exceeds the exchange interval"
            );
          }

          simulation_->setInputs(message.signals);
          simulation_->advance(message.time);
          outgoing_.signals = simulation_->outputs();
          break;
        case CoSim::MessageType::UPDATE:
        case CoSim::MessageType::COMMIT:
          if (!sameTime(message.time, simulation_->time()))
          {
            throw std::runtime_error(
              "interface update has the wrong simulation time"
            );
          }

          if (message.type == CoSim::MessageType::UPDATE)
          {
            simulation_->setInputs(message.signals);
            simulation_->refresh();
            outgoing_.signals = simulation_->outputs();
          }
          else
          {
            simulation_->commit();
          }
          break;
        case CoSim::MessageType::END:
          return finishStudy(message);
        default:
          throw std::runtime_error("unexpected co-simulation request");
        }
        return true;
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
#endif
