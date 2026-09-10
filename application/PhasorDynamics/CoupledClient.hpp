#pragma once

#include <GridKit/Definitions.hpp>

#ifdef GRIDKIT_ENABLE_ZMQ

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <GridKit/Solver/SteadyState/Kinsol.hpp>

#include "CoSimulationTypes.hpp"
#include "SimulationOptions.hpp"
#include "SimulationRunner.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Implementation of a client coupled through co-simulation.
    class CoupledClient
    {
      using InterfaceSolver =
        AnalysisManager::Sundials::Kinsol<double, std::size_t>;

      const SimulationOptions&           options_;
      NextCoupledStudy                   next_;
      CompleteCoupledStudy               complete_;
      bool                               shutdown_peer_;
      bool                               print_results_;
      std::optional<StudyData>           data_;
      std::unique_ptr<CoupledSimulation> simulation_;
      std::unique_ptr<InterfaceSolver>   interface_solver_;
      CoSimulationMessage                outgoing_;
      std::uint64_t                      session_{0};
      std::vector<double>                boundaries_;
      std::size_t                        boundary_{0};
      std::size_t                        step_{1};

      bool beginStudy()
      {
        data_ = next_();
        if (!data_)
        {
          return false;
        }

        session_++;
        simulation_ = std::make_unique<CoupledSimulation>(
          *data_, options_.inputs, options_.outputs
        );

        outgoing_ = CoSimulationMessage(CoSim::MessageType::INIT, session_);

        outgoing_.signals    = simulation_->outputs();
        outgoing_.final_time = data_->tmax;
        outgoing_.interval   = options_.interval;
        boundary_            = 0;
        step_                = 1;
        return true;
      }

      void startStudy(const CoSimulationMessage& reply)
      {
        if (!sameTime(reply.final_time, data_->tmax)
            || !sameTime(reply.interval, options_.interval))
        {
          throw std::runtime_error("incompatible co-simulation study timing");
        }

        boundaries_ = reply.event_times;
        for (const auto& event : data_->events)
        {
          boundaries_.push_back(event.time);
        }

        for (auto time : boundaries_)
        {
          if (time < 0.0 || time > data_->tmax)
          {
            throw std::runtime_error(
              "peer event time is outside the study interval"
            );
          }
        }

        boundaries_.push_back(data_->tmax);
        std::sort(boundaries_.begin(), boundaries_.end());
        simulation_->setInputs(reply.signals);
        simulation_->start();
      }

      bool finishStudy()
      {
        complete_(simulation_->finish(print_results_));
        simulation_.reset();
        if (beginStudy())
        {
          return true;
        }
        if (!shutdown_peer_)
        {
          return false;
        }

        outgoing_ = CoSimulationMessage(CoSim::MessageType::SHUTDOWN, session_);
        return true;
      }

      /// Solve $F(inputs) - inputs = 0$ at the current time, holding
      /// differential states fixed in both partitions.
      void solveInterface(
        const Signals&                                initial,
        const CoSim::MessageExchangeFunction<double>& exchange
      )
      {
        simulation_->setInputs(initial);
        Signals             trial = initial;
        std::vector<double> values, scales, residual(initial.size());
        bool                evaluated = false;

        for (const auto& [name, value] : initial)
        {
          values.push_back(value);
          scales.push_back(std::max(1.0, std::abs(value)));
        }

        const auto evaluate =
          [&](std::span<const double> guess, std::span<double> result)
        {
          if (evaluated
              && std::equal(
                guess.begin(),
                guess.end(),
                trial.begin(),
                [](double value, const auto& input)
                { return value == input.second; }
              ))
          {
            if (result.data() != residual.data())
            {
              std::copy(residual.begin(), residual.end(), result.begin());
            }
            return 0;
          }
          std::size_t i = 0;

          for (auto& [name, value] : trial)
          {
            value = guess[i++];
          }

          simulation_->setInputs(trial);
          simulation_->refresh();

          auto message = CoSimulationMessage(
            CoSim::MessageType::UPDATE, session_, simulation_->time()
          );

          message.signals  = simulation_->outputs();
          const auto reply = exchange(message);
          if (reply.signals.size() != trial.size())
          {
            throw std::runtime_error(
              "Peer signal names do not match the configured inputs"
            );
          }

          i = 0;
          for (const auto& [name, value] : trial)
          {
            residual[i] = (reply.signals.at(name) - value) / scales[i];
            i++;
          }

          if (result.data() != residual.data())
          {
            std::copy(residual.begin(), residual.end(), result.begin());
          }

          evaluated = true;
          return 0;
        };
        const auto converged = [&]
        {
          for (std::size_t i = 0; i < values.size(); ++i)
          {
            if (std::abs(residual[i]) * scales[i]
                  / std::max(1.0, std::abs(values[i]))
                >= data_->cosim_rel_tol)
            {
              return false;
            }
          }

          return true;
        };

        // check the initial guess and cache its residual for kinsol
        evaluate(values, residual);

        if (!converged())
        {
          if (!interface_solver_)
          {
            interface_solver_ =
              std::make_unique<InterfaceSolver>(values.size());
            interface_solver_->configureSimulation();
            interface_solver_->setJacobianUpdateInterval(1);
          }

          interface_solver_->setFunctionError(data_->cosim_function_error);
          interface_solver_->setTolerance(
            data_->cosim_rel_tol, data_->cosim_step_tol
          );
          interface_solver_->runSimulation(values, evaluate, false);

          evaluate(values, residual);

          // try tightening the tolerance further
          if (!converged())
          {
            interface_solver_->setTolerance(
              std::max(
                0.01 * data_->cosim_rel_tol,
                std::numeric_limits<double>::denorm_min()
              ),
              data_->cosim_step_tol
            );
            interface_solver_->runSimulation(values, evaluate, false);

            evaluate(values, residual);
            if (!converged())
            {
              throw std::runtime_error(
                std::format(
                  "co-simulation interface did not converge at t={}",
                  simulation_->time()
                )
              );
            }
          }
        }

        outgoing_ = CoSimulationMessage(
          CoSim::MessageType::COMMIT, session_, simulation_->time()
        );
      }

      void advance()
      {
        const auto time = simulation_->time();
        while (boundary_ < boundaries_.size())
        {
          const auto boundary = boundaries_[boundary_];
          if (!(boundary < time || sameTime(boundary, time)))
          {
            break;
          }
          ++boundary_;
        }

        const double exchange_time =
          static_cast<double>(step_) * options_.interval;
        const double target =
          std::min(exchange_time, boundaries_.at(boundary_));
        simulation_->advance(target);
        outgoing_ =
          CoSimulationMessage(CoSim::MessageType::STEP, session_, target);

        outgoing_.signals = simulation_->outputs();
        if (sameTime(target, exchange_time))
        {
          step_++;
        }
      }

    public:
      CoupledClient(
        const SimulationOptions& options,
        NextCoupledStudy         next,
        CompleteCoupledStudy     complete,
        bool                     shutdown_peer,
        bool                     print_results
      )
        : options_(options),
          next_(std::move(next)),
          complete_(std::move(complete)),
          shutdown_peer_(shutdown_peer),
          print_results_(print_results)
      {
      }

      CoSimulationMessage extract()
      {
        if (session_ == 0 && !beginStudy())
        {
          outgoing_ = CoSimulationMessage(CoSim::MessageType::SHUTDOWN, 0);
        }
        return outgoing_;
      }

      bool finished() const noexcept
      {
        return simulation_ == nullptr;
      }

      bool handleReply(
        const CoSimulationMessage&                    reply,
        const CoSim::MessageExchangeFunction<double>& exchange
      )
      {
        switch (outgoing_.type)
        {
        case CoSim::MessageType::INIT:
          startStudy(reply);
          solveInterface(reply.signals, exchange);
          return true;
        case CoSim::MessageType::STEP:
          solveInterface(reply.signals, exchange);
          return true;
        case CoSim::MessageType::COMMIT:
          simulation_->commit();
          if (sameTime(simulation_->time(), data_->tmax))
          {
            outgoing_ = CoSimulationMessage(
              CoSim::MessageType::END, session_, data_->tmax
            );
          }
          else
          {
            advance();
          }

          return true;
        case CoSim::MessageType::END:
          return finishStudy();
        case CoSim::MessageType::SHUTDOWN:
          return false;
        default:
          throw std::logic_error("unexpected co-simulation client state");
        }
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
#endif
