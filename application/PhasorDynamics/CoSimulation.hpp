#pragma once

#include "SimulationOptions.hpp"
#include "SimulationRunner.hpp"

#ifdef GRIDKIT_ENABLE_ZMQ

#include <optional>
#include <stdexcept>
#include <utility>

#include <GridKit/CoSim/CoSim.hpp>

#include "CoSimulationTypes.hpp"
#include "CoupledClient.hpp"
#include "CoupledServer.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    inline CoSim::RunResult runCoupledClientStudies(
      const SimulationOptions& options,
      NextCoupledStudy         next,
      CompleteCoupledStudy     complete,
      bool                     shutdown_peer,
      bool                     print_results = false
    )
    {
      CoupledClient client(
        options,
        std::move(next),
        std::move(complete),
        shutdown_peer,
        print_results
      );

      CoSim::CoSim<double, CoSim::Disposition::CLIENT> adapter(
        options.endpoint,
        CoSim::MessageStepFunction<double>(
          [&](
            const CoSimulationMessage&                    reply,
            const CoSim::MessageExchangeFunction<double>& exchange
          )
          { return client.handleReply(reply, exchange); }
        ),
        CoSim::MessageExtractFunction<double>([&](const auto&)
                                              { return client.extract(); }),
        options.timeout_ms
      );

      const auto result = adapter.run();

      if (result == CoSim::RunResult::COMPLETED && !client.finished())
      {
        throw std::runtime_error(
          "co-simulation peer ended before the study completed"
        );
      }

      return result;
    }

    inline Testing::TestStatus runCoupledServer(
      const StudyData& base, const SimulationOptions& options
    )
    {
      CoupledServer server(base, options);

      CoSim::CoSim<double, CoSim::Disposition::SERVER> adapter(
        options.endpoint,
        CoSim::MessageStepFunction<double>(
          [&](const CoSimulationMessage& message, const auto&)
          { return server.handleRequest(message); }
        ),
        CoSim::MessageExtractFunction<double>([&](const auto&)
                                              { return server.extract(); }),
        options.timeout_ms
      );

      if (adapter.run() == CoSim::RunResult::STOPPED)
      {
        return false;
      }

      if (!server.finished())
      {
        throw std::runtime_error(
          "co-simulation peer closed before the study completed"
        );
      }

      return server.status();
    }

    inline Testing::TestStatus runCoupledStudy(
      const StudyData& data, const SimulationOptions& options
    )
    {
      if (options.role == "server")
      {
        return runCoupledServer(data, options);
      }

      bool                first  = true;
      Testing::TestStatus status = true;
      const auto          next   = [&]() -> std::optional<StudyData>
      {
        if (!first)
        {
          return std::nullopt;
        }

        first = false;
        return data;
      };

      if (runCoupledClientStudies(
            options,
            next,
            [&](Testing::TestStatus result)
            { status *= result; },
            false,
            true
          )
          == CoSim::RunResult::STOPPED)
      {
        status = false;
      }

      return status;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
#endif
