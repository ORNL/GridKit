#pragma once

#include <cmath>
#include <format>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#include <GridKit/Definitions.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/String.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Options to a PhasorDynamics simulation application.
    struct SimulationOptions
    {
      /// Study file.
      std::string study_file;

      /// Co-simulation role.
      std::string role;

      /// Co-simulation endpoint.
      std::string endpoint;

      /// Co-simulation exchange interval.
      double interval{0.0};

      /// Co-simulation send/receive timeout.
      int timeout_ms{30000};

      /// Whether to persistently serve when running as a co-simulation server.
      bool serve{false};

      /// Networked signal inputs.
      std::map<std::string, std::size_t> inputs;

      /// Networked signal outputs.
      std::map<std::string, std::size_t> outputs;

      /// Whether the options specified a co-simulation.
      bool coSimulated() const
      {
        return !role.empty();
      }
    };

    inline std::size_t parseSignalId(const std::string& text)
    {
      if (text.empty()
          || text.find_first_not_of("0123456789") != std::string::npos)
      {
        throw std::invalid_argument(
          std::format("expected a nonnegative signal ID: {}", text)
        );
      }

      const auto value = Utilities::parse<std::size_t>(text);

      if (value >= std::numeric_limits<std::size_t>::max())
      {
        throw std::invalid_argument(
          std::format("signal ID is out of range: {}", text)
        );
      }

      return value;
    }

    /// Parse command line options to a simulation application.
    inline SimulationOptions parseSimulationOptions(
      int argc, const char* argv[], bool contingency = false
    )
    {
      using Utilities::ArgType;
      Utilities::CliArgs args{
        {.name       = {"study"},
         .help       = "Path to the study JSON file",
         .required   = true,
         .type       = ArgType::String,
         .positional = true},
        {.name = {"--cosim-role"},
         .help = "Co-simulation role (client, server)",
         .type = ArgType::String},
        {.name = {"--cosim-endpoint"},
         .help = "ZeroMQ endpoint, e.g. tcp://127.0.0.1:8800",
         .type = ArgType::String},
        {.name = {"--cosim-dt"},
         .help = "Positive exchange interval in seconds",
         .type = ArgType::Real},
        {.name       = {"--cosim-input"},
         .help       = "Receive JSON key NAME into signal ID: "
                       "NAME=ID (repeatable)",
         .type       = ArgType::String,
         .repeatable = true},
        {.name       = {"--cosim-output"},
         .help       = "Send signal ID as JSON key NAME: NAME=ID (repeatable)",
         .type       = ArgType::String,
         .repeatable = true},
        {.name = {"--cosim-timeout-ms"},
         .help = "Send/receive timeout in milliseconds (default: 30000)",
         .type = ArgType::Integer},
        {.name = {"--cosim-serve"},
         .help = "Server continues to accept sessions until shutdown",
         .flag = true}
      };
      args.parseArgs(argc, argv);

      SimulationOptions options;
      options.study_file = args.get("study");
      options.serve      = args.get<bool>("cosim-serve");
      if (!args["cosim-role"].empty())
      {
        options.role = args.get("cosim-role");
      }
      if (!args["cosim-endpoint"].empty())
      {
        options.endpoint = args.get("cosim-endpoint");
      }
      if (!args["cosim-dt"].empty())
      {
        options.interval = args.get<double>("cosim-dt");
        if (!std::isfinite(options.interval) || options.interval <= 0.0)
        {
          throw std::invalid_argument(
            "co-simulation exchange interval must be positive and finite"
          );
        }
      }
      if (!args["cosim-timeout-ms"].empty())
      {
        options.timeout_ms = args.get<int>("cosim-timeout-ms");
        if (options.timeout_ms <= 0)
        {
          throw std::invalid_argument("co-simulation timeout must be positive");
        }
      }

      const auto bindings = [&](const std::string& name, auto& mappings)
      {
        const auto& values = args[name];
        for (std::size_t i = 0; i < values.size(); ++i)
        {
          const auto& value = values[i]();
          const auto  equal = value.find('=');
          if (equal == std::string::npos || equal == 0)
          {
            throw std::invalid_argument(
              std::format("signal mapping must be NAME=ID: {}", value)
            );
          }
          if (!mappings
                 .emplace(
                   value.substr(0, equal),
                   parseSignalId(value.substr(equal + 1))
                 )
                 .second)
          {
            throw std::invalid_argument(
              std::format("duplicate signal name: {}", value)
            );
          }
        }
      };
      bindings("cosim-input", options.inputs);
      bindings("cosim-output", options.outputs);

      if (options.study_file.empty())
      {
        throw std::invalid_argument(
          "no study file provided. use --help for usage"
        );
      }
      const bool is_cosim =
        options.serve || !args["cosim-role"].empty()
        || !args["cosim-endpoint"].empty() || !args["cosim-dt"].empty()
        || !args["cosim-input"].empty() || !args["cosim-output"].empty()
        || !args["cosim-timeout-ms"].empty();
      if (is_cosim)
      {
#ifndef GRIDKIT_ENABLE_ZMQ
        throw std::invalid_argument(
          "co-simulation requires a build with GRIDKIT_ENABLE_ZMQ=ON"
        );
#endif
        if ((options.role != "client" && options.role != "server")
            || options.endpoint.empty() || options.interval <= 0.0
            || options.inputs.empty() || options.outputs.empty())
        {
          throw std::invalid_argument(
            "co-simulation requires role, endpoint, dt, input and output "
            "mappings"
          );
        }
        if (contingency && options.role != "client")
        {
          throw std::invalid_argument(
            "ContingencyAnalysis requires the co-simulation client role. this "
            "is a temporary limitation"
          );
        }
        if (options.serve && options.role != "server")
        {
          throw std::invalid_argument("--cosim-serve requires the server role");
        }
        std::set<std::size_t> ids;
        for (const auto* mappings : {&options.inputs, &options.outputs})
        {
          for (const auto& [name, id] : *mappings)
          {
            if (!ids.insert(id).second)
            {
              throw std::invalid_argument(
                std::format("signal ID mapped more than once: {}", id)
              );
            }
          }
        }
      }
      return options;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
