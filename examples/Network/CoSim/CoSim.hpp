#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSimOptions.hpp"
#include <zmq.hpp>
#ifdef GRIDKIT_COSIM_ARCADE
#include "ArcadeClient.hpp"
#endif

namespace GridKit
{
  namespace CoSim
  {
    using json = nlohmann::json;

    inline std::pair<double, double> decode(const std::string& message,
                                            const std::string& real_tag,
                                            const std::string& imag_tag)
    {
      const auto data  = json::parse(message);
      auto       value = [&](const std::string& tag)
      {
        if (!data.is_object() || !data.contains(tag) || !data.at(tag).is_number())
        {
          throw std::invalid_argument("Expected numeric JSON input: " + tag);
        }
        const auto result = data.at(tag).get<double>();
        if (!std::isfinite(result))
        {
          throw std::invalid_argument("Non-finite JSON input: " + tag);
        }
        return result;
      };
      return {value(real_tag), value(imag_tag)};
    }

    inline std::string encode(const std::pair<double, double>& values,
                              const std::string&               real_tag,
                              const std::string&               imag_tag)
    {
      if (!std::isfinite(values.first) || !std::isfinite(values.second))
      {
        throw std::runtime_error("Cannot send non-finite co-simulation output");
      }
      return json{{real_tag, values.first}, {imag_tag, values.second}}.dump();
    }

    inline int run(int argc, const char* argv[], bool server)
    {
      using ScalarT = double;
      using IdxT    = std::size_t;
      using Log     = Utilities::Logger;
      try
      {
        auto args = commandLine();
        args.parseArgs(argc, argv);
        const auto dt      = args.get<double>("time-step");
        const auto tf      = args.get<double>("final-time");
        const auto timeout = args.get<int>("timeout-ms");
        if (!std::isfinite(dt) || dt <= 0.0
            || !std::isfinite(tf) || tf <= 0.0
            || timeout <= 0)
        {
          throw std::invalid_argument(
              "Time step, final time, and timeout must be finite and positive");
        }
        if (tf / dt >= std::numeric_limits<int>::max())
        {
          throw std::invalid_argument("Too many co-simulation steps");
        }
        const auto transport = args["transport"]();
        if (transport != "zmq" && transport != "arcade")
        {
          throw std::invalid_argument("Expected --transport zmq or arcade");
        }
        if (server && transport == "arcade")
        {
          throw std::invalid_argument("ARCADE uses CoSimClient for shared memory and CoSimServer for ZeroMQ");
        }
#ifndef GRIDKIT_COSIM_ARCADE
        if (transport == "arcade")
        {
          throw std::invalid_argument("ARCADE transport requires Unix");
        }
#endif
        const auto [in_real, in_imag]   = tagNames(args, "input-tags", server);
        const auto [out_real, out_imag] = tagNames(args, "output-tags", !server);
        if (in_real.empty() || in_imag.empty() || out_real.empty() || out_imag.empty()
            || in_real == in_imag || out_real == out_imag)
        {
          throw std::invalid_argument(
              "Each JSON tag pair must contain distinct, nonempty names");
        }

        const auto                                 data = PhasorDynamics::parseSystemModelData(args["case-file"]());
        PhasorDynamics::SystemModel<ScalarT, IdxT> sys(data);
        auto                                       output = outputReader<ScalarT, IdxT>(sys, args, !server);
        const auto [real_id, imag_id]                     = inputSignalIds(data, args, server);
        auto* real_input                                  = sys.getSignalNode(real_id);
        auto* imag_input                                  = sys.getSignalNode(imag_id);
        if (real_input == imag_input)
        {
          throw std::invalid_argument("Input signals must have distinct IDs");
        }
        ScalarT real{}, imag{};
        IdxT    real_index{INVALID_INDEX<IdxT>}, imag_index{INVALID_INDEX<IdxT>};
        real_input->link(&real, &real_index);
        imag_input->link(&imag, &imag_index);
#ifdef GRIDKIT_COSIM_ARCADE
        std::optional<ArcadeClient> arcade;
        if (transport == "arcade")
        {
          arcade.emplace(std::array{in_real, in_imag},
                         std::array{out_real, out_imag},
                         dt,
                         timeout);
        }
#endif
        std::optional<zmq::context_t> context;
        std::optional<zmq::socket_t>  socket;
        if (transport == "zmq")
        {
          context.emplace();
          socket.emplace(*context, server ? zmq::socket_type::rep : zmq::socket_type::req);
          socket->set(zmq::sockopt::linger, server ? timeout : 0);
          socket->set(zmq::sockopt::rcvtimeo, timeout);
          socket->set(zmq::sockopt::sndtimeo, timeout);
          if (server)
          {
            socket->bind(args["endpoint"]());
          }
          else
          {
            socket->connect(args["endpoint"]());
          }
        }
        auto send = [&](double time)
        {
#ifdef GRIDKIT_COSIM_ARCADE
          if (arcade)
          {
            arcade->send(output(), time);
            return;
          }
#else
          (void) time;
#endif
          const auto message = encode(output(), out_real, out_imag);
          if (!socket->send(zmq::buffer(message), zmq::send_flags::none))
          {
            throw std::runtime_error("Co-simulation send timed out");
          }
        };
        auto receive = [&]([[maybe_unused]] double time)
        {
#ifdef GRIDKIT_COSIM_ARCADE
          if (arcade)
          {
            std::tie(real, imag) = arcade->receive(time);
            return;
          }
#endif
          zmq::message_t message;
          if (!socket->recv(message, zmq::recv_flags::none))
          {
            throw std::runtime_error("Co-simulation receive timed out");
          }
          std::tie(real, imag) = decode(message.to_string(), in_real, in_imag);
        };
        if (server)
        {
          receive(0.0);
        }
        sys.allocate();
        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setTolerance(server ? 1.0e-5 : 1.0e-7, server ? 1.0e-7 : 1.0e-9);
        ida.configureSimulation();
        const auto steps = std::max(1, ida.getStepCount(tf, dt));

        ida.initializeSimulation(0.0, server);
        sys.evaluateResidual();
#ifdef GRIDKIT_COSIM_ARCADE
        if (arcade)
        {
          arcade->start(output());
        }
#endif
        send(0.0);
        if (!server)
        {
          receive(0.0);
        }
        for (int step = 0; step < steps; ++step)
        {
          const auto t      = step * dt;
          const auto target = step + 1 == steps ? tf : (step + 1) * dt;
          if (server)
          {
            receive(target);
          }
          ida.initializeSimulation(t, true, false);
          ida.runSimulationStep(target, target - t, 1, 1, [&](double time)
                                {
                                 sys.evaluateResidual();
                                 send(time); });
          if (!server)
          {
            receive(target);
          }
        }
        return 0;
      }
      catch (const std::exception& error)
      {
        Log::error() << "Co-simulation: " << error.what() << '\n';
        return 1;
      }
    }
  } // namespace CoSim
} // namespace GridKit
