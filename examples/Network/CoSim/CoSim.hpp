#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSimOptions.hpp"
#include <zmq.hpp>

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

        PhasorDynamics::BusBase<ScalarT, IdxT>* fault_bus      = nullptr;
        int                                     fault_on_step  = -1;
        int                                     fault_off_step = -1;
        ScalarT                                 fault_R{}, fault_X{};
        if (!args["fault"].empty())
        {
          const auto& fault = args["fault"];
          fault_bus         = sys.getBus(parseId<IdxT>(fault[0]()));
          const auto on     = fault[1].as<double>();
          const auto off    = fault[2].as<double>();
          fault_R           = fault[3].as<double>();
          fault_X           = fault[4].as<double>();
          if (!std::isfinite(on) || !std::isfinite(off)
              || on < 0.0 || off <= on || off > tf
              || !std::isfinite(fault_R) || !std::isfinite(fault_X)
              || (fault_R == 0.0 && fault_X == 0.0))
          {
            throw std::invalid_argument(
                "Fault requires 0 <= ON_TIME < OFF_TIME <= final time and finite, nonzero impedance");
          }
          if (std::abs(on / dt - std::round(on / dt)) > 1.0e-9
              || std::abs(off / dt - std::round(off / dt)) > 1.0e-9)
          {
            throw std::invalid_argument("Fault times must lie on communication boundaries");
          }
          fault_on_step  = static_cast<int>(std::round(on / dt));
          fault_off_step = static_cast<int>(std::round(off / dt));
        }

        auto output                   = outputReader<ScalarT, IdxT>(sys, args, !server);
        const auto [real_id, imag_id] = inputSignalIds(data, args, server);
        auto* real_input              = sys.getSignalNode(real_id);
        auto* imag_input              = sys.getSignalNode(imag_id);
        if (real_input == imag_input)
        {
          throw std::invalid_argument("Input signals must have distinct IDs");
        }
        ScalarT real{}, imag{};
        IdxT    real_index{INVALID_INDEX<IdxT>}, imag_index{INVALID_INDEX<IdxT>};
        real_input->link(&real, &real_index);
        imag_input->link(&imag, &imag_index);
        zmq::context_t context;
        zmq::socket_t  socket(context, server ? zmq::socket_type::rep : zmq::socket_type::req);
        socket.set(zmq::sockopt::linger, 0);
        socket.set(zmq::sockopt::rcvtimeo, timeout);
        socket.set(zmq::sockopt::sndtimeo, timeout);
        if (server)
        {
          socket.bind(args["endpoint"]());
        }
        else
        {
          socket.connect(args["endpoint"]());
        }
        auto send = [&]()
        {
          const auto message = encode(output(), out_real, out_imag);
          if (!socket.send(zmq::buffer(message), zmq::send_flags::none))
          {
            throw std::runtime_error("Co-simulation send timed out");
          }
        };
        auto receive = [&]()
        {
          zmq::message_t message;
          if (!socket.recv(message, zmq::recv_flags::none))
          {
            throw std::runtime_error("Co-simulation receive timed out");
          }

          const auto values = decode(message.to_string(), in_real, in_imag);
          real              = values.first;
          imag              = values.second;
        };
        if (server)
        {
          receive();
        }
        sys.allocate();
        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setTolerance(server ? 1.0e-5 : 1.0e-7, server ? 1.0e-7 : 1.0e-9);
        ida.configureSimulation();
        const auto steps = ida.getStepCount(tf, dt);

        if (!server)
        {
          ida.initializeSimulation(0.0, false);
        }
        for (int step = 0; step <= steps; ++step)
        {
          if (server && step > 0)
          {
            receive();
          }
          if (step == 0)
          {
            if (server)
            {
              ida.initializeSimulation(0.0);
            }
          }
          else
          {
            const auto t      = std::fma(static_cast<double>(step - 1), dt, 0.0);
            const auto target = step == steps ? tf : std::fma(static_cast<double>(step), dt, 0.0);
            ida.initializeSimulation(t);
            ida.runSimulationStep(target, dt, 1, 1);
          }
          if (step == fault_on_step || step == fault_off_step)
          {
            const int status = step == fault_on_step
                                   ? fault_bus->setFault(true, fault_R, fault_X)
                                   : fault_bus->clearFault();
            if (status != 0)
            {
              throw std::runtime_error("Failed to apply or clear co-simulation fault");
            }
            const auto t = step == steps ? tf : std::fma(static_cast<double>(step), dt, 0.0);
            ida.initializeSimulation(t);
          }
          sys.evaluateResidual();
          send();
          if (!server)
          {
            receive();
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
