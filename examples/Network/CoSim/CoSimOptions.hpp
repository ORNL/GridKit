#pragma once

#include <array>
#include <charconv>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>

namespace GridKit
{
  namespace CoSim
  {
    inline Utilities::CliArgs commandLine()
    {
      return Utilities::CliArgs{
          {.name       = {"case-file"},
           .help       = "Path to the partition case JSON file",
           .required   = true,
           .positional = true},
          {.name  = {"--input-signals", "-i"},
           .help  = "Input signal IDs for --output-signals mode",
           .nargs = 2},
          {.name = {"--output-bus", "-b"},
           .help = "Interface bus ID. Infer input signals and send bus voltage (client) or current (server)"},
          {.name  = {"--output-signals", "-o"},
           .help  = "IDs of the output signals, instead of a bus",
           .nargs = 2},
          {.name     = {"--time-step"},
           .help     = "Communication step in seconds (these must match across client/server)",
           .defaults = 1.0 / 16.0 / 60.0},
          {.name     = {"--final-time"},
           .help     = "Final simulation time in seconds (default: 1s)",
           .defaults = 1.0},
          {.name     = {"--endpoint"},
           .help     = "ZeroMQ address to bind (server) or connect (client)",
           .defaults = "tcp://127.0.0.1:5556"},
          {.name  = {"--input-tags"},
           .help  = "JSON input keys (default: BUS_vr/BUS_vi on server, BUS_ir/BUS_ii on client)",
           .nargs = 2},
          {.name  = {"--output-tags"},
           .help  = "JSON output keys (default: BUS_ir/BUS_ii on server, BUS_vr/BUS_vi on client)",
           .nargs = 2},
          {.name     = {"--timeout-ms"},
           .help     = "Maximum wait for a peer message in milliseconds",
           .defaults = 10000}};
    }

    template <typename IdxT>
    IdxT parseId(const std::string& text)
    {
      IdxT value{};
      const auto [end, error] = std::from_chars(text.data(),
                                                text.data() + text.size(),
                                                value);
      if (error != std::errc{} || end != text.data() + text.size())
      {
        throw std::invalid_argument("Invalid bus or signal ID: " + text);
      }
      return value;
    }

    template <typename RealT, typename IdxT>
    const PhasorDynamics::BusData<RealT, IdxT>* findBusData(
        const PhasorDynamics::SystemModelData<RealT, IdxT>& data, IdxT id)
    {
      for (const auto& bus : data.bus)
      {
        if (bus.bus_id == id)
        {
          return &bus;
        }
      }
      for (const auto& subsystem : data.subsystems)
      {
        if (const auto* bus = findBusData(subsystem, id))
        {
          return bus;
        }
      }
      return nullptr;
    }

    template <typename RealT, typename IdxT>
    std::array<IdxT, 2> inputSignalIds(
        const PhasorDynamics::SystemModelData<RealT, IdxT>& data,
        const Utilities::CliArgs&                           args,
        bool                                                server)
    {
      if (args["output-bus"].empty())
      {
        if (args["input-signals"].empty())
        {
          throw std::invalid_argument("--output-signals requires --input-signals");
        }
        return {parseId<IdxT>(args["input-signals"][0]()),
                parseId<IdxT>(args["input-signals"][1]())};
      }
      if (!args["input-signals"].empty())
      {
        throw std::invalid_argument(
            "--output-bus infers input signals. omit --input-signals");
      }
      const auto* bus = findBusData(data, parseId<IdxT>(args["output-bus"]()));
      using BusType   = typename PhasorDynamics::BusData<RealT, IdxT>::BusType;
      if (!bus || bus->bus_type != (server ? BusType::SIGNAL_VOLTAGE_IN : BusType::SIGNAL_VOLTAGE_OUT))
      {
        throw std::invalid_argument(
            "Interface bus must be SignalVoltageIn on the server or SignalVoltageOut on the client");
      }
      using Input     = PhasorDynamics::BusSignalInputs;
      const auto real = server ? Input::vr : Input::ir;
      const auto imag = server ? Input::vi : Input::ii;
      if (!bus->signal_inputs.contains(real) || !bus->signal_inputs.contains(imag))
      {
        throw std::invalid_argument("Interface bus is missing a required input signal port");
      }
      return {bus->signal_inputs.at(real), bus->signal_inputs.at(imag)};
    }

    inline std::array<std::string, 2> tagNames(const Utilities::CliArgs& args,
                                               const std::string&        option,
                                               bool                      voltage)
    {
      if (!args[option].empty())
      {
        return args[option].as<2>();
      }
      if (args["output-bus"].empty())
      {
        throw std::invalid_argument("--output-signals requires explicit --" + option);
      }
      const auto prefix = std::to_string(parseId<std::size_t>(args["output-bus"]())) + "_";
      return {prefix + (voltage ? "vr" : "ir"), prefix + (voltage ? "vi" : "ii")};
    }

    template <typename ScalarT, typename IdxT>
    std::function<std::pair<ScalarT, ScalarT>()> outputReader(
        PhasorDynamics::SystemModel<ScalarT, IdxT>& sys,
        const Utilities::CliArgs&                   args,
        bool                                        voltage)
    {
      const bool bus     = !args["output-bus"].empty();
      const bool signals = !args["output-signals"].empty();
      if (bus == signals)
      {
        throw std::invalid_argument(
            "Specify exactly one of --output-bus and --output-signals");
      }
      if (bus)
      {
        auto* output = sys.getBus(parseId<IdxT>(args["output-bus"]()));
        return [output, voltage]() -> std::pair<ScalarT, ScalarT>
        {
          return voltage ? std::pair{output->Vr(), output->Vi()}
                         : std::pair{output->Ir(), output->Ii()};
        };
      }
      auto* real = sys.getSignalNode(parseId<IdxT>(args["output-signals"][0]()));
      auto* imag = sys.getSignalNode(parseId<IdxT>(args["output-signals"][1]()));
      return [real, imag]() -> std::pair<ScalarT, ScalarT>
      { return {real->read(), imag->read()}; };
    }
  } // namespace CoSim
} // namespace GridKit
