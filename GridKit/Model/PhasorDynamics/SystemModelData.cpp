#include "SystemModelData.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace
    {
      /// Whether the state sets `flag` of device `id` to `value`
      bool hasFlag(const Model::StateData& state, const std::string& id, const std::string& flag, bool value)
      {
        const Model::StateRecord* record = state.device(id);
        return record != nullptr && record->flags.contains(flag) && record->flags.at(flag) == value;
      }

      /// Remove the devices whose `flag` the state sets to `value`
      template <typename DeviceDataT>
      void removeDevices(std::vector<DeviceDataT>& devices,
                         const Model::StateData&   state,
                         const std::string&        flag,
                         bool                      value)
      {
        std::vector<DeviceDataT> kept;
        for (auto& device : devices)
        {
          if (!hasFlag(state, device.disambiguation_string, flag, value))
          {
            kept.push_back(std::move(device));
          }
        }
        devices = std::move(kept);
      }

      /// Machine `p0` and `q0` from the power the state injects
      template <typename MachineDataT>
      void applyInjections(std::vector<MachineDataT>& machines, const Model::StateData& state)
      {
        using Parameters = typename MachineDataT::Parameters;
        using Buses      = typename MachineDataT::Buses;

        for (auto& machine : machines)
        {
          const auto& id = machine.disambiguation_string;
          if (hasFlag(state, id, "online", false))
          {
            throw std::invalid_argument("Offline machine " + id + " is not supported yet");
          }

          double p = 0.0;
          double q = 0.0;
          if (Model::terminalPower(state, id, machine.buses.at(Buses::bus), 0, 1, p, q))
          {
            machine.parameters[Parameters::p0] = p;
            machine.parameters[Parameters::q0] = q;
          }
        }
      }

      /// Machine currents that inject `p0` and `q0`
      template <typename MachineDataT>
      void extractInjections(Model::StateData& state, const std::vector<MachineDataT>& machines)
      {
        using Parameters = typename MachineDataT::Parameters;
        using Buses      = typename MachineDataT::Buses;

        for (const auto& machine : machines)
        {
          Model::setTerminalCurrent(state,
                                    machine.disambiguation_string,
                                    machine.buses.at(Buses::bus),
                                    0,
                                    1,
                                    realParameter(machine, Parameters::p0, 0.0),
                                    realParameter(machine, Parameters::q0, 0.0));
        }
      }
    } // namespace

    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream)
    {
      SystemModelData<double, size_t> data(json::parse(stream));
      return data;
    }

    SystemModelData<double, size_t> parseSystemModelData(std::istream&& stream)
    {
      SystemModelData<double, size_t> data(json::parse(stream));
      return data;
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& filePath)
    {
      auto stream = std::ifstream(filePath);
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << filePath;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      return parseSystemModelData(stream);
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::string& fileName)
    {
      auto stream = std::ifstream(fileName);
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << fileName;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      return parseSystemModelData(stream);
    }

    void applyState(SystemModelData<double, size_t>& data, const Model::StateData& state)
    {
      for (auto& bus : data.bus)
      {
        const Model::StateRecord* record = state.bus(bus.bus_id);
        if (record != nullptr)
        {
          bus.Vr0 = record->value("vr", bus.Vr0);
          bus.Vi0 = record->value("vi", bus.Vi0);
        }
      }

      applyInjections(data.genrou, state);
      applyInjections(data.gensal, state);
      applyInjections(data.genclassical, state);
      applyInjections(data.regca, state);

      // LoadZIP anchors at its initial voltage, so it draws Pnom + jQnom at t = 0
      for (auto& load : data.loadzip)
      {
        const auto& id = load.disambiguation_string;

        double p         = 0.0;
        double q         = 0.0;
        bool   has_power = Model::terminalPower(state, id, load.buses.at(LoadZIPBuses::bus), 0, 1, p, q);
        if (hasFlag(state, id, "online", false))
        {
          p         = 0.0;
          q         = 0.0;
          has_power = true;
        }

        if (has_power)
        {
          load.parameters[LoadZIPParameters::Pnom] = -p;
          load.parameters[LoadZIPParameters::Qnom] = -q;
        }
      }

      removeDevices(data.loadz, state, "online", false);
      removeDevices(data.branch, state, "open", true);

      for (auto& branch : data.branch)
      {
        const Model::StateRecord* record = state.device(branch.disambiguation_string);
        if (record != nullptr)
        {
          if (record->values.contains("tap"))
          {
            branch.parameters[BranchParameters::tap] = record->values.at("tap");
          }
          if (record->values.contains("phase"))
          {
            branch.parameters[BranchParameters::phase] = record->values.at("phase");
          }
        }
      }
    }

    Model::StateData extractState(const SystemModelData<double, size_t>& data)
    {
      Model::StateData state;

      for (const auto& bus : data.bus)
      {
        auto& voltage = state.buses[Model::busKey(bus.bus_id)].values;
        voltage["vr"] = bus.Vr0;
        voltage["vi"] = bus.Vi0;
      }

      extractInjections(state, data.genrou);
      extractInjections(state, data.gensal);
      extractInjections(state, data.genclassical);
      extractInjections(state, data.regca);

      for (const auto& load : data.loadzip)
      {
        Model::setTerminalCurrent(state,
                                  load.disambiguation_string,
                                  load.buses.at(LoadZIPBuses::bus),
                                  0,
                                  1,
                                  -realParameter(load, LoadZIPParameters::Pnom, 0.0),
                                  -realParameter(load, LoadZIPParameters::Qnom, 0.0));
      }

      // Power into the bus is -conj(Y) V^2 with Y = 1 / (R + jX)
      for (const auto& load : data.loadz)
      {
        const size_t number = load.buses.at(LoadZBuses::bus);
        const double r      = realParameter(load, LoadZParameters::R);
        const double x      = realParameter(load, LoadZParameters::X);
        const double g      = r / (r * r + x * x);
        const double b      = -x / (r * r + x * x);

        const auto&  voltage = state.buses.at(Model::busKey(number)).values;
        const double v2      = voltage.at("vr") * voltage.at("vr") + voltage.at("vi") * voltage.at("vi");

        Model::setTerminalCurrent(state, load.disambiguation_string, number, 0, 1, -g * v2, b * v2);
      }

      for (const auto& branch : data.branch)
      {
        auto& settings    = state.devices[branch.disambiguation_string].values;
        settings["tap"]   = realParameter(branch, BranchParameters::tap, 1.0);
        settings["phase"] = realParameter(branch, BranchParameters::phase, 0.0);
      }

      return state;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
