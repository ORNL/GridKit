/**
 * @file StateUtilities.hpp
 * @brief Apply operating states to case data in PhasorDynamics applications.
 */

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace detail
    {
      inline bool hasFlag(const Model::StateData& state, const std::string& id, const std::string& flag, bool value)
      {
        const Model::StateRecord* record = state.device(id);
        return record != nullptr && record->flags.contains(flag) && record->flags.at(flag) == value;
      }

      template <typename DeviceDataT>
      void removeDevices(std::vector<DeviceDataT>& devices,
                         const Model::StateData&   state,
                         const std::string&        flag,
                         bool                      value)
      {
        std::erase_if(devices, [&](const auto& device)
                      { return hasFlag(state, device.disambiguation_string, flag, value); });
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
    } // namespace detail

    /**
     * @brief Apply voltages, terminal injections, and device settings to a case
     */
    inline void applyState(SystemModelData<double, size_t>& data, const Model::StateData& state)
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

      detail::applyInjections(data.genrou, state);
      detail::applyInjections(data.gensal, state);
      detail::applyInjections(data.genclassical, state);
      detail::applyInjections(data.regca, state);

      // LoadZIP anchors at its initial voltage, so it draws Pnom + jQnom at t = 0
      for (auto& load : data.loadzip)
      {
        const auto& id = load.disambiguation_string;

        double p         = 0.0;
        double q         = 0.0;
        bool   has_power = Model::terminalPower(state, id, load.buses.at(LoadZIPBuses::bus), 0, 1, p, q);
        if (detail::hasFlag(state, id, "online", false))
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

      detail::removeDevices(data.loadz, state, "online", false);
      detail::removeDevices(data.branch, state, "open", true);

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
  } // namespace PhasorDynamics
} // namespace GridKit
