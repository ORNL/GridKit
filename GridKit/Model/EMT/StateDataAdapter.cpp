/**
 * @file StateDataAdapter.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Application of portable model state to EMT system model data
 *
 */
#include "StateDataAdapter.hpp"

#include <string>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Apply a parsed operating point to parsed system model data.
     *
     * The adapter writes only data: bus instantaneous voltages land in the
     * bus initial conditions, machine dispatch lands in the machine `p0` and
     * `q0` parameters, and switch commands land in the switch `open`
     * parameter. Model constructors remain the single ingestion path.
     */
    void applyState(SystemModelData<double, size_t>& model_data,
                    const Model::StateData&          state_data)
    {
      for (auto& bus : model_data.bus)
      {
        const auto key   = "bus_id_" + std::to_string(bus.bus_id);
        const auto entry = state_data.buses.find(key);
        if (entry == state_data.buses.end())
        {
          continue;
        }

        const auto& state = entry->second;
        if (state.va.has_value())
        {
          bus.va0 = *state.va;
        }
        if (state.vb.has_value())
        {
          bus.vb0 = *state.vb;
        }
        if (state.vc.has_value())
        {
          bus.vc0 = *state.vc;
        }
      }

      for (auto& machine : model_data.synchronous_machine)
      {
        const auto entry = state_data.devices.find(machine.disambiguation_string);
        if (entry == state_data.devices.end())
        {
          continue;
        }

        const auto& state = entry->second;
        using Parameter   = SynchronousMachineParameters;
        if (state.p.has_value())
        {
          machine.parameters[Parameter::p0] = *state.p;
        }
        if (state.q.has_value())
        {
          machine.parameters[Parameter::q0] = *state.q;
        }
      }

      for (auto& sw : model_data.sw)
      {
        const auto entry = state_data.devices.find(sw.disambiguation_string);
        if (entry == state_data.devices.end())
        {
          continue;
        }

        const auto& state = entry->second;
        if (state.open.has_value())
        {
          sw.parameters[SwitchParameters::open] = *state.open;
        }
      }
    }
  } // namespace EMT
} // namespace GridKit
