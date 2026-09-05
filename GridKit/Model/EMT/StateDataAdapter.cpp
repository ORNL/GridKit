/**
 * @file StateDataAdapter.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Application of portable model state to EMT system model data
 *
 */
#include "StateDataAdapter.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace
    {
      void applyContainerState(ContainerData<double, size_t>& data,
                               const Model::StateData&        state_data,
                               const std::string&             prefix)
      {
        for (auto& bus : data.bus)
        {
          const auto entry = state_data.buses.find(prefix + bus.id);
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

        for (auto& machine : data.machine)
        {
          const auto entry = state_data.devices.find(prefix + machine.id);
          if (entry == state_data.devices.end())
          {
            continue;
          }

          const auto& state = entry->second;
          using Parameter   = MachineParameters;
          if (state.p.has_value())
          {
            machine.parameters[Parameter::p0] = *state.p;
          }
          if (state.q.has_value())
          {
            machine.parameters[Parameter::q0] = *state.q;
          }
        }

        for (auto& sw : data.sw)
        {
          const auto entry = state_data.devices.find(prefix + sw.id);
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

        for (auto& child : data.container)
        {
          applyContainerState(child, state_data, prefix + child.id + ".");
        }
      }
    } // namespace

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
      applyContainerState(model_data, state_data, "");
    }
  } // namespace EMT
} // namespace GridKit
