#include "StateDataAdapter.hpp"

#include <map>
#include <string>
#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace
    {
      template <typename DeviceDataT>
      void attachDeviceStates(std::vector<DeviceDataT>&                        devices,
                              const std::map<std::string, Model::DeviceState>& states)
      {
        for (auto& device : devices)
        {
          auto state = states.find(device.disambiguation_string);
          if (state != states.end())
          {
            device.initial_state = state->second;
          }
        }
      }
    } // namespace

    void applyState(SystemModelData<>& model_data, const Model::StateData& state_data)
    {
      for (auto& bus : model_data.bus)
      {
        const auto state_id = std::string{"bus_id_"} + std::to_string(bus.bus_id);
        auto       state    = state_data.buses.find(state_id);
        if (state != state_data.buses.end())
        {
          bus.initial_state = state->second;
        }
      }

      attachDeviceStates(model_data.branch, state_data.devices);
      attachDeviceStates(model_data.genrou, state_data.devices);
      attachDeviceStates(model_data.gensal, state_data.devices);
      attachDeviceStates(model_data.genclassical, state_data.devices);
      attachDeviceStates(model_data.loadz, state_data.devices);
      attachDeviceStates(model_data.loadzip, state_data.devices);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
