#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// Parse one weighted signal-node junction input.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeJunctionInputData<RealT, IdxT>& input)
    {
      j.at("signal_id").get_to(input.signal_id);
      input.gain = j.value("gain", RealT{1});
    }

    /// Parse an algebraic signal-node junction definition.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeJunctionData<RealT, IdxT>& junction)
    {
      junction.bias = j.value("bias", RealT{0});
      j.at("initialization_input").get_to(junction.initialization_input);
      j.at("inputs").get_to(junction.inputs);
    }

    /// JSON parser function implementation for the `SignalNodeData` type.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeData<RealT, IdxT>& sd)
    {
      j.at("name").get_to(sd.name);
      j.at("signal_id").get_to(sd.signal_id);

      if (j.contains("junction"))
      {
        sd.junction.emplace();
        j.at("junction").get_to(*sd.junction);
      }
      else
      {
        sd.junction.reset();
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
