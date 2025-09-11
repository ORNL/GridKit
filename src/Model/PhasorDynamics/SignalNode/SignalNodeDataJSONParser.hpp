#pragma once

#include <stdexcept>

#include <Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `BusData` type
    ///
    /// See the `README.md` in `src/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeData<RealT, IdxT>& sd)
    {
      j.at("name").get_to(sd.name);
      j.at("signal_id").get_to(sd.signal_id);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
