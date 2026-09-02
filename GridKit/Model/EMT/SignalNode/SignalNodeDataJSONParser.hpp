#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `SignalNodeData` type
    ///
    /// See the `README.md` in `GridKit/Model/EMT` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeData<RealT, IdxT>& sd)
    {
      j.at("name").get_to(sd.name);
      j.at("signal_id").get_to(sd.signal_id);
    }
  } // namespace EMT
} // namespace GridKit
