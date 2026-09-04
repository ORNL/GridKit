#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Signal/SignalData.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `SignalData` type
    ///
    /// See the `README.md` in `GridKit/Model/EMT` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalData<RealT, IdxT>& sd)
    {
      j.at("name").get_to(sd.name);
      j.at("signal_id").get_to(sd.signal_id);
    }
  } // namespace EMT
} // namespace GridKit
