#pragma once

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
      j.at("id").get_to(sd.id);
    }
  } // namespace EMT
} // namespace GridKit
