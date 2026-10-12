#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `BusData` type
    ///
    /// See the `README.md` in `GridKit/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SignalNodeData<RealT, IdxT>& sd)
    {
      j.at("name").get_to(sd.name);
      j.at("signal_id").get_to(sd.signal_id);
      sd.source.reset();
      sd.delay   = 0;
      sd.monitor = false;
      if (j.contains("source") != j.contains("delay"))
      {
        throw std::invalid_argument("A delayed signal requires both source and delay");
      }
      if (j.contains("source"))
      {
        sd.source = j.at("source").template get<IdxT>();
        j.at("delay").get_to(sd.delay);
        if (!std::isfinite(sd.delay) || sd.delay < 0)
        {
          throw std::invalid_argument("Signal delay must be finite and nonnegative");
        }
      }
      if (j.contains("mon"))
      {
        for (const auto& variable : j.at("mon").template get<std::vector<std::string>>())
        {
          if (variable != "z" || !sd.source)
          {
            throw std::invalid_argument("Signal monitor z requires a positive delay");
          }
          sd.monitor = true;
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
