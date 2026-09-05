#pragma once

#include <cmath>
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
      j.at("id").get_to(sd.id);
      sd.value.reset();
      if (j.contains("value"))
      {
        if (!j.at("value").is_number())
        {
          throw std::invalid_argument("A constant signal value must be numeric");
        }
        const auto value = j.at("value").template get<RealT>();
        if (!std::isfinite(value))
        {
          throw std::invalid_argument("A constant signal value must be finite");
        }
        sd.value = value;
      }
    }
  } // namespace EMT
} // namespace GridKit
