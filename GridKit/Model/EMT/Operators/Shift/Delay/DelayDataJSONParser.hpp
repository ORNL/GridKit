#pragma once

#include <GridKit/Model/EMT/Operators/Rational/RationalDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Delay/DelayData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const nlohmann::json& j, DelayData<RealT, IdxT>& data)
    {
      validateJsonFields(j, "Delay", {"M", "tau", "limit_step"});
      data.M = parseRationalDimension<IdxT>(j, "M");
      if (!j.at("tau").is_array())
        throw std::invalid_argument("Delay: tau must be an array");
      data.tau.clear();
      for (const auto& value : j.at("tau"))
        data.tau.push_back(parseFiniteReal<RealT>(value, "Delay tau"));
      data.limit_step = j.value("limit_step", false);
      if (data.validate())
        throw std::invalid_argument("Delay: invalid channel delays");
    }
  } // namespace EMT
} // namespace GridKit
