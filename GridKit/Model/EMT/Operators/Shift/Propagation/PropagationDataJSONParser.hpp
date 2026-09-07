#pragma once

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Propagation/PropagationData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const nlohmann::json& j, PropagationData<RealT, IdxT>& data)
    {
      validateJsonFields(j, "Propagation", {"K", "modes", "limit_step"});
      data.K          = parseRationalDimension<IdxT>(j, "K");
      data.limit_step = j.value("limit_step", false);
      if (!j.at("modes").is_array())
        throw std::invalid_argument("Propagation: modes must be an array");
      data.modes.clear();
      for (const auto& mode : j.at("modes"))
      {
        validateJsonFields(mode, "Propagation mode", {"tau", "H"});
        data.modes.push_back({parseFiniteReal<RealT>(mode.at("tau"), "Propagation tau"),
                              parseVectorFitOperand<RealT, IdxT>(mode.at("H"))});
      }
      if (data.validate())
        throw std::invalid_argument("Propagation: invalid mode coefficients");
    }

    template <typename RealT, typename IdxT>
    PropagationData<RealT, IdxT> parsePropagationOperand(const nlohmann::json& value)
    {
      if (value.is_string())
      {
        const auto    path = value.template get<std::string>();
        std::ifstream stream(path);
        if (!stream)
          throw std::runtime_error("Could not open propagation file: " + path);
        return nlohmann::json::parse(stream).template get<PropagationData<RealT, IdxT>>();
      }
      return value.template get<PropagationData<RealT, IdxT>>();
    }
  } // namespace EMT
} // namespace GridKit
