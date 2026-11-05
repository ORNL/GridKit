/**
 * @file ComponentDataJSONParser.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief JSON parsing of EMT component data.
 */

#pragma once

#include <tuple>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    template <typename KeyT, typename DataT, typename MemberT>
    void readParameter(const json& value, KeyT key, const Parameter<KeyT, DataT, MemberT>& parameter, DataT& data)
    {
      if (parameter.key == key)
      {
        value.get_to(data.*parameter.member);
      }
    }

    /// An omitted parameter keeps its default
    template <typename DataT>
      requires requires { DataT::PARAMETERS; }
    void from_json(const json& j, DataT& data)
    {
      for (const auto& [name, value] : j.items())
      {
        const auto key = Utilities::enumKey<typename DataT::Parameters>(name, "parameter", "params");
        std::apply([&](const auto&... parameter)
                   { (readParameter(value, key, parameter, data), ...); },
                   DataT::PARAMETERS);
      }
    }
  } // namespace EMT
} // namespace GridKit
