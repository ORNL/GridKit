#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>

#include <nlohmann/json.hpp>

namespace GridKit::EMT
{
  inline void validateJsonFields(const nlohmann::json&                   object,
                                 const std::string&                      context,
                                 std::initializer_list<std::string_view> allowed,
                                 std::initializer_list<std::string_view> extra = {})
  {
    if (!object.is_object())
      throw std::invalid_argument(context + " must be an object");
    for (const auto& [key, value] : object.items())
    {
      static_cast<void>(value);
      if (std::find(allowed.begin(), allowed.end(), key) == allowed.end()
          && std::find(extra.begin(), extra.end(), key) == extra.end())
        throw std::invalid_argument(context + " has unknown field \"" + key + "\"");
    }
  }

  template <typename RealT>
  RealT parseFiniteReal(const nlohmann::json& value, const std::string& context)
  {
    if (!value.is_number())
      throw std::invalid_argument(context + " must be numeric");
    const auto result = value.template get<RealT>();
    if (!std::isfinite(result))
      throw std::invalid_argument(context + " must be finite");
    return result;
  }
} // namespace GridKit::EMT
