#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <magic_enum/magic_enum.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Three-phase vector parameter value
    template <typename T>
    using ABCVector = std::array<T, 3>;

    /// Three-phase matrix parameter value
    template <typename T>
    using ABCMatrix = std::array<std::array<T, 3>, 3>;

    /// Value type held by a component parameter map entry
    template <typename real_type, typename index_type>
    using ParameterValue = std::variant<bool,
                                        real_type,
                                        index_type,
                                        ABCVector<real_type>,
                                        ABCVector<index_type>,
                                        ABCMatrix<real_type>>;

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     */
    template <typename real_type,
              typename index_type,
              typename Parameters,
              typename Inputs,
              typename Outputs,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Inputs>
               && std::is_enum_v<Outputs>
               && std::is_enum_v<MonitorableVariables>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;
      /// Index type
      using IdxT  = index_type;

      /// Class of device this is for
      std::string device_class;

      /// Unique component identifier within the system
      std::string id;

      /// Mapping of parameters to parameter values
      std::map<Parameters, ParameterValue<RealT, IdxT>> parameters;

      /// Mapping of model inputs to component or signal identifiers
      std::map<Inputs, std::string> inputs;

      /// Mapping of model outputs to signal identifiers
      std::map<Outputs, std::string> outputs;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

    protected:
      ComponentData() = default;
    };

    namespace detail
    {
      /// Convert parameter values without treating booleans as numbers or wrapping indices.
      template <typename T, typename U>
      T parameterValue(const U& value)
      {
        if constexpr (std::is_same_v<T, bool> && std::is_same_v<U, bool>)
        {
          return value;
        }
        else if constexpr (std::is_arithmetic_v<T> && std::is_arithmetic_v<U>
                           && !std::is_same_v<T, bool> && !std::is_same_v<U, bool>)
        {
          if constexpr (std::is_integral_v<T>)
          {
            if constexpr (std::is_integral_v<U>)
            {
              if (value < 0 || !std::in_range<T>(value))
                throw std::invalid_argument("requires a nonnegative integer within the index range");
            }
            else if (!std::isfinite(value) || value < 0 || std::trunc(value) != value
                     || value >= std::ldexp(U{1}, std::numeric_limits<T>::digits))
            {
              throw std::invalid_argument("requires a nonnegative integer within the index range");
            }
          }
          const auto result = static_cast<T>(value);
          if (!std::isfinite(result))
            throw std::invalid_argument("requires finite numeric values");
          return result;
        }
        else if constexpr (requires { typename T::value_type; typename U::value_type; })
        {
          T result{};
          for (size_t n = 0; n < result.size(); ++n)
            result[n] = parameterValue<typename T::value_type>(value[n]);
          return result;
        }
        throw std::invalid_argument("has an incompatible value type");
      }
    } // namespace detail

    /// Read a required, typed parameter with component and parameter context on errors.
    template <typename T, typename Data, typename Parameter>
    T parameter(const Data& data, Parameter key)
    {
      try
      {
        const auto entry = data.parameters.find(key);
        if (entry == data.parameters.end())
          throw std::invalid_argument("is required");
        return std::visit([](const auto& value)
                          { return detail::parameterValue<T>(value); },
                          entry->second);
      }
      catch (const std::invalid_argument& error)
      {
        throw std::invalid_argument(data.device_class + " \"" + data.id + "\" parameter \""
                                    + std::string(magic_enum::enum_name(key)) + "\" " + error.what());
      }
    }

    /// Read an optional parameter, preserving the model default when omitted.
    template <typename T, typename Data, typename Parameter>
    T parameter(const Data& data, Parameter key, const T& fallback)
    {
      return data.parameters.contains(key) ? parameter<T>(data, key) : fallback;
    }
  } // namespace EMT
} // namespace GridKit
