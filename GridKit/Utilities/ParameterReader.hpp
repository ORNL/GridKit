/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter loading from model data containers.
 */

#pragma once

#include <cmath>
#include <map>
#include <variant>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Utilities/ConfigurationChecks.hpp>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Reads typed parameters out of a model data container.
     *
     * Every method leaves the target untouched and reports one error through
     * the shared checks object when a provided value has the wrong type or is
     * not finite; an omitted optional parameter keeps the model default
     * silently. Real parameters accept integer values. Parameter names in
     * messages come from the parameter enumeration.
     *
     * @tparam ModelDataT A model data container exposing `RealT`, `IdxT`,
     *         the `Parameters` enumeration, and the `parameters` map.
     */
    template <typename ModelDataT>
    class ParameterReader
    {
    public:
      using RealT       = typename ModelDataT::RealT;
      using IdxT        = typename ModelDataT::IdxT;
      using ParametersT = typename ModelDataT::Parameters;

      ParameterReader(const ModelDataT& data, ConfigurationChecks& checks)
        : parameters_(data.parameters),
          checks_(checks)
      {
      }

      /**
       * @brief Load an optional real parameter. Integer values are accepted.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the finite numeric value when provided.
       * @return true when the parameter was provided and stored.
       */
      bool loadReal(ParametersT key, RealT& target)
      {
        if (!parameters_.contains(key))
        {
          return false;
        }

        const auto& value = parameters_.at(key);
        RealT       parsed_value{};
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          parsed_value = *real_value;
        }
        else if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          parsed_value = static_cast<RealT>(*index_value);
        }
        else
        {
          checks_.fail() << "parameter '" << magic_enum::enum_name(key)
                         << "' must be numeric\n";
          return false;
        }

        if (!std::isfinite(parsed_value))
        {
          checks_.fail() << "parameter '" << magic_enum::enum_name(key)
                         << "' must be finite\n";
          return false;
        }

        target = parsed_value;
        return true;
      }

      /**
       * @brief Load a real parameter that must be provided.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the finite numeric value when provided.
       * @return true when the parameter was provided and stored.
       */
      bool requireReal(ParametersT key, RealT& target)
      {
        if (!parameters_.contains(key))
        {
          checks_.fail() << "missing required parameter '"
                         << magic_enum::enum_name(key) << "'\n";
          return false;
        }
        return loadReal(key, target);
      }

      /**
       * @brief Load an optional boolean switch parameter.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the boolean value when provided.
       * @return true when the parameter was provided and stored.
       */
      bool loadSwitch(ParametersT key, bool& target)
      {
        if (!parameters_.contains(key))
        {
          return false;
        }

        const auto& value = parameters_.at(key);
        if (const auto* bool_value = std::get_if<bool>(&value))
        {
          target = *bool_value;
          return true;
        }

        checks_.fail() << "parameter '" << magic_enum::enum_name(key)
                       << "' must be boolean\n";
        return false;
      }

      /**
       * @brief Load a switch parameter that must be provided. A boolean or an
       *        integer 0/1 value is accepted.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the switch value when provided.
       * @return true when the parameter was provided and stored.
       */
      bool requireSwitch(ParametersT key, bool& target)
      {
        if (!parameters_.contains(key))
        {
          checks_.fail() << "missing required parameter '"
                         << magic_enum::enum_name(key) << "'\n";
          return false;
        }

        const auto& value = parameters_.at(key);
        if (const auto* bool_value = std::get_if<bool>(&value))
        {
          target = *bool_value;
          return true;
        }
        if (const auto* index_value = std::get_if<IdxT>(&value);
            index_value && (*index_value == 0 || *index_value == 1))
        {
          target = (*index_value == 1);
          return true;
        }

        checks_.fail() << "parameter '" << magic_enum::enum_name(key)
                       << "' must be bool or 0/1\n";
        return false;
      }

      /**
       * @brief Load an optional integer selector parameter.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the integer value when provided.
       * @return true when the parameter was provided and stored.
       */
      bool loadSelector(ParametersT key, IdxT& target)
      {
        if (!parameters_.contains(key))
        {
          return false;
        }

        const auto& value = parameters_.at(key);
        if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          target = *index_value;
          return true;
        }

        checks_.fail() << "parameter '" << magic_enum::enum_name(key)
                       << "' must be an integer selector\n";
        return false;
      }

    private:
      const std::map<ParametersT, std::variant<bool, RealT, IdxT>>& parameters_;
      ConfigurationChecks&                                          checks_;
    };
  } // namespace Utilities
} // namespace GridKit
