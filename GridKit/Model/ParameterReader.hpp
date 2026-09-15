/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter loading from model data containers.
 */

#pragma once

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/PhasorDynamics/ModelData.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Reads typed parameters out of a model data container.
     *
     * A provided value with the wrong type or a non-finite value cannot
     * produce a valid model, so every method throws std::invalid_argument
     * at the first such value; the message names the model and the
     * parameter. An omitted optional parameter leaves the target untouched.
     * Real parameters accept integer values.
     *
     * @tparam ModelDataT A model data container satisfying
     *         PhasorDynamics::ModelData.
     */
    template <PhasorDynamics::ModelData ModelDataT>
    class ParameterReader
    {
    public:
      using RealT       = typename ModelDataT::RealT;
      using IdxT        = typename ModelDataT::IdxT;
      using ParametersT = typename ModelDataT::Parameters;

      /**
       * @param[in] data Model data container to read from.
       * @param[in] model Model name used in rejection messages.
       */
      ParameterReader(const ModelDataT& data, const char* model)
        : parameters_(data.parameters),
          model_(model)
      {
      }

      /**
       * @brief Load an optional real parameter. Integer values are accepted.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the finite numeric value when provided.
       * @return true when the parameter was provided.
       */
      bool loadReal(ParametersT key, RealT& target) const
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
          reject(key, "must be numeric");
        }

        if (!std::isfinite(parsed_value))
        {
          reject(key, "must be finite");
        }

        target = parsed_value;
        return true;
      }

      /**
       * @brief Load a real parameter that must be provided.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the finite numeric value.
       * @return true, since a missing parameter is rejected.
       */
      bool requireReal(ParametersT key, RealT& target) const
      {
        if (!parameters_.contains(key))
        {
          reject(key, "is required");
        }
        return loadReal(key, target);
      }

      /**
       * @brief Load an optional boolean switch parameter.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the boolean value when provided.
       * @return true when the parameter was provided.
       */
      bool loadSwitch(ParametersT key, bool& target) const
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

        reject(key, "must be boolean");
      }

      /**
       * @brief Load a switch parameter that must be provided. A boolean or an
       *        integer 0/1 value is accepted.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the switch value.
       * @return true, since a missing parameter is rejected.
       */
      bool requireSwitch(ParametersT key, bool& target) const
      {
        if (!parameters_.contains(key))
        {
          reject(key, "is required");
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

        reject(key, "must be bool or 0/1");
      }

      /**
       * @brief Load an optional integer selector parameter.
       *
       * @param[in] key Parameter to look up.
       * @param[out] target Stores the integer value when provided.
       * @return true when the parameter was provided.
       */
      bool loadSelector(ParametersT key, IdxT& target) const
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

        reject(key, "must be an integer selector");
      }

    private:
      /// Reject the model data with a message naming the model and parameter.
      [[noreturn]] void reject(ParametersT key, const char* reason) const
      {
        throw std::invalid_argument(std::string(model_) + ": parameter '"
                                    + std::string(magic_enum::enum_name(key))
                                    + "' " + reason);
      }

      const std::map<ParametersT, std::variant<bool, RealT, IdxT>>& parameters_;
      const char*                                                   model_;
    };
  } // namespace Model
} // namespace GridKit
