/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter loading from model data containers.
 */

#pragma once

#include <cmath>
#include <format>
#include <stdexcept>
#include <string_view>
#include <variant>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/ModelData.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Reads typed parameters out of a model data container.
     *
     * A load method stores a provided value in its target and returns true.
     * An omitted parameter leaves the target untouched, so the target's
     * initial value is the default. A require method rejects an omitted
     * parameter. A value of the wrong type cannot produce a valid model, so
     * it is rejected with std::invalid_argument naming the model and the
     * parameter. Real parameters accept integer values and must be finite.
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

      /// Load an optional real parameter.
      bool loadReal(ParametersT key, RealT& target) const
      {
        const auto* value = find(key);
        if (value == nullptr)
        {
          return false;
        }

        RealT real{};
        if (const auto* real_value = std::get_if<RealT>(value))
        {
          real = *real_value;
        }
        else if (const auto* integer_value = std::get_if<IdxT>(value))
        {
          real = static_cast<RealT>(*integer_value);
        }
        else
        {
          reject(key, "must be numeric");
        }

        if (!std::isfinite(real))
        {
          reject(key, "must be finite");
        }

        target = real;
        return true;
      }

      /// Load an optional non-negative real parameter.
      bool loadNonNegative(ParametersT key, RealT& target) const
      {
        if (!loadReal(key, target))
        {
          return false;
        }

        if (target < ZERO<RealT>)
        {
          reject(key, "must be non-negative");
        }
        return true;
      }

      /// Load an optional boolean parameter.
      bool loadBool(ParametersT key, bool& target) const
      {
        return loadExact(key, target, "must be boolean");
      }

      /// Load an optional integer parameter.
      bool loadInteger(ParametersT key, IdxT& target) const
      {
        return loadExact(key, target, "must be an integer");
      }

      /// Read a real parameter that must be provided.
      RealT requireReal(ParametersT key) const
      {
        RealT real{};
        if (!loadReal(key, real))
        {
          reject(key, "is required");
        }
        return real;
      }

      /// Read a boolean parameter that must be provided.
      bool requireBool(ParametersT key) const
      {
        bool flag{};
        if (!loadBool(key, flag))
        {
          reject(key, "is required");
        }
        return flag;
      }

    private:
      using ParameterMapT = decltype(ModelDataT::parameters);
      using ValueT        = typename ParameterMapT::mapped_type;

      /// The provided value, or nullptr when the parameter is omitted.
      const ValueT* find(ParametersT key) const
      {
        const auto entry = parameters_.find(key);
        if (entry == parameters_.end())
        {
          return nullptr;
        }
        return &entry->second;
      }

      /// Load an optional parameter that must hold exactly a T.
      template <typename T>
      bool loadExact(ParametersT key, T& target, std::string_view reason) const
      {
        const auto* value = find(key);
        if (value == nullptr)
        {
          return false;
        }

        const auto* exact_value = std::get_if<T>(value);
        if (exact_value == nullptr)
        {
          reject(key, reason);
        }

        target = *exact_value;
        return true;
      }

      /// Reject the model data with a message naming the model and parameter.
      [[noreturn]] void reject(ParametersT key, std::string_view reason) const
      {
        throw std::invalid_argument(
            std::format("{}: parameter '{}' {}", model_, magic_enum::enum_name(key), reason));
      }

      const ParameterMapT& parameters_;
      const char*          model_;
    };
  } // namespace Model
} // namespace GridKit
