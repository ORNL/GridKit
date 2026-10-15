/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter loading from model data containers.
 */

#pragma once

#include <cmath>
#include <format>
#include <optional>
#include <stdexcept>
#include <string>
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
     * A read method returns the provided value or std::nullopt when the
     * parameter is omitted. The model chooses a default or requires a value.
     * Invalid values are rejected with std::invalid_argument naming the
     * model and parameter. Real parameters accept integer values and must
     * be finite.
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

      /// Read an optional real parameter.
      [[nodiscard]] std::optional<RealT> readReal(ParametersT key) const
      {
        const auto* value = find(key);
        if (value == nullptr)
        {
          return std::nullopt;
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

        return real;
      }

      /// Read an optional non-negative real parameter.
      [[nodiscard]] std::optional<RealT> readNonNegative(ParametersT key) const
      {
        auto value = readReal(key);
        if (value && *value < ZERO<RealT>)
        {
          reject(key, "must be non-negative");
        }
        return value;
      }

      /// Read an optional boolean parameter.
      [[nodiscard]] std::optional<bool> readBool(ParametersT key) const
      {
        return readExact<bool>(key, "must be boolean");
      }

      /// Read an optional integer parameter.
      [[nodiscard]] std::optional<IdxT> readInteger(ParametersT key) const
      {
        return readExact<IdxT>(key, "must be an integer");
      }

      /// Read an optional string parameter.
      [[nodiscard]] std::optional<std::string> readString(ParametersT key) const
      {
        return readExact<std::string>(key, "must be a string");
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

      /// Read an optional parameter that must hold exactly a T.
      template <typename T>
      [[nodiscard]] std::optional<T> readExact(ParametersT key, std::string_view reason) const
      {
        const auto* value = find(key);
        if (value == nullptr)
        {
          return std::nullopt;
        }

        const auto* exact_value = std::get_if<T>(value);
        if (exact_value == nullptr)
        {
          reject(key, reason);
        }

        return *exact_value;
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
