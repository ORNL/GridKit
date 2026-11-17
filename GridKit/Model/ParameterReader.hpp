/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter loading from model data containers.
 */

#pragma once

#include <cmath>
#include <concepts>
#include <format>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <variant>

#include <GridKit/Constants.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Model
  {
    /// Model data container with a map from each parameter to its value
    template <typename T>
    concept ParameterData =
        requires(const T& data) {
          typename T::RealT;
          typename T::IdxT;
          typename T::Parameters;
          data.parameters;
        }
        && Utilities::SizedEnum<typename T::Parameters>;

    /**
     * @brief Reads typed parameters out of a model data container.
     *
     * A read method returns the provided value or std::nullopt when the
     * parameter is omitted. The model chooses a default or requires a value.
     * Invalid values are rejected with std::invalid_argument naming the
     * model and parameter. Real parameters accept integer values and must
     * be finite.
     *
     * @tparam ModelDataT A model data container satisfying ParameterData.
     */
    template <ParameterData ModelDataT>
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
        : data_(data),
          model_(model)
      {
      }

      /// Read an optional real parameter or initial value.
      // TODO: Unlike `ParametersT`, the initial values enum is not a class-level alias because the
      // Optimal Power Flow model data does not define initial values yet. Once it does, add
      // `InitialValues` to `ParameterData`, introduce an `InitialValuesT` alias, and replace this
      // constrained template with plain `ParametersT` and `InitialValuesT` overloads.
      template <typename KeyT>
        requires std::same_as<KeyT, ParametersT> || std::same_as<KeyT, typename ModelDataT::InitialValues>
      [[nodiscard]] std::optional<RealT> readReal(KeyT key) const
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

    private:
      using ParameterMapT = decltype(ModelDataT::parameters);
      using ValueT        = typename ParameterMapT::mapped_type;

      /// The provided value, or nullptr when the parameter is omitted.
      const ValueT* find(ParametersT key) const
      {
        const auto entry = data_.parameters.find(key);
        if (entry == data_.parameters.end())
        {
          return nullptr;
        }
        return &entry->second;
      }

      /// The provided value, or nullptr when the initial value is omitted.
      // TODO: `InitialValuesT` is a constrained template parameter, not a class-level alias like
      // `ParametersT`, until the Optimal Power Flow model data supports initial values. Then make
      // this a plain `find(InitialValuesT)` overload.
      template <typename InitialValuesT>
        requires std::same_as<InitialValuesT, typename ModelDataT::InitialValues>
      const ValueT* find(InitialValuesT key) const
      {
        const auto entry = data_.initial_values.find(key);
        if (entry == data_.initial_values.end())
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
        reject("parameter", Utilities::enum_name(key), reason);
      }

      /// Reject the model data with a message naming the model and initial value.
      template <typename InitialValuesT>
        requires std::same_as<InitialValuesT, typename ModelDataT::InitialValues>
      [[noreturn]] void reject(InitialValuesT key, std::string_view reason) const
      {
        reject("initial value", Utilities::enum_name(key), reason);
      }

      /// Reject the model data with a message naming the model, the kind of entry, and its name.
      [[noreturn]] void reject(std::string_view kind, std::string_view name, std::string_view reason) const
      {
        throw std::invalid_argument(std::format("{}: {} '{}' {}", model_, kind, name, reason));
      }

      const ModelDataT& data_;
      const char*       model_;
    };
  } // namespace Model
} // namespace GridKit
