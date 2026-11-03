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

#include <magic_enum/magic_enum.hpp>

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

    /// Key into a model data container: a parameter or an initial value
    template <typename KeyT, typename T>
    concept DataKey = std::same_as<KeyT, typename T::Parameters>
                      || std::same_as<KeyT, typename T::InitialValues>;

    /**
     * @brief Reads typed parameters and initial values out of a model data container.
     *
     * A read method returns the provided value or std::nullopt when the
     * value is omitted. The model chooses a default or requires a value.
     * Invalid values are rejected with std::invalid_argument naming the
     * model and key. Real values accept integers and must be finite.
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
      template <DataKey<ModelDataT> KeyT>
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

      /// The provided value, or nullptr when the key is omitted.
      template <DataKey<ModelDataT> KeyT>
      const ValueT* find(KeyT key) const
      {
        if constexpr (std::same_as<KeyT, ParametersT>)
        {
          return lookup(data_.parameters, key);
        }
        else
        {
          return lookup(data_.initial_values, key);
        }
      }

      /// The value of `key` in `values`, or nullptr when it is omitted.
      template <typename MapT, typename KeyT>
      static const ValueT* lookup(const MapT& values, KeyT key)
      {
        const auto entry = values.find(key);
        if (entry == values.end())
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

      /// Reject the model data with a message naming the model and key.
      template <DataKey<ModelDataT> KeyT>
      [[noreturn]] void reject(KeyT key, std::string_view reason) const
      {
        std::string_view kind = "parameter";
        if constexpr (!std::same_as<KeyT, ParametersT>)
        {
          kind = "initial value";
        }
        throw std::invalid_argument(
            std::format("{}: {} '{}' {}", model_, kind, magic_enum::enum_name(key), reason));
      }

      const ModelDataT& data_;
      const char*       model_;
    };
  } // namespace Model
} // namespace GridKit
