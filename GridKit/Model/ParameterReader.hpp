/**
 * @file ParameterReader.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Typed parameter and initial value loading from model data containers.
 */

#pragma once

#include <cmath>
#include <format>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
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

    /// Model data container with a map from each initial value to its value
    template <typename T>
    concept InitialValueData =
        requires(const T& data) {
          typename T::RealT;
          typename T::IdxT;
          typename T::InitialValues;
          data.initial_values;
        }
        && Utilities::SizedEnum<typename T::InitialValues>;

    /**
     * @brief Reads typed values out of one map of a model data container.
     *
     * A read method returns the provided value or std::nullopt when the
     * entry is omitted. The model chooses a default or requires a value.
     * Invalid values are rejected with std::invalid_argument naming the
     * model and entry. Real entries accept integer values and must
     * be finite.
     *
     * @tparam KeyT Enumeration over the entries of the map.
     * @tparam RealT Real value type.
     * @tparam IdxT Integer value type.
     */
    template <Utilities::SizedEnum KeyT, typename RealT, typename IdxT>
    class DataReader
    {
    private:
      using ValueT = std::variant<std::string, bool, RealT, IdxT>;
      using MapT   = std::map<KeyT, ValueT>;

    public:
      /**
       * @param[in] values Map to read from.
       * @param[in] model Model name used in rejection messages.
       * @param[in] kind Kind of entry used in rejection messages.
       */
      DataReader(const MapT& values, const char* model, const char* kind)
        : values_(values),
          model_(model),
          kind_(kind)
      {
      }

      /// Read an optional real parameter or initial value.
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

      /// Read an optional non-negative real value.
      [[nodiscard]] std::optional<RealT> readNonNegative(KeyT key) const
      {
        auto value = readReal(key);
        if (value && *value < ZERO<RealT>)
        {
          reject(key, "must be non-negative");
        }
        return value;
      }

      /// Read an optional boolean value.
      [[nodiscard]] std::optional<bool> readBool(KeyT key) const
      {
        return readExact<bool>(key, "must be boolean");
      }

      /// Read an optional integer value.
      [[nodiscard]] std::optional<IdxT> readInteger(KeyT key) const
      {
        return readExact<IdxT>(key, "must be an integer");
      }

    private:
      /// The provided value, or nullptr when the entry is omitted.
      const ValueT* find(KeyT key) const
      {
        const auto entry = values_.find(key);
        if (entry == values_.end())
        {
          return nullptr;
        }
        return &entry->second;
      }

      /// Read an optional entry that must hold exactly a T.
      template <typename T>
      [[nodiscard]] std::optional<T> readExact(KeyT key, std::string_view reason) const
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

      /// Reject the model data with a message naming the model and the entry.
      [[noreturn]] void reject(KeyT key, std::string_view reason) const
      {
        throw std::invalid_argument(
            std::format("{}: {} '{}' {}", model_, kind_, Utilities::enum_name(key), reason));
      }

      const MapT& values_;
      const char* model_;
      const char* kind_;
    };

    /**
     * @brief Reads typed parameters out of a model data container.
     *
     * @tparam ModelDataT A model data container satisfying ParameterData.
     */
    template <ParameterData ModelDataT>
    class ParameterReader
      : public DataReader<typename ModelDataT::Parameters, typename ModelDataT::RealT, typename ModelDataT::IdxT>
    {
    public:
      /**
       * @param[in] data Model data container to read from.
       * @param[in] model Model name used in rejection messages.
       */
      ParameterReader(const ModelDataT& data, const char* model)
        : ParameterReader::DataReader(data.parameters, model, "parameter")
      {
      }
    };

    /**
     * @brief Reads typed initial values out of a model data container.
     *
     * @tparam ModelDataT A model data container satisfying InitialValueData.
     */
    template <InitialValueData ModelDataT>
    class InitialValueReader
      : public DataReader<typename ModelDataT::InitialValues, typename ModelDataT::RealT, typename ModelDataT::IdxT>
    {
    public:
      /**
       * @param[in] data Model data container to read from.
       * @param[in] model Model name used in rejection messages.
       */
      InitialValueReader(const ModelDataT& data, const char* model)
        : InitialValueReader::DataReader(data.initial_values, model, "initial value")
      {
      }
    };
  } // namespace Model
} // namespace GridKit
