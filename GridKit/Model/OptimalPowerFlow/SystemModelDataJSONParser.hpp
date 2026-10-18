/**
 * @file SystemModelDataJSONParser.hpp
 * @brief JSON parser functions for optimal power flow system model data.
 */

#pragma once

#include <format>
#include <stdexcept>
#include <string>
#include <string_view>

#include <nlohmann/json.hpp>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using json = nlohmann::json;

    /// Member `key` of `object`, which the input format requires
    inline const json& required(const json& object, const char* key, std::string_view owner)
    {
      if (!object.contains(key))
      {
        throw std::invalid_argument(std::format("{}: missing '{}'", owner, key));
      }
      return object.at(key);
    }

    /// Reject member `key` of `object`, which this module does not support
    inline void rejectUnsupported(const json& object, const char* key, std::string_view owner)
    {
      if (object.contains(key))
      {
        throw std::invalid_argument(std::format("{}: '{}' is not supported", owner, key));
      }
    }

    /// Nonnegative bus number below the reserved invalid index
    template <typename IdxT>
    IdxT busNumber(const json& value, std::string_view owner)
    {
      if (!value.is_number_integer()
          || (!value.is_number_unsigned() && value.get<json::number_integer_t>() < 0)
          || value.get<json::number_unsigned_t>() >= INVALID_INDEX<IdxT>)
      {
        throw std::invalid_argument(
            std::format("{}: bus number {} must be a nonnegative integer in the supported range", owner, value.dump()));
      }
      return value.get<IdxT>();
    }

    /// Enumerator named `key`, or throw naming its `kind` and `owner`
    template <Utilities::SizedEnum EnumT>
    EnumT enumKey(const std::string& key, std::string_view kind, std::string_view owner)
    {
      const auto value = Utilities::enum_cast<EnumT>(key);
      if (!value)
      {
        throw std::invalid_argument(std::format("{}: unknown {} '{}'", owner, kind, key));
      }
      return *value;
    }

    /// Value of each entry of `params`, typed as in the input file
    template <typename Parameters, typename RealT, typename IdxT>
    ParameterMap<Parameters, RealT, IdxT> readParameters(const json& params, std::string_view owner)
    {
      ParameterMap<Parameters, RealT, IdxT> parameters;
      for (const auto& [key, value] : params.items())
      {
        const auto parameter = enumKey<Parameters>(key, "parameter", owner);
        if (value.is_boolean())
        {
          parameters[parameter] = value.get<bool>();
        }
        else if (value.is_number_unsigned())
        {
          parameters[parameter] = value.get<IdxT>();
        }
        else if (value.is_number())
        {
          parameters[parameter] = value.get<RealT>();
        }
        else if (value.is_string())
        {
          parameters[parameter] = value.get<std::string>();
        }
        else
        {
          throw std::invalid_argument(
              std::format("{}: parameter '{}' must be a number, a boolean, or a string", owner, key));
        }
      }
      return parameters;
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bus)
    {
      bus.number = busNumber<IdxT>(required(j, "number", "bus"), "bus");
      rejectUnsupported(j, "init", bus.id());
      rejectUnsupported(j, "mon", bus.id());
      bus.parameters = readParameters<BusParameters, RealT, IdxT>(required(j, "params", bus.id()), bus.id());
    }

    template <typename RealT, typename IdxT, typename Parameters, typename Ports>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports>& device)
    {
      const json& id = required(j, "id", "device");
      if (!id.is_string())
      {
        throw std::invalid_argument(std::format("device {}: 'id' must be a string", id.dump()));
      }
      device.id = id.get<std::string>();

      rejectUnsupported(j, "init", device.id);
      rejectUnsupported(j, "mon", device.id);
      device.parameters = readParameters<Parameters, RealT, IdxT>(required(j, "params", device.id), device.id);
      for (const auto& [key, value] : required(j, "ports", device.id).items())
      {
        device.ports[enumKey<Ports>(key, "port", device.id)] = busNumber<IdxT>(value, device.id);
      }
    }

    /// See the `README.md` in `GridKit/Model/OptimalPowerFlow`
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& data)
    {
      using DataT = SystemModelData<RealT, IdxT>;

      constexpr std::string_view owner = "OptimalPowerFlow";

      rejectUnsupported(j, "signals", owner);
      rejectUnsupported(j, "monitors", owner);

      const json& va_base = required(required(j, "params", owner), "va_base", owner);
      if (!va_base.is_number())
      {
        throw std::invalid_argument(std::format("{}: 'va_base' must be a number", owner));
      }
      data.va_base = va_base.get<RealT>();

      for (const auto& bus : required(j, "buses", owner))
      {
        if (required(bus, "class", "bus") != "Bus")
        {
          throw std::invalid_argument(std::format("bus {}: unknown class {}", bus.at("number").dump(), bus.at("class").dump()));
        }
        data.bus.push_back(bus.get<typename DataT::BusDataT>());
      }

      for (const auto& device : required(j, "devices", owner))
      {
        const json& kind = required(device, "class", "device");
        if (kind == "Branch")
        {
          data.branch.push_back(device.get<typename DataT::BranchDataT>());
        }
        else if (kind == "Generator")
        {
          data.generator.push_back(device.get<typename DataT::GeneratorDataT>());
        }
        else if (kind == "Load")
        {
          data.load.push_back(device.get<typename DataT::LoadDataT>());
        }
        else
        {
          throw std::invalid_argument(std::format("{}: unknown class {}", device.value("id", std::string{}), kind.dump()));
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
