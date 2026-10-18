/**
 * @file SystemModelDataJSONParser.hpp
 * @brief JSON parser functions for optimal power flow system model data.
 */

#pragma once

#include <format>
#include <map>
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

    /// Nonnegative bus number below the reserved invalid index
    template <typename IdxT>
    IdxT busNumber(const json& value)
    {
      if (!value.is_number_integer()
          || (!value.is_number_unsigned() && value.get<json::number_integer_t>() < 0)
          || value.get<json::number_unsigned_t>() >= INVALID_INDEX<IdxT>)
      {
        throw std::invalid_argument("Bus number must be a nonnegative integer in the supported range");
      }
      return value.get<IdxT>();
    }

    /// Enumerator named `key`, or throw naming its `owner`
    template <Utilities::SizedEnum EnumT>
    EnumT enumKey(const std::string& key, const std::string& owner)
    {
      const auto value = Utilities::enum_cast<EnumT>(key);
      if (!value)
      {
        throw std::invalid_argument(std::format("{}: unknown entry \"{}\"", owner, key));
      }
      return *value;
    }

    /// Numeric `params` entries into `parameters`, and the boolean `flag_key` into `flag`
    template <typename Parameters, typename RealT>
    void readParameters(const json&                  params,
                        std::map<Parameters, RealT>& parameters,
                        std::string_view             flag_key,
                        bool&                        flag,
                        const std::string&           owner)
    {
      for (const auto& [key, value] : params.items())
      {
        if (key == flag_key)
        {
          value.get_to(flag);
        }
        else
        {
          parameters[enumKey<Parameters>(key, owner)] = value.get<RealT>();
        }
      }
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bus)
    {
      bus.number = busNumber<IdxT>(j.at("number"));
      readParameters(j.at("params"), bus.parameters, "reference", bus.reference, std::format("bus {}", bus.number));
    }

    template <typename RealT, typename IdxT, typename Parameters, typename Ports>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports>& device)
    {
      j.at("id").get_to(device.id);
      readParameters(j.at("params"), device.parameters, "in_service", device.in_service, device.id);
      for (const auto& [key, value] : j.at("ports").items())
      {
        device.ports[enumKey<Ports>(key, device.id)] = busNumber<IdxT>(value);
      }
    }

    /// See the `README.md` in `GridKit/Model/OptimalPowerFlow`
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& data)
    {
      using DataT = SystemModelData<RealT, IdxT>;

      j.at("params").at("va_base").get_to(data.va_base);

      for (const auto& bus : j.at("buses"))
      {
        if (bus.at("class") != "Bus")
        {
          throw std::invalid_argument(std::format("bus {}: unknown class {}", bus.at("number").dump(), bus.at("class").dump()));
        }
        data.bus.push_back(bus.get<typename DataT::BusDataT>());
      }

      for (const auto& device : j.at("devices"))
      {
        const auto kind = device.at("class").get<std::string>();
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
          throw std::invalid_argument(std::format("{}: unknown class \"{}\"", device.at("id").get<std::string>(), kind));
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
