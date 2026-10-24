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

#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using json = nlohmann::json;

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
        else
        {
          parameters[parameter] = value.get<std::string>();
        }
      }
      return parameters;
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bus)
    {
      j.at("number").get_to(bus.number);
      bus.parameters = readParameters<BusParameters, RealT, IdxT>(j.at("params"), std::format("bus {}", bus.number));
    }

    template <typename RealT, typename IdxT, typename Parameters, typename Ports>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports>& device)
    {
      j.at("id").get_to(device.id);
      device.parameters = readParameters<Parameters, RealT, IdxT>(j.at("params"), device.id);
      for (const auto& [key, value] : j.at("ports").items())
      {
        value.get_to(device.ports[enumKey<Ports>(key, "port", device.id)]);
      }
    }

    /// See the `README.md` in `GridKit/Model/OptimalPowerFlow`
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& data)
    {
      using DataT = SystemModelData<RealT, IdxT>;

      j.at("params").at("va_base").get_to(data.va_base);
      j.at("buses").get_to(data.bus);

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
          throw std::invalid_argument(std::format("OptimalPowerFlow: unknown device class '{}'", kind));
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
