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

    template <typename RealT, typename IdxT, typename Parameters, typename Ports, typename Signals>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports, Signals>& device)
    {
      j.at("id").get_to(device.id);
      device.parameters = readParameters<Parameters, RealT, IdxT>(j.at("params"), device.id);
      for (const auto& [key, value] : j.at("ports").items())
      {
        if (const auto port = Utilities::enum_cast<Ports>(key))
        {
          value.get_to(device.ports[*port]);
        }
        else
        {
          value.get_to(device.signals[enumKey<Signals>(key, "port", device.id)]);
        }
      }
    }

    /// See the `README.md` in `GridKit/Model/OptimalPowerFlow`
    template <typename RealT, typename IdxT>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& data)
    {
      using DataT = SystemModelData<RealT, IdxT>;

      j.at("params").at("va_base").get_to(data.va_base);
      j.at("buses").get_to(data.bus);

      for (const auto& signal : j.value("signals", json::array()))
      {
        data.signal.push_back(signal.at("signal_id").get<IdxT>());
      }

      for (const auto& device : j.at("devices"))
      {
        const auto kind = device.at("class").get<std::string>();
        if (kind == "Branch")
        {
          data.branch.push_back(device.get<typename DataT::BranchDataT>());
        }
        else if (kind == "Esdc1a")
        {
          data.esdc1a.push_back(device.get<typename DataT::Esdc1aDataT>());
        }
        else if (kind == "GastPti")
        {
          data.gast_pti.push_back(device.get<typename DataT::GastPtiDataT>());
        }
        else if (kind == "GenClassical")
        {
          data.gen_classical.push_back(device.get<typename DataT::GenClassicalDataT>());
        }
        else if (kind == "Generator")
        {
          data.generator.push_back(device.get<typename DataT::GeneratorDataT>());
        }
        else if (kind == "Genrou")
        {
          data.genrou.push_back(device.get<typename DataT::GenrouDataT>());
        }
        else if (kind == "Gensal")
        {
          data.gensal.push_back(device.get<typename DataT::GensalDataT>());
        }
        else if (kind == "Hygov")
        {
          data.hygov.push_back(device.get<typename DataT::HygovDataT>());
        }
        else if (kind == "Ieeet1")
        {
          data.ieeet1.push_back(device.get<typename DataT::Ieeet1DataT>());
        }
        else if (kind == "Load")
        {
          data.load.push_back(device.get<typename DataT::LoadDataT>());
        }
        else if (kind == "Reecb")
        {
          data.reecb.push_back(device.get<typename DataT::ReecbDataT>());
        }
        else if (kind == "Regca")
        {
          data.regca.push_back(device.get<typename DataT::RegcaDataT>());
        }
        else if (kind == "Repca")
        {
          data.repca.push_back(device.get<typename DataT::RepcaDataT>());
        }
        else if (kind == "SexsPti")
        {
          data.sexs_pti.push_back(device.get<typename DataT::SexsPtiDataT>());
        }
        else if (kind == "Tgov1")
        {
          data.tgov1.push_back(device.get<typename DataT::Tgov1DataT>());
        }
        else
        {
          throw std::invalid_argument(std::format("OptimalPowerFlow: unknown device class '{}'", kind));
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
