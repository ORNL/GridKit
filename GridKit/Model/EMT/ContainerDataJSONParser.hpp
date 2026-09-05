#pragma once

#include <set>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Bus/BusDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceDataJSONParser.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/ContainerData.hpp>
#include <GridKit/Model/EMT/Signal/SignalDataJSONParser.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    inline void validateLocalName(const std::string& name,
                                  const std::string& kind,
                                  const std::string& scope)
    {
      if (name.empty())
      {
        throw std::runtime_error(kind + " name must not be empty in \"" + scope + "\"");
      }
      if (name.find('.') != std::string::npos)
      {
        throw std::runtime_error(kind + " name \"" + name
                                 + "\" must not contain '.' in \"" + scope + "\"");
      }
    }

    /**
     * @brief Parse the signals, devices, and public boundary of one scope.
     */
    template <typename RealT, typename IdxT>
    void parseContainerData(const json&                 j,
                            ContainerData<RealT, IdxT>& data,
                            const std::string&          scope)
    {
      if (j.contains("buses"))
      {
        throw std::runtime_error(
            "Legacy 'buses' is not supported; list Bus entries in 'devices'");
      }

      if (j.contains("inputs"))
      {
        j.at("inputs").get_to(data.inputs);
      }
      if (j.contains("outputs"))
      {
        j.at("outputs").get_to(data.outputs);
      }

      for (const auto& [name, reference] : data.inputs)
      {
        validateLocalName(name, "Input", scope);
        if (reference.empty())
        {
          throw std::runtime_error("Input \"" + name + "\" has an empty reference in \""
                                   + scope + "\"");
        }
        if (data.outputs.contains(name))
        {
          throw std::runtime_error("Boundary name \"" + name
                                   + "\" is both an input and output in \"" + scope + "\"");
        }
      }
      for (const auto& [name, reference] : data.outputs)
      {
        validateLocalName(name, "Output", scope);
        if (reference.empty())
        {
          throw std::runtime_error("Output \"" + name + "\" has an empty reference in \""
                                   + scope + "\"");
        }
      }

      std::set<std::string> local_names;
      for (const auto& [name, reference] : data.inputs)
      {
        static_cast<void>(reference);
        local_names.insert(name);
      }
      if (j.contains("signals"))
      {
        for (const auto& raw_signal : j.at("signals"))
        {
          const auto id = raw_signal.at("id").template get<std::string>();
          validateLocalName(id, "Signal", scope);
          if (!local_names.insert(id).second)
          {
            throw std::runtime_error("Duplicate local name \"" + id + "\" in \"" + scope + "\"");
          }
        }
        j.at("signals").get_to(data.signal);
      }

      const auto& devices = j.at("devices");
      for (const auto& raw_device : devices)
      {
        const auto id = raw_device.at("id").template get<std::string>();
        validateLocalName(id, "Device", scope);
        if (!local_names.insert(id).second)
        {
          throw std::runtime_error("Duplicate local name \"" + id + "\" in \"" + scope + "\"");
        }
      }

      for (const auto& raw_device : devices)
      {
        const auto kind = raw_device.at("class").template get<std::string>();
        if (kind == "Container")
        {
          static const std::set<std::string> allowed_keys{
              "class", "devices", "id", "inputs", "outputs", "signals"};
          for (const auto& [key, value] : raw_device.items())
          {
            static_cast<void>(value);
            if (!allowed_keys.contains(key))
            {
              throw std::runtime_error("Invalid Container field \"" + key + "\" in \""
                                       + scope + "."
                                       + raw_device.at("id").template get<std::string>() + "\"");
            }
          }

          auto& child = data.container.emplace_back();
          raw_device.at("id").get_to(child.id);
          parseContainerData(raw_device, child, scope + "." + child.id);
        }
        else if (kind == "Bus")
        {
          raw_device.get_to(data.bus.emplace_back());
        }
        else if (kind == "DependentVoltageSource")
        {
          raw_device.get_to(data.dependent_voltage_source.emplace_back());
        }
        else if (kind == "VoltageSource")
        {
          raw_device.get_to(data.voltage_source.emplace_back());
        }
        else if (kind == "Machine")
        {
          raw_device.get_to(data.machine.emplace_back());
        }
        else if (kind == "LineLumped")
        {
          raw_device.get_to(data.line_lumped.emplace_back());
        }
        else if (kind == "LoadZ")
        {
          raw_device.get_to(data.loadz.emplace_back());
        }
        else if (kind == "Switch")
        {
          raw_device.get_to(data.sw.emplace_back());
        }
        else if (kind == "Tgov1")
        {
          raw_device.get_to(data.gov.emplace_back());
        }
        else
        {
          Log::error() << "\n\tInvalid device class: \"" << kind << "\". "
                       << "\n\tSee the \"devices\" list in scope \"" << scope << "\"."
                       << std::endl;
          throw std::runtime_error("JSON parser failed");
        }
      }
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, ContainerData<RealT, IdxT>& data)
    {
      j.at("id").get_to(data.id);
      validateLocalName(data.id, "Container", data.id);
      parseContainerData(j, data, data.id);
    }
  } // namespace EMT
} // namespace GridKit
