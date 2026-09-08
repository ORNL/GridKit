#pragma once

#include <set>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Bus/BusDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControlDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControlDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/PwmDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributedDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceDataJSONParser.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/ContainerData.hpp>
#include <GridKit/Model/EMT/Operators/Converter/ConverterDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/ModulationDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/ParkDataJSONParser.hpp>
#include <GridKit/Model/EMT/Signal/SignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

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
        if (!j.at("signals").is_array())
          throw std::invalid_argument("Signals must be an array in \"" + scope + "\"");
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
      if (!devices.is_array())
        throw std::invalid_argument("Devices must be an array in \"" + scope + "\"");
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
          validateJsonFields(raw_device, "Container \"" + scope + "." + raw_device.at("id").template get<std::string>() + "\"", {"class", "devices", "id", "inputs", "outputs", "signals"});

          auto& child = data.container.emplace_back();
          raw_device.at("id").get_to(child.id);
          parseContainerData(raw_device, child, scope + "." + child.id);
        }
        else if (kind == "InnerCurrentControl")
        {
          raw_device.get_to(data.inner_current_control.emplace_back());
        }
        else if (kind == "OuterVoltageControl")
        {
          raw_device.get_to(data.outer_voltage_control.emplace_back());
        }
        else if (kind == "Park")
        {
          raw_device.get_to(data.park.emplace_back());
        }
        else if (kind == "Angle")
        {
          raw_device.get_to(data.angle.emplace_back());
        }
        else if (kind == "Modulation")
        {
          raw_device.get_to(data.modulation.emplace_back());
        }
        else if (kind == "PWM")
        {
          raw_device.get_to(data.pwm.emplace_back());
        }
        else if (kind == "DCLink")
        {
          raw_device.get_to(data.dc_link.emplace_back());
        }
        else if (kind == "Converter")
        {
          raw_device.get_to(data.converter.emplace_back());
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
        else if (kind == "LineDistributed")
        {
          raw_device.get_to(data.line_distributed.emplace_back());
        }
        else if (kind == "LoadZ")
        {
          raw_device.get_to(data.loadz.emplace_back());
        }
        else if (kind == "Switch")
        {
          raw_device.get_to(data.sw.emplace_back());
        }
        else if (kind == "Ieeest" || kind == "IEEEST")
        {
          raw_device.get_to(data.ieeest.emplace_back());
        }
        else if (kind == "GastPti" || kind == "GASTPTI" || kind == "GAST")
        {
          raw_device.get_to(data.gastpti.emplace_back());
        }
        else if (kind == "Tgov1")
        {
          raw_device.get_to(data.gov.emplace_back());
        }
        else if (kind == "SexsPti" || kind == "SEXS-PTI" || kind == "SEXS")
        {
          raw_device.get_to(data.sexs_pti.emplace_back());
        }
        else if (kind == "Ieeet1" || kind == "IEEET1")
        {
          raw_device.get_to(data.exciter.emplace_back());
        }
        else
        {
          throw std::invalid_argument("Invalid device class \"" + kind + "\" in \"" + scope + "\"");
        }
      }
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, ContainerData<RealT, IdxT>& data)
    {
      validateJsonFields(j, "Container", {"class", "devices", "id", "inputs", "outputs", "signals"});
      j.at("id").get_to(data.id);
      validateLocalName(data.id, "Container", data.id);
      parseContainerData(j, data, data.id);
    }
  } // namespace EMT
} // namespace GridKit
