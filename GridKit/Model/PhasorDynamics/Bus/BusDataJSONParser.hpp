#pragma once

#include <sstream>
#include <stdexcept>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function implementation for the `BusData` type
    ///
    /// See the `README.md` in `GridKit/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bd)
    {
      j.at("name").get_to(bd.name);
      bd.disambiguation_string = bd.name;

      std::stringstream error_context;
      error_context << "\n\tSee bus number " << bd.bus_id
                    << " (\"name\": \"" << bd.name << "\") "
                    << "in the \"buses\" list of your JSON file.";

      if (j.contains("init"))
      {
        for (auto& raw_parameter : j.at("init").items())
        {
          if (raw_parameter.key() == "Vi")
          {
            raw_parameter.value().get_to(bd.Vi0);
          }
          else if (raw_parameter.key() == "Vr")
          {
            raw_parameter.value().get_to(bd.Vr0);
          }
          else
          {
            Log::error() << "\n\tInvalid initial parameter \""
                         << raw_parameter.key() << "\" in \"init\" section."
                         << error_context.str() << std::endl;
          }
        }
      }

      j.at("number").get_to(bd.bus_id);

      auto string_class = j.at("class").get<std::string>();
      if (string_class == "bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
      }
      else if (string_class == "infinite_bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
      }
      else
      {
        Log::error() << "\n\tInvalid bus class: \"" << string_class << "\"."
                     << error_context.str() << std::endl;
      }

      j.at("v_base").get_to(bd.v_base);

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(bd.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(bd.va_base);
      }

      if (j.contains("mon"))
      {
        using magic_enum::case_insensitive;
        using magic_enum::enum_cast;
        using MonitorableVariables = typename BusData<RealT, IdxT>::MonitorableVariables;
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto var_name  = raw_monitored_variable.get<std::string>();
          auto monitored = enum_cast<MonitorableVariables>(var_name, case_insensitive);
          if (monitored.has_value())
          {
            bd.monitored_variables.insert(monitored.value());
          }
          else
          {
            Log::error() << "\n\tInvalid monitored variable: \"" << var_name
                         << "\" in \"mon\" list." << error_context.str()
                         << std::endl;
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
