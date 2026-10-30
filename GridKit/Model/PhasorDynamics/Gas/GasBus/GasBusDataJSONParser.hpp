#pragma once

#include <sstream>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBusData.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/String.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      using json = nlohmann::json;
      using Log  = ::GridKit::Utilities::Logger;

      /// JSON parser function implementation for the `GasBusData` type
      ///
      /// See `INPUT_FORMAT.md` in `GridKit/Model/PhasorDynamics`
      template <typename RealT, typename IdxT>
      void from_json(const json& j, GasBusData<RealT, IdxT>& gd)
      {
        j.at("number").get_to(gd.bus_id);
        j.at("class").get_to(gd.device_class);
        j.at("name").get_to(gd.name);

        std::stringstream error_context;
        error_context << "\n\tSee bus number " << gd.bus_id
                      << " (\"name\": \"" << gd.name << "\") "
                      << "in the \"buses\" list of your JSON file.";

        if (!Utilities::isValidIdentifier(gd.name))
        {
          Log::error()
              << "\n\tBus name may contain only ASCII letters, digits, "
                 "spaces, '_', '.' and '-'."
              << error_context.str() << std::endl;
          throw std::runtime_error("JSON parser failed");
        }

        if (j.contains("init"))
        {
          for (auto& raw_value : j.at("init").items())
          {
            if (raw_value.key() == "p")
            {
              raw_value.value().get_to(gd.p0);
            }
            else
            {
              Log::error() << "\n\tInvalid initial value \""
                           << raw_value.key() << "\" in \"init\" section."
                           << error_context.str() << std::endl;
            }
          }
        }

        // Only a GasBus takes an injection; a GasBusInfinite supplies any flow
        if (j.contains("params"))
        {
          using Parameters = typename GasBusData<RealT, IdxT>::Parameters;
          for (auto& raw_parameter : j.at("params").items())
          {
            auto key = Utilities::enum_cast<Parameters>(raw_parameter.key());
            if (key.has_value() && gd.device_class == "GasBus" && raw_parameter.value().is_number())
            {
              gd.parameters[key.value()] = raw_parameter.value().template get<RealT>();
            }
            else
            {
              Log::error() << "\n\tInvalid parameter \""
                           << raw_parameter.key() << "\" in \"params\" section."
                           << error_context.str() << std::endl;
            }
          }
        }

        if (j.contains("mon"))
        {
          using MonitorableVariables = typename GasBusData<RealT, IdxT>::MonitorableVariables;
          for (auto& raw_monitored_variable : j.at("mon"))
          {
            auto var_name  = raw_monitored_variable.get<std::string>();
            auto monitored = Utilities::enum_cast<MonitorableVariables>(var_name, magic_enum::case_insensitive);
            if (monitored.has_value())
            {
              gd.monitored_variables.insert(monitored.value());
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
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
