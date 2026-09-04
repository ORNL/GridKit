#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function for the `ComponentData` class and descendants
    ///
    /// Scalar values follow the PhasorDynamics convention: booleans, floats,
    /// and integers map to the matching variant alternatives, so real-valued
    /// parameters must be written with a decimal point. A length-3 array maps
    /// to a three-phase vector, integer-valued when every element is an
    /// integer, and a 3x3 nested array maps to a real three-phase matrix.
    template <typename RealT,
              typename IdxT,
              typename Parameters,
              typename Inputs,
              typename Outputs,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Inputs>
               && std::is_enum_v<Outputs>
               && std::is_enum_v<MonitorableVariables>
    void from_json(const json&                          j,
                   ComponentData<RealT,
                                 IdxT,
                                 Parameters,
                                 Inputs,
                                 Outputs,
                                 MonitorableVariables>& c)
    {
      j.at("class").get_to(c.device_class);

      j.at("id").get_to(c.id);

      std::stringstream error_context;
      error_context << "\n\tSee the \"" << c.device_class
                    << "\" device with \"id\": \"" << c.id
                    << "\" in the \"devices\" list of your JSON file.";

      if (j.contains("ports"))
      {
        throw std::runtime_error(
            "Legacy 'ports' is not supported; use 'inputs' and 'outputs'");
      }

      auto parse_vector = [&error_context](const json& value, auto& parameter_slot, const std::string& key) -> bool
      {
        if (value.size() != 3)
        {
          Log::error() << "\n\tInvalid three-phase vector length for \""
                       << key << "\": " << value.size() << "."
                       << error_context.str() << std::endl;
          return false;
        }

        bool all_integers = true;
        for (const auto& element : value)
        {
          if (!element.is_number_integer())
          {
            all_integers = false;
          }
        }

        if (all_integers)
        {
          ABCVector<IdxT> vec{};
          for (size_t n = 0; n < 3; ++n)
          {
            value.at(n).get_to(vec[n]);
          }
          parameter_slot = vec;
          return true;
        }

        ABCVector<RealT> vec{};
        for (size_t n = 0; n < 3; ++n)
        {
          value.at(n).get_to(vec[n]);
        }
        parameter_slot = vec;
        return true;
      };

      auto parse_matrix = [&error_context](const json& value, auto& parameter_slot, const std::string& key) -> bool
      {
        if (value.size() != 3)
        {
          Log::error() << "\n\tInvalid three-phase matrix row count for \""
                       << key << "\": " << value.size() << "."
                       << error_context.str() << std::endl;
          return false;
        }

        ABCMatrix<RealT> mat{};
        for (size_t n = 0; n < 3; ++n)
        {
          const auto& row = value.at(n);
          if (!row.is_array() || row.size() != 3)
          {
            Log::error() << "\n\tInvalid three-phase matrix row for \""
                         << key << "\"." << error_context.str() << std::endl;
            return false;
          }
          for (size_t k = 0; k < 3; ++k)
          {
            row.at(k).get_to(mat[n][k]);
          }
        }
        parameter_slot = mat;
        return true;
      };

      if (j.contains("params"))
      {
        for (auto& raw_parameter : j.at("params").items())
        {
          auto key = magic_enum::enum_cast<Parameters>(raw_parameter.key());
          if (key.has_value())
          {
            // NOTE: this is necessary because it doesn't seem like nlohmann/json
            //       handles std::variant out of the box
            if (raw_parameter.value().is_boolean())
            {
              c.parameters[key.value()] = raw_parameter.value().template get<bool>();
            }
            else if (raw_parameter.value().is_number_float())
            {
              c.parameters[key.value()] = raw_parameter.value().template get<RealT>();
            }
            else if (raw_parameter.value().is_number_integer())
            {
              c.parameters[key.value()] = raw_parameter.value().template get<IdxT>();
            }
            else if (raw_parameter.value().is_array() && !raw_parameter.value().empty()
                     && raw_parameter.value().at(0).is_array())
            {
              parse_matrix(raw_parameter.value(), c.parameters[key.value()], raw_parameter.key());
            }
            else if (raw_parameter.value().is_array())
            {
              parse_vector(raw_parameter.value(), c.parameters[key.value()], raw_parameter.key());
            }
            else
            {
              Log::error() << "\n\tInvalid initial parameter value type: "
                           << "\"" << raw_parameter.key() << "\": "
                           << raw_parameter.value()
                           << " (typed as \"" << raw_parameter.value().type_name()
                           << "\")." << error_context.str() << std::endl;
            }
          }
          else
          {
            Log::error() << "\n\tInitial parameter \"" << raw_parameter.key()
                         << "\" has no value." << error_context.str()
                         << std::endl;
          }
        }
      }

      if (j.contains("inputs"))
      {
        for (auto& raw_input : j.at("inputs").items())
        {
          auto input = magic_enum::enum_cast<Inputs>(raw_input.key());
          if (!input.has_value() || input.value() == Inputs::SIZE)
          {
            Log::error() << "\n\tInvalid input mapping: \"" << raw_input.key()
                         << "\" has no value." << error_context.str()
                         << std::endl;
            throw std::runtime_error("JSON parser failed");
          }
          raw_input.value().get_to(c.inputs[input.value()]);
        }
      }

      if (j.contains("outputs"))
      {
        for (auto& raw_output : j.at("outputs").items())
        {
          auto output = magic_enum::enum_cast<Outputs>(raw_output.key());
          if (!output.has_value() || output.value() == Outputs::SIZE)
          {
            Log::error() << "\n\tInvalid output mapping: \"" << raw_output.key()
                         << "\" has no value." << error_context.str()
                         << std::endl;
            throw std::runtime_error("JSON parser failed");
          }
          raw_output.value().get_to(c.outputs[output.value()]);
        }
      }

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto var_name  = raw_monitored_variable.get<std::string>();
          auto monitored = magic_enum::enum_cast<MonitorableVariables>(var_name);
          if (monitored.has_value())
          {
            c.monitored_variables.insert(monitored.value());
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
  } // namespace EMT
} // namespace GridKit
