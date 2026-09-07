#pragma once

#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/JsonValidation.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    /// JSON parser function for the `ComponentData` class and descendants
    ///
    /// Numeric syntax does not determine the model parameter type: typed reads
    /// validate scalar and three-phase values against the requested type.
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
    void from_json(const json&                             j,
                   ComponentData<RealT,
                                 IdxT,
                                 Parameters,
                                 Inputs,
                                 Outputs,
                                 MonitorableVariables>&    c,
                   std::initializer_list<std::string_view> extra_fields = {})
    {
      j.at("class").get_to(c.device_class);
      j.at("id").get_to(c.id);
      const auto context = c.device_class + " \"" + c.id + "\"";
      validateJsonFields(j, context, {"class", "id", "params", "inputs", "outputs", "mon"}, extra_fields);
      c.parameters.clear();
      c.inputs.clear();
      c.outputs.clear();
      c.monitored_variables.clear();

      auto is_index = [](const json& value)
      {
        if (!value.is_number_integer())
          return false;
        if (!value.is_number_unsigned() && value.template get<int64_t>() < 0)
          return false;
        return value.template get<uint64_t>() <= static_cast<uint64_t>(std::numeric_limits<IdxT>::max());
      };
      if (j.contains("params"))
      {
        const auto& params = j.at("params");
        if (!params.is_object())
          throw std::invalid_argument(context + " params must be an object");
        for (const auto& [name, value] : params.items())
        {
          std::optional<Parameters> key;
          if constexpr (magic_enum::enum_count<Parameters>() != 0)
            key = magic_enum::enum_cast<Parameters>(name);
          const auto field = context + " parameter \"" + name + "\"";
          if (!key || name == "SIZE")
            throw std::invalid_argument(field + " is unknown");
          auto& slot = c.parameters[*key];
          if (value.is_boolean())
            slot = value.template get<bool>();
          else if (is_index(value))
            slot = value.template get<IdxT>();
          else if (value.is_number())
            slot = parseFiniteReal<RealT>(value, field);
          else if (value.is_array())
          {
            if (value.size() != 3)
              throw std::invalid_argument(field + " requires a length-3 vector or 3x3 matrix");
            if (value.at(0).is_array())
            {
              ABCMatrix<RealT> matrix{};
              for (size_t n = 0; n < 3; ++n)
              {
                if (!value.at(n).is_array() || value.at(n).size() != 3)
                  throw std::invalid_argument(field + " requires a 3x3 matrix");
                for (size_t k = 0; k < 3; ++k)
                  matrix[n][k] = parseFiniteReal<RealT>(value.at(n).at(k), field);
              }
              slot = matrix;
            }
            else if (std::all_of(value.begin(), value.end(), is_index))
            {
              ABCVector<IdxT> vector{};
              for (size_t n = 0; n < 3; ++n)
                vector[n] = value.at(n).template get<IdxT>();
              slot = vector;
            }
            else
            {
              ABCVector<RealT> vector{};
              for (size_t n = 0; n < 3; ++n)
                vector[n] = parseFiniteReal<RealT>(value.at(n), field);
              slot = vector;
            }
          }
          else
            throw std::invalid_argument(field + " has an invalid value type");
        }
      }

      auto reference = [&context](const json& value, const std::string& name)
      {
        if (!value.is_string() || value.template get<std::string>().empty())
          throw std::invalid_argument(context + " mapping \"" + name + "\" requires a nonempty signal ID");
        return value.template get<std::string>();
      };
      if (j.contains("inputs"))
      {
        if (!j.at("inputs").is_object())
          throw std::invalid_argument(context + " inputs must be an object");
        for (const auto& [name, value] : j.at("inputs").items())
        {
          if (name == "bus" || name == "bus1" || name == "bus2")
          {
            const std::string prefix = name == "bus" ? "v" : "v" + name.substr(3);
            for (const char phase : {'a', 'b', 'c'})
            {
              const auto input = magic_enum::enum_cast<Inputs>(prefix + phase);
              if (!input || j.at("inputs").contains(prefix + phase))
                throw std::invalid_argument(context + " has an invalid or duplicate bus shortcut \"" + name + "\"");
              c.inputs[*input] = reference(value, name) + ".v" + phase;
            }
            continue;
          }
          const auto input = magic_enum::enum_cast<Inputs>(name);
          if (!input || *input == Inputs::SIZE)
            throw std::invalid_argument(context + " has unknown input \"" + name + "\"");
          c.inputs[*input] = reference(value, name);
        }
      }
      if (j.contains("outputs"))
      {
        if (!j.at("outputs").is_object())
          throw std::invalid_argument(context + " outputs must be an object");
        for (const auto& [name, value] : j.at("outputs").items())
        {
          const auto output = magic_enum::enum_cast<Outputs>(name);
          if (!output || *output == Outputs::SIZE)
            throw std::invalid_argument(context + " has unknown output \"" + name + "\"");
          c.outputs[*output] = reference(value, name);
        }
      }
      if (j.contains("mon"))
      {
        if (!j.at("mon").is_array())
          throw std::invalid_argument(context + " mon must be an array");
        for (const auto& value : j.at("mon"))
        {
          const auto name      = value.template get<std::string>();
          const auto monitored = magic_enum::enum_cast<MonitorableVariables>(name);
          if (!monitored || name == "SIZE")
            throw std::invalid_argument(context + " has unknown monitored variable \"" + name + "\"");
          c.monitored_variables.insert(*monitored);
        }
      }
    }
  } // namespace EMT
} // namespace GridKit
