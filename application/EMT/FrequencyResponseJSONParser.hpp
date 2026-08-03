/**
 * @file FrequencyResponseJSONParser.hpp
 *
 * @brief JSON parser for the EMT FrequencyResponse application.
 *
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "FrequencyResponseData.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Application
    {
      namespace Detail
      {
        namespace fs = std::filesystem;

        using json = ::nlohmann::json;

        inline std::ifstream openJsonFile(const fs::path& file_path)
        {
          if (!fs::exists(file_path))
          {
            throw std::runtime_error("FrequencyResponse solver file not found: " + file_path.string());
          }

          auto stream = std::ifstream(file_path);
          if (!stream)
          {
            throw std::runtime_error("Failed to open FrequencyResponse solver file: " + file_path.string());
          }
          return stream;
        }

        inline std::string lower(std::string value)
        {
          std::ranges::transform(value, value.begin(), [](unsigned char c)
                                 { return static_cast<char>(std::tolower(c)); });
          return value;
        }

        inline double requirePositive(const json& j, const char* key, const std::string& context)
        {
          const auto value = j.at(key).template get<double>();
          if (!std::isfinite(value) || value <= 0.0)
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be positive");
          }
          return value;
        }

        inline size_t requirePointCount(const json& j, const char* key, const std::string& context)
        {
          const auto value = j.at(key).template get<long long>();
          if (value < 2)
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be at least 2");
          }
          if (static_cast<unsigned long long>(value) > std::numeric_limits<size_t>::max())
          {
            throw std::runtime_error(context + ": \"" + key + "\" is too large");
          }
          return static_cast<size_t>(value);
        }

        inline void rejectUnsupportedVocabulary(const json& j)
        {
          for (const auto* key : {"unit", "omega", "coordinate", "sweep", "plots", "frequency_hz"})
          {
            if (j.contains(key))
            {
              throw std::runtime_error("FrequencyResponse solver JSON does not support \"" + std::string(key) + "\"");
            }
          }
        }

        inline std::optional<FrequencyResponseData::MonitorVariable>
        parseMonitorVariable(const std::string& name)
        {
          using Variable = FrequencyResponseData::MonitorVariable;

          const auto key = lower(name);
          if (key == "r")
          {
            return Variable::R;
          }
          if (key == "l")
          {
            return Variable::L;
          }
          if (key == "g")
          {
            return Variable::G;
          }
          if (key == "c")
          {
            return Variable::C;
          }
          if (key == "tv")
          {
            return Variable::Tv;
          }
          if (key == "ti")
          {
            return Variable::Ti;
          }
          if (key == "alpha")
          {
            return Variable::Alpha;
          }
          if (key == "beta")
          {
            return Variable::Beta;
          }
          if (key == "tau")
          {
            return Variable::Tau;
          }
          if (key == "h")
          {
            return Variable::H;
          }
          if (key == "yc")
          {
            return Variable::Yc;
          }
          if (key == "zc")
          {
            return Variable::Zc;
          }
          return std::nullopt;
        }
      } // namespace Detail

      inline FrequencyResponseData
      parseFrequencyResponseData(const std::filesystem::path& file_path)
      {
        using json = Detail::json;

        const auto j = json::parse(Detail::openJsonFile(file_path));
        Detail::rejectUnsupportedVocabulary(j);

        FrequencyResponseData data;

        data.model = j.at("model").template get<std::filesystem::path>();
        if (data.model.empty())
        {
          throw std::runtime_error("FrequencyResponse solver JSON requires nonempty \"model\"");
        }

        if (j.contains("ida"))
        {
          const auto& ida = j.at("ida");
          if (!ida.is_object())
          {
            throw std::runtime_error("FrequencyResponse \"ida\" must be an object");
          }

          data.ida.suppress_algebraic_error =
              ida.value("suppress_algebraic_error",
                        data.ida.suppress_algebraic_error);
        }

        const auto& frequency = j.at("frequency");
        Detail::rejectUnsupportedVocabulary(frequency);
        data.frequency.start  = Detail::requirePositive(frequency, "start", "frequency");
        data.frequency.stop   = Detail::requirePositive(frequency, "stop", "frequency");
        data.frequency.points = Detail::requirePointCount(frequency, "points", "frequency");
        data.frequency.scale  = Detail::lower(frequency.value("scale", std::string{"log"}));
        if (data.frequency.scale != "log")
        {
          throw std::runtime_error("FrequencyResponse supports only log frequency scale");
        }
        if (data.frequency.stop <= data.frequency.start)
        {
          throw std::runtime_error("frequency: \"stop\" must be greater than \"start\"");
        }

        const auto& output = j.at("output");
        Detail::rejectUnsupportedVocabulary(output);
        data.output_file = output.at("file").template get<std::filesystem::path>();
        if (data.output_file.empty())
        {
          throw std::runtime_error("FrequencyResponse output requires nonempty \"file\"");
        }

        const auto& variables = output.at("variables");
        if (!variables.is_array() || variables.empty())
        {
          throw std::runtime_error("FrequencyResponse output requires a nonempty \"variables\" array");
        }
        for (const auto& raw_variable : variables)
        {
          const auto variable_name = raw_variable.template get<std::string>();
          const auto variable      = Detail::parseMonitorVariable(variable_name);
          if (!variable.has_value())
          {
            throw std::runtime_error("Invalid FrequencyResponse output variable: " + variable_name);
          }
          data.variables.insert(variable.value());
        }

        const auto base_path = file_path.parent_path();
        if (!data.model.is_absolute())
        {
          data.model = base_path / data.model;
        }
        if (!data.output_file.is_absolute())
        {
          data.output_file = base_path / data.output_file;
        }

        return data;
      }
    } // namespace Application
  } // namespace EMT
} // namespace GridKit
