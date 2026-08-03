/**
 * @file OverheadDataJSONParser.hpp
 *
 * @brief JSON parser for physical overhead transmission line data.
 *
 */

#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/OverheadData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      namespace Detail
      {
        namespace fs = std::filesystem;

        using json = ::nlohmann::json;

        inline std::ifstream openJsonFile(const fs::path& file_path)
        {
          if (!fs::exists(file_path))
          {
            throw std::runtime_error("Line model file not found: " + file_path.string());
          }

          auto stream = std::ifstream(file_path);
          if (!stream)
          {
            throw std::runtime_error("Failed to open line model file: " + file_path.string());
          }
          return stream;
        }

        template <typename ScalarT>
        ScalarT requireFinite(const json& j, const char* key, const std::string& context)
        {
          const auto value = j.at(key).template get<ScalarT>();
          if (!std::isfinite(static_cast<double>(value)))
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be finite");
          }
          return value;
        }

        template <typename ScalarT>
        ScalarT requirePositive(const json& j, const char* key, const std::string& context)
        {
          const auto value = requireFinite<ScalarT>(j, key, context);
          if (value <= static_cast<ScalarT>(0.0))
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be positive");
          }
          return value;
        }

        template <typename ScalarT>
        ScalarT requireNonnegative(const json& j, const char* key, const std::string& context)
        {
          const auto value = requireFinite<ScalarT>(j, key, context);
          if (value < static_cast<ScalarT>(0.0))
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be nonnegative");
          }
          return value;
        }

        template <typename ScalarT>
        ScalarT
        requireAtLeast(const json& j, const char* key, const std::string& context, ScalarT lower)
        {
          const auto value = requireFinite<ScalarT>(j, key, context);
          if (value < lower)
          {
            throw std::runtime_error(context + ": \"" + key + "\" is below the allowed minimum");
          }
          return value;
        }

        inline bool validConductorPhase(const std::string& phase)
        {
          return phase == "a" || phase == "b" || phase == "c" || phase == "n" || phase == "g";
        }

        template <typename ScalarT>
        ScalarT requireBounded(const json&        j,
                               const char*        key,
                               const std::string& context,
                               ScalarT            lower,
                               ScalarT            upper)
        {
          const auto value = requireFinite<ScalarT>(j, key, context);
          if (value < lower || value > upper)
          {
            throw std::runtime_error(context + ": \"" + key + "\" is out of range");
          }
          return value;
        }

        template <typename ScalarT>
        std::vector<PathPoint<ScalarT>> parsePath(const json& raw_path)
        {
          if (!raw_path.is_array() || raw_path.size() < 2)
          {
            throw std::runtime_error("Overhead \"path\" must be an array with at least two points");
          }

          std::vector<PathPoint<ScalarT>> path;
          path.reserve(raw_path.size());
          for (size_t i = 0; i < raw_path.size(); ++i)
          {
            const auto& raw_point = raw_path.at(i);
            const auto  context   = "path[" + std::to_string(i) + "]";
            if (!raw_point.is_object())
            {
              throw std::runtime_error(context + " must be an object");
            }

            PathPoint<ScalarT> point;
            point.latitude  = requireBounded<ScalarT>(raw_point,
                                                     "latitude",
                                                     context,
                                                     static_cast<ScalarT>(-90.0),
                                                     static_cast<ScalarT>(90.0));
            point.longitude = requireBounded<ScalarT>(raw_point,
                                                      "longitude",
                                                      context,
                                                      static_cast<ScalarT>(-180.0),
                                                      static_cast<ScalarT>(180.0));
            path.push_back(point);
          }
          return path;
        }

        template <typename ScalarT>
        std::vector<ScalarT>
        requireFiniteVector(const json&        j,
                            const char*        key,
                            size_t             expected_size,
                            const std::string& context)
        {
          const auto& raw_values = j.at(key);
          if (!raw_values.is_array() || raw_values.size() != expected_size)
          {
            throw std::runtime_error(context + ": \"" + key
                                     + "\" size must match conductor count");
          }

          std::vector<ScalarT> values;
          values.reserve(expected_size);
          for (size_t i = 0; i < raw_values.size(); ++i)
          {
            const auto value = raw_values.at(i).template get<ScalarT>();
            if (!std::isfinite(static_cast<double>(value)))
            {
              throw std::runtime_error(context + ": \"" + key + "\" entries must be finite");
            }
            values.push_back(value);
          }
          return values;
        }

        template <typename ScalarT>
        std::vector<ScalarT>
        requirePositiveVector(const json&        j,
                              const char*        key,
                              size_t             expected_size,
                              const std::string& context)
        {
          auto values = requireFiniteVector<ScalarT>(j, key, expected_size, context);
          for (const auto value : values)
          {
            if (value <= static_cast<ScalarT>(0.0))
            {
              throw std::runtime_error(context + ": \"" + key + "\" entries must be positive");
            }
          }
          return values;
        }

      } // namespace Detail

      template <typename RealT = double, typename IdxT = size_t>
      OverheadData<RealT, IdxT>
      parseOverheadData(const std::filesystem::path& file_path)
      {
        using json = Detail::json;

        const auto j = json::parse(Detail::openJsonFile(file_path));

        const auto model_class = j.at("class").template get<std::string>();
        if (model_class != "Overhead")
        {
          throw std::runtime_error("Unsupported line model class: " + model_class);
        }
        if (j.contains("skin_effect"))
        {
          throw std::runtime_error("Overhead line JSON does not support \"skin_effect\"");
        }
        if (j.contains("span"))
        {
          throw std::runtime_error("Overhead line JSON requires \"span\" under \"tower\"");
        }
        if (j.contains("segment"))
        {
          throw std::runtime_error("Overhead line JSON does not support \"segment\"");
        }
        if (j.contains("earth"))
        {
          throw std::runtime_error(
              "Overhead line JSON uses top-level \"earth_conductivity\" and \"earth_permittivity\"");
        }

        const auto& conductors = j.at("conductors");
        if (!conductors.is_array() || conductors.empty())
        {
          throw std::runtime_error("Overhead line JSON requires a nonempty \"conductors\" array");
        }

        OverheadData<RealT, IdxT> data;
        if (j.contains("length"))
        {
          data.path.length = Detail::requirePositive<RealT>(j, "length", "Overhead");
        }
        if (j.contains("path"))
        {
          data.path.path = Detail::parsePath<RealT>(j.at("path"));
        }
        if (!data.path.length.has_value() && data.path.path.empty())
        {
          throw std::runtime_error("Overhead line JSON requires \"length\" or \"path\"");
        }
        if (conductors.size() > std::numeric_limits<IdxT>::max())
        {
          throw std::runtime_error("Overhead line JSON has too many conductors");
        }

        const auto conductor_count = static_cast<IdxT>(conductors.size());
        data.tower.K               = conductor_count;
        data.conductor.K           = conductor_count;

        const auto& raw_tower = j.at("tower");
        if (!raw_tower.is_object())
        {
          throw std::runtime_error("Overhead \"tower\" must be an object");
        }
        for (const auto& entry : raw_tower.items())
        {
          const auto& key = entry.key();
          if (key != "x" && key != "height" && key != "span" && key != "tension")
          {
            throw std::runtime_error("Overhead \"tower\" contains unsupported field: " + key);
          }
        }
        data.tower.position =
            Detail::requireFiniteVector<RealT>(raw_tower, "x", conductors.size(), "tower");
        data.tower.height =
            Detail::requirePositiveVector<RealT>(raw_tower, "height", conductors.size(), "tower");
        data.tower.span = Detail::requirePositive<RealT>(raw_tower, "span", "tower");
        if (raw_tower.contains("tension"))
        {
          data.tower.tension =
              Detail::requirePositiveVector<RealT>(raw_tower,
                                                   "tension",
                                                   conductors.size(),
                                                   "tower");
        }
        data.conductor.radius.reserve(conductors.size());
        data.conductor.inner_radius.reserve(conductors.size());
        data.conductor.sigma.reserve(conductors.size());
        data.conductor.mu.reserve(conductors.size());
        const bool has_weight = conductors.at(0).contains("weight");
        const bool has_phase  = conductors.at(0).contains("phase");
        if (!has_weight)
        {
          throw std::runtime_error("Overhead line JSON requires conductor weights");
        }
        data.conductor.weight.reserve(conductors.size());
        if (has_phase)
        {
          data.conductor.phase.reserve(conductors.size());
        }

        for (size_t i = 0; i < conductors.size(); ++i)
        {
          const auto& raw     = conductors.at(i);
          const auto  context = "conductors[" + std::to_string(i) + "]";
          if (raw.contains("x") || raw.contains("height"))
          {
            throw std::runtime_error(context + ": tower placement belongs in \"tower\"");
          }
          const auto radius       = Detail::requirePositive<RealT>(raw, "radius", context);
          const auto inner_radius = Detail::requireNonnegative<RealT>(raw, "inner_radius", context);
          if (inner_radius >= radius)
          {
            throw std::runtime_error(context + ": \"inner_radius\" must be less than \"radius\"");
          }

          data.conductor.radius.push_back(radius);
          data.conductor.inner_radius.push_back(inner_radius);
          data.conductor.sigma.push_back(Detail::requirePositive<RealT>(raw, "conductivity", context));
          data.conductor.mu.push_back(Detail::requirePositive<RealT>(raw, "permeability", context));
          if (raw.contains("weight") != has_weight)
          {
            throw std::runtime_error(
                "Overhead line JSON conductor weights must be specified for all conductors or none");
          }
          data.conductor.weight.push_back(Detail::requirePositive<RealT>(raw, "weight", context));
          if (raw.contains("phase") != has_phase)
          {
            throw std::runtime_error(
                "Overhead line JSON conductor phases must be specified for all conductors or none");
          }
          if (has_phase)
          {
            auto phase = raw.at("phase").template get<std::string>();
            if (!Detail::validConductorPhase(phase))
            {
              throw std::runtime_error(context
                                       + ": \"phase\" must be one of \"a\", \"b\", \"c\", \"n\", or \"g\"");
            }
            data.conductor.phase.push_back(std::move(phase));
          }
        }

        Conductor<RealT, IdxT> conductor(data.conductor);
        conductor.initialize();
        Tower<RealT, IdxT> tower(data.tower, conductor);
        tower.initialize();
        Path<RealT, IdxT> parsed_path(data.path, tower);
        parsed_path.initialize();

        data.carson.earth_sigma =
            Detail::requireNonnegative<RealT>(j, "earth_conductivity", "Overhead");
        data.carson.earth_eps =
            Detail::requireAtLeast<RealT>(j,
                                          "earth_permittivity",
                                          "Overhead",
                                          Constants::epsilon0<RealT>());

        return data;
      }
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
