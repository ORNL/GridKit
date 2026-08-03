/**
 * @file OverheadDataJSONParser.hpp
 *
 * @brief Parser for overhead transmission line documents in the line 1.0
 * format.
 *
 * The format is described by the JSON Schema authored in
 * docs/_ext/gridkit/line_schema.py; tests/Schema/validate_line_inputs.py
 * holds the reference resolver semantics and golden fixtures this parser
 * mirrors. A line document maps physical conductors onto the attachment
 * points of a tower type; conductor and tower types come from the
 * document's own catalog section or from included JSON catalog files.
 * Type and attachment names do not survive resolution.
 */

#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
          const auto& raw = j.at(key);
          if (!raw.is_number())
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be a number");
          }
          const auto value = raw.template get<ScalarT>();
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

        inline std::string
        requireString(const json& j, const char* key, const std::string& context)
        {
          const auto& raw = j.at(key);
          if (!raw.is_string())
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be a string");
          }
          return raw.template get<std::string>();
        }

        inline const json&
        requireObject(const json& j, const char* key, const std::string& context)
        {
          const auto& raw = j.at(key);
          if (!raw.is_object())
          {
            throw std::runtime_error(context + ": \"" + key + "\" must be an object");
          }
          return raw;
        }

        inline void rejectUnknownKeys(const json&                        j,
                                      std::initializer_list<const char*> allowed,
                                      const std::string&                 context)
        {
          for (const auto& entry : j.items())
          {
            const auto& key   = entry.key();
            bool        found = false;
            for (const auto* candidate : allowed)
            {
              if (key == candidate)
              {
                found = true;
                break;
              }
            }
            if (!found)
            {
              throw std::runtime_error(context + " contains unsupported field: " + key);
            }
          }
        }

        /// Catalog and attachment names are lower-case kebab-case.
        inline bool validName(const std::string& name)
        {
          if (name.empty() || name.front() == '-')
          {
            return false;
          }
          for (const char c : name)
          {
            if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '-'))
            {
              return false;
            }
          }
          return true;
        }

        inline std::string
        requireName(const json& j, const char* key, const std::string& context)
        {
          auto name = requireString(j, key, context);
          if (!validName(name))
          {
            throw std::runtime_error(context + ": \"" + key
                                     + "\" must be a lower-case kebab-case name");
          }
          return name;
        }

        inline bool validConductorPhase(const std::string& phase)
        {
          return phase == "a" || phase == "b" || phase == "c" || phase == "n" || phase == "g";
        }

        template <typename ScalarT>
        void validateConductorType(const json& j, const std::string& context)
        {
          rejectUnknownKeys(j,
                            {"description", "source", "radius", "conductivity", "permeability", "weight"},
                            context);
          const auto& radius = requireObject(j, "radius", context);
          rejectUnknownKeys(radius, {"outer", "inner"}, context + " \"radius\"");
          const auto outer = requirePositive<ScalarT>(radius, "outer", context + " \"radius\"");
          if (radius.contains("inner"))
          {
            const auto inner =
                requireNonnegative<ScalarT>(radius, "inner", context + " \"radius\"");
            if (inner >= outer)
            {
              throw std::runtime_error(context + ": inner radius must be below the outer radius");
            }
          }
          requirePositive<ScalarT>(j, "conductivity", context);
          if (j.contains("permeability"))
          {
            // Relative permeability; the lower bound rejects absolute values
            // entered by mistake.
            requireAtLeast<ScalarT>(j, "permeability", context, static_cast<ScalarT>(0.1));
          }
          requirePositive<ScalarT>(j, "weight", context);
        }

        template <typename ScalarT>
        void validateTowerType(const json& j, const std::string& context)
        {
          rejectUnknownKeys(j, {"description", "source", "attachments"}, context);
          const auto& attachments = requireObject(j, "attachments", context);
          if (attachments.empty())
          {
            throw std::runtime_error(context + ": \"attachments\" must not be empty");
          }
          for (const auto& entry : attachments.items())
          {
            const auto point_context = context + " attachment \"" + entry.key() + "\"";
            if (!validName(entry.key()))
            {
              throw std::runtime_error(point_context + ": name must be lower-case kebab-case");
            }
            const auto& point = entry.value();
            if (!point.is_object())
            {
              throw std::runtime_error(point_context + " must be an object");
            }
            rejectUnknownKeys(point, {"x", "h"}, point_context);
            requireFinite<ScalarT>(point, "x", point_context);
            requirePositive<ScalarT>(point, "h", point_context);
          }
        }

        /// Named conductor and tower types merged from a document's catalog
        /// section and its included catalog files.
        struct CatalogTables
        {
          std::map<std::string, json> conductors;
          std::map<std::string, json> towers;
        };

        inline void mergeCatalogSections(CatalogTables&     tables,
                                         const json&        sections,
                                         const std::string& origin)
        {
          const auto merge = [&](const char* kind, std::map<std::string, json>& table)
          {
            if (!sections.contains(kind))
            {
              return;
            }
            const auto& section = requireObject(sections, kind, origin);
            if (section.empty())
            {
              throw std::runtime_error(origin + ": \"" + std::string(kind)
                                       + "\" must not be empty");
            }
            for (const auto& entry : section.items())
            {
              const auto& name = entry.key();
              if (!validName(name))
              {
                throw std::runtime_error(origin + " " + kind + " \"" + name
                                         + "\": name must be lower-case kebab-case");
              }
              if (table.count(name) != 0)
              {
                throw std::runtime_error(std::string(kind) + " name collision: " + name + " ("
                                         + origin + ")");
              }
              table.emplace(name, entry.value());
            }
          };
          merge("conductors", tables.conductors);
          merge("towers", tables.towers);
        }

        inline void loadIncludedCatalogs(CatalogTables&  tables,
                                         const json&     doc,
                                         const fs::path& doc_path)
        {
          if (!doc.contains("include"))
          {
            return;
          }
          const auto& includes = doc.at("include");
          if (!includes.is_array() || includes.empty())
          {
            throw std::runtime_error("Line \"include\" must be a nonempty array");
          }
          for (const auto& raw : includes)
          {
            if (!raw.is_string())
            {
              throw std::runtime_error("Line \"include\" entries must be strings");
            }
            const auto name     = raw.template get<std::string>();
            const auto inc_path = doc_path.parent_path() / name;
            if (!fs::is_regular_file(inc_path))
            {
              throw std::runtime_error("include not found: " + name);
            }
            auto       stream  = openJsonFile(inc_path);
            const auto catalog = json::parse(stream);
            const auto context = "include \"" + name + "\"";
            if (!catalog.is_object())
            {
              throw std::runtime_error(context + " is not a valid catalog");
            }
            rejectUnknownKeys(catalog,
                              {"catalog", "$schema", "name", "description", "source", "conductors", "towers"},
                              context);
            if (requireString(catalog, "catalog", context) != "1.0")
            {
              throw std::runtime_error(context + ": unsupported catalog format version");
            }
            for (const auto* key : {"$schema", "name", "description", "source"})
            {
              if (catalog.contains(key))
              {
                requireString(catalog, key, context);
              }
            }
            if (!catalog.contains("conductors") && !catalog.contains("towers"))
            {
              throw std::runtime_error(context + " defines no conductor or tower types");
            }
            mergeCatalogSections(tables, catalog, context);
          }
        }

        template <typename ScalarT>
        std::vector<PathPoint<ScalarT>> parsePoints(const json& raw_points)
        {
          if (!raw_points.is_array() || raw_points.size() < 2)
          {
            throw std::runtime_error("path \"points\" must be an array with at least two points");
          }

          std::vector<PathPoint<ScalarT>> points;
          points.reserve(raw_points.size());
          for (size_t i = 0; i < raw_points.size(); ++i)
          {
            const auto& raw_point = raw_points.at(i);
            const auto  context   = "points[" + std::to_string(i) + "]";
            if (!raw_point.is_object())
            {
              throw std::runtime_error(context + " must be an object");
            }
            rejectUnknownKeys(raw_point, {"latitude", "longitude"}, context);

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
            points.push_back(point);
          }
          return points;
        }
      } // namespace Detail

      template <typename RealT = double, typename IdxT = size_t>
      OverheadData<RealT, IdxT>
      parseOverheadData(const std::filesystem::path& file_path)
      {
        using json = Detail::json;

        auto       stream = Detail::openJsonFile(file_path);
        const auto j      = json::parse(stream);
        if (!j.is_object())
        {
          throw std::runtime_error("Line document must be an object");
        }

        Detail::rejectUnknownKeys(j,
                                  {"line", "$schema", "name", "description", "source", "include", "catalog", "tower", "conductors", "path", "earth"},
                                  "Line");
        if (Detail::requireString(j, "line", "Line") != "1.0")
        {
          throw std::runtime_error("Unsupported line format version");
        }
        for (const auto* key : {"$schema", "name", "description", "source"})
        {
          if (j.contains(key))
          {
            Detail::requireString(j, key, "Line");
          }
        }

        Detail::CatalogTables tables;
        if (j.contains("catalog"))
        {
          const auto& catalog = Detail::requireObject(j, "catalog", "Line");
          Detail::rejectUnknownKeys(catalog, {"conductors", "towers"}, "catalog");
          if (catalog.empty())
          {
            throw std::runtime_error("Line \"catalog\" must not be empty");
          }
          Detail::mergeCatalogSections(tables, catalog, "catalog");
        }
        Detail::loadIncludedCatalogs(tables, j, file_path);

        for (const auto& [name, ctype] : tables.conductors)
        {
          Detail::validateConductorType<RealT>(ctype, "conductor type \"" + name + "\"");
        }
        for (const auto& [name, ttype] : tables.towers)
        {
          Detail::validateTowerType<RealT>(ttype, "tower type \"" + name + "\"");
        }

        const auto tower_name = Detail::requireName(j, "tower", "Line");
        const auto tower_type = tables.towers.find(tower_name);
        if (tower_type == tables.towers.end())
        {
          throw std::runtime_error("unknown tower type: " + tower_name);
        }
        const auto& attachments = tower_type->second.at("attachments");

        const auto& conductors = j.at("conductors");
        if (!conductors.is_array() || conductors.empty())
        {
          throw std::runtime_error("Line requires a nonempty \"conductors\" array");
        }
        if (conductors.size() > std::numeric_limits<IdxT>::max())
        {
          throw std::runtime_error("Line has too many conductors");
        }

        const auto& path = Detail::requireObject(j, "path", "Line");
        Detail::rejectUnknownKeys(path, {"span", "length", "points"}, "path");
        const auto span = Detail::requirePositive<RealT>(path, "span", "path");
        if (path.contains("length") == path.contains("points"))
        {
          throw std::runtime_error("path requires exactly one of \"length\" or \"points\"");
        }

        OverheadData<RealT, IdxT> data;
        const auto                conductor_count = static_cast<IdxT>(conductors.size());
        data.tower.K                              = conductor_count;
        data.tower.span                           = span;
        data.conductor.K                          = conductor_count;

        if (path.contains("length"))
        {
          data.path.length = Detail::requirePositive<RealT>(path, "length", "path");
        }
        else
        {
          data.path.path = Detail::parsePoints<RealT>(path.at("points"));
        }

        data.tower.position.reserve(conductors.size());
        data.tower.height.reserve(conductors.size());
        data.conductor.radius.reserve(conductors.size());
        data.conductor.inner_radius.reserve(conductors.size());
        data.conductor.sigma.reserve(conductors.size());
        data.conductor.mu.reserve(conductors.size());
        data.conductor.weight.reserve(conductors.size());
        data.conductor.phase.reserve(conductors.size());
        data.conductor.circuit.reserve(conductors.size());

        std::vector<std::optional<RealT>> tension(conductors.size());
        std::set<std::string>             used;
        bool                              has_tension = false;

        for (size_t i = 0; i < conductors.size(); ++i)
        {
          const auto& raw     = conductors.at(i);
          const auto  context = "conductors[" + std::to_string(i) + "]";
          if (!raw.is_object())
          {
            throw std::runtime_error(context + " must be an object");
          }
          Detail::rejectUnknownKeys(raw, {"at", "phase", "circuit", "type", "tension"}, context);

          const auto at         = Detail::requireName(raw, "at", context);
          const auto attachment = attachments.find(at);
          if (attachment == attachments.end())
          {
            throw std::runtime_error(context + " references unknown attachment: " + at);
          }
          if (!used.insert(at).second)
          {
            throw std::runtime_error("attachment used more than once: " + at);
          }

          const auto type_name = Detail::requireName(raw, "type", context);
          const auto type      = tables.conductors.find(type_name);
          if (type == tables.conductors.end())
          {
            throw std::runtime_error(context + " references unknown conductor type: "
                                     + type_name);
          }
          const auto& ctype  = type->second;
          const auto& radius = ctype.at("radius");

          const auto phase = Detail::requireString(raw, "phase", context);
          if (!Detail::validConductorPhase(phase))
          {
            throw std::runtime_error(
                context + ": \"phase\" must be one of \"a\", \"b\", \"c\", \"n\", or \"g\"");
          }

          IdxT circuit{1};
          if (raw.contains("circuit"))
          {
            const auto& raw_circuit = raw.at("circuit");
            if (!raw_circuit.is_number_integer() || raw_circuit.template get<long long>() < 1)
            {
              throw std::runtime_error(context
                                       + ": \"circuit\" must be an integer of at least one");
            }
            circuit = raw_circuit.template get<IdxT>();
          }

          if (raw.contains("tension"))
          {
            tension[i]  = Detail::requirePositive<RealT>(raw, "tension", context);
            has_tension = true;
          }

          data.tower.position.push_back(Detail::requireFinite<RealT>(*attachment, "x", context));
          data.tower.height.push_back(Detail::requirePositive<RealT>(*attachment, "h", context));
          data.conductor.radius.push_back(radius.at("outer").template get<RealT>());
          data.conductor.inner_radius.push_back(
              radius.contains("inner") ? radius.at("inner").template get<RealT>()
                                       : static_cast<RealT>(0.0));
          data.conductor.sigma.push_back(ctype.at("conductivity").template get<RealT>());
          data.conductor.mu.push_back(
              (ctype.contains("permeability") ? ctype.at("permeability").template get<RealT>()
                                              : static_cast<RealT>(1.0))
              * Constants::mu0<RealT>());
          data.conductor.weight.push_back(ctype.at("weight").template get<RealT>());
          data.conductor.phase.push_back(phase);
          data.conductor.circuit.push_back(circuit);
        }

        if (has_tension)
        {
          // Supplied tensions must leave a positive minimum conductor height;
          // rejecting here reports the offending position instead of a later
          // geometry failure inside the parameter models.
          for (size_t i = 0; i < tension.size(); ++i)
          {
            if (!tension[i].has_value())
            {
              continue;
            }
            const RealT a   = tension[i].value() / data.conductor.weight[i];
            const RealT sag = a
                              * (std::cosh(span / (static_cast<RealT>(2.0) * a))
                                 - static_cast<RealT>(1.0));
            if (!(data.tower.height[i] - sag > static_cast<RealT>(0.0)))
            {
              throw std::runtime_error(
                  "tension leaves nonpositive minimum height at x="
                  + std::to_string(static_cast<double>(data.tower.position[i])));
            }
          }
          data.tower.tension = std::move(tension);
        }

        const auto& earth = Detail::requireObject(j, "earth", "Line");
        Detail::rejectUnknownKeys(earth, {"conductivity", "permittivity"}, "earth");
        data.carson.earth_sigma =
            Detail::requireNonnegative<RealT>(earth, "conductivity", "earth");
        // Relative permittivity; the lower bound of one rejects absolute
        // values entered by mistake.
        data.carson.earth_eps =
            (earth.contains("permittivity")
                 ? Detail::requireAtLeast<RealT>(earth,
                                                 "permittivity",
                                                 "earth",
                                                 static_cast<RealT>(1.0))
                 : static_cast<RealT>(1.0))
            * Constants::epsilon0<RealT>();

        Conductor<RealT, IdxT> conductor(data.conductor);
        conductor.initialize();
        Tower<RealT, IdxT> tower(data.tower, conductor);
        tower.initialize();
        Path<RealT, IdxT> parsed_path(data.path, tower);
        parsed_path.initialize();

        return data;
      }
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
