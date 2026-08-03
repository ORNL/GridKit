#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>
#include <GridKit/Model/EMT/Parameters/OverheadDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <application/EMT/FrequencyResponse/FrequencyResponseJSONParser.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;
  namespace fs      = std::filesystem;

  constexpr scalar_type pi   = GridKit::EMT::Constants::pi<scalar_type>();
  constexpr scalar_type mu0  = GridKit::EMT::Constants::mu0<scalar_type>();
  constexpr scalar_type eps0 = GridKit::EMT::Constants::epsilon0<scalar_type>();

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-12)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  fs::path testPath(const std::string& file_name)
  {
    return fs::temp_directory_path() / file_name;
  }

  void writeFile(const fs::path& file_path, const std::string& contents)
  {
    std::ofstream stream(file_path);
    stream << contents;
  }

  scalar_type parsedPathLength(
      const GridKit::EMT::Parameters::OverheadData<scalar_type, index_type>& data)
  {
    GridKit::EMT::Parameters::Conductor<scalar_type, index_type> conductor(data.conductor);
    conductor.initialize();
    GridKit::EMT::Parameters::Tower<scalar_type, index_type> tower(data.tower, conductor);
    tower.initialize();
    GridKit::EMT::Parameters::Path<scalar_type, index_type> path(data.path, tower);
    path.initialize();
    return path.length();
  }

  scalar_type parsedTowerSaggedToSpanRatio(
      const GridKit::EMT::Parameters::OverheadData<scalar_type, index_type>& data)
  {
    GridKit::EMT::Parameters::Conductor<scalar_type, index_type> conductor(data.conductor);
    conductor.initialize();
    GridKit::EMT::Parameters::Tower<scalar_type, index_type> tower(data.tower, conductor);
    tower.initialize();
    return tower.saggedToSpanRatio();
  }

  template <typename FuncT>
  bool throws(FuncT&& func)
  {
    try
    {
      func();
    }
    catch (const std::exception&)
    {
      return true;
    }
    return false;
  }

  /// A line 1.0 document assembled section by section, so a test can replace
  /// or drop one member. An empty member is omitted from the rendered
  /// document. Mirrors the mutation pattern of
  /// tests/Schema/validate_line_inputs.py.
  struct LineDoc
  {
    std::string line    = R"("1.0")";
    std::string catalog = R"({
      "conductors": {
        "wire": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 }
      },
      "towers": {
        "flat": {
          "attachments": {
            "left":  { "x":  0.0, "h": 10.0 },
            "right": { "x": 30.0, "h": 10.0 }
          }
        }
      }
    })";
    std::string include;
    std::string tower      = R"("flat")";
    std::string conductors = R"([{ "at": "left", "phase": "a", "type": "wire" }])";
    std::string path       = R"({ "span": 300.0, "length": 100000.0 })";
    std::string earth      = R"({ "conductivity": 1.0e-6 })";
    std::string extra;
  };

  std::string render(const LineDoc& doc)
  {
    std::string rendered = "{";
    const auto  member   = [&](const char* key, const std::string& value)
    {
      if (value.empty())
      {
        return;
      }
      if (rendered.size() > 1)
      {
        rendered += ", ";
      }
      rendered += "\"" + std::string(key) + "\": " + value;
    };
    member("line", doc.line);
    member("catalog", doc.catalog);
    member("include", doc.include);
    member("tower", doc.tower);
    member("conductors", doc.conductors);
    member("path", doc.path);
    member("earth", doc.earth);
    if (!doc.extra.empty())
    {
      rendered += ", " + doc.extra;
    }
    rendered += "}";
    return rendered;
  }

  bool rejects(const std::string& file_name, const LineDoc& doc)
  {
    using GridKit::EMT::Parameters::parseOverheadData;

    const auto path = testPath(file_name);
    writeFile(path, render(doc));
    const bool thrown = throws([&]
                               { parseOverheadData<scalar_type, index_type>(path); });
    fs::remove(path);
    return thrown;
  }

  GridKit::Testing::TestOutcome line_parser_valid()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    LineDoc doc;
    doc.catalog    = R"({
      "conductors": {
        "first":  { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 },
        "second": { "radius": { "outer": 0.02, "inner": 0.005 }, "conductivity": 5.8e7,
                    "permeability": 2.0, "weight": 12.0 }
      },
      "towers": {
        "flat": {
          "attachments": {
            "left":  { "x":  0.0, "h": 10.0 },
            "right": { "x": 30.0, "h": 10.0 }
          }
        }
      }
    })";
    doc.conductors = R"([
      { "at": "left",  "phase": "a", "type": "first", "tension": 200000.0 },
      { "at": "right", "phase": "b", "circuit": 2, "type": "second" }
    ])";
    doc.earth      = R"({ "conductivity": 1.0e-6, "permittivity": 2.5 })";

    const auto path = testPath("gridkit_line_valid.line.json");
    writeFile(path, render(doc));

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    TestStatus success  = true;
    success            *= (data.tower.K == 2);
    success            *= (data.conductor.K == 2);
    success            *= close(data.tower.position[0], 0.0);
    success            *= close(data.tower.position[1], 30.0);
    success            *= close(data.tower.height[0], 10.0);
    success            *= close(data.tower.span, 300.0);
    success            *= data.path.length.has_value();
    success            *= close(data.path.length.value(), 100000.0);
    success            *= data.path.path.empty();
    success            *= (data.tower.tension.size() == 2);
    success            *= data.tower.tension[0].has_value();
    success            *= close(data.tower.tension[0].value(), 200000.0);
    success            *= !data.tower.tension[1].has_value();
    success            *= close(data.conductor.radius[0], 0.01);
    success            *= close(data.conductor.radius[1], 0.02);
    success            *= close(data.conductor.inner_radius[0], 0.0);
    success            *= close(data.conductor.inner_radius[1], 0.005);
    success            *= close(data.conductor.sigma[0], 3.7e7);
    success            *= close(data.conductor.sigma[1], 5.8e7);
    success            *= close(data.conductor.mu[0], mu0);
    success            *= close(data.conductor.mu[1], 2.0 * mu0);
    success            *= close(data.conductor.weight[0], 11.0);
    success            *= close(data.conductor.weight[1], 12.0);
    success            *= (data.conductor.phase[0] == "a");
    success            *= (data.conductor.phase[1] == "b");
    success            *= (data.conductor.circuit[0] == 1);
    success            *= (data.conductor.circuit[1] == 2);
    success            *= close(data.carson.earth_sigma, 1.0e-6);
    success            *= close(data.carson.earth_eps, 2.5 * eps0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome line_parser_gis_path()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    LineDoc doc;
    doc.path = R"({
      "span": 300.0,
      "points": [
        { "latitude": 0.0, "longitude": 0.0 },
        { "latitude": 0.0, "longitude": 1.0 }
      ]
    })";

    const auto path = testPath("gridkit_line_gis.line.json");
    writeFile(path, render(doc));

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    const scalar_type ref =
        GridKit::EMT::Parameters::Path<scalar_type, index_type>::earthRadius() * pi / 180.0
        * parsedTowerSaggedToSpanRatio(data);

    TestStatus success  = true;
    success            *= !data.path.length.has_value();
    success            *= (data.path.path.size() == 2);
    success            *= data.tower.tension.empty();
    success            *= close(parsedPathLength(data), ref, 1.0e-12);

    return success.report(__func__);
  }

  // The golden coordinates below duplicate the GOLDEN table of
  // tests/Schema/validate_line_inputs.py, so the C++ loader and the Python
  // reference resolver are pinned to the same fixtures.
  GridKit::Testing::TestOutcome line_parser_golden_wood_pole()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    // examples/EMT/Lines/69kv-wood-pole.line.json, document catalog only.
    const auto path = testPath("gridkit_line_wood_pole.line.json");
    writeFile(path, R"({
      "line": "1.0",
      "catalog": {
        "conductors": {
          "phase-acsr":   { "radius": { "outer": 0.00915 }, "conductivity": 3.5e7, "weight": 6.7 },
          "neutral-acsr": { "radius": { "outer": 0.00501 }, "conductivity": 3.5e7, "weight": 2.1 }
        },
        "towers": {
          "69kv-wood-pole": {
            "attachments": {
              "top":     { "x": -1.0668, "h": 14.9352 },
              "middle":  { "x":  1.0668, "h": 13.7160 },
              "bottom":  { "x": -1.0668, "h": 12.4968 },
              "neutral": { "x":  0.0,    "h": 10.3632 }
            }
          }
        }
      },
      "tower": "69kv-wood-pole",
      "conductors": [
        { "at": "top",     "phase": "a", "type": "phase-acsr" },
        { "at": "middle",  "phase": "b", "type": "phase-acsr" },
        { "at": "bottom",  "phase": "c", "type": "phase-acsr" },
        { "at": "neutral", "phase": "n", "type": "neutral-acsr" }
      ],
      "path": { "span": 90.0, "length": 16000.0 },
      "earth": { "conductivity": 0.01 }
    })");

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    const std::vector<scalar_type> x = {-1.0668, 1.0668, -1.0668, 0.0};
    const std::vector<scalar_type> h = {14.9352, 13.7160, 12.4968, 10.3632};

    TestStatus success  = true;
    success            *= (data.tower.K == 4);
    for (size_t i = 0; i < x.size(); ++i)
    {
      success *= (data.tower.position[i] == x[i]);
      success *= (data.tower.height[i] == h[i]);
    }
    success *= (data.conductor.phase[3] == "n");
    success *= data.tower.tension.empty();

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome line_parser_included_catalog()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    // examples/EMT/Lines/345kv-horizontal.line.json and the shipped catalog,
    // reduced to the types the document references.
    const auto catalog_path = testPath("gridkit_line_included.catalog.json");
    writeFile(catalog_path, R"({
      "catalog": "1.0",
      "name": "Included test catalog",
      "conductors": {
        "drake-acsr": {
          "radius": { "outer": 0.01407, "inner": 0.00514 },
          "conductivity": 3.5e7,
          "weight": 16.0
        },
        "ehs-3-8-steel": {
          "radius": { "outer": 0.00457 },
          "conductivity": 5.8e6,
          "permeability": 1.0,
          "weight": 4.0
        }
      },
      "towers": {
        "345kv-h-frame": {
          "attachments": {
            "left-1":       { "x": -7.4285, "h": 24.0 },
            "left-2":       { "x": -6.9715, "h": 24.0 },
            "center-1":     { "x": -0.2285, "h": 24.0 },
            "center-2":     { "x":  0.2285, "h": 24.0 },
            "right-1":      { "x":  6.9715, "h": 24.0 },
            "right-2":      { "x":  7.4285, "h": 24.0 },
            "shield-left":  { "x": -4.3,    "h": 29.5 },
            "shield-right": { "x":  4.3,    "h": 29.5 }
          }
        }
      }
    })");

    const auto path = testPath("gridkit_line_included.line.json");
    writeFile(path, R"({
      "line": "1.0",
      "include": ["gridkit_line_included.catalog.json"],
      "tower": "345kv-h-frame",
      "conductors": [
        { "at": "left-1",       "phase": "a", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "left-2",       "phase": "a", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "center-1",     "phase": "b", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "center-2",     "phase": "b", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "right-1",      "phase": "c", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "right-2",      "phase": "c", "type": "drake-acsr", "tension": 28000.0 },
        { "at": "shield-left",  "phase": "g", "type": "ehs-3-8-steel", "tension": 12000.0 },
        { "at": "shield-right", "phase": "g", "type": "ehs-3-8-steel", "tension": 12000.0 }
      ],
      "path": { "span": 350.0, "length": 100000.0 },
      "earth": { "conductivity": 0.01 }
    })");

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);
    fs::remove(catalog_path);

    const std::vector<scalar_type> x = {-7.4285, -6.9715, -0.2285, 0.2285, 6.9715, 7.4285, -4.3, 4.3};
    const std::vector<scalar_type> h = {24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 29.5, 29.5};

    TestStatus success  = true;
    success            *= (data.tower.K == 8);
    for (size_t i = 0; i < x.size(); ++i)
    {
      success *= (data.tower.position[i] == x[i]);
      success *= (data.tower.height[i] == h[i]);
    }
    success *= (data.tower.tension.size() == 8);
    success *= close(data.tower.tension[0].value(), 28000.0);
    success *= close(data.tower.tension[7].value(), 12000.0);
    success *= close(data.conductor.radius[0], 0.01407);
    success *= close(data.conductor.inner_radius[0], 0.00514);
    success *= close(data.conductor.mu[6], mu0);
    success *= (data.conductor.phase[6] == "g");

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome line_parser_invalid_documents()
  {
    using GridKit::Testing::TestStatus;

    TestStatus success = true;

    // Ports of the schema-level negative cases of
    // tests/Schema/validate_line_inputs.py.
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0, "length": 100000.0,
                      "points": [{ "latitude": 0.0, "longitude": 0.0 },
                                 { "latitude": 1.0, "longitude": 1.0 }] })";
      success  *= rejects("gridkit_line_both_length_points.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0 })";
      success  *= rejects("gridkit_line_neither_length_points.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "length": 100000.0 })";
      success  *= rejects("gridkit_line_no_span.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([{ "at": "left", "phase": "s", "type": "wire" }])";
      success        *= rejects("gridkit_line_bad_phase.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors =
          R"([{ "at": "left", "phase": "a", "type": "wire", "radius": { "outer": 0.01 } }])";
      success *= rejects("gridkit_line_inline_conductor_data.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([{ "phase": "a", "type": "wire" }])";
      success        *= rejects("gridkit_line_no_attachment.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([{ "at": "left", "phase": "a" }])";
      success        *= rejects("gridkit_line_no_type.line.json", doc);
    }
    {
      LineDoc doc;
      doc.tower  = R"({ "x": [0.0], "height": [10.0] })";
      success   *= rejects("gridkit_line_inline_tower.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([{ "at": "left", "phase": "a", "type": "wire", "tension": 0.0 }])";
      success        *= rejects("gridkit_line_zero_tension.line.json", doc);
    }
    {
      LineDoc doc;
      doc.catalog  = R"({
        "conductors": {
          "Drake": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 }
        },
        "towers": {
          "flat": { "attachments": { "left": { "x": 0.0, "h": 10.0 } } }
        }
      })";
      success     *= rejects("gridkit_line_uppercase_type.line.json", doc);
    }
    {
      LineDoc doc;
      doc.catalog  = R"({
        "conductors": {
          "wire": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 }
        },
        "towers": {}
      })";
      success     *= rejects("gridkit_line_empty_catalog_section.line.json", doc);
    }
    {
      LineDoc doc;
      doc.catalog  = R"({
        "conductors": {
          "wire": { "radius": { "inner": 0.005 }, "conductivity": 3.7e7, "weight": 11.0 }
        },
        "towers": {
          "flat": { "attachments": { "left": { "x": 0.0, "h": 10.0 } } }
        }
      })";
      success     *= rejects("gridkit_line_no_outer_radius.line.json", doc);
    }
    {
      LineDoc doc;
      doc.catalog  = R"({
        "conductors": {
          "wire": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7,
                    "permeability": 1.2566370614e-6, "weight": 11.0 }
        },
        "towers": {
          "flat": { "attachments": { "left": { "x": 0.0, "h": 10.0 } } }
        }
      })";
      success     *= rejects("gridkit_line_absolute_permeability.line.json", doc);
    }
    {
      LineDoc doc;
      doc.earth  = R"({ "conductivity": 1.0e-6, "permittivity": 8.8541878128e-12 })";
      success   *= rejects("gridkit_line_absolute_permittivity.line.json", doc);
    }
    {
      LineDoc doc;
      doc.catalog  = R"({
        "conductors": {
          "wire": { "radius": { "outer": 0.01, "inner": 0.02 }, "conductivity": 3.7e7,
                    "weight": 11.0 }
        },
        "towers": {
          "flat": { "attachments": { "left": { "x": 0.0, "h": 10.0 } } }
        }
      })";
      success     *= rejects("gridkit_line_inner_above_outer.line.json", doc);
    }
    {
      LineDoc doc;
      doc.line  = "";
      success  *= rejects("gridkit_line_missing_version.line.json", doc);
    }
    {
      LineDoc doc;
      doc.line  = R"("2.0")";
      success  *= rejects("gridkit_line_unsupported_version.line.json", doc);
    }
    {
      LineDoc doc;
      doc.extra  = R"("segment": { "span": 300.0 })";
      success   *= rejects("gridkit_line_unknown_field.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = "[]";
      success        *= rejects("gridkit_line_empty_conductors.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0,
                      "points": [{ "latitude": 91.0, "longitude": 0.0 },
                                 { "latitude": 0.0, "longitude": 0.0 }] })";
      success  *= rejects("gridkit_line_bad_latitude.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0,
                      "points": [{ "latitude": 0.0, "longitude": 181.0 },
                                 { "latitude": 0.0, "longitude": 0.0 }] })";
      success  *= rejects("gridkit_line_bad_longitude.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0, "points": [{ "latitude": 0.0, "longitude": 0.0 }] })";
      success  *= rejects("gridkit_line_short_points.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 300.0,
                      "points": [{ "latitude": 0.0, "longitude": 0.0 },
                                 { "latitude": 0.0, "longitude": 0.0 }] })";
      success  *= rejects("gridkit_line_zero_length_path.line.json", doc);
    }
    {
      LineDoc doc;
      doc.path  = R"({ "span": 0.0, "length": 100000.0 })";
      success  *= rejects("gridkit_line_zero_span.line.json", doc);
    }

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome line_parser_invalid_references()
  {
    using GridKit::Testing::TestStatus;

    TestStatus success = true;

    // Ports of the loader-level negative cases of
    // tests/Schema/validate_line_inputs.py.
    {
      LineDoc doc;
      doc.conductors  = R"([{ "at": "left", "phase": "a", "type": "nope" }])";
      success        *= rejects("gridkit_line_unknown_type.line.json", doc);
    }
    {
      LineDoc doc;
      doc.tower  = R"("nope")";
      success   *= rejects("gridkit_line_unknown_tower.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([{ "at": "nowhere", "phase": "a", "type": "wire" }])";
      success        *= rejects("gridkit_line_unknown_attachment.line.json", doc);
    }
    {
      LineDoc doc;
      doc.conductors  = R"([
        { "at": "left", "phase": "a", "type": "wire" },
        { "at": "left", "phase": "b", "type": "wire" }
      ])";
      success        *= rejects("gridkit_line_attachment_reuse.line.json", doc);
    }
    {
      LineDoc doc;
      doc.include  = R"(["missing.catalog.json"])";
      success     *= rejects("gridkit_line_missing_include.line.json", doc);
    }
    {
      // Tension whose catenary sag exceeds the attachment height.
      LineDoc doc;
      doc.conductors  = R"([{ "at": "left", "phase": "a", "type": "wire", "tension": 2000.0 }])";
      success        *= rejects("gridkit_line_sagged_underground.line.json", doc);
    }

    // Include-based negatives share one catalog file on disk.
    const auto catalog_path = testPath("gridkit_line_negative.catalog.json");
    writeFile(catalog_path, R"({
      "catalog": "1.0",
      "conductors": {
        "wire": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 }
      }
    })");
    {
      // Local catalog and include both define the type name.
      LineDoc doc;
      doc.include  = R"(["gridkit_line_negative.catalog.json"])";
      success     *= rejects("gridkit_line_name_collision.line.json", doc);
    }
    fs::remove(catalog_path);

    const auto headerless_path = testPath("gridkit_line_headerless.catalog.json");
    writeFile(headerless_path, R"({
      "conductors": {
        "spare": { "radius": { "outer": 0.01 }, "conductivity": 3.7e7, "weight": 11.0 }
      }
    })");
    {
      LineDoc doc;
      doc.include  = R"(["gridkit_line_headerless.catalog.json"])";
      success     *= rejects("gridkit_line_headerless_include.line.json", doc);
    }
    fs::remove(headerless_path);

    const auto stray_path = testPath("gridkit_line_stray.catalog.json");
    writeFile(stray_path, R"({ "catalog": "1.0", "earths": {} })");
    {
      LineDoc doc;
      doc.include  = R"(["gridkit_line_stray.catalog.json"])";
      success     *= rejects("gridkit_line_stray_include.line.json", doc);
    }
    fs::remove(stray_path);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome frequency_response_parser_valid()
  {
    using GridKit::EMT::Application::parseFrequencyResponseData;
    using GridKit::Testing::TestStatus;
    using MonitorVariable = GridKit::EMT::Application::FrequencyResponseData::MonitorVariable;

    const auto path = testPath("gridkit_frequency_response_valid.solver.json");
    writeFile(path, R"({
      "model": "overhead.line.json",
      "frequency": {
        "start": 0.15915494309189535,
        "stop": 159154.94309189535,
        "points": 21
      },
      "output": {
        "file": "overhead.response.csv",
        "variables": ["R", "Yc", "zc", "Tv", "ti", "Alpha", "beta", "Tau", "H"]
      }
    })");

    const auto data = parseFrequencyResponseData(path);
    fs::remove(path);

    const auto dir = fs::temp_directory_path();

    TestStatus success  = true;
    success            *= (data.model == dir / "overhead.line.json");
    success            *= close(data.frequency.start, 0.15915494309189535);
    success            *= close(data.frequency.stop, 159154.94309189535);
    success            *= (data.frequency.points == 21);
    success            *= (data.frequency.scale == "log");
    success            *= (data.ida.suppress_algebraic_error == true);
    success            *= (data.output_file == dir / "overhead.response.csv");
    success            *= data.variables.contains(MonitorVariable::R);
    success            *= data.variables.contains(MonitorVariable::Tv);
    success            *= data.variables.contains(MonitorVariable::Ti);
    success            *= data.variables.contains(MonitorVariable::Alpha);
    success            *= data.variables.contains(MonitorVariable::Beta);
    success            *= data.variables.contains(MonitorVariable::Tau);
    success            *= data.variables.contains(MonitorVariable::H);
    success            *= data.variables.contains(MonitorVariable::Yc);
    success            *= data.variables.contains(MonitorVariable::Zc);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome frequency_response_parser_ida_options()
  {
    using GridKit::EMT::Application::parseFrequencyResponseData;
    using GridKit::Testing::TestStatus;

    const auto path = testPath("gridkit_frequency_response_ida.solver.json");
    writeFile(path, R"({
      "model": "overhead.line.json",
      "ida": {
        "suppress_algebraic_error": false
      },
      "frequency": {
        "start": 1.0,
        "stop": 2.0,
        "points": 2
      },
      "output": {
        "file": "overhead.response.csv",
        "variables": ["R"]
      }
    })");

    const auto data = parseFrequencyResponseData(path);
    fs::remove(path);

    TestStatus success  = true;
    success            *= (data.ida.suppress_algebraic_error == false);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome frequency_response_parser_invalid_inputs()
  {
    using GridKit::EMT::Application::parseFrequencyResponseData;
    using GridKit::Testing::TestStatus;

    const auto bad_points = testPath("gridkit_frequency_response_bad_points.solver.json");
    writeFile(bad_points, R"({
      "model":"overhead.line.json",
      "frequency":{"start":1.0,"stop":2.0,"points":1},
      "output":{"file":"overhead.response.csv","variables":["R"]}
    })");

    const auto bad_scale = testPath("gridkit_frequency_response_bad_scale.solver.json");
    writeFile(bad_scale, R"({
      "model":"overhead.line.json",
      "frequency":{"start":1.0,"stop":2.0,"points":2,"scale":"linear"},
      "output":{"file":"overhead.response.csv","variables":["R"]}
    })");

    const auto bad_variable = testPath("gridkit_frequency_response_bad_variable.solver.json");
    writeFile(bad_variable, R"({
      "model":"overhead.line.json",
      "frequency":{"start":1.0,"stop":2.0,"points":2},
      "output":{"file":"overhead.response.csv","variables":["not_a_variable"]}
    })");

    const auto bad_derived_variable = testPath("gridkit_frequency_response_bad_derived_variable.solver.json");
    writeFile(bad_derived_variable, R"({
      "model":"overhead.line.json",
      "frequency":{"start":1.0,"stop":2.0,"points":2},
      "output":{"file":"overhead.response.csv","variables":["X","B"]}
    })");

    const auto bad_vocabulary = testPath("gridkit_frequency_response_bad_vocabulary.solver.json");
    writeFile(bad_vocabulary, R"({
      "model":"overhead.line.json",
      "frequency_hz":{"start":1.0,"stop":2.0,"points":2},
      "frequency":{"start":1.0,"stop":2.0,"points":2},
      "output":{"file":"overhead.response.csv","variables":["R"]}
    })");

    const auto bad_ida_type = testPath("gridkit_frequency_response_bad_ida_type.solver.json");
    writeFile(bad_ida_type, R"({
      "model":"overhead.line.json",
      "ida":false,
      "frequency":{"start":1.0,"stop":2.0,"points":2},
      "output":{"file":"overhead.response.csv","variables":["R"]}
    })");

    const auto bad_ida_value = testPath("gridkit_frequency_response_bad_ida_value.solver.json");
    writeFile(bad_ida_value, R"({
      "model":"overhead.line.json",
      "ida":{"suppress_algebraic_error":"false"},
      "frequency":{"start":1.0,"stop":2.0,"points":2},
      "output":{"file":"overhead.response.csv","variables":["R"]}
    })");

    TestStatus success  = true;
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_points); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_scale); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_variable); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_derived_variable); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_vocabulary); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_ida_type); });
    success            *= throws([&]
                      { parseFrequencyResponseData(bad_ida_value); });

    fs::remove(bad_points);
    fs::remove(bad_scale);
    fs::remove(bad_variable);
    fs::remove(bad_derived_variable);
    fs::remove(bad_vocabulary);
    fs::remove(bad_ida_type);
    fs::remove(bad_ida_value);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += line_parser_valid();
  result += line_parser_gis_path();
  result += line_parser_golden_wood_pole();
  result += line_parser_included_catalog();
  result += line_parser_invalid_documents();
  result += line_parser_invalid_references();
  result += frequency_response_parser_valid();
  result += frequency_response_parser_ida_options();
  result += frequency_response_parser_invalid_inputs();
  return result.summary();
}
