#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>
#include <GridKit/Model/EMT/Parameters/OverheadDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <application/EMT/FrequencyResponseJSONParser.hpp>

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

  std::string validLineJson()
  {
    return R"({
      "class": "Overhead",
      "length": 100000.0,
      "tower": {
        "x": [0.0, 30.0],
        "height": [10.0, 10.0],
        "span": 300.0,
        "tension": [200000.0, 200000.0]
      },
      "conductors": [
        {
          "radius": 0.01,
          "inner_radius": 0.0,
          "conductivity": 37000000.0,
          "permeability": 0.0000012566370614359173,
          "weight": 11.0,
          "phase": "a"
        },
        {
          "radius": 0.02,
          "inner_radius": 0.0,
          "conductivity": 58000000.0,
          "permeability": 0.0000012566370614359173,
          "weight": 12.0,
          "phase": "b"
        }
      ],
      "earth_conductivity": 0.000001,
      "earth_permittivity": 0.000000000022135469532
    })";
  }

  std::string gisLineJson()
  {
    return R"({
      "class": "Overhead",
      "tower": {
        "x": [0.0],
        "height": [10.0],
        "span": 300.0
      },
      "path": [
        {"latitude": 0.0, "longitude": 0.0},
        {"latitude": 0.0, "longitude": 1.0}
      ],
      "conductors": [
        {
          "radius": 0.01,
          "inner_radius": 0.0,
          "conductivity": 37000000.0,
          "permeability": 0.0000012566370614359173,
          "weight": 11.0
        }
      ],
      "earth_conductivity": 0.000001,
      "earth_permittivity": 0.000000000022135469532
    })";
  }

  GridKit::Testing::TestOutcome overhead_line_parser_valid()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    const auto path = testPath("gridkit_frequency_response_valid.line.json");
    writeFile(path, validLineJson());

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    TestStatus success  = true;
    success            *= (data.tower.K == 2);
    success            *= (data.conductor.K == 2);
    success            *= close(data.tower.position[0], 0.0);
    success            *= data.path.length.has_value();
    success            *= close(data.path.length.value(), 100000.0);
    success            *= data.path.path.empty();
    success            *= close(data.tower.span, 300.0);
    success            *= data.tower.tension.has_value();
    success            *= close(data.tower.tension.value()[0], 200000.0);
    success            *= close(data.tower.position[1], 30.0);
    success            *= close(data.tower.height[0], 10.0);
    success            *= close(data.conductor.radius[0], 0.01);
    success            *= close(data.conductor.radius[1], 0.02);
    success            *= close(data.conductor.inner_radius[0], 0.0);
    success            *= close(data.conductor.sigma[0], 3.7e7);
    success            *= close(data.conductor.mu[0], mu0);
    success            *= close(data.conductor.weight[0], 11.0);
    success            *= close(data.conductor.weight[1], 12.0);
    success            *= (data.conductor.phase[0] == "a");
    success            *= (data.conductor.phase[1] == "b");
    success            *= close(data.carson.earth_sigma, 1.0e-6);
    success            *= close(data.carson.earth_eps, 2.5 * eps0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_line_parser_gis_path()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    const auto path = testPath("gridkit_frequency_response_gis.line.json");
    writeFile(path, gisLineJson());

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    const scalar_type ref =
        GridKit::EMT::Parameters::Path<scalar_type, index_type>::earthRadius() * pi / 180.0
        * parsedTowerSaggedToSpanRatio(data);

    TestStatus success  = true;
    success            *= !data.path.length.has_value();
    success            *= (data.path.path.size() == 2);
    success            *= close(parsedPathLength(data), ref, 1.0e-12);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_line_parser_length_override()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    const auto path = testPath("gridkit_frequency_response_override.line.json");
    writeFile(path, R"({
      "class": "Overhead",
      "length": 100000.0,
      "tower": {
        "x": [0.0],
        "height": [10.0],
        "span": 300.0
      },
      "path": [
        {"latitude": 0.0, "longitude": 0.0},
        {"latitude": 0.0, "longitude": 1.0}
      ],
      "conductors": [{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto data = parseOverheadData<scalar_type, index_type>(path);
    fs::remove(path);

    TestStatus success  = true;
    success            *= data.path.length.has_value();
    success            *= (data.path.path.size() == 2);
    success            *= close(parsedPathLength(data), 100000.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_line_parser_invalid_inputs()
  {
    using GridKit::EMT::Parameters::parseOverheadData;
    using GridKit::Testing::TestStatus;

    const auto invalid_class = testPath("gridkit_frequency_response_invalid_class.line.json");
    writeFile(invalid_class, R"({"class":"UndergroundTransmission","conductors":[]})");

    const auto empty_conductors = testPath("gridkit_frequency_response_empty.line.json");
    writeFile(empty_conductors, R"({
      "class":"Overhead",
      "conductors":[],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_radius = testPath("gridkit_frequency_response_bad_radius.line.json");
    writeFile(bad_radius, R"({
      "class":"Overhead",
      "length":100000.0,
      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
      "conductors":[{
        "radius":0.0,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_length = testPath("gridkit_frequency_response_bad_length.line.json");
    writeFile(bad_length, R"({
      "class":"Overhead",
      "length":0.0,
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto missing_length_path = testPath("gridkit_frequency_response_missing_path.line.json");
    writeFile(missing_length_path, R"({
      "class":"Overhead",
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_path_type = testPath("gridkit_frequency_response_bad_path_type.line.json");
    writeFile(bad_path_type, R"({
      "class":"Overhead",
      "path":{"latitude":0.0,"longitude":0.0},
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto short_path = testPath("gridkit_frequency_response_short_path.line.json");
    writeFile(short_path, R"({
      "class":"Overhead",
      "path":[{"latitude":0.0,"longitude":0.0}],
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_latitude = testPath("gridkit_frequency_response_bad_latitude.line.json");
    writeFile(bad_latitude, R"({
      "class":"Overhead",
      "path":[{"latitude":91.0,"longitude":0.0},{"latitude":0.0,"longitude":0.0}],
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_longitude = testPath("gridkit_frequency_response_bad_longitude.line.json");
    writeFile(bad_longitude, R"({
      "class":"Overhead",
      "path":[{"latitude":0.0,"longitude":181.0},{"latitude":0.0,"longitude":0.0}],
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto zero_path = testPath("gridkit_frequency_response_zero_path.line.json");
    writeFile(zero_path, R"({
      "class":"Overhead",
      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
      "path":[{"latitude":0.0,"longitude":0.0},{"latitude":0.0,"longitude":0.0}],
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto stale_skin_effect = testPath("gridkit_frequency_response_stale_skin_effect.line.json");
    writeFile(stale_skin_effect, R"({
      "class":"Overhead",
      "length":100000.0,
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0,
      "skin_effect":{"root_count":8}
    })");

    const auto bad_span = testPath("gridkit_frequency_response_bad_span.line.json");
    writeFile(bad_span, R"({
      "class":"Overhead",
      "length":100000.0,
      "tower":{"x":[0.0],"height":[10.0],"span":0.0},
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto stale_segment = testPath("gridkit_frequency_response_stale_segment.line.json");
    writeFile(stale_segment, R"({
      "class":"Overhead",
      "length":100000.0,
      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
      "segment":{"span":300.0},
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto unsupported_tower_field =
        testPath("gridkit_frequency_response_unsupported_tower_field.line.json");
    writeFile(unsupported_tower_field, R"({
      "class":"Overhead",
      "length":100000.0,
      "tower":{
        "x":[0.0],
        "height":[10.0],
        "span":300.0,
        "unsupported":true
      },
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto missing_weight = testPath("gridkit_frequency_response_missing_weight.line.json");
    writeFile(missing_weight, R"({
      "class":"Overhead",
      "length":100000.0,
      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0
      }],
      "earth_conductivity":1.0,
      "earth_permittivity":1.0
    })");

    const auto bad_phase = testPath("gridkit_frequency_response_bad_phase.line.json");
    writeFile(bad_phase, R"({
	      "class":"Overhead",
	      "length":100000.0,
	      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
      "conductors":[{
        "radius":0.01,
        "inner_radius":0.0,
        "conductivity":1.0,
        "permeability":1.0,
        "weight":1.0,
        "phase":"x"
      }],
	      "earth_conductivity":1.0,
	      "earth_permittivity":1.0
	    })");

    const auto stale_earth = testPath("gridkit_frequency_response_stale_earth.line.json");
    writeFile(stale_earth, R"({
	      "class":"Overhead",
	      "length":100000.0,
	      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
	      "conductors":[{
	        "radius":0.01,
	        "inner_radius":0.0,
	        "conductivity":1.0,
	        "permeability":1.0,
	        "weight":1.0
	      }],
	      "earth":{"conductivity":1.0,"permittivity":1.0}
	    })");

    const auto bad_earth_permittivity =
        testPath("gridkit_frequency_response_bad_earth_permittivity.line.json");
    writeFile(bad_earth_permittivity, R"({
	      "class":"Overhead",
	      "length":100000.0,
	      "tower":{"x":[0.0],"height":[10.0],"span":300.0},
	      "conductors":[{
	        "radius":0.01,
	        "inner_radius":0.0,
	        "conductivity":1.0,
	        "permeability":1.0,
	        "weight":1.0
	      }],
	      "earth_conductivity":1.0,
	      "earth_permittivity":0.0
	    })");

    TestStatus success  = true;
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(invalid_class); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(empty_conductors); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_radius); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_length); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(missing_length_path); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_path_type); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(short_path); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_latitude); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_longitude); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(zero_path); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(stale_skin_effect); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_span); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(stale_segment); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(unsupported_tower_field); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(missing_weight); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_phase); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(stale_earth); });
    success            *= throws([&]
                      { parseOverheadData<scalar_type, index_type>(bad_earth_permittivity); });

    fs::remove(invalid_class);
    fs::remove(empty_conductors);
    fs::remove(bad_radius);
    fs::remove(bad_length);
    fs::remove(missing_length_path);
    fs::remove(bad_path_type);
    fs::remove(short_path);
    fs::remove(bad_latitude);
    fs::remove(bad_longitude);
    fs::remove(zero_path);
    fs::remove(stale_skin_effect);
    fs::remove(bad_span);
    fs::remove(stale_segment);
    fs::remove(unsupported_tower_field);
    fs::remove(missing_weight);
    fs::remove(bad_phase);
    fs::remove(stale_earth);
    fs::remove(bad_earth_permittivity);

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
  result += overhead_line_parser_valid();
  result += overhead_line_parser_gis_path();
  result += overhead_line_parser_length_override();
  result += overhead_line_parser_invalid_inputs();
  result += frequency_response_parser_valid();
  result += frequency_response_parser_ida_options();
  result += frequency_response_parser_invalid_inputs();
  return result.summary();
}
