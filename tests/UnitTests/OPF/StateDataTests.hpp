#pragma once

#include <chrono>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class StateDataTests
    {
    public:
      TestOutcome parseCompleteState()
      {
        const auto state = parse(R"({
          "header": {
            "version": 3,
            "time": 12.5,
            "created": "2026-07-31T12:00:00-05:00",
            "description": "complete state",
            "ignored": true
          },
          "buses": {
            "bus_id_1": {
              "vr": 0.9,
              "vi": -0.1,
              "injections": {
                "generator": {"ir": 1.25, "ii": -0.75, "ignored": 1},
                "null_injection": null
              },
              "ignored": "value"
            },
            "null_bus": null
          },
          "devices": {
            "branch": {
              "active": true,
              "online": false,
              "open": true,
              "p": 2.0,
              "q": -3.0,
              "tap": 1.05,
              "phase": 0.2,
              "ignored": []
            },
            "null_device": null
          },
          "ignored": {"field": true}
        })");

        TestStatus success  = true;
        success            *= state.header.has_value();
        success            *= state.header->version == 3U;
        success            *= state.header->time == 12.5;
        success            *= state.header->created == "2026-07-31T12:00:00-05:00";
        success            *= state.header->description == "complete state";
        success            *= state.buses.size() == 1;
        success            *= state.buses.at("bus_id_1").vr == 0.9;
        success            *= state.buses.at("bus_id_1").vi == -0.1;
        success            *= state.buses.at("bus_id_1").injections.size() == 1;
        success            *= state.buses.at("bus_id_1")
                       .injections.at("generator")
                       .ir
                   == 1.25;
        success *= state.buses.at("bus_id_1")
                       .injections.at("generator")
                       .ii
                   == -0.75;
        success *= state.devices.size() == 1;
        success *= state.devices.at("branch").active == true;
        success *= state.devices.at("branch").online == false;
        success *= state.devices.at("branch").open == true;
        success *= state.devices.at("branch").p == 2.0;
        success *= state.devices.at("branch").q == -3.0;
        success *= state.devices.at("branch").tap == 1.05;
        success *= state.devices.at("branch").phase == 0.2;

        return success.report(__func__);
      }

      TestOutcome parseOmittedAndNullValues()
      {
        const auto omitted       = parse("{}");
        const auto null_sections = parse(
            R"({"header": null, "buses": null, "devices": null})");
        const auto null_fields = parse(R"({
          "header": {
            "version": null,
            "time": null,
            "created": null,
            "description": null
          },
          "buses": {
            "bus": {"vr": null, "vi": null, "injections": null}
          },
          "devices": {
            "device": {
              "active": null,
              "online": null,
              "open": null,
              "p": null,
              "q": null,
              "tap": null,
              "phase": null
            }
          }
        })");

        TestStatus success  = true;
        success            *= !omitted.header.has_value();
        success            *= omitted.buses.empty();
        success            *= omitted.devices.empty();
        success            *= !null_sections.header.has_value();
        success            *= null_sections.buses.empty();
        success            *= null_sections.devices.empty();
        success            *= null_fields.header.has_value();
        success            *= !null_fields.header->version.has_value();
        success            *= !null_fields.header->time.has_value();
        success            *= !null_fields.header->created.has_value();
        success            *= !null_fields.header->description.has_value();
        success            *= !null_fields.buses.at("bus").vr.has_value();
        success            *= !null_fields.buses.at("bus").vi.has_value();
        success            *= null_fields.buses.at("bus").injections.empty();
        success            *= !null_fields.devices.at("device").active.has_value();
        success            *= !null_fields.devices.at("device").online.has_value();
        success            *= !null_fields.devices.at("device").open.has_value();
        success            *= !null_fields.devices.at("device").p.has_value();
        success            *= !null_fields.devices.at("device").q.has_value();
        success            *= !null_fields.devices.at("device").tap.has_value();
        success            *= !null_fields.devices.at("device").phase.has_value();

        return success.report(__func__);
      }

      TestOutcome parseLegacyGoldenState()
      {
        const auto state = parse(R"({
    "header": {
        "version": 1,
        "time": null,
        "created": "2026-07-31T05:55:01-05:00",
        "description": "Steady state solution at nominal dispatch"
    },
    "buses": {
        "bus_id_0": {
            "vr": 0.9949877346411762,
            "vi": 0.09999703952427966,
            "injections": {
                "DV1": {
                    "ir": 1.0000005862325285,
                    "ii": 0.0501183043867175
                }
            }
        },
        "bus_id_1": {
            "vr": 1.0,
            "vi": 0.0,
            "injections": {}
        }
    },
    "devices": {
        "BR1": {
            "open": false,
            "tap": 1.0,
            "phase": 0.0
        },
        "DV1": {
            "online": true,
            "p": 1.0,
            "q": 0.05013
        }
    }
}
)");

        TestStatus success  = true;
        success            *= state.header.has_value();
        success            *= state.header->version == 1U;
        success            *= !state.header->time.has_value();
        success            *= state.buses.at("bus_id_0").vr == 0.9949877346411762;
        success            *= state.buses.at("bus_id_0").vi == 0.09999703952427966;
        success            *= state.buses.at("bus_id_0").injections.at("DV1").ir
                   == 1.0000005862325285;
        success *= state.buses.at("bus_id_0").injections.at("DV1").ii
                   == 0.0501183043867175;
        success *= state.devices.at("BR1").open == false;
        success *= state.devices.at("DV1").online == true;
        success *= state.devices.at("DV1").p == 1.0;
        success *= state.devices.at("DV1").q == 0.05013;
        return success.report(__func__);
      }

      TestOutcome rejectInvalidInput()
      {
        TestStatus success  = true;
        success            *= rejects("[]");
        success            *= rejects(R"({"header": []})");
        success            *= rejects(R"({"buses": []})");
        success            *= rejects(R"({"buses": {"bus": []}})");
        success            *= rejects(
            R"({"buses": {"bus": {"injections": []}}})");
        success *= rejects(
            R"({"buses": {"bus": {"injections": {"device": []}}}})");
        success *= rejects(R"({"devices": []})");
        success *= rejects(R"({"devices": {"device": []}})");
        success *= rejects(R"({"header": {"version": "three"}})");
        success *= rejects(R"({"header": {"time": "now"}})");
        success *= rejects(R"({"buses": {"bus": {"vr": "one"}}})");
        success *= rejects(
            R"({"buses": {"bus": {"injections": {"device": {"ir": "one"}}}}})");
        success *= rejects(R"({"devices": {"device": {"online": 1}}})");
        success *= rejects(R"({"devices": {"device": {"p": "one"}}})");
        success *= rejects("{");

        return success.report(__func__);
      }

      TestOutcome writeDeterministicState()
      {
        Model::StateData state;
        state.header.emplace();
        state.buses["bus_id_1"];
        state.devices["load"];

        std::ostringstream output;
        Model::writeStateData(state, output);

        const std::string expected = R"({
    "buses": {
        "bus_id_1": {
            "injections": {}
        }
    },
    "devices": {
        "load": {}
    },
    "header": {}
}
)";

        std::ostringstream second_output;
        Model::writeStateData(state, second_output);

        std::istringstream round_trip_input(output.str());
        const auto         round_trip = Model::parseStateData(round_trip_input);

        TestStatus success  = true;
        success            *= output.str() == expected;
        success            *= output.str() == second_output.str();
        success            *= output.str().find("null") == std::string::npos;
        success            *= round_trip.header.has_value();
        success            *= round_trip.buses.size() == 1;
        success            *= round_trip.buses.at("bus_id_1").injections.empty();
        success            *= round_trip.devices.size() == 1;

        Model::StateData   no_header;
        std::ostringstream no_header_output;
        Model::writeStateData(no_header, no_header_output);
        success *= no_header_output.str()
                   == R"({
    "buses": {},
    "devices": {}
}
)";

        return success.report(__func__);
      }

      TestOutcome writeCompleteState()
      {
        Model::StateData state;
        state.header.emplace();
        state.header->version     = 3U;
        state.header->time        = 12.5;
        state.header->created     = "2026-07-31T12:00:00-05:00";
        state.header->description = "complete state";

        auto& first_bus = state.buses["z_bus"];
        first_bus.vr    = 0.9;
        first_bus.vi    = -0.1;
        auto& injection = first_bus.injections["z_injection"];
        injection.ir    = 1.25;
        injection.ii    = -0.75;
        state.buses["a_bus"].injections["a_injection"];

        auto& device  = state.devices["z_device"];
        device.active = true;
        device.online = false;
        device.open   = true;
        device.p      = 2.0;
        device.q      = -3.0;
        device.tap    = 1.05;
        device.phase  = 0.2;
        state.devices["a_device"];

        std::ostringstream output;
        Model::writeStateData(state, output);
        std::istringstream input(output.str());
        const auto         round_trip = Model::parseStateData(input);

        TestStatus success  = true;
        success            *= output.str().find("null") == std::string::npos;
        success            *= output.str().find("\"a_bus\"")
                   < output.str().find("\"z_bus\"");
        success *= output.str().find("\"a_device\"")
                   < output.str().find("\"z_device\"");
        success *= round_trip.header->version == state.header->version;
        success *= round_trip.header->time == state.header->time;
        success *= round_trip.header->created == state.header->created;
        success *= round_trip.header->description == state.header->description;
        success *= round_trip.buses.at("z_bus").vr == first_bus.vr;
        success *= round_trip.buses.at("z_bus").vi == first_bus.vi;
        success *= round_trip.buses.at("z_bus").injections.at("z_injection").ir
                   == injection.ir;
        success *= round_trip.buses.at("z_bus").injections.at("z_injection").ii
                   == injection.ii;
        success *= round_trip.devices.at("z_device").active == device.active;
        success *= round_trip.devices.at("z_device").online == device.online;
        success *= round_trip.devices.at("z_device").open == device.open;
        success *= round_trip.devices.at("z_device").p == device.p;
        success *= round_trip.devices.at("z_device").q == device.q;
        success *= round_trip.devices.at("z_device").tap == device.tap;
        success *= round_trip.devices.at("z_device").phase == device.phase;

        return success.report(__func__);
      }

      TestOutcome fileAndOutputFailures()
      {
        Model::StateData state;
        state.buses["bus_id_1"].vr = 1.0;

        const auto nonce = std::chrono::steady_clock::now()
                               .time_since_epoch()
                               .count();
        const auto file_path = std::filesystem::temp_directory_path()
                               / ("gridkit_state_data_" + std::to_string(nonce)
                                  + ".json");

        TestStatus success  = true;
        success            *= rejectsFile(file_path);

        try
        {
          Model::writeStateData(state, file_path);
          const auto parsed  = Model::parseStateData(file_path);
          success           *= parsed.buses.at("bus_id_1").vr == 1.0;
        }
        catch (const std::exception&)
        {
          success = false;
        }

        Model::StateData invalid_replacement;
        invalid_replacement.header.emplace();
        invalid_replacement.header->time =
            std::numeric_limits<double>::infinity();
        bool rejected_invalid_replacement = false;
        try
        {
          Model::writeStateData(invalid_replacement, file_path);
        }
        catch (const std::invalid_argument&)
        {
          rejected_invalid_replacement = true;
        }
        success *= rejected_invalid_replacement;

        try
        {
          const auto preserved  = Model::parseStateData(file_path);
          success              *= preserved.buses.at("bus_id_1").vr == 1.0;
        }
        catch (const std::exception&)
        {
          success = false;
        }

        bool rejected_bad_path = false;
        try
        {
          Model::writeStateData(state, file_path / "state.json");
        }
        catch (const std::runtime_error&)
        {
          rejected_bad_path = true;
        }
        success *= rejected_bad_path;

        std::ostringstream failed_stream;
        failed_stream.setstate(std::ios::badbit);
        bool rejected_bad_stream = false;
        try
        {
          Model::writeStateData(state, failed_stream);
        }
        catch (const std::runtime_error&)
        {
          rejected_bad_stream = true;
        }
        success *= rejected_bad_stream;

        std::error_code error;
        success *= std::filesystem::remove(file_path, error);
        success *= !error;

        return success.report(__func__);
      }

      TestOutcome rejectNonfiniteOutput()
      {
        Model::StateData header;
        header.header.emplace();
        header.header->time = std::numeric_limits<double>::infinity();

        Model::StateData bus;
        bus.buses["bus"].vi = std::numeric_limits<double>::quiet_NaN();

        Model::StateData injection;
        injection.buses["bus"].injections["device"].ir =
            -std::numeric_limits<double>::infinity();

        Model::StateData device;
        device.devices["device"].phase =
            std::numeric_limits<double>::quiet_NaN();

        TestStatus success  = true;
        success            *= rejectsNonfinite(header);
        success            *= rejectsNonfinite(bus);
        success            *= rejectsNonfinite(injection);
        success            *= rejectsNonfinite(device);

        return success.report(__func__);
      }

    private:
      static Model::StateData parse(const std::string& input)
      {
        std::istringstream stream(input);
        return Model::parseStateData(stream);
      }

      static bool rejects(const std::string& input)
      {
        try
        {
          parse(input);
        }
        catch (const std::exception&)
        {
          return true;
        }
        return false;
      }

      static bool rejectsFile(const std::filesystem::path& file_path)
      {
        try
        {
          Model::parseStateData(file_path);
        }
        catch (const std::runtime_error&)
        {
          return true;
        }
        return false;
      }

      static bool rejectsNonfinite(const Model::StateData& state)
      {
        std::ostringstream output;
        try
        {
          Model::writeStateData(state, output);
        }
        catch (const std::invalid_argument&)
        {
          return output.str().empty();
        }
        return false;
      }
    };
  } // namespace Testing
} // namespace GridKit
