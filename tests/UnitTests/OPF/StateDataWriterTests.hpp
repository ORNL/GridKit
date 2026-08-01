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
    class StateDataWriterTests
    {
    public:
      TestOutcome writeCompleteState()
      {
        Model::StateData state;
        state.header.emplace();
        state.header->version     = 3U;
        state.header->time        = 12.5;
        state.header->created     = "2026-07-31T12:00:00-05:00";
        state.header->description = "complete state";

        auto& bus       = state.buses["z_bus"];
        bus.vr          = 0.9;
        bus.vi          = -0.1;
        auto& injection = bus.injections["z_injection"];
        injection.ir    = 1.25;
        injection.ii    = -0.75;

        state.buses["a_bus"];
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

        std::ostringstream first;
        std::ostringstream second;
        Model::writeStateData(state, first);
        Model::writeStateData(state, second);

        const std::string serialized  = first.str();
        TestStatus        success     = true;
        success                      *= serialized == second.str();
        success                      *= serialized.find("null") == std::string::npos;
        success                      *= serialized.find("\"a_bus\"")
                   < serialized.find("\"z_bus\"");
        success *= serialized.find("\"a_device\"")
                   < serialized.find("\"z_device\"");

        std::istringstream input(serialized);
        const auto         round_trip  = Model::parseStateData(input);
        success                       *= round_trip.header->version == state.header->version;
        success                       *= round_trip.header->time == state.header->time;
        success                       *= round_trip.header->created == state.header->created;
        success                       *= round_trip.header->description == state.header->description;
        success                       *= round_trip.buses.at("z_bus").vr == bus.vr;
        success                       *= round_trip.buses.at("z_bus").vi == bus.vi;
        success                       *= round_trip.buses.at("z_bus").injections.at("z_injection").ir
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
        success *= !round_trip.buses.at("a_bus").vr.has_value();
        success *= !round_trip.buses.at("a_bus").vi.has_value();
        success *= !round_trip.buses.at("a_bus")
                        .injections.at("a_injection")
                        .ir.has_value();
        success *= !round_trip.buses.at("a_bus")
                        .injections.at("a_injection")
                        .ii.has_value();
        success *= !round_trip.devices.at("a_device").active.has_value();
        success *= !round_trip.devices.at("a_device").online.has_value();
        success *= !round_trip.devices.at("a_device").open.has_value();
        success *= !round_trip.devices.at("a_device").p.has_value();
        success *= !round_trip.devices.at("a_device").q.has_value();
        success *= !round_trip.devices.at("a_device").tap.has_value();
        success *= !round_trip.devices.at("a_device").phase.has_value();

        Model::StateData omitted_header_fields;
        omitted_header_fields.header.emplace();
        std::ostringstream omitted_output;
        Model::writeStateData(omitted_header_fields, omitted_output);
        success *= omitted_output.str().find("null") == std::string::npos;
        std::istringstream omitted_input(omitted_output.str());
        const auto         omitted_round_trip  = Model::parseStateData(omitted_input);
        success                               *= omitted_round_trip.header.has_value();
        success                               *= !omitted_round_trip.header->version.has_value();
        success                               *= !omitted_round_trip.header->time.has_value();
        success                               *= !omitted_round_trip.header->created.has_value();
        success                               *= !omitted_round_trip.header->description.has_value();

        return success.report(__func__);
      }

      TestOutcome writeFileAndRejectFailures()
      {
        Model::StateData state;
        state.buses["bus_id_1"].vr = 1.0;

        const auto nonce = std::chrono::steady_clock::now()
                               .time_since_epoch()
                               .count();
        const auto file_path = std::filesystem::temp_directory_path()
                               / ("gridkit_opf_state_" + std::to_string(nonce)
                                  + ".json");

        TestStatus success = true;
        try
        {
          Model::writeStateData(state, file_path);
          const auto parsed  = Model::parseStateData(file_path);
          success           *= !parsed.header.has_value();
          success           *= parsed.buses.at("bus_id_1").vr == 1.0;
        }
        catch (const std::exception&)
        {
          success = false;
        }

        Model::StateData invalid_replacement;
        invalid_replacement.header.emplace();
        invalid_replacement.header->time  = std::numeric_limits<double>::infinity();
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

        Model::StateData nonfinite_header;
        nonfinite_header.header.emplace();
        nonfinite_header.header->time  = std::numeric_limits<double>::infinity();
        success                       *= rejectsNonfinite(nonfinite_header);

        Model::StateData nonfinite_bus;
        nonfinite_bus.buses["bus"].vr  = std::numeric_limits<double>::quiet_NaN();
        success                       *= rejectsNonfinite(nonfinite_bus);

        Model::StateData nonfinite_injection;
        nonfinite_injection.buses["bus"].injections["device"].ii  = std::numeric_limits<double>::infinity();
        success                                                  *= rejectsNonfinite(nonfinite_injection);

        Model::StateData nonfinite_device;
        nonfinite_device.devices["device"].p  = std::numeric_limits<double>::quiet_NaN();
        success                              *= rejectsNonfinite(nonfinite_device);

        std::error_code error;
        success *= std::filesystem::remove(file_path, error);
        success *= !error;

        return success.report(__func__);
      }

    private:
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
