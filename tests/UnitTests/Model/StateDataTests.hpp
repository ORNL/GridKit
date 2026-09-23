#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class StateDataTests
    {
    public:
      /// Bus and current keys follow `STATE.md`
      TestOutcome keys()
      {
        TestStatus success = true;

        success *= Model::busKey(7) == "bus_id_7";
        success *= Model::currentKeys(0, 1) == std::pair<std::string, std::string>{"ir", "ii"};
        success *= Model::currentKeys(1, 2) == std::pair<std::string, std::string>{"ir2", "ii2"};

        return success.report(__func__);
      }

      /// Numbers and flags are read into separate maps and survive writing
      TestOutcome roundTrip()
      {
        TestStatus success = true;

        std::istringstream     input(R"({
          "header": {"version": 1, "description": "round trip"},
          "buses": {"bus_id_2": {"vr": 0.98, "vi": -0.12}},
          "devices": {
            "gen_2": {"online": true, "ir": 0.8, "ii": -0.2},
            "branch_1_2": {"open": false, "tap": 1, "phase": 0.05, "note": null}
          }
        })");
        const Model::StateData state = Model::parseStateData(input);

        success *= state.header.version == 1U;
        success *= state.header.description == "round trip";
        success *= !state.header.created.has_value();
        success *= !state.header.time.has_value();

        const Model::StateRecord* bus  = state.bus(2);
        success                       *= bus != nullptr;
        success                       *= state.bus(3) == nullptr;

        const Model::StateRecord* generator  = state.device("gen_2");
        success                             *= generator != nullptr;
        success                             *= generator->flag("online", false);
        success                             *= !generator->values.contains("online");

        const Model::StateRecord* branch  = state.device("branch_1_2");
        success                          *= branch != nullptr;
        success                          *= !branch->flag("open", true);
        success                          *= isEqual(branch->value("tap", 0.0), 1.0);
        success                          *= !branch->values.contains("note");
        success                          *= isEqual(branch->value("missing", 0.25), 0.25);

        std::stringstream buffer;
        Model::writeStateData(state, buffer);
        const Model::StateData copy = Model::parseStateData(buffer);

        success *= copy.header.version == state.header.version;
        success *= copy.header.description == state.header.description;
        success *= !copy.header.created.has_value();
        success *= isEqual(copy.buses.at("bus_id_2").values, bus->values);
        success *= isEqual(copy.devices.at("gen_2").values, generator->values);
        success *= copy.devices.at("gen_2").flags == generator->flags;
        success *= isEqual(copy.devices.at("branch_1_2").values, branch->values);
        success *= copy.devices.at("branch_1_2").flags == branch->flags;

        return success.report(__func__);
      }

      /// Entries that are neither numbers nor booleans are rejected
      TestOutcome rejectsText()
      {
        TestStatus success = true;

        success *= throws<std::invalid_argument>(parseText);

        return success.report(__func__);
      }

      /// Terminal current and terminal power are inverse conversions
      TestOutcome terminalConversion()
      {
        TestStatus success = true;

        Model::StateData state;
        state.buses[Model::busKey(4)].values = {{"vr", 1.02}, {"vi", -0.15}};
        state.buses[Model::busKey(5)].values = {{"vr", 0.97}, {"vi", 0.08}};

        double p = 0.0;
        double q = 0.0;

        Model::setTerminalCurrent(state, "gen_4", 4, 0, 1, 0.7, -0.25);
        success *= Model::terminalPower(state, "gen_4", 4, 0, 1, p, q);
        success *= isEqual(p, 0.7);
        success *= isEqual(q, -0.25);

        Model::setTerminalCurrent(state, "branch_4_5", 5, 1, 2, -0.4, 0.1);
        success *= state.devices.at("branch_4_5").values.contains("ir2");
        success *= state.devices.at("branch_4_5").values.contains("ii2");
        success *= Model::terminalPower(state, "branch_4_5", 5, 1, 2, p, q);
        success *= isEqual(p, -0.4);
        success *= isEqual(q, 0.1);

        // No current for this terminal, and no voltage for bus 6
        success *= !Model::terminalPower(state, "branch_4_5", 4, 0, 2, p, q);
        success *= !Model::terminalPower(state, "gen_4", 6, 0, 1, p, q);

        return success.report(__func__);
      }

    private:
      static void parseText()
      {
        std::istringstream input(R"({"devices": {"gen_2": {"ir": "0.8"}}})");
        Model::parseStateData(input);
      }
    };
  } // namespace Testing
} // namespace GridKit
