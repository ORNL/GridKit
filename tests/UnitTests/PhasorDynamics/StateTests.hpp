#pragma once

#include <stdexcept>
#include <string>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "StateUtilities.hpp"

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Application state initialization of PhasorDynamics case data
     */
    class StateTests
    {
      using DataT = PhasorDynamics::SystemModelData<double, size_t>;

    public:
      /// Supplied voltage, current, and tap values set the case operating point
      TestOutcome applyState()
      {
        TestStatus success = true;

        DataT            data = caseData();
        Model::StateData state;
        state.buses[Model::busKey(1)].values = {{"vr", 1.0}, {"vi", 0.0}};
        state.buses[Model::busKey(2)].values = {{"vr", 1.0}, {"vi", -0.1}};
        state.devices["genrou_1"].values     = {{"ir", 0.8}, {"ii", -0.2}};
        state.devices["loadzip_2"].values    = {{"ir", -0.7}, {"ii", 0.1}};
        state.devices["branch_1_2"].values   = {{"tap", 1.05}, {"phase", 0.02}};
        PhasorDynamics::applyState(data, state);

        success *= isEqual(data.bus[1].Vr0, 1.0) && isEqual(data.bus[1].Vi0, -0.1);
        success *= isEqual(real(data.genrou[0].parameters.at(PhasorDynamics::GenrouParameters::p0)), 0.8);
        success *= isEqual(real(data.genrou[0].parameters.at(PhasorDynamics::GenrouParameters::q0)), 0.2);
        success *= isEqual(real(data.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Pnom)), 0.71);
        success *= isEqual(real(data.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Qnom)), 0.03);
        success *= isEqual(real(data.branch[0].parameters.at(PhasorDynamics::BranchParameters::tap)), 1.05);
        success *= isEqual(real(data.branch[0].parameters.at(PhasorDynamics::BranchParameters::phase)), 0.02);

        return success.report(__func__);
      }

      /// Open branches and offline shunt loads are removed, offline loads draw nothing
      TestOutcome removals()
      {
        TestStatus success = true;

        DataT            data = caseData();
        Model::StateData state;
        state.devices["branch_1_2"].flags["open"]  = true;
        state.devices["loadz_2"].flags["online"]   = false;
        state.devices["loadzip_2"].flags["online"] = false;
        PhasorDynamics::applyState(data, state);

        success *= data.branch.empty();
        success *= data.loadz.empty();
        success *= isEqual(real(data.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Pnom)), 0.0);
        success *= isEqual(real(data.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Qnom)), 0.0);

        return success.report(__func__);
      }

      /// Machines cannot be offline yet
      TestOutcome offlineMachine()
      {
        TestStatus success = true;

        success *= throws<std::invalid_argument>(applyOfflineMachine);

        return success.report(__func__);
      }

    private:
      static constexpr double P0   = 0.8;
      static constexpr double Q0   = 0.2;
      static constexpr double PNOM = 0.7;
      static constexpr double QNOM = 0.1;
      static constexpr double R    = 1.0;
      static constexpr double X    = 0.5;
      static constexpr double TAP  = 1.02;

      static double real(const std::variant<std::string, bool, double, size_t>& value)
      {
        return std::get<double>(value);
      }

      static void applyOfflineMachine()
      {
        DataT            data = caseData();
        Model::StateData state;
        state.devices["genrou_1"].flags["online"] = false;
        PhasorDynamics::applyState(data, state);
      }

      static DataT caseData()
      {
        using namespace PhasorDynamics;

        DataT data;
        data.bus.resize(2);
        data.bus[0].bus_id = 1;
        data.bus[0].Vr0    = 1.02;
        data.bus[0].Vi0    = 0.0;
        data.bus[1].bus_id = 2;
        data.bus[1].Vr0    = 0.98;
        data.bus[1].Vi0    = -0.09;

        auto& branch                             = data.branch.emplace_back();
        branch.disambiguation_string             = "branch_1_2";
        branch.buses[BranchBuses::bus1]          = 1;
        branch.buses[BranchBuses::bus2]          = 2;
        branch.parameters[BranchParameters::X]   = 0.1;
        branch.parameters[BranchParameters::tap] = TAP;

        auto& genrou                            = data.genrou.emplace_back();
        genrou.disambiguation_string            = "genrou_1";
        genrou.buses[GenrouBuses::bus]          = 1;
        genrou.parameters[GenrouParameters::p0] = P0;
        genrou.parameters[GenrouParameters::q0] = Q0;

        auto& loadzip                               = data.loadzip.emplace_back();
        loadzip.disambiguation_string               = "loadzip_2";
        loadzip.buses[LoadZIPBuses::bus]            = 2;
        loadzip.parameters[LoadZIPParameters::Pnom] = PNOM;
        loadzip.parameters[LoadZIPParameters::Qnom] = QNOM;

        auto& loadz                          = data.loadz.emplace_back();
        loadz.disambiguation_string          = "loadz_2";
        loadz.buses[LoadZBuses::bus]         = 2;
        loadz.parameters[LoadZParameters::R] = R;
        loadz.parameters[LoadZParameters::X] = X;

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
