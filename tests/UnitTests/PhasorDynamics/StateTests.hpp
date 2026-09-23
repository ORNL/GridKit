#pragma once

#include <cmath>
#include <stdexcept>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief `applyState` and `extractState` on system model data
     */
    class StateTests
    {
      using DataT = PhasorDynamics::SystemModelData<double, size_t>;

    public:
      /// Applying the extracted state keeps the case operating point
      TestOutcome roundTrip()
      {
        TestStatus success = true;

        const DataT      data  = caseData();
        Model::StateData state = PhasorDynamics::extractState(data);

        DataT copy                                                          = caseData();
        copy.genrou[0].parameters[PhasorDynamics::GenrouParameters::p0]     = 0.0;
        copy.loadzip[0].parameters[PhasorDynamics::LoadZIPParameters::Pnom] = 0.0;
        PhasorDynamics::applyState(copy, state);

        success *= isEqual(copy.bus[1].Vr0, data.bus[1].Vr0) && isEqual(copy.bus[1].Vi0, data.bus[1].Vi0);
        success *= isEqual(real(copy.genrou[0].parameters.at(PhasorDynamics::GenrouParameters::p0)), P0);
        success *= isEqual(real(copy.genrou[0].parameters.at(PhasorDynamics::GenrouParameters::q0)), Q0);
        success *= isEqual(real(copy.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Pnom)), PNOM);
        success *= isEqual(real(copy.loadzip[0].parameters.at(PhasorDynamics::LoadZIPParameters::Qnom)), QNOM);
        success *= isEqual(real(copy.branch[0].parameters.at(PhasorDynamics::BranchParameters::tap)), TAP);

        // Power of the LoadZ at its bus voltage
        const Model::StateRecord* bus  = state.bus(2);
        const double              v2   = bus->value("vr", 0.0) * bus->value("vr", 0.0) + bus->value("vi", 0.0) * bus->value("vi", 0.0);
        double                    p    = 0.0;
        double                    q    = 0.0;
        success                       *= Model::terminalPower(state, "loadz_2", 2, 0, 1, p, q);
        success                       *= isEqual(p, -v2 * R / (R * R + X * X));
        success                       *= isEqual(q, -v2 * X / (R * R + X * X));

        return success.report(__func__);
      }

      /// Open branches and offline shunt loads are removed, offline loads draw nothing
      TestOutcome removals()
      {
        TestStatus success = true;

        DataT            data                      = caseData();
        Model::StateData state                     = PhasorDynamics::extractState(data);
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

      static double real(const std::variant<bool, double, size_t>& value)
      {
        return std::get<double>(value);
      }

      static void applyOfflineMachine()
      {
        DataT            data                     = caseData();
        Model::StateData state                    = PhasorDynamics::extractState(data);
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
