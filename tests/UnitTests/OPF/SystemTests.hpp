#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <GridKit/Model/OPF/System.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class OPFSystemTests
    {
    public:
      using IdxT        = std::size_t;
      using RealT       = double;
      using SystemDataT = OPF::SystemData<RealT, IdxT>;
      using SystemT     = OPF::System<RealT, IdxT>;

      TestOutcome rejectsInvalidSystems()
      {
        TestStatus success = true;

        {
          auto [data, state]   = validProblem();
          data.params.va_base  = 0.0;
          success             *= rejects(data, state);
        }
        {
          auto [data, state]     = validProblem();
          data.params.freq_base  = std::numeric_limits<RealT>::infinity();
          success               *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          data.buses.clear();
          success *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          data.buses[0].id.clear();
          success *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          data.buses[1].id    = data.buses[0].id;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]    = validProblem();
          data.buses[1].number  = data.buses[0].number;
          success              *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          data.buses[0].kv    = 0.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          data.buses[0].vmin  = -0.01;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          data.buses[0].vmax  = 0.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          data.buses[0].vmin  = 1.1;
          data.buses[0].vmax  = 1.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          state.buses.erase("bus_id_1");
          success *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          state.buses["bus_id_0"].vi.reset();
          success *= rejects(data, state);
        }
        {
          auto [data, state]          = validProblem();
          state.buses["bus_id_0"].vr  = 0.0;
          state.buses["bus_id_0"].vi  = 0.0;
          success                    *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& branch        = std::get<SystemDataT::BranchDataT>(data.devices[0]);
          branch.to           = branch.from;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& branch        = std::get<SystemDataT::BranchDataT>(data.devices[0]);
          branch.R            = 0.0;
          branch.X            = 0.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& branch        = std::get<SystemDataT::BranchDataT>(data.devices[0]);
          branch.R            = std::numeric_limits<RealT>::denorm_min();
          branch.X            = 0.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& branch        = std::get<SystemDataT::BranchDataT>(data.devices[0]);
          branch.smax         = -1.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]       = validProblem();
          state.devices["BR"].tap  = 0.0;
          success                 *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          state.devices["BR"].phase =
              std::numeric_limits<RealT>::infinity();
          success *= rejects(data, state);
        }
        {
          auto [data, state]        = validProblem();
          state.devices["BR"].open  = true;
          success                  *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& generator     = std::get<SystemDataT::GeneratorDataT>(data.devices[1]);
          generator.pmin      = 2.0;
          generator.pmax      = 1.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& generator     = std::get<SystemDataT::GeneratorDataT>(data.devices[1]);
          generator.mva       = 0.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& generator     = std::get<SystemDataT::GeneratorDataT>(data.devices[1]);
          generator.c2        = std::numeric_limits<RealT>::quiet_NaN();
          success            *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          state.devices.erase("LD");
          success *= rejects(data, state);
        }
        {
          auto [data, state]     = validProblem();
          state.devices["LD"].p  = -2.0;
          success               *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          state.devices["LD"].q.reset();
          success *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& load          = std::get<SystemDataT::LoadDataT>(data.devices[2]);
          load.qmin           = -0.1;
          load.qmax           = -0.2;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& shunt         = std::get<SystemDataT::ShuntDataT>(data.devices[3]);
          shunt.B             = std::numeric_limits<RealT>::infinity();
          success            *= rejects(data, state);
        }
        {
          auto [data, state]                                    = validProblem();
          std::get<SystemDataT::LoadDataT>(data.devices[2]).id  = "BR";
          success                                              *= rejects(data, state);
        }
        {
          auto [data, state]                                      = validProblem();
          std::get<SystemDataT::ShuntDataT>(data.devices[3]).bus  = 9;
          success                                                *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          auto isolated      = data.buses.back();
          isolated.id        = "B2";
          isolated.number    = 2;
          data.buses.push_back(isolated);
          state.buses["bus_id_2"].vr  = 1.0;
          state.buses["bus_id_2"].vi  = 0.0;
          success                    *= rejects(data, state);
        }
        {
          auto [data, state] = validProblem();
          auto& branch       = std::get<SystemDataT::BranchDataT>(data.devices[0]);
          branch.R           = std::numeric_limits<RealT>::max();
          branch.X           = std::numeric_limits<RealT>::max();
          SystemT model(data, state);
          success *= model.allocate() == 0;
        }

        return success.report(__func__);
      }

      TestOutcome handlesCompositeTopologyAndStateCorners()
      {
        auto [data, state] = validProblem();

        auto& first_branch =
            std::get<SystemDataT::BranchDataT>(data.devices[0]);
        first_branch.smax.reset();

        auto parallel_branch = first_branch;
        parallel_branch.id   = "BR2";
        parallel_branch.R    = 0.03;
        parallel_branch.X    = 0.12;
        parallel_branch.smax = 1.5;
        data.devices.emplace_back(parallel_branch);

        SystemDataT::GeneratorDataT offline_generator;
        offline_generator.id   = "G1";
        offline_generator.bus  = 1;
        offline_generator.pmin = 0.25;
        offline_generator.pmax = 1.0;
        offline_generator.qmin = -0.5;
        offline_generator.qmax = 0.5;
        offline_generator.mva  = 50.0;
        offline_generator.c0   = 7.0;
        offline_generator.c2   = 2.0;
        data.devices.emplace_back(offline_generator);

        SystemDataT::LoadDataT offline_load;
        offline_load.id   = "LD2";
        offline_load.bus  = 0;
        offline_load.pmin = -0.2;
        offline_load.pmax = 0.0;
        offline_load.qmin = -0.1;
        offline_load.qmax = 0.0;
        data.devices.emplace_back(offline_load);

        SystemDataT::ShuntDataT absent_shunt;
        absent_shunt.id  = "SH2";
        absent_shunt.bus = 0;
        absent_shunt.G   = 0.01;
        absent_shunt.B   = 0.02;
        data.devices.emplace_back(absent_shunt);

        state.devices["G0"].p       = 3.0;
        state.devices["G0"].q       = 4.0;
        state.devices["G1"].online  = false;
        state.devices["G1"].p       = 0.5;
        state.devices["G1"].q       = 0.25;
        state.devices["LD2"].online = false;
        state.devices["LD2"].p      = -0.1;
        state.devices["LD2"].q      = -0.05;
        state.devices["SH"].online  = false;

        SystemT    model(data, state);
        TestStatus success = model.allocate() == 0 && model.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        success        *= model.size() == 8;
        success        *= model.sizeConstraints() == 6;
        auto* jacobian  = model.getCsrJacobian();
        auto* hessian   = model.getCsrHessian();
        success        *= jacobian->getNnz() == 28;
        success        *= hessian->getNnz() == 12;
        success        *= validCsr(*jacobian, false);
        success        *= validCsr(*hessian, true);

        const auto* variables  = model.variables().getData(memory::HOST);
        success               *= near(variables[4], 3.0);
        success               *= near(variables[5], 4.0);
        success               *= near(variables[6], 0.0);
        success               *= near(variables[7], 0.0);

        const auto* lower  = model.variableLowerBounds().getData(memory::HOST);
        const auto* upper  = model.variableUpperBounds().getData(memory::HOST);
        success           *= near(lower[6], 0.0) && near(upper[6], 0.0);
        success           *= near(lower[7], 0.0) && near(upper[7], 0.0);

        const std::array<RealT, 6> multipliers{{0.2, -0.1, 0.3, 0.4, 0.5, -0.6}};
        success *= model.evaluateConstraints() == 0;
        success *= model.evaluateJacobian() == 0;
        success *= model.evaluateHessian(
                       1.25, multipliers.data(), multipliers.size())
                   == 0;
        success *= allFinite(matrixValues(*jacobian));
        success *= allFinite(matrixValues(*hessian));

        const auto solution  = model.solutionState();
        success             *= !solution.devices.contains("BR2");
        success             *= !solution.devices.contains("SH2");
        success             *= solution.buses.at("bus_id_0").injections.contains("SH2");
        success             *= solution.buses.at("bus_id_0").injections.contains("LD2");
        success             *= solution.devices.at("G1").online == false;
        success             *= solution.devices.at("G1").p == 0.0;
        success             *= solution.devices.at("G1").q == 0.0;
        success             *= solution.devices.at("LD2").p == -0.1;
        success             *= solution.devices.at("LD2").q == -0.05;
        success             *= solution.buses.at("bus_id_0").injections.at("LD2").ir == 0.0;
        success             *= solution.buses.at("bus_id_0").injections.at("LD2").ii == 0.0;
        success             *= solution.buses.at("bus_id_1").injections.at("SH").ir == 0.0;
        success             *= solution.buses.at("bus_id_1").injections.at("SH").ii == 0.0;

        return success.report(__func__);
      }

      TestOutcome allocatesExactStructures()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = true;

        success *= throws<std::logic_error>([&]()
                                            { model.initialize(); });
        success *= model.allocate() == 0;
        success *= model.allocate() == 0;
        success *= model.size() == 6;
        success *= model.sizeConstraints() == 6;
        success *= model.hasJacobian();
        success *= model.hasHessian();

        auto* jacobian  = model.getCsrJacobian();
        auto* hessian   = model.getCsrHessian();
        success        *= jacobian != nullptr;
        success        *= hessian != nullptr;
        if (jacobian != nullptr)
        {
          success *= jacobian->getNumRows() == 6;
          success *= jacobian->getNumColumns() == 6;
          success *= jacobian->getNnz() == 26;
          success *= validCsr(*jacobian, false);
        }
        if (hessian != nullptr)
        {
          success *= hessian->getNumRows() == 6;
          success *= hessian->getNumColumns() == 6;
          success *= hessian->getNnz() == 11;
          success *= validCsr(*hessian, true);
        }

        success *= throws<std::logic_error>([&]()
                                            { model.solutionState(); });
        success *= model.initialize() == 0;

        const auto* values  = model.variables().getData(memory::HOST);
        success            *= near(values[0], 1.0);
        success            *= near(values[1], 0.0);
        success            *= near(values[2], std::hypot(0.98, -0.10));
        success            *= near(values[3], std::atan2(-0.10, 0.98));
        success            *= near(values[4], 0.0);
        success            *= near(values[5], 0.0);

        const auto* lower  = model.variableLowerBounds().getData(memory::HOST);
        const auto* upper  = model.variableUpperBounds().getData(memory::HOST);
        success           *= near(lower[0], 0.90) && near(upper[0], 1.10);
        success           *= near(lower[1], 0.0) && near(upper[1], 0.0);
        success           *= near(lower[2], 0.80) && near(upper[2], 1.20);
        success           *= std::isinf(lower[3]) && lower[3] < 0.0;
        success           *= std::isinf(upper[3]) && upper[3] > 0.0;
        success           *= near(lower[4], -2.0) && near(upper[4], 2.0);
        success           *= near(lower[5], -3.0) && near(upper[5], 3.0);

        const auto* constraint_lower =
            model.constraintLowerBounds().getData(memory::HOST);
        const auto* constraint_upper =
            model.constraintUpperBounds().getData(memory::HOST);
        for (IdxT row = 0; row < 4; ++row)
        {
          success *= near(constraint_lower[row], 0.0);
          success *= near(constraint_upper[row], 0.0);
        }
        for (IdxT row = 4; row < 6; ++row)
        {
          success *= near(constraint_lower[row], 0.0);
          success *= near(constraint_upper[row], 4.0);
        }

        auto [reversed_data, reversed_state] = validProblem();
        std::reverse(reversed_data.buses.begin(), reversed_data.buses.end());
        const RealT gauge_angle             = 0.25;
        reversed_state.buses["bus_id_0"].vr = std::cos(gauge_angle);
        reversed_state.buses["bus_id_0"].vi = std::sin(gauge_angle);
        SystemT reversed(reversed_data, reversed_state);
        success *= reversed.allocate() == 0;
        lower    = reversed.variableLowerBounds().getData(memory::HOST);
        upper    = reversed.variableUpperBounds().getData(memory::HOST);
        success *= std::isinf(lower[1]) && lower[1] < 0.0;
        success *= std::isinf(upper[1]) && upper[1] > 0.0;
        success *= near(lower[3], gauge_angle);
        success *= near(upper[3], gauge_angle);

        return success.report(__func__);
      }

      TestOutcome evaluatesExactDerivatives()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = model.allocate() == 0 && model.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        auto*                      variables = model.variables().getData(memory::HOST);
        const std::array<RealT, 6> point{{1.01, 0.02, 0.99, -0.08, 0.80, 0.15}};
        std::copy(point.begin(), point.end(), variables);
        model.variables().setDataUpdated(memory::HOST);

        auto*       jacobian         = model.getCsrJacobian();
        auto*       hessian          = model.getCsrHessian();
        const auto* jacobian_rows    = jacobian->getRowData(memory::HOST);
        const auto* jacobian_columns = jacobian->getColData(memory::HOST);
        const auto* hessian_rows     = hessian->getRowData(memory::HOST);
        const auto* hessian_columns  = hessian->getColData(memory::HOST);

        success                   *= model.evaluateObjective() == 0;
        success                   *= near(model.objective(), 2.92, 1.0e-12);
        success                   *= model.evaluateObjectiveGradient() == 0;
        const auto gradient_first  = copy(model.objectiveGradient());
        success                   *= near(gradient_first[4], 2.8, 1.0e-12);
        for (IdxT index = 0; index < model.size(); ++index)
        {
          if (index != 4)
          {
            success *= near(gradient_first[index], 0.0, 1.0e-12);
          }
        }

        success                      *= model.evaluateConstraints() == 0;
        const auto constraints_first  = copy(model.constraints());
        success                      *= model.evaluateConstraints() == 0;
        success                      *= same(constraints_first, copy(model.constraints()), 1.0e-13);

        success                   *= model.evaluateJacobian() == 0;
        const auto jacobian_first  = matrixValues(*jacobian);
        success                   *= model.evaluateJacobian() == 0;
        success                   *= same(jacobian_first, matrixValues(*jacobian), 1.0e-13);

        const std::array<RealT, 6> multipliers{{0.3, -0.2, 0.1, 0.4, 0.5, -0.7}};
        success *= model.evaluateHessian(
                       2.0, multipliers.data(), multipliers.size())
                   == 0;
        const auto hessian_first  = matrixValues(*hessian);
        success                  *= model.evaluateHessian(
                       2.0, multipliers.data(), multipliers.size())
                   == 0;
        success *= same(hessian_first, matrixValues(*hessian), 1.0e-13);
        success *= near(matrixValue(*hessian, 4, 4), 2.0, 1.0e-12);

        success *= model.evaluateHessian(
                       1.0, multipliers.data(), multipliers.size() - 1)
                   != 0;
        success *= model.evaluateHessian(1.0, nullptr, multipliers.size()) != 0;

        std::fill(variables, variables + model.size(), RealT{0});
        variables[0] = 1.0;
        variables[2] = 1.0;
        model.variables().setDataUpdated(memory::HOST);
        const std::array<RealT, 6> zeros{};
        success *= model.evaluateJacobian() == 0;
        success *= model.evaluateHessian(0.0, zeros.data(), zeros.size()) == 0;
        success *= allFinite(matrixValues(*jacobian));
        success *= allFinite(matrixValues(*hessian));

        success *= jacobian_rows == jacobian->getRowData(memory::HOST);
        success *= jacobian_columns == jacobian->getColData(memory::HOST);
        success *= hessian_rows == hessian->getRowData(memory::HOST);
        success *= hessian_columns == hessian->getColData(memory::HOST);

        return success.report(__func__);
      }

      TestOutcome writesCompatibleSolutionState()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = model.allocate() == 0 && model.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        auto*                      values = model.variables().getData(memory::HOST);
        const std::array<RealT, 6> point{{1.02, 0.10, 0.97, -0.20, 0.75, 0.12}};
        std::copy(point.begin(), point.end(), values);
        model.variables().setDataUpdated(memory::HOST);

        const auto solution  = model.solutionState();
        success             *= solution.header.has_value();
        success             *= solution.header->description == "preserve me";
        success             *= solution.devices.contains("OTHER");
        success             *= solution.devices.at("OTHER").online == false;
        success             *= solution.devices.at("BR").open == false;
        success             *= solution.devices.at("BR").tap == 1.0;
        success             *= solution.devices.at("BR").phase == 0.0;
        success             *= solution.devices.at("LD").p == -1.0;
        success             *= solution.devices.at("LD").q == -0.2;
        success             *= solution.devices.at("G0").p == 0.75;
        success             *= solution.devices.at("G0").q == 0.12;
        success             *= !solution.devices.at("G0").online.has_value();

        const auto& bus0  = solution.buses.at("bus_id_0");
        const auto& bus1  = solution.buses.at("bus_id_1");
        success          *= near(*bus0.vr, 1.02 * std::cos(0.10), 1.0e-12);
        success          *= near(*bus0.vi, 1.02 * std::sin(0.10), 1.0e-12);
        success          *= near(*bus1.vr, 0.97 * std::cos(-0.20), 1.0e-12);
        success          *= near(*bus1.vi, 0.97 * std::sin(-0.20), 1.0e-12);
        success          *= bus0.injections.contains("G0");
        success          *= bus0.injections.contains("OTHER_INJECTION");
        success          *= bus1.injections.contains("LD");
        success          *= bus1.injections.contains("SH");

        std::ostringstream serialized;
        Model::writeStateData(solution, serialized);
        std::istringstream serialized_input(serialized.str());
        const auto         round_trip  = Model::parseStateData(serialized_input);
        success                       *= round_trip.header->description == "preserve me";
        success                       *= round_trip.devices.at("G0").p == 0.75;
        success                       *= round_trip.devices.at("LD").p == -1.0;
        success                       *= round_trip.buses.at("bus_id_0").injections.contains("G0");
        success                       *= round_trip.buses.at("bus_id_1").injections.contains("SH");

        return success.report(__func__);
      }

    private:
      static std::pair<SystemDataT, Model::StateData> validProblem()
      {
        SystemDataT data;
        data.header.format_version  = 0;
        data.header.format_revision = 1;
        data.header.case_name       = "OPF System unit test";
        data.params.freq_base       = 60.0;
        data.params.va_base         = 100.0e6;

        SystemDataT::BusDataT bus0;
        bus0.id     = "B0";
        bus0.number = 0;
        bus0.kv     = 230.0;
        bus0.vmin   = 0.90;
        bus0.vmax   = 1.10;
        data.buses.push_back(bus0);

        SystemDataT::BusDataT bus1;
        bus1.id     = "B1";
        bus1.number = 1;
        bus1.kv     = 230.0;
        bus1.vmin   = 0.80;
        bus1.vmax   = 1.20;
        data.buses.push_back(bus1);

        SystemDataT::BranchDataT branch;
        branch.id   = "BR";
        branch.from = 0;
        branch.to   = 1;
        branch.R    = 0.02;
        branch.X    = 0.10;
        branch.G    = 0.01;
        branch.B    = 0.04;
        branch.smax = 2.0;
        data.devices.emplace_back(branch);

        SystemDataT::GeneratorDataT generator;
        generator.id   = "G0";
        generator.bus  = 0;
        generator.pmin = -2.0;
        generator.pmax = 2.0;
        generator.qmin = -3.0;
        generator.qmax = 3.0;
        generator.mva  = 100.0;
        generator.c0   = 1.0;
        generator.c1   = 2.0;
        generator.c2   = 0.5;
        data.devices.emplace_back(generator);

        SystemDataT::LoadDataT load;
        load.id   = "LD";
        load.bus  = 1;
        load.pmin = -1.1;
        load.pmax = -0.9;
        load.qmin = -0.3;
        load.qmax = -0.1;
        data.devices.emplace_back(load);

        SystemDataT::ShuntDataT shunt;
        shunt.id  = "SH";
        shunt.bus = 1;
        shunt.G   = 0.02;
        shunt.B   = -0.08;
        data.devices.emplace_back(shunt);

        Model::StateData state;
        state.header.emplace();
        state.header->description                                = "preserve me";
        state.buses["bus_id_0"].vr                               = 1.0;
        state.buses["bus_id_0"].vi                               = 0.0;
        state.buses["bus_id_0"].injections["OTHER_INJECTION"].ir = 4.0;
        state.buses["bus_id_1"].vr                               = 0.98;
        state.buses["bus_id_1"].vi                               = -0.10;
        state.devices["BR"].open                                 = false;
        state.devices["BR"].tap                                  = 1.0;
        state.devices["BR"].phase                                = 0.0;
        state.devices["LD"].online                               = true;
        state.devices["LD"].p                                    = -1.0;
        state.devices["LD"].q                                    = -0.2;
        state.devices["SH"].online                               = true;
        state.devices["OTHER"].online                            = false;
        return {data, state};
      }

      static bool rejects(const SystemDataT& data, const Model::StateData& state)
      {
        try
        {
          SystemT model(data, state);
          static_cast<void>(model.allocate());
        }
        catch (const std::invalid_argument&)
        {
          return true;
        }
        return false;
      }

      static bool near(RealT value,
                       RealT reference,
                       RealT tolerance = 1.0e-10)
      {
        return std::abs(value - reference)
               <= tolerance * (1.0 + std::abs(reference));
      }

      template <typename VectorT>
      static std::vector<RealT> copy(const VectorT& vector)
      {
        const auto* values = vector.getData(memory::HOST);
        return {values, values + vector.getSize()};
      }

      template <typename MatrixT>
      static std::vector<RealT> matrixValues(MatrixT& matrix)
      {
        const auto* values = matrix.getValues(memory::HOST);
        return {values, values + matrix.getNnz()};
      }

      static bool same(const std::vector<RealT>& left,
                       const std::vector<RealT>& right,
                       RealT                     tolerance)
      {
        if (left.size() != right.size())
        {
          return false;
        }
        for (std::size_t index = 0; index < left.size(); ++index)
        {
          if (!near(left[index], right[index], tolerance))
          {
            return false;
          }
        }
        return true;
      }

      static bool allFinite(const std::vector<RealT>& values)
      {
        return std::all_of(values.begin(), values.end(), [](RealT value)
                           { return std::isfinite(value); });
      }

      template <typename MatrixT>
      static bool validCsr(MatrixT& matrix, bool lower_triangle)
      {
        const auto* rows    = matrix.getRowData(memory::HOST);
        const auto* columns = matrix.getColData(memory::HOST);
        if (rows[0] != 0 || rows[matrix.getNumRows()] != matrix.getNnz())
        {
          return false;
        }
        for (IdxT row = 0; row < matrix.getNumRows(); ++row)
        {
          if (rows[row] > rows[row + 1])
          {
            return false;
          }
          for (IdxT entry = rows[row]; entry < rows[row + 1]; ++entry)
          {
            if (entry > rows[row] && columns[entry - 1] >= columns[entry])
            {
              return false;
            }
            if (lower_triangle && columns[entry] > row)
            {
              return false;
            }
          }
        }
        return true;
      }

      template <typename MatrixT>
      static RealT matrixValue(MatrixT& matrix, IdxT row, IdxT column)
      {
        const auto* rows    = matrix.getRowData(memory::HOST);
        const auto* columns = matrix.getColData(memory::HOST);
        const auto* values  = matrix.getValues(memory::HOST);
        for (IdxT entry = rows[row]; entry < rows[row + 1]; ++entry)
        {
          if (columns[entry] == column)
          {
            return values[entry];
          }
        }
        return std::numeric_limits<RealT>::quiet_NaN();
      }
    };
  } // namespace Testing
} // namespace GridKit
