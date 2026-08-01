#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <GridKit/Model/OPF/System.hpp>
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

      TestOutcome rejectInvalidSystems()
      {
        TestStatus success = true;

        {
          auto [data, state]   = validProblem();
          data.params.va_base  = 0.0;
          success             *= rejects(data, state);
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
          auto [data, state]       = validProblem();
          data.buses[0].bus_class  = OPF::BusClass::BUS;
          success                 *= rejects(data, state);
        }
        {
          auto [data, state]       = validProblem();
          data.buses[1].bus_class  = OPF::BusClass::SLACK;
          success                 *= rejects(data, state);
        }
        {
          auto [data, state]       = validProblem();
          data.buses[1].bus_class  = static_cast<OPF::BusClass>(99);
          success                 *= rejects(data, state);
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
          branch.to           = 99;
          success            *= rejects(data, state);
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
          branch.smax         = -1.0;
          success            *= rejects(data, state);
        }
        {
          auto [data, state]       = validProblem();
          state.devices["BR"].tap  = 0.0;
          success                 *= rejects(data, state);
        }
        {
          auto [data, state]       = validProblem();
          state.devices["BR"].tap  = std::numeric_limits<RealT>::denorm_min();
          success                 *= rejects(data, state);
        }
        {
          auto [data, state]        = validProblem();
          state.devices["BR"].open  = true;
          success                  *= rejects(data, state);
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
          generator.pmin      = 2.0;
          generator.pmax      = 1.0;
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
          auto [data, state]  = validProblem();
          auto& shunt         = std::get<SystemDataT::ShuntDataT>(data.devices[3]);
          shunt.B             = std::numeric_limits<RealT>::infinity();
          success            *= rejects(data, state);
        }
        {
          auto [data, state]  = validProblem();
          auto& shunt         = std::get<SystemDataT::ShuntDataT>(data.devices[3]);
          shunt.id            = "LD";
          success            *= rejects(data, state);
        }

        return success.report(__func__);
      }

      TestOutcome allocateInitializeAndBounds()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = true;

        success *= throws<std::logic_error>([&]()
                                            { model.initialize(); });
        success *= model.allocate() == 0;
        success *= model.allocate() == 0;
        success *= model.size() == 6;
        success *= model.sizeResidual() == 6;
        success *= model.nnz() == 26;
        success *= model.hasJacobian();
        success *= throws<std::logic_error>([&]()
                                            { model.solutionState(); });
        success *= model.initialize() == 0;

        const RealT* values  = model.y().getData();
        const RealT  vm1     = std::hypot(0.98, -0.10);
        const RealT  va1     = std::atan2(-0.10, 0.98);
        success             *= near(values[0], 1.0);
        success             *= near(values[1], 0.0);
        success             *= near(values[2], vm1);
        success             *= near(values[3], va1);
        success             *= near(values[4], 0.0);
        success             *= near(values[5], 0.0);

        SystemT::RealVectorT lower;
        SystemT::RealVectorT upper;
        success                     *= model.getVariableBounds(lower, upper) == 0;
        const RealT* variable_lower  = lower.getData();
        const RealT* variable_upper  = upper.getData();
        success                     *= near(variable_lower[0], 0.90);
        success                     *= near(variable_upper[0], 1.10);
        success                     *= near(variable_lower[1], 0.0);
        success                     *= near(variable_upper[1], 0.0);
        success                     *= near(variable_lower[2], 0.80);
        success                     *= near(variable_upper[2], 1.20);
        success                     *= std::isinf(variable_lower[3]) && variable_lower[3] < 0.0;
        success                     *= std::isinf(variable_upper[3]) && variable_upper[3] > 0.0;
        success                     *= near(variable_lower[4], -2.0);
        success                     *= near(variable_upper[4], 2.0);
        success                     *= near(variable_lower[5], -3.0);
        success                     *= near(variable_upper[5], 3.0);

        success                     *= model.getResidualBounds(lower, upper) == 0;
        const RealT* residual_lower  = lower.getData();
        const RealT* residual_upper  = upper.getData();
        for (IdxT i = 0; i < 4; ++i)
        {
          success *= near(residual_lower[i], 0.0);
          success *= near(residual_upper[i], 0.0);
        }
        for (IdxT i = 4; i < 6; ++i)
        {
          success *= near(residual_lower[i], 0.0);
          success *= near(residual_upper[i], 4.0);
        }

        auto* jacobian  = model.getCsrJacobian();
        success        *= jacobian != nullptr;
        if (jacobian != nullptr)
        {
          success          *= jacobian->getNumRows() == model.sizeResidual();
          success          *= jacobian->getNumColumns() == model.size();
          success          *= jacobian->getNnz() == model.nnz();
          const IdxT* rows  = jacobian->getRowData();
          const IdxT* cols  = jacobian->getColData();
          success          *= rows[0] == 0;
          success          *= rows[model.sizeResidual()] == model.nnz();
          for (IdxT row = 0; row < model.sizeResidual(); ++row)
          {
            success *= rows[row] <= rows[row + 1];
            for (IdxT slot = rows[row] + 1; slot < rows[row + 1]; ++slot)
            {
              success *= cols[slot - 1] < cols[slot];
            }
          }
        }

        auto [outside_data, outside_state] = validProblem();
        outside_state.buses["bus_id_1"].vr = 1.30;
        outside_state.buses["bus_id_1"].vi = 0.0;
        outside_state.devices["G0"].p      = 3.0;
        outside_state.devices["G0"].q      = -4.0;
        SystemT outside_model(outside_data, outside_state);
        success                     *= outside_model.allocate() == 0;
        success                     *= outside_model.initialize() == 0;
        const RealT* outside_values  = outside_model.y().getData();
        success                     *= near(outside_values[2], 1.30);
        success                     *= near(outside_values[4], 3.0);
        success                     *= near(outside_values[5], -4.0);

        return success.report(__func__);
      }

      TestOutcome rejectNonfiniteEvaluations()
      {
        SystemDataT data;
        data.header.format_version  = 0;
        data.header.format_revision = 1;
        data.header.case_name       = "OPF finite evaluation test";
        data.params.freq_base       = 60.0;
        data.params.va_base         = 100.0e6;

        SystemDataT::BusDataT bus;
        bus.bus_class = OPF::BusClass::SLACK;
        bus.id        = "B0";
        bus.number    = 0;
        bus.kv        = 230.0;
        data.buses.push_back(bus);

        const RealT half_max = std::numeric_limits<RealT>::max() / 2.0;
        for (int i = 0; i < 3; ++i)
        {
          SystemDataT::ShuntDataT shunt;
          shunt.id  = "SH" + std::to_string(i);
          shunt.bus = 0;
          shunt.G   = half_max;
          shunt.B   = 0.0;
          data.devices.emplace_back(shunt);
        }

        Model::StateData state;
        state.buses["bus_id_0"].vr = 1.0;
        state.buses["bus_id_0"].vi = 0.0;

        SystemT    model(data, state);
        TestStatus success  = model.allocate() == 0 && model.initialize() == 0;
        success            *= model.evaluateResidual() != 0;
        success            *= model.evaluateJacobian() != 0;

        SystemDataT objective_data;
        objective_data.header = data.header;
        objective_data.params = data.params;
        objective_data.buses.push_back(bus);
        for (int i = 0; i < 2; ++i)
        {
          SystemDataT::GeneratorDataT generator;
          generator.id  = "G" + std::to_string(i);
          generator.bus = 0;
          generator.mva = 100.0;
          generator.c0  = 0.75 * std::numeric_limits<RealT>::max();
          objective_data.devices.emplace_back(generator);
        }

        SystemT objective_model(objective_data, state);
        success *= objective_model.allocate() == 0;
        success *= objective_model.initialize() == 0;
        success *= objective_model.evaluateObjective() != 0;

        return success.report(__func__);
      }

      TestOutcome evaluateAndDifferentiate()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = model.allocate() == 0 && model.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        RealT* variables = model.y().getData();
        variables[0]     = 1.01;
        variables[1]     = 0.02;
        variables[2]     = 0.99;
        variables[3]     = -0.08;
        variables[4]     = 0.80;
        variables[5]     = 0.15;
        model.y().setDataUpdated();

        success                   *= model.evaluateObjective() == 0;
        success                   *= near(model.objective(), 2.92, 1.0e-12);
        success                   *= model.evaluateResidual() == 0;
        const auto residual_first  = copy(model.getResidual());
        success                   *= model.evaluateResidual() == 0;
        success                   *= same(residual_first, copy(model.getResidual()), 1.0e-13);

        success                   *= model.evaluateJacobian() == 0;
        const auto gradient_first  = copy(model.getObjectiveGradient());
        const auto jacobian_first  = jacobianValues(model);
        success                   *= model.evaluateJacobian() == 0;
        success                   *= same(gradient_first,
                        copy(model.getObjectiveGradient()),
                        1.0e-13);
        success                   *= same(jacobian_first, jacobianValues(model), 1.0e-13);
        success                   *= near(gradient_first[4], 2.8, 1.0e-12);
        for (IdxT i = 0; i < model.size(); ++i)
        {
          if (i != 4)
          {
            success *= near(gradient_first[i], 0.0, 1.0e-12);
          }
        }

        const RealT step           = 1.0e-6;
        const auto  base           = copy(model.y());
        const auto  dense_jacobian = denseJacobian(model);
        for (IdxT column = 0; column < model.size(); ++column)
        {
          variables[column] = base[column] + step;
          model.y().setDataUpdated();
          success         *= model.evaluateResidual() == 0;
          const auto plus  = copy(model.getResidual());

          variables[column] = base[column] - step;
          model.y().setDataUpdated();
          success          *= model.evaluateResidual() == 0;
          const auto minus  = copy(model.getResidual());

          for (IdxT row = 0; row < model.sizeResidual(); ++row)
          {
            const RealT difference  = (plus[row] - minus[row]) / (2.0 * step);
            success                *= near(dense_jacobian[row * model.size() + column],
                            difference,
                            2.0e-5);
          }
          variables[column] = base[column];
        }
        model.y().setDataUpdated();

        for (IdxT column = 0; column < model.size(); ++column)
        {
          variables[column] = base[column] + step;
          model.y().setDataUpdated();
          success          *= model.evaluateObjective() == 0;
          const RealT plus  = model.objective();

          variables[column] = base[column] - step;
          model.y().setDataUpdated();
          success           *= model.evaluateObjective() == 0;
          const RealT minus  = model.objective();

          const RealT difference  = (plus - minus) / (2.0 * step);
          success                *= near(gradient_first[column], difference, 1.0e-7);
          variables[column]       = base[column];
        }
        model.y().setDataUpdated();

        return success.report(__func__);
      }

      TestOutcome writeSolutionState()
      {
        auto [data, state] = validProblem();
        SystemT    model(data, state);
        TestStatus success = model.allocate() == 0 && model.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        RealT* values = model.y().getData();
        values[0]     = 1.02;
        values[1]     = 0.10;
        values[2]     = 0.97;
        values[3]     = -0.20;
        values[4]     = 0.75;
        values[5]     = 0.12;
        model.y().setDataUpdated();

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
        success             *= solution.devices.contains("G0");
        success             *= solution.devices.at("G0").p == 0.75;
        success             *= solution.devices.at("G0").q == 0.12;

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
        success                       *= round_trip.header.has_value();
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

        SystemDataT::BusDataT slack;
        slack.bus_class = OPF::BusClass::SLACK;
        slack.id        = "B0";
        slack.number    = 0;
        slack.kv        = 230.0;
        slack.vmin      = 0.90;
        slack.vmax      = 1.10;
        data.buses.push_back(slack);

        SystemDataT::BusDataT bus;
        bus.bus_class = OPF::BusClass::BUS;
        bus.id        = "B1";
        bus.number    = 1;
        bus.kv        = 230.0;
        bus.vmin      = 0.80;
        bus.vmax      = 1.20;
        data.buses.push_back(bus);

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

      static bool rejects(const SystemDataT&      data,
                          const Model::StateData& state)
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

      static bool near(RealT value, RealT reference, RealT tolerance = 1.0e-10)
      {
        return std::abs(value - reference)
               <= tolerance * (1.0 + std::abs(reference));
      }

      template <typename VectorT>
      static std::vector<RealT> copy(const VectorT& vector)
      {
        const RealT* values = vector.getData();
        return {values, values + vector.getSize()};
      }

      static bool same(const std::vector<RealT>& left,
                       const std::vector<RealT>& right,
                       RealT                     tolerance)
      {
        if (left.size() != right.size())
        {
          return false;
        }
        for (std::size_t i = 0; i < left.size(); ++i)
        {
          if (!near(left[i], right[i], tolerance))
          {
            return false;
          }
        }
        return true;
      }

      static std::vector<RealT> jacobianValues(SystemT& model)
      {
        auto*        jacobian = model.getCsrJacobian();
        const RealT* values   = jacobian->getValues();
        return {values, values + jacobian->getNnz()};
      }

      static std::vector<RealT> denseJacobian(SystemT& model)
      {
        std::vector<RealT> dense(model.sizeResidual() * model.size(), 0.0);
        auto*              jacobian = model.getCsrJacobian();
        const IdxT*        rows     = jacobian->getRowData();
        const IdxT*        columns  = jacobian->getColData();
        const RealT*       values   = jacobian->getValues();
        for (IdxT row = 0; row < model.sizeResidual(); ++row)
        {
          for (IdxT slot = rows[row]; slot < rows[row + 1]; ++slot)
          {
            dense[row * model.size() + columns[slot]] = values[slot];
          }
        }
        return dense;
      }
    };

  } // namespace Testing
} // namespace GridKit
