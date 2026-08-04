#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include <GridKit/Model/OPF/Bus/Bus.hpp>
#include <GridKit/Model/OPF/Generator/Generator.hpp>
#include <GridKit/Model/OPF/Load/Load.hpp>
#include <GridKit/Model/OPF/Shunt/Shunt.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class OPFComponentTests
    {
    public:
      using IdxT  = std::size_t;
      using RealT = double;

      TestOutcome busAndLoadExactEmptyStructures()
      {
        OPF::BusData<RealT, IdxT> bus_data;
        bus_data.id     = "B7";
        bus_data.number = 7;
        bus_data.kv     = 230.0;
        bus_data.vmin   = 0.9;
        bus_data.vmax   = 1.1;

        Model::BusState bus_state;
        bus_state.vr = 0.6;
        bus_state.vi = 0.8;

        OPF::Bus<RealT, IdxT> bus(bus_data, bus_state, {3, 1}, {7, 8});
        std::array<RealT, 9>  variables{};
        std::array<RealT, 9>  lower{};
        std::array<RealT, 9>  upper{};
        std::array<RealT, 0>  no_multipliers{};

        TestStatus success  = bus.hasJacobian() && bus.hasHessian();
        success            *= bus.variableIndices().size() == 2;
        success            *= bus.constraintIndices().empty();
        success            *= bus.balanceIndices().size() == 2;
        success            *= bus.voltageMagnitudeIndex() == 3;
        success            *= bus.voltageAngleIndex() == 1;
        success            *= bus.activeBalanceIndex() == 7;
        success            *= bus.reactiveBalanceIndex() == 8;
        success            *= bus.jacobianEntries().empty();
        success            *= bus.hessianEntries().empty();
        success            *= bus.initialize(variables.data()) == 0;
        success            *= near(variables[3], 1.0);
        success            *= near(variables[1], std::atan2(0.8, 0.6));
        success            *= bus.setVariableBounds(lower.data(), upper.data()) == 0;
        success            *= near(lower[3], 0.9) && near(upper[3], 1.1);
        success            *= std::isinf(lower[1]) && lower[1] < 0.0;
        success            *= std::isinf(upper[1]) && upper[1] > 0.0;
        success            *= bus.gatherVariables(variables.data()) == 0;
        success            *= bus.evaluateObjectiveGradient() == 0;
        success            *= bus.objectiveGradientValues().size() == 2;
        success            *= allNear(bus.objectiveGradientValues(), 0.0);
        success            *= bus.evaluateJacobian() == 0;
        success            *= bus.evaluateHessian(1.0, no_multipliers) == 0;
        success            *= bus.jacobianValues().empty();
        success            *= bus.hessianValues().empty();

        OPF::LoadData<RealT, IdxT> load_data;
        load_data.id  = "LD";
        load_data.bus = 7;
        Model::DeviceState load_state;
        load_state.online = true;
        load_state.p      = -0.9;
        load_state.q      = -0.25;

        OPF::Load<RealT, IdxT>     load(load_data, load_state, {7, 8}, {3, 1});
        const std::array<RealT, 2> load_multipliers{{0.4, -0.7}};
        success *= load.hasJacobian() && load.hasHessian();
        success *= load.variableIndices().empty();
        success *= load.jacobianEntries().empty();
        success *= load.hessianEntries().empty();
        success *= load.gatherVariables(variables.data()) == 0;
        success *= load.evaluateConstraints() == 0;
        success *= near(load.constraintValues()[0], -0.9);
        success *= near(load.constraintValues()[1], -0.25);
        success *= load.evaluateJacobian() == 0;
        success *= load.evaluateHessian(1.0, load_multipliers) == 0;
        success *= load.jacobianValues().empty();
        success *= load.hessianValues().empty();

        load_state.online = false;
        OPF::Load<RealT, IdxT> offline_load(
            load_data, load_state, {7, 8}, {3, 1});
        success *= offline_load.gatherVariables(variables.data()) == 0;
        success *= offline_load.evaluateConstraints() == 0;
        success *= allNear(offline_load.constraintValues(), 0.0);
        success *= offline_load.hasJacobian() && offline_load.hasHessian();

        return success.report(__func__);
      }

      TestOutcome generatorExactDerivatives()
      {
        OPF::GeneratorData<RealT, IdxT> data;
        data.id   = "G";
        data.bus  = 2;
        data.pmin = -1.0;
        data.pmax = 2.0;
        data.qmin = -0.5;
        data.qmax = 0.8;
        data.mva  = 100.0;
        data.c0   = 1.2;
        data.c1   = -0.5;
        data.c2   = 2.0;

        Model::DeviceState state;
        state.online = true;
        state.p      = 0.7;
        state.q      = -0.2;

        OPF::Generator<RealT, IdxT> generator(
            data, state, {2, 4}, {3, 5}, {0, 1});
        std::array<RealT, 6>       variables{{1.03, -0.16, 0.7, 0.0, -0.2, 0.0}};
        const std::array<RealT, 2> multipliers{{0.3, -0.4}};

        TestStatus success  = generator.hasJacobian() && generator.hasHessian();
        success            *= generator.jacobianEntries().size() == 2;
        success            *= generator.jacobianEntries()[0].variable == 0;
        success            *= generator.jacobianEntries()[0].constraint == 0;
        success            *= generator.jacobianEntries()[1].variable == 1;
        success            *= generator.jacobianEntries()[1].constraint == 1;
        success            *= generator.hessianEntries().size() == 1;
        success            *= generator.hessianEntries()[0].row == 0;
        success            *= generator.hessianEntries()[0].column == 0;

        success *= generator.gatherVariables(variables.data()) == 0;
        success *= generator.evaluateObjective() == 0;
        success *= near(generator.objective(), 1.83);
        success *= generator.evaluateObjectiveGradient() == 0;
        success *= near(generator.objectiveGradientValues()[0], 2.3);
        success *= near(generator.objectiveGradientValues()[1], 0.0);
        success *= generator.evaluateConstraints() == 0;
        success *= near(generator.constraintValues()[0], 0.7);
        success *= near(generator.constraintValues()[1], -0.2);
        success *= generator.evaluateJacobian() == 0;
        success *= near(generator.jacobianValues()[0], 1.0);
        success *= near(generator.jacobianValues()[1], 1.0);
        success *= generator.evaluateHessian(1.5, multipliers) == 0;
        success *= near(generator.hessianValues()[0], 6.0);

        state.online = false;
        OPF::Generator<RealT, IdxT> offline(
            data, state, {2, 4}, {3, 5}, {0, 1});
        variables[2]  = 0.0;
        variables[4]  = 0.0;
        success      *= offline.gatherVariables(variables.data()) == 0;
        success      *= offline.evaluateObjective() == 0;
        success      *= near(offline.objective(), 0.0);
        success      *= offline.evaluateObjectiveGradient() == 0;
        success      *= allNear(offline.objectiveGradientValues(), 0.0);
        success      *= offline.evaluateConstraints() == 0;
        success      *= allNear(offline.constraintValues(), 0.0);
        success      *= offline.evaluateJacobian() == 0;
        success      *= near(offline.jacobianValues()[0], 1.0);
        success      *= near(offline.jacobianValues()[1], 1.0);
        success      *= offline.evaluateHessian(1.5, multipliers) == 0;
        success      *= allNear(offline.hessianValues(), 0.0);
        success      *= offline.jacobianEntries().size()
                   == generator.jacobianEntries().size();
        success *= offline.hessianEntries().size()
                   == generator.hessianEntries().size();

        return success.report(__func__);
      }

      TestOutcome shuntExactDerivatives()
      {
        OPF::ShuntData<RealT, IdxT> data;
        data.id  = "SH";
        data.bus = 4;
        data.G   = 0.04;
        data.B   = -0.12;

        Model::DeviceState state;
        state.online = true;
        OPF::Shunt<RealT, IdxT>    shunt(data, state, 2, {5, 6}, {2, 1});
        std::array<RealT, 7>       variables{{0.0, 0.31, 1.07, 0.0, 0.0, 0.0, 0.0}};
        const std::array<RealT, 2> multipliers{{0.7, -0.4}};

        TestStatus success  = shunt.hasJacobian() && shunt.hasHessian();
        success            *= shunt.jacobianEntries().size() == 2;
        success            *= shunt.jacobianEntries()[0].variable == 0;
        success            *= shunt.jacobianEntries()[0].constraint == 0;
        success            *= shunt.jacobianEntries()[1].variable == 0;
        success            *= shunt.jacobianEntries()[1].constraint == 1;
        success            *= shunt.hessianEntries().size() == 1;
        success            *= shunt.gatherVariables(variables.data()) == 0;
        success            *= shunt.evaluateConstraints() == 0;
        success            *= near(shunt.constraintValues()[0], -0.045796);
        success            *= near(shunt.constraintValues()[1], -0.137388);
        success            *= shunt.evaluateJacobian() == 0;
        success            *= near(shunt.jacobianValues()[0], -0.0856);
        success            *= near(shunt.jacobianValues()[1], -0.2568);
        success            *= shunt.evaluateHessian(3.0, multipliers) == 0;
        success            *= near(shunt.hessianValues()[0], 0.04);

        state.online = false;
        OPF::Shunt<RealT, IdxT> offline(data, state, 2, {5, 6}, {2, 1});
        success *= offline.gatherVariables(variables.data()) == 0;
        success *= offline.evaluateConstraints() == 0;
        success *= offline.evaluateJacobian() == 0;
        success *= offline.evaluateHessian(3.0, multipliers) == 0;
        success *= allNear(offline.constraintValues(), 0.0);
        success *= allNear(offline.jacobianValues(), 0.0);
        success *= allNear(offline.hessianValues(), 0.0);
        success *= offline.jacobianEntries().size() == shunt.jacobianEntries().size();
        success *= offline.hessianEntries().size() == shunt.hessianEntries().size();

        return success.report(__func__);
      }

    private:
      static bool near(RealT value, RealT expected, RealT tolerance = 1.0e-11)
      {
        return std::abs(value - expected)
               <= tolerance * (RealT{1} + std::abs(expected));
      }

      template <typename SpanT>
      static bool allNear(const SpanT& values,
                          RealT        expected,
                          RealT        tolerance = 1.0e-11)
      {
        for (const auto value : values)
        {
          if (!near(static_cast<RealT>(value), expected, tolerance))
          {
            return false;
          }
        }
        return true;
      }
    };
  } // namespace Testing
} // namespace GridKit
