#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <GridKit/Model/OPF/Bus/Bus.hpp>
#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/Generator/Generator.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class BusGeneratorTests
    {
    public:
      using ScalarT    = double;
      using IdxT       = std::size_t;
      using ComponentT = OPF::Component<ScalarT, IdxT>;
      using BusT       = OPF::Bus<ScalarT, IdxT>;
      using GeneratorT = OPF::Generator<ScalarT, IdxT>;

      TestOutcome polymorphicSizesAndBindings()
      {
        OPF::BusData<ScalarT, IdxT> bus_data;
        bus_data.number = 4;
        bus_data.kv     = 115.0;
        Model::BusState bus_state;
        bus_state.vr = 1.0;
        bus_state.vi = 0.0;
        BusT bus(bus_data, bus_state);

        ComponentT* bus_component = &bus;

        TestStatus success  = true;
        success            *= bus_component->sizeInternalVariables() == 2;
        success            *= bus_component->sizeInternalConstraints() == 2;
        success            *= bus_component->nnz() == 0;

        bus.setVariableOffset(3);
        bus.setConstraintOffset(7);
        success *= bus.variables().template internalIndex<OPF::BusInternalVariables::VM>() == 3;
        success *= bus.variables().template internalIndex<OPF::BusInternalVariables::VA>() == 4;
        success *= bus.constraints().template internalIndex<OPF::BusInternalConstraints::DIVP>() == 7;
        success *= bus.constraints().template internalIndex<OPF::BusInternalConstraints::DIVQ>() == 8;

        OPF::GeneratorData<ScalarT, IdxT> generator_data;
        generator_data.id  = "G1";
        generator_data.bus = 4;
        generator_data.mva = 100.0;
        Model::DeviceState generator_state;
        GeneratorT         generator(generator_data, generator_state);

        ComponentT* generator_component  = &generator;
        success                         *= generator_component->sizeInternalVariables() == 2;
        success                         *= generator_component->sizeInternalConstraints() == 0;
        success                         *= generator_component->nnz() == 2;

        generator.setVariableOffset(11);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VM>(3);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VA>(4);
        generator.constraints().template bindExternal<OPF::GeneratorExternalConstraints::DIVP>(7);
        generator.constraints().template bindExternal<OPF::GeneratorExternalConstraints::DIVQ>(8);

        success *= generator.variables().template internalIndex<OPF::GeneratorInternalVariables::PG>() == 11;
        success *= generator.variables().template internalIndex<OPF::GeneratorInternalVariables::QG>() == 12;
        success *= generator.variables().template externalIndex<OPF::GeneratorExternalVariables::VM>() == 3;
        success *= generator.variables().template externalIndex<OPF::GeneratorExternalVariables::VA>() == 4;
        success *= generator.constraints().template externalIndex<OPF::GeneratorExternalConstraints::DIVP>() == 7;
        success *= generator.constraints().template externalIndex<OPF::GeneratorExternalConstraints::DIVQ>() == 8;

        return success.report(__func__);
      }

      TestOutcome busInitializationAndBounds()
      {
        OPF::BusData<ScalarT, IdxT> data;
        data.number = 3;
        data.kv     = 230.0;

        Model::BusState state;
        state.vr = 0.8;
        state.vi = 0.6;

        BusT bus(data, state);
        bus.setVariableOffset(2);

        std::array<ScalarT, 6> values{};
        TestStatus             success  = true;
        success                        *= bus.initialize(values.data()) == 0;
        success                        *= isEqual(values[2], 1.0, 1.0e-14);
        success                        *= isEqual(values[3], std::atan2(0.6, 0.8), 1.0e-14);

        std::array<ScalarT, 6> lower{};
        std::array<ScalarT, 6> upper{};
        success *= bus.setVariableBounds(lower.data(), upper.data()) == 0;
        success *= std::isinf(lower[2]) && lower[2] < 0.0;
        success *= std::isinf(upper[2]) && upper[2] > 0.0;
        success *= std::isinf(lower[3]) && lower[3] < 0.0;
        success *= std::isinf(upper[3]) && upper[3] > 0.0;

        Model::BusState missing_state;
        BusT            missing_bus(data, missing_state);
        missing_bus.setVariableOffset(0);
        success *= missing_bus.initialize(values.data()) != 0;

        return success.report(__func__);
      }

      TestOutcome boundsAndStateOutput()
      {
        OPF::BusData<ScalarT, IdxT> data;
        data.number = 9;
        data.kv     = 115.0;
        data.vmin   = 0.9;
        data.vmax   = 1.1;

        Model::BusState initial_state;
        initial_state.vr = 0.0;
        initial_state.vi = 1.0;

        BusT bus(data, initial_state);
        bus.setVariableOffset(1);

        std::array<ScalarT, 4> lower{};
        std::array<ScalarT, 4> upper{};
        TestStatus             success  = true;
        success                        *= bus.setVariableBounds(lower.data(), upper.data()) == 0;
        success                        *= isEqual(lower[1], 0.9);
        success                        *= isEqual(upper[1], 1.1);
        success                        *= std::isinf(lower[2]) && lower[2] < 0.0;
        success                        *= std::isinf(upper[2]) && upper[2] > 0.0;

        std::array<ScalarT, 4> values{};
        values[1] = 2.0;
        values[2] = std::acos(-1.0) / 6.0;

        Model::StateData output;
        output.header.emplace();
        output.header->description    = "preserve me";
        output.devices["other"].p     = 4.0;
        auto& bus_state               = output.buses["bus_id_9"];
        bus_state.injections["L1"].ir = -0.25;
        bus_state.injections["L1"].ii = 0.5;

        success *= bus.updateSolutionState(values.data(), output) == 0;
        success *= output.buses["bus_id_9"].vr.has_value();
        success *= output.buses["bus_id_9"].vi.has_value();
        if (output.buses["bus_id_9"].vr.has_value()
            && output.buses["bus_id_9"].vi.has_value())
        {
          success *= isEqual(*output.buses["bus_id_9"].vr,
                             std::sqrt(3.0),
                             1.0e-13);
          success *= isEqual(*output.buses["bus_id_9"].vi, 1.0, 1.0e-13);
        }
        success *= output.buses["bus_id_9"].injections["L1"].ir == -0.25;
        success *= output.buses["bus_id_9"].injections["L1"].ii == 0.5;
        success *= output.header->description == "preserve me";
        success *= output.devices["other"].p == 4.0;

        return success.report(__func__);
      }

      TestOutcome generatorInitializationBoundsAndObjective()
      {
        OPF::GeneratorData<ScalarT, IdxT> data;
        data.id   = "G1";
        data.bus  = 2;
        data.pmin = -0.2;
        data.pmax = 1.4;
        data.qmin = -0.5;
        data.qmax = 0.6;
        data.mva  = 100.0;
        data.c0   = 3.0;
        data.c1   = 4.0;
        data.c2   = 5.0;

        Model::DeviceState state;
        state.p = 0.8;
        state.q = -0.1;

        GeneratorT generator(data, state);
        generator.setVariableOffset(3);

        std::array<ScalarT, 7> values{};
        TestStatus             success  = true;
        success                        *= generator.initialize(values.data()) == 0;
        success                        *= isEqual(values[3], 0.8);
        success                        *= isEqual(values[4], -0.1);

        std::array<ScalarT, 7> lower{};
        std::array<ScalarT, 7> upper{};
        success *= generator.setVariableBounds(lower.data(), upper.data()) == 0;
        success *= isEqual(lower[3], -0.2);
        success *= isEqual(upper[3], 1.4);
        success *= isEqual(lower[4], -0.5);
        success *= isEqual(upper[4], 0.6);

        ScalarT objective  = 7.0;
        success           *= generator.addObjective(values.data(), objective) == 0;
        success           *= isEqual(objective, 16.4, 1.0e-14);

        std::array<ScalarT, 7> gradient{};
        gradient.fill(1.0);
        success *= generator.addObjectiveGradient(values.data(), gradient.data()) == 0;
        success *= isEqual(gradient[3], 13.0, 1.0e-14);
        success *= isEqual(gradient[4], 1.0, 1.0e-14);

        const ScalarT step       = 1.0e-6;
        auto          plus       = values;
        auto          minus      = values;
        plus[3]                 += step;
        minus[3]                -= step;
        ScalarT objective_plus   = 0.0;
        ScalarT objective_minus  = 0.0;
        success                 *= generator.addObjective(plus.data(), objective_plus) == 0;
        success                 *= generator.addObjective(minus.data(), objective_minus) == 0;
        const ScalarT finite_difference =
            (objective_plus - objective_minus) / (2.0 * step);
        success *= isEqual(finite_difference, gradient[3] - 1.0, 1.0e-9);

        OPF::GeneratorData<ScalarT, IdxT> default_data;
        default_data.id  = "G2";
        default_data.bus = 2;
        default_data.mva = 50.0;
        Model::DeviceState default_state;
        GeneratorT         default_generator(default_data, default_state);
        default_generator.setVariableOffset(0);
        std::array<ScalarT, 2> default_values{4.0, 5.0};
        std::array<ScalarT, 2> default_lower{};
        std::array<ScalarT, 2> default_upper{};
        success *= default_generator.initialize(default_values.data()) == 0;
        success *= isEqual(default_values[0], 0.0);
        success *= isEqual(default_values[1], 0.0);
        success *= default_generator.setVariableBounds(default_lower.data(),
                                                       default_upper.data())
                   == 0;
        success *= std::isinf(default_lower[0]) && default_lower[0] < 0.0;
        success *= std::isinf(default_upper[0]) && default_upper[0] > 0.0;
        success *= std::isinf(default_lower[1]) && default_lower[1] < 0.0;
        success *= std::isinf(default_upper[1]) && default_upper[1] > 0.0;

        auto overflow_data = data;
        overflow_data.c0   = 0.0;
        overflow_data.c1   = 0.0;
        overflow_data.c2   = std::numeric_limits<ScalarT>::max();
        GeneratorT overflow_generator(overflow_data, default_state);
        overflow_generator.setVariableOffset(0);
        std::array<ScalarT, 2> overflow_values{2.0, 0.0};
        ScalarT                overflow_objective = 7.0;
        std::array<ScalarT, 2> overflow_gradient{3.0, 4.0};
        success *= overflow_generator.addObjective(
                       overflow_values.data(), overflow_objective)
                   != 0;
        success *= overflow_objective == 7.0;
        success *= overflow_generator.addObjectiveGradient(
                       overflow_values.data(), overflow_gradient.data())
                   != 0;
        success *= overflow_gradient[0] == 3.0;
        success *= overflow_gradient[1] == 4.0;

        return success.report(__func__);
      }

      TestOutcome generatorConstraintsAndJacobian()
      {
        OPF::GeneratorData<ScalarT, IdxT> data;
        data.id  = "G1";
        data.bus = 1;
        data.mva = 100.0;
        Model::DeviceState state;

        GeneratorT generator(data, state);
        generator.setVariableOffset(2);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VM>(9);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VA>(10);
        generator.constraints().template bindExternal<OPF::GeneratorExternalConstraints::DIVP>(5);
        generator.constraints().template bindExternal<OPF::GeneratorExternalConstraints::DIVQ>(8);

        std::array<ScalarT, 12> values{};
        values[2] = 1.25;
        values[3] = -0.35;
        std::array<ScalarT, 12> constraints{};
        constraints[5] = 3.0;
        constraints[8] = -2.0;

        TestStatus success  = true;
        success            *= generator.addConstraints(values.data(), constraints.data()) == 0;
        success            *= isEqual(constraints[5], 4.25);
        success            *= isEqual(constraints[8], -2.35);

        std::vector<ComponentT::JacobianEntry> entries{{99, 99}};
        generator.addJacobianSparsity(entries);
        success *= entries.size() == 3;
        if (entries.size() == 3)
        {
          success *= entries[1] == ComponentT::JacobianEntry{5, 2};
          success *= entries[2] == ComponentT::JacobianEntry{8, 3};
        }

        success *= generator.setJacobianSlots({4}) != 0;
        std::array<ScalarT, 6> jacobian_values{};
        success            *= generator.addJacobian(values.data(), jacobian_values.data()) != 0;
        success            *= generator.setJacobianSlots({4, 1}) == 0;
        jacobian_values[4]  = 2.0;
        jacobian_values[1]  = -3.0;
        success            *= generator.addJacobian(values.data(), jacobian_values.data()) == 0;
        success            *= isEqual(jacobian_values[4], 3.0);
        success            *= isEqual(jacobian_values[1], -2.0);
        success            *= isEqual(jacobian_values[0], 0.0);

        return success.report(__func__);
      }

      TestOutcome generatorOfflineBehavior()
      {
        OPF::GeneratorData<ScalarT, IdxT> data;
        data.id   = "G1";
        data.bus  = 1;
        data.pmin = -4.0;
        data.pmax = 5.0;
        data.qmin = -6.0;
        data.qmax = 7.0;
        data.mva  = 100.0;
        data.c0   = 10.0;
        data.c1   = 20.0;
        data.c2   = 30.0;

        Model::DeviceState state;
        state.online = false;
        state.p      = 1.2;
        state.q      = -0.4;

        GeneratorT generator(data, state);
        generator.setVariableOffset(0);

        std::array<ScalarT, 2> values{9.0, 8.0};
        std::array<ScalarT, 2> lower{-1.0, -1.0};
        std::array<ScalarT, 2> upper{1.0, 1.0};
        TestStatus             success  = true;
        success                        *= generator.initialize(values.data()) == 0;
        success                        *= isEqual(values[0], 0.0);
        success                        *= isEqual(values[1], 0.0);
        success                        *= generator.setVariableBounds(lower.data(), upper.data()) == 0;
        success                        *= isEqual(lower[0], 0.0);
        success                        *= isEqual(upper[0], 0.0);
        success                        *= isEqual(lower[1], 0.0);
        success                        *= isEqual(upper[1], 0.0);

        values[0]                        = 4.0;
        values[1]                        = 5.0;
        ScalarT                objective = 7.0;
        std::array<ScalarT, 2> gradient{3.0, 4.0};
        success *= generator.addObjective(values.data(), objective) == 0;
        success *= generator.addObjectiveGradient(values.data(), gradient.data()) == 0;
        success *= isEqual(objective, 7.0);
        success *= isEqual(gradient[0], 3.0);
        success *= isEqual(gradient[1], 4.0);

        values[0] = 0.0;
        values[1] = 0.0;
        Model::StateData output;
        output.buses["bus_id_1"];
        success *= generator.updateSolutionState(values.data(), output) == 0;
        success *= output.devices.at("G1").p == 0.0;
        success *= output.devices.at("G1").q == 0.0;
        success *= output.buses.at("bus_id_1").injections.at("G1").ir == 0.0;
        success *= output.buses.at("bus_id_1").injections.at("G1").ii == 0.0;

        return success.report(__func__);
      }

      TestOutcome generatorStateOutput()
      {
        OPF::GeneratorData<ScalarT, IdxT> data;
        data.id  = "G1";
        data.bus = 17;
        data.mva = 100.0;
        Model::DeviceState initial_state;

        GeneratorT generator(data, initial_state);
        generator.setVariableOffset(0);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VM>(2);
        generator.variables().template bindExternal<OPF::GeneratorExternalVariables::VA>(3);

        std::array<ScalarT, 4> values{0.8, 0.3, 1.25, 0.4};
        Model::StateData       output;
        output.header.emplace();
        output.header->description                       = "preserve me";
        output.buses["bus_id_17"].injections["other"].ir = -2.0;
        output.buses["bus_id_17"].injections["other"].ii = 3.0;
        output.devices["G1"].active                      = true;
        output.devices["G1"].online                      = true;
        output.devices["G1"].tap                         = 1.1;

        TestStatus success  = true;
        success            *= generator.updateSolutionState(values.data(), output) == 0;

        const auto& device  = output.devices["G1"];
        success            *= device.p.has_value();
        success            *= device.q.has_value();
        if (device.p.has_value() && device.q.has_value())
        {
          success *= isEqual(*device.p, 0.8);
          success *= isEqual(*device.q, 0.3);
        }
        success *= device.active == true;
        success *= device.online == true;
        success *= device.tap == 1.1;

        const ScalarT vr          = values[2] * std::cos(values[3]);
        const ScalarT vi          = values[2] * std::sin(values[3]);
        const ScalarT expected_ir = (values[0] * vr + values[1] * vi)
                                    / (values[2] * values[2]);
        const ScalarT expected_ii = (values[0] * vi - values[1] * vr)
                                    / (values[2] * values[2]);
        const auto& injection  = output.buses["bus_id_17"].injections["G1"];
        success               *= injection.ir.has_value();
        success               *= injection.ii.has_value();
        if (injection.ir.has_value() && injection.ii.has_value())
        {
          success *= isEqual(*injection.ir, expected_ir, 1.0e-14);
          success *= isEqual(*injection.ii, expected_ii, 1.0e-14);
        }
        success *= output.buses["bus_id_17"].injections["other"].ir == -2.0;
        success *= output.buses["bus_id_17"].injections["other"].ii == 3.0;
        success *= output.header->description == "preserve me";

        const Model::StateData before_failure  = output;
        values[2]                              = 0.0;
        success                               *= generator.updateSolutionState(values.data(), output) != 0;
        success                               *= output.devices["G1"].p == before_failure.devices.at("G1").p;
        success                               *= output.devices["G1"].q == before_failure.devices.at("G1").q;
        success                               *= output.buses["bus_id_17"].injections["G1"].ir
                   == before_failure.buses.at("bus_id_17").injections.at("G1").ir;
        success *= output.buses["bus_id_17"].injections["G1"].ii
                   == before_failure.buses.at("bus_id_17").injections.at("G1").ii;

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
