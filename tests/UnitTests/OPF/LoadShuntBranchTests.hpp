#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <vector>

#include <GridKit/Model/OPF/Branch/Branch.hpp>
#include <GridKit/Model/OPF/Load/Load.hpp>
#include <GridKit/Model/OPF/Shunt/Shunt.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class LoadShuntBranchTests
    {
    public:
      using ScalarT    = double;
      using RealT      = double;
      using IdxT       = std::size_t;
      using ComponentT = OPF::Component<ScalarT, IdxT>;
      using LoadT      = OPF::Load<ScalarT, IdxT>;
      using ShuntT     = OPF::Shunt<ScalarT, IdxT>;
      using BranchT    = OPF::Branch<ScalarT, IdxT>;

      TestOutcome loadBehavior()
      {
        TestStatus success = true;

        typename LoadT::DataT data;
        data.id  = "L1";
        data.bus = 7;

        Model::DeviceState state;
        state.active = false; // Load status is controlled only by online.
        state.p      = -1.75;
        state.q      = -0.8;

        LoadT       load(data, state);
        ComponentT* component = &load;
        bindLoad(load, 0, 1, 2, 3);

        success *= component->sizeInternalVariables() == 0;
        success *= component->sizeInternalConstraints() == 0;
        success *= component->nnz() == 0;
        success *= load.variables().externalIndex<OPF::LoadExternalVariables::VM>() == 0;
        success *= load.variables().externalIndex<OPF::LoadExternalVariables::VA>() == 1;
        success *= load.constraints().externalIndex<OPF::LoadExternalConstraints::DIVP>() == 2;
        success *= load.constraints().externalIndex<OPF::LoadExternalConstraints::DIVQ>() == 3;

        std::vector<RealT>   lower(4, -6.0);
        std::vector<RealT>   upper(4, 6.0);
        std::vector<ScalarT> gradient(4, 0.3);
        RealT                objective  = 2.0;
        success                        *= component->initialize(nullptr) == 0;
        success                        *= component->setVariableBounds(lower.data(), upper.data()) == 0;
        success                        *= component->setConstraintBounds(lower.data(), upper.data()) == 0;
        success                        *= component->addObjective(nullptr, objective) == 0;
        success                        *= component->addObjectiveGradient(nullptr, gradient.data()) == 0;
        success                        *= lower[0] == -6.0 && upper[0] == 6.0;
        success                        *= objective == 2.0 && gradient[0] == 0.3;

        std::vector<ScalarT> values{0.97, -0.1356993493425591};
        std::vector<ScalarT> constraints(4, 0.0);
        constraints[2]  = 0.4;
        constraints[3]  = -0.2;
        success        *= load.addConstraints(values.data(), constraints.data()) == 0;
        success        *= isEqual(constraints[2], 0.4 - 1.75, 1.0e-13);
        success        *= isEqual(constraints[3], -0.2 - 0.8, 1.0e-13);

        Model::StateData output;
        output.buses["bus_id_7"].injections["keep"].ir  = 4.0;
        output.devices["L1"]                            = state;
        success                                        *= load.updateSolutionState(values.data(), output) == 0;

        const auto& injection   = output.buses.at("bus_id_7").injections.at("L1");
        const RealT expected_ir = (*state.p * std::cos(values[1])
                                   + *state.q * std::sin(values[1]))
                                  / values[0];
        const RealT expected_ii = (*state.p * std::sin(values[1])
                                   - *state.q * std::cos(values[1]))
                                  / values[0];
        success *= isEqual(*injection.ir, expected_ir, 1.0e-13);
        success *= isEqual(*injection.ii, expected_ii, 1.0e-13);
        success *= output.devices.at("L1").p == state.p;
        success *= output.devices.at("L1").q == state.q;
        success *= output.devices.at("L1").active == state.active;
        success *= output.buses.at("bus_id_7").injections.at("keep").ir == 4.0;

        auto tiny_voltage = values;
        tiny_voltage[0]   = std::numeric_limits<RealT>::denorm_min();
        Model::StateData overflow_output;
        overflow_output.buses["bus_id_7"];
        success *= load.updateSolutionState(tiny_voltage.data(), overflow_output) != 0;
        success *= !overflow_output.buses.at("bus_id_7").injections.contains("L1");

        Model::DeviceState offline_state = state;
        offline_state.active             = true;
        offline_state.online             = false;
        LoadT offline_load(data, offline_state);
        bindLoad(offline_load, 0, 1, 2, 3);

        std::vector<ScalarT> offline_constraints(4, 0.0);
        success *= offline_load.addConstraints(values.data(), offline_constraints.data()) == 0;
        success *= offline_constraints[2] == 0.0;
        success *= offline_constraints[3] == 0.0;

        Model::StateData offline_output;
        offline_output.buses["bus_id_7"];
        offline_output.devices["L1"]   = offline_state;
        success                       *= offline_load.updateSolutionState(values.data(), offline_output) == 0;
        const auto& offline_injection  = offline_output.buses.at("bus_id_7").injections.at("L1");
        success                       *= offline_injection.ir == 0.0;
        success                       *= offline_injection.ii == 0.0;
        success                       *= offline_output.devices.at("L1").p == state.p;
        success                       *= offline_output.devices.at("L1").q == state.q;

        return success.report(__func__);
      }

      TestOutcome shuntBehavior()
      {
        TestStatus success = true;

        typename ShuntT::DataT data;
        data.id  = "S1";
        data.bus = 3;
        data.G   = 0.03;
        data.B   = 0.2;

        Model::DeviceState state;
        state.active = false; // Shunt status is controlled only by online.

        ShuntT      shunt(data, state);
        ComponentT* component = &shunt;
        bindShunt(shunt, 0, 1, 2, 3);

        success *= component->sizeInternalVariables() == 0;
        success *= component->sizeInternalConstraints() == 0;
        success *= component->nnz() == 2;
        success *= shunt.variables().externalIndex<OPF::ShuntExternalVariables::VM>() == 0;
        success *= shunt.variables().externalIndex<OPF::ShuntExternalVariables::VA>() == 1;

        std::vector<RealT> lower(4, -6.0);
        std::vector<RealT> upper(4, 6.0);
        RealT              objective  = 2.0;
        success                      *= component->initialize(nullptr) == 0;
        success                      *= component->setVariableBounds(lower.data(), upper.data()) == 0;
        success                      *= component->setConstraintBounds(lower.data(), upper.data()) == 0;
        success                      *= component->addObjective(nullptr, objective) == 0;
        success                      *= lower[0] == -6.0 && upper[0] == 6.0;
        success                      *= objective == 2.0;

        std::vector<typename ComponentT::JacobianEntry> sparsity;
        shunt.addJacobianSparsity(sparsity);
        const std::vector<typename ComponentT::JacobianEntry> expected_sparsity{
            {2, 0},
            {3, 0}};
        success *= sparsity == expected_sparsity;
        success *= shunt.setJacobianSlots({0}) != 0;
        success *= shunt.setJacobianSlots({4, 1}) == 0;

        std::vector<ScalarT> values{1.04, -0.35};
        std::vector<ScalarT> constraints(4, 0.0);
        constraints[2]  = 0.5;
        constraints[3]  = -0.25;
        success        *= shunt.addConstraints(values.data(), constraints.data()) == 0;
        success        *= isEqual(constraints[2], 0.5 - data.G * values[0] * values[0], 1.0e-13);
        success        *= isEqual(constraints[3], -0.25 + data.B * values[0] * values[0], 1.0e-13);

        std::vector<RealT> jacobian(5, 0.25);
        success *= shunt.addJacobian(values.data(), jacobian.data()) == 0;
        success *= isEqual(jacobian[4], 0.25 - 2.0 * data.G * values[0], 1.0e-13);
        success *= isEqual(jacobian[1], 0.25 + 2.0 * data.B * values[0], 1.0e-13);

        Model::StateData output;
        output.buses["bus_id_3"];
        output.devices["S1"].online  = true;
        output.devices["S1"].p       = 9.0;
        output.devices["S1"].q       = -4.0;
        success                     *= shunt.updateSolutionState(values.data(), output) == 0;

        const RealT vr         = values[0] * std::cos(values[1]);
        const RealT vi         = values[0] * std::sin(values[1]);
        const auto& injection  = output.buses.at("bus_id_3").injections.at("S1");
        success               *= isEqual(*injection.ir, -(data.G * vr - data.B * vi), 1.0e-13);
        success               *= isEqual(*injection.ii, -(data.B * vr + data.G * vi), 1.0e-13);
        success               *= output.devices.at("S1").p == 9.0;
        success               *= output.devices.at("S1").q == -4.0;

        Model::DeviceState offline_state;
        offline_state.active = true;
        offline_state.online = false;
        ShuntT offline_shunt(data, offline_state);
        bindShunt(offline_shunt, 0, 1, 2, 3);
        std::vector<typename ComponentT::JacobianEntry> offline_sparsity;
        offline_shunt.addJacobianSparsity(offline_sparsity);
        success *= offline_sparsity == expected_sparsity;
        success *= offline_shunt.nnz() == 2;
        success *= offline_shunt.setJacobianSlots({0, 1}) == 0;

        std::vector<ScalarT> offline_constraints(4, 0.0);
        std::vector<RealT>   offline_jacobian(2, 0.75);
        success *= offline_shunt.addConstraints(values.data(), offline_constraints.data()) == 0;
        success *= offline_shunt.addJacobian(values.data(), offline_jacobian.data()) == 0;
        success *= offline_constraints[2] == 0.0;
        success *= offline_constraints[3] == 0.0;
        success *= offline_jacobian[0] == 0.75;
        success *= offline_jacobian[1] == 0.75;

        std::vector<ScalarT> nonfinite_offline_values{
            std::numeric_limits<RealT>::infinity(),
            std::numeric_limits<RealT>::quiet_NaN()};
        success *= offline_shunt.addConstraints(
                       nonfinite_offline_values.data(), offline_constraints.data())
                   == 0;
        success *= offline_shunt.addJacobian(
                       nonfinite_offline_values.data(), offline_jacobian.data())
                   == 0;
        success *= offline_constraints[2] == 0.0;
        success *= offline_constraints[3] == 0.0;
        success *= offline_jacobian[0] == 0.75;
        success *= offline_jacobian[1] == 0.75;

        Model::StateData offline_output;
        offline_output.buses["bus_id_3"];
        success *= offline_shunt.updateSolutionState(values.data(), offline_output) == 0;
        success *= offline_output.buses.at("bus_id_3").injections.at("S1").ir == 0.0;
        success *= offline_output.buses.at("bus_id_3").injections.at("S1").ii == 0.0;

        auto overflow_data = data;
        overflow_data.G    = std::numeric_limits<RealT>::max();
        ShuntT overflow_shunt(overflow_data, state);
        bindShunt(overflow_shunt, 0, 1, 2, 3);
        success *= overflow_shunt.setJacobianSlots({0, 1}) == 0;
        std::vector<ScalarT> overflow_values{2.0, 0.0};
        std::vector<ScalarT> overflow_constraints(4, 0.0);
        std::vector<RealT>   overflow_jacobian(2, 0.0);
        success *= overflow_shunt.addConstraints(
                       overflow_values.data(), overflow_constraints.data())
                   != 0;
        success *= overflow_shunt.addJacobian(
                       overflow_values.data(), overflow_jacobian.data())
                   != 0;
        Model::StateData overflow_output;
        overflow_output.buses["bus_id_3"];
        success *= overflow_shunt.updateSolutionState(
                       overflow_values.data(), overflow_output)
                   != 0;
        success *= !overflow_output.buses.at("bus_id_3").injections.contains("S1");

        return success.report(__func__);
      }

      TestOutcome shuntFiniteDifferenceJacobian()
      {
        TestStatus success = true;

        typename ShuntT::DataT data;
        data.id  = "S1";
        data.bus = 3;
        data.G   = 0.07;
        data.B   = -0.16;

        Model::DeviceState state;
        ShuntT             shunt(data, state);
        bindShunt(shunt, 0, 1, 0, 1);
        success *= shunt.setJacobianSlots({0, 1}) == 0;

        std::vector<ScalarT> values{0.93, 0.27};
        std::vector<RealT>   jacobian(2, 0.4);
        success *= shunt.addJacobian(values.data(), jacobian.data()) == 0;

        const RealT h = 1.0e-6;
        for (std::size_t column = 0; column < 2; ++column)
        {
          auto plus      = values;
          auto minus     = values;
          plus[column]  += h;
          minus[column] -= h;
          std::vector<ScalarT> plus_constraints(2, 0.0);
          std::vector<ScalarT> minus_constraints(2, 0.0);
          success *= shunt.addConstraints(plus.data(), plus_constraints.data()) == 0;
          success *= shunt.addConstraints(minus.data(), minus_constraints.data()) == 0;

          for (std::size_t row = 0; row < 2; ++row)
          {
            const RealT finite_difference =
                (plus_constraints[row] - minus_constraints[row]) / (2.0 * h);
            const RealT analytic  = column == 0 ? jacobian[row] - 0.4 : 0.0;
            success              *= isEqual(analytic, finite_difference, 1.0e-8);
          }
        }

        return success.report(__func__);
      }

      TestOutcome branchBehavior()
      {
        TestStatus success = true;

        auto               data = branchData();
        Model::DeviceState state;
        state.online = false; // Branch status is controlled only by open.
        state.tap    = 1.07;
        state.phase  = 0.11;

        BranchT     branch(data, state);
        ComponentT* component = &branch;
        bindBranch(branch);

        success *= component->sizeInternalVariables() == 0;
        success *= component->sizeInternalConstraints() == 2;
        success *= component->nnz() == 24;
        success *= branch.variables().externalIndex<OPF::BranchExternalVariables::VMF>() == 0;
        success *= branch.variables().externalIndex<OPF::BranchExternalVariables::VAF>() == 1;
        success *= branch.variables().externalIndex<OPF::BranchExternalVariables::VMT>() == 2;
        success *= branch.variables().externalIndex<OPF::BranchExternalVariables::VAT>() == 3;
        success *= branch.constraints().internalIndex<OPF::BranchInternalConstraints::SF2>() == 4;
        success *= branch.constraints().internalIndex<OPF::BranchInternalConstraints::ST2>() == 5;

        std::vector<RealT> variable_lower(4, -7.0);
        std::vector<RealT> variable_upper(4, 7.0);
        RealT              objective  = 3.0;
        success                      *= component->initialize(nullptr) == 0;
        success                      *= component->setVariableBounds(
                       variable_lower.data(), variable_upper.data())
                   == 0;
        success *= component->addObjective(nullptr, objective) == 0;
        success *= variable_lower[0] == -7.0 && variable_upper[0] == 7.0;
        success *= objective == 3.0;

        const auto                                      expected_sparsity = branchSparsity(true);
        std::vector<typename ComponentT::JacobianEntry> sparsity;
        branch.addJacobianSparsity(sparsity);
        success *= sparsity == expected_sparsity;
        success *= branch.setJacobianSlots(std::vector<IdxT>(23, 0)) != 0;

        std::vector<IdxT> slots(24);
        std::iota(slots.begin(), slots.end(), IdxT{0});
        success *= branch.setJacobianSlots(slots) == 0;

        std::vector<RealT> lower(6, -9.0);
        std::vector<RealT> upper(6, 9.0);
        success *= branch.setConstraintBounds(lower.data(), upper.data()) == 0;
        success *= lower[0] == -9.0 && upper[0] == 9.0;
        success *= lower[3] == -9.0 && upper[3] == 9.0;
        success *= lower[4] == 0.0 && lower[5] == 0.0;
        success *= upper[4] == 6.25 && upper[5] == 6.25;

        std::vector<ScalarT> values{1.03, 0.08, 0.97, -0.04};
        const auto           powers = branchReference(data, state, values);
        std::vector<ScalarT> constraints(6, 0.75);
        success *= branch.addConstraints(values.data(), constraints.data()) == 0;
        for (std::size_t row = 0; row < 4; ++row)
        {
          success *= isEqual(constraints[row], 0.75 + powers[row], 1.0e-12);
        }
        success *= isEqual(constraints[4], 0.75 + powers[0] * powers[0] + powers[1] * powers[1], 1.0e-12);
        success *= isEqual(constraints[5], 0.75 + powers[2] * powers[2] + powers[3] * powers[3], 1.0e-12);

        Model::StateData output;
        output.devices[data.id].tap    = state.tap;
        output.devices[data.id].phase  = state.phase;
        output.devices[data.id].open   = false;
        success                       *= branch.updateSolutionState(values.data(), output) == 0;
        success                       *= output.devices.at(data.id).tap == state.tap;
        success                       *= output.devices.at(data.id).phase == state.phase;
        success                       *= output.devices.at(data.id).open == false;

        auto unlimited_data = data;
        unlimited_data.smax.reset();
        BranchT unlimited_branch(unlimited_data, state);
        bindBranch(unlimited_branch);
        success *= unlimited_branch.sizeInternalConstraints() == 0;
        success *= unlimited_branch.nnz() == 16;
        std::vector<typename ComponentT::JacobianEntry> unlimited_sparsity;
        unlimited_branch.addJacobianSparsity(unlimited_sparsity);
        success *= unlimited_sparsity == branchSparsity(false);
        std::vector<RealT> unlimited_lower(6, -8.0);
        std::vector<RealT> unlimited_upper(6, 8.0);
        success *= unlimited_branch.setConstraintBounds(
                       unlimited_lower.data(), unlimited_upper.data())
                   == 0;
        success *= unlimited_lower[4] == -8.0 && unlimited_upper[4] == 8.0;

        Model::DeviceState open_state = state;
        open_state.open               = true;
        BranchT open_branch(data, open_state);
        bindBranch(open_branch);
        std::vector<typename ComponentT::JacobianEntry> open_sparsity;
        open_branch.addJacobianSparsity(open_sparsity);
        success *= open_branch.nnz() == 24;
        success *= open_sparsity == expected_sparsity;
        success *= open_branch.setJacobianSlots(slots) == 0;

        std::vector<ScalarT> open_constraints(6, 0.0);
        std::vector<RealT>   open_jacobian(24, 0.6);
        success *= open_branch.addConstraints(values.data(), open_constraints.data()) == 0;
        success *= open_branch.addJacobian(values.data(), open_jacobian.data()) == 0;
        for (const ScalarT value : open_constraints)
        {
          success *= value == 0.0;
        }
        for (const RealT value : open_jacobian)
        {
          success *= value == 0.6;
        }

        Model::DeviceState default_state;
        Model::DeviceState explicit_default_state;
        explicit_default_state.open  = false;
        explicit_default_state.tap   = 1.0;
        explicit_default_state.phase = 0.0;
        BranchT default_branch(data, default_state);
        BranchT explicit_default_branch(data, explicit_default_state);
        bindBranch(default_branch);
        bindBranch(explicit_default_branch);
        success *= default_branch.setJacobianSlots(slots) == 0;
        success *= explicit_default_branch.setJacobianSlots(slots) == 0;
        std::vector<ScalarT> default_constraints(6, 0.0);
        std::vector<ScalarT> explicit_default_constraints(6, 0.0);
        std::vector<RealT>   default_jacobian(24, 0.0);
        std::vector<RealT>   explicit_default_jacobian(24, 0.0);
        success *= default_branch.addConstraints(
                       values.data(), default_constraints.data())
                   == 0;
        success *= explicit_default_branch.addConstraints(
                       values.data(), explicit_default_constraints.data())
                   == 0;
        success *= default_branch.addJacobian(
                       values.data(), default_jacobian.data())
                   == 0;
        success *= explicit_default_branch.addJacobian(
                       values.data(), explicit_default_jacobian.data())
                   == 0;
        for (std::size_t i = 0; i < default_constraints.size(); ++i)
        {
          success *= isEqual(default_constraints[i],
                             explicit_default_constraints[i],
                             1.0e-13);
        }
        for (std::size_t i = 0; i < default_jacobian.size(); ++i)
        {
          success *= isEqual(default_jacobian[i],
                             explicit_default_jacobian[i],
                             1.0e-13);
        }

        Model::DeviceState invalid_tap_state = state;
        invalid_tap_state.tap                = std::numeric_limits<RealT>::denorm_min();
        BranchT invalid_tap_branch(data, invalid_tap_state);
        bindBranch(invalid_tap_branch);
        success *= invalid_tap_branch.setJacobianSlots(slots) == 0;
        std::vector<ScalarT> invalid_tap_constraints(6, 0.0);
        std::vector<RealT>   invalid_tap_jacobian(24, 0.0);
        success *= invalid_tap_branch.addConstraints(
                       values.data(), invalid_tap_constraints.data())
                   != 0;
        success *= invalid_tap_branch.addJacobian(
                       values.data(), invalid_tap_jacobian.data())
                   != 0;

        auto overflow_data = data;
        overflow_data.R    = 1.0;
        overflow_data.X    = 0.0;
        overflow_data.G    = 0.0;
        overflow_data.B    = 0.0;
        Model::DeviceState overflow_state;
        BranchT            overflow_branch(overflow_data, overflow_state);
        bindBranch(overflow_branch);
        success *= overflow_branch.setJacobianSlots(slots) == 0;
        std::vector<ScalarT> overflow_values{1.0e103, 0.0, 0.0, 0.0};
        std::vector<ScalarT> overflow_constraints(6, 0.25);
        std::vector<RealT>   overflow_jacobian(24, 0.5);
        success *= overflow_branch.addConstraints(
                       overflow_values.data(), overflow_constraints.data())
                   != 0;
        success *= overflow_branch.addJacobian(
                       overflow_values.data(), overflow_jacobian.data())
                   != 0;
        success *= std::all_of(overflow_constraints.begin(),
                               overflow_constraints.end(),
                               [](ScalarT value)
                               { return value == 0.25; });
        success *= std::all_of(overflow_jacobian.begin(),
                               overflow_jacobian.end(),
                               [](RealT value)
                               { return value == 0.5; });

        return success.report(__func__);
      }

      TestOutcome branchFiniteDifferenceJacobian()
      {
        TestStatus success = true;

        const auto         data = branchData();
        Model::DeviceState state;
        state.tap   = 1.07;
        state.phase = 0.11;

        BranchT branch(data, state);
        bindBranch(branch);
        std::vector<IdxT> slots(24);
        std::iota(slots.begin(), slots.end(), IdxT{0});
        success *= branch.setJacobianSlots(slots) == 0;

        std::vector<ScalarT> values{1.03, 0.08, 0.97, -0.04};
        std::vector<RealT>   jacobian(24, 0.5);
        success *= branch.addJacobian(values.data(), jacobian.data()) == 0;

        const RealT h = 1.0e-6;
        for (std::size_t column = 0; column < 4; ++column)
        {
          auto plus      = values;
          auto minus     = values;
          plus[column]  += h;
          minus[column] -= h;
          std::vector<ScalarT> plus_constraints(6, 0.0);
          std::vector<ScalarT> minus_constraints(6, 0.0);
          success *= branch.addConstraints(plus.data(), plus_constraints.data()) == 0;
          success *= branch.addConstraints(minus.data(), minus_constraints.data()) == 0;

          for (std::size_t row = 0; row < 6; ++row)
          {
            const RealT finite_difference =
                (plus_constraints[row] - minus_constraints[row]) / (2.0 * h);
            const RealT analytic  = jacobian[4 * row + column] - 0.5;
            success              *= isEqual(analytic, finite_difference, 2.0e-7);
          }
        }

        return success.report(__func__);
      }

    private:
      static void bindLoad(LoadT& load,
                           IdxT   vm,
                           IdxT   va,
                           IdxT   divp,
                           IdxT   divq)
      {
        load.variables().bindExternal<OPF::LoadExternalVariables::VM>(vm);
        load.variables().bindExternal<OPF::LoadExternalVariables::VA>(va);
        load.constraints().bindExternal<OPF::LoadExternalConstraints::DIVP>(divp);
        load.constraints().bindExternal<OPF::LoadExternalConstraints::DIVQ>(divq);
      }

      static void bindShunt(ShuntT& shunt,
                            IdxT    vm,
                            IdxT    va,
                            IdxT    divp,
                            IdxT    divq)
      {
        shunt.variables().bindExternal<OPF::ShuntExternalVariables::VM>(vm);
        shunt.variables().bindExternal<OPF::ShuntExternalVariables::VA>(va);
        shunt.constraints().bindExternal<OPF::ShuntExternalConstraints::DIVP>(divp);
        shunt.constraints().bindExternal<OPF::ShuntExternalConstraints::DIVQ>(divq);
      }

      static void bindBranch(BranchT& branch)
      {
        branch.variables().bindExternal<OPF::BranchExternalVariables::VMF>(0);
        branch.variables().bindExternal<OPF::BranchExternalVariables::VAF>(1);
        branch.variables().bindExternal<OPF::BranchExternalVariables::VMT>(2);
        branch.variables().bindExternal<OPF::BranchExternalVariables::VAT>(3);
        branch.constraints().bindExternal<OPF::BranchExternalConstraints::DIVPF>(0);
        branch.constraints().bindExternal<OPF::BranchExternalConstraints::DIVQF>(1);
        branch.constraints().bindExternal<OPF::BranchExternalConstraints::DIVPT>(2);
        branch.constraints().bindExternal<OPF::BranchExternalConstraints::DIVQT>(3);
        branch.setConstraintOffset(4);
      }

      static typename BranchT::DataT branchData()
      {
        typename BranchT::DataT data;
        data.id   = "BR1";
        data.from = 1;
        data.to   = 2;
        data.R    = 0.02;
        data.X    = 0.17;
        data.G    = 0.01;
        data.B    = 0.04;
        data.smax = 2.5;
        return data;
      }

      static std::vector<typename ComponentT::JacobianEntry> branchSparsity(bool limited)
      {
        std::vector<typename ComponentT::JacobianEntry> entries;
        const std::size_t                               rows = limited ? 6 : 4;
        for (std::size_t row = 0; row < rows; ++row)
        {
          for (std::size_t column = 0; column < 4; ++column)
          {
            entries.emplace_back(row, column);
          }
        }
        return entries;
      }

      static std::array<RealT, 4> branchReference(
          const typename BranchT::DataT& data,
          const Model::DeviceState&      state,
          const std::vector<ScalarT>&    values)
      {
        if (state.open.value_or(false))
        {
          return {};
        }

        using Complex       = std::complex<RealT>;
        const RealT   tap   = state.tap.value_or(1.0);
        const RealT   phase = state.phase.value_or(0.0);
        const Complex y     = Complex{1.0, 0.0} / Complex{data.R, data.X};
        const Complex y_shunt{data.G, data.B};
        const Complex rotation = std::polar(1.0, phase);

        const Complex yff = -(y + 0.5 * y_shunt) / (tap * tap);
        const Complex yft = y * rotation / tap;
        const Complex ytf = y * std::conj(rotation) / tap;
        const Complex ytt = -(y + 0.5 * y_shunt);

        const Complex vf = std::polar(values[0], values[1]);
        const Complex vt = std::polar(values[2], values[3]);
        const Complex sf = vf * std::conj(yff * vf + yft * vt);
        const Complex st = vt * std::conj(ytf * vf + ytt * vt);
        return {sf.real(), sf.imag(), st.real(), st.imag()};
      }
    };

  } // namespace Testing
} // namespace GridKit
