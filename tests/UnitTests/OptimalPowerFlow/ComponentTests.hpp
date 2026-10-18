#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

#include "TestHelpers.hpp"

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    /**
     * @brief Component derivatives against dependency tracking and closed forms
     *
     * Each device is composed ahead of its buses, so its local variables have
     * their global indices.
     */
    class ComponentTests
    {
      using RealT         = double;
      using IdxT          = size_t;
      using VariableT     = DependencyTracking::Variable;
      using DependencyMap = VariableT::DependencyMap;
      using CooBufferT    = OPF::Component<RealT, IdxT>::CooBuffer;
      using BranchT       = OPF::Branch<RealT, IdxT>;

    public:
      // Branch Hessian rows reach 5.2 eps, and 8.2 eps on x86-64-v4 with ASan.
      static constexpr RealT kTol =
          static_cast<RealT>(10.0) * std::numeric_limits<RealT>::epsilon();

      /// Power into each bus at V1 = 1 and V2 = j, from the admittances by hand
      TestOutcome branchPower()
      {
        TestStatus success = true;

        auto branch_owner  = std::make_unique<BranchT>(branchData());
        auto system        = branchSystem(std::move(branch_owner));
        success           *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, POWER_X);
        success *= system.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT* g = system.g().getData();

        // S1 = -3.9 - j3.45 and S2 = -7.6 - j7.3
        success *= isEqual(g[BRANCH_BALANCE_ROWS[0]], -3.9, kTol);
        success *= isEqual(g[BRANCH_BALANCE_ROWS[1]], -3.45, kTol);
        success *= isEqual(g[BRANCH_BALANCE_ROWS[2]], -7.6, kTol);
        success *= isEqual(g[BRANCH_BALANCE_ROWS[3]], -7.3, kTol);
        success *= isEqual(g[BRANCH_FLOW_ROWS[0]], 3.9 * 3.9 + 3.45 * 3.45, kTol);
        success *= isEqual(g[BRANCH_FLOW_ROWS[1]], 7.6 * 7.6 + 7.3 * 7.3, kTol);

        return success.report(__func__);
      }

      /// Enzyme Jacobian entries equal the dependency-tracking derivatives
      TestOutcome branchJacobian()
      {
        TestStatus success = true;

        auto  branch_owner  = std::make_unique<BranchT>(branchData());
        auto& branch        = *branch_owner;
        auto  system        = branchSystem(std::move(branch_owner));
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, BRANCH_X);
        success *= branch.evaluateJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        auto tracked  = branchSystem(std::make_unique<OPF::Branch<VariableT, IdxT>>(branchData()));
        success      *= tracked.allocate() == 0 && tracked.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(tracked, BRANCH_X);
        success *= tracked.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        const VariableT* g = tracked.g().getData();

        const auto enzyme = rowsOf(sum(branch.jacobianEntries()));
        for (const IdxT row : BRANCH_FLOW_ROWS)
        {
          success *= isEqual(enzyme.at(row), g[row].getDependencies(), kTol);
        }
        for (const IdxT row : BRANCH_BALANCE_ROWS)
        {
          success *= isEqual(enzyme.at(row), g[row].getDependencies(), kTol);
        }
        success *= enzyme.size() == BRANCH_FLOW_ROWS.size() + BRANCH_BALANCE_ROWS.size();
        success *= branch.jacobianEntries().nnz == BRANCH_JACOBIAN_NNZ;

        return success.report(__func__);
      }

      /// Enzyme Hessian entries equal the closed-form Cartesian second derivatives
      TestOutcome branchHessian()
      {
        TestStatus success = true;

        auto  branch_owner  = std::make_unique<BranchT>(branchData());
        auto& branch        = *branch_owner;
        auto  system        = branchSystem(std::move(branch_owner));
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, BRANCH_X);
        setMultipliers(system, BRANCH_LAMBDA);
        success *= branch.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        // Balances: P1, Q1, P2, Q2
        EntryMap expected;
        addPowerHessian(expected, 0, BUS1_ADMITTANCE, BRANCH_LAMBDA[BRANCH_BALANCE_ROWS[0]], BRANCH_LAMBDA[BRANCH_BALANCE_ROWS[1]]);
        addPowerHessian(expected, 1, BUS2_ADMITTANCE, BRANCH_LAMBDA[BRANCH_BALANCE_ROWS[2]], BRANCH_LAMBDA[BRANCH_BALANCE_ROWS[3]]);

        // Flows: |S|^2 = P^2 + Q^2 by the chain rule over tracked gradients
        auto tracked  = branchSystem(std::make_unique<OPF::Branch<VariableT, IdxT>>(branchData()));
        success      *= tracked.allocate() == 0 && tracked.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(tracked, BRANCH_X);
        success *= tracked.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        const VariableT* g = tracked.g().getData();
        for (IdxT port = 0; port < 2; ++port)
        {
          const RealT     weight = TWO<RealT> * BRANCH_LAMBDA[BRANCH_FLOW_ROWS[port]];
          const VariableT p      = g[BRANCH_BALANCE_ROWS[2 * port]];
          const VariableT q      = g[BRANCH_BALANCE_ROWS[1 + 2 * port]];
          addOuterProduct(expected, p.getDependencies(), weight);
          addOuterProduct(expected, q.getDependencies(), weight);
        }
        const RealT weight1 = TWO<RealT> * BRANCH_LAMBDA[BRANCH_FLOW_ROWS[0]];
        const RealT weight2 = TWO<RealT> * BRANCH_LAMBDA[BRANCH_FLOW_ROWS[1]];
        addPowerHessian(expected, 0, BUS1_ADMITTANCE, weight1 * g[BRANCH_BALANCE_ROWS[0]].getValue(), weight1 * g[BRANCH_BALANCE_ROWS[1]].getValue());
        addPowerHessian(expected, 1, BUS2_ADMITTANCE, weight2 * g[BRANCH_BALANCE_ROWS[2]].getValue(), weight2 * g[BRANCH_BALANCE_ROWS[3]].getValue());

        success *= matches(sum(branch.hessianEntries()), expected, kTol);
        success *= branch.hessianEntries().nnz == BRANCH_HESSIAN_NNZ;

        return success.report(__func__);
      }

      /// Entry counts are structural, so an out-of-service branch at a flat
      /// start with zero multipliers has the same entries, all zero
      TestOutcome branchPattern()
      {
        TestStatus success = true;

        auto  active_owner   = std::make_unique<BranchT>(branchData());
        auto& active_branch  = *active_owner;
        auto  active         = branchSystem(std::move(active_owner));
        success             *= active.allocate() == 0 && active.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(active, BRANCH_X);
        setMultipliers(active, BRANCH_LAMBDA);
        success *= active_branch.evaluateJacobian() == 0;
        success *= active_branch.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        auto data                                          = branchData();
        data.parameters[OPF::BranchParameters::in_service] = false;

        auto  branch_owner  = std::make_unique<BranchT>(data);
        auto& branch        = *branch_owner;
        auto  system        = branchSystem(std::move(branch_owner));
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        // The multipliers start at zero
        setVariables(system, FLAT_START);

        success *= branch.evaluateJacobian() == 0;
        success *= branch.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        success *= pattern(branch.jacobianEntries()) == pattern(active_branch.jacobianEntries());
        success *= pattern(branch.hessianEntries()) == pattern(active_branch.hessianEntries());
        success *= branch.jacobianEntries().nnz == BRANCH_JACOBIAN_NNZ;
        success *= branch.hessianEntries().nnz == BRANCH_HESSIAN_NNZ;
        for (const RealT value : values(branch.jacobianEntries()))
        {
          success *= isEqual(value, ZERO<RealT>);
        }
        for (const RealT value : values(branch.hessianEntries()))
        {
          success *= isEqual(value, ZERO<RealT>);
        }

        return success.report(__func__);
      }

      /// Generator cost, bounds, and derivatives in closed form
      TestOutcome generator()
      {
        TestStatus success = true;

        auto  bus_owner        = std::make_unique<OPF::Bus<RealT, IdxT>>(busData(1));
        auto& bus              = *bus_owner;
        auto  generator_owner  = std::make_unique<OPF::Generator<RealT, IdxT>>(generatorData());
        auto& generator        = *generator_owner;
        success               *= generator.connect(0, bus) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        auto system  = makeSystem<RealT>(std::move(generator_owner), std::move(bus_owner));
        success     *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* x = system.x().getData();

        setVariables(system, GENERATOR_X);
        success *= system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(x[0], ZERO<RealT>) && isEqual(x[1], ZERO<RealT>);

        const RealT* lower = system.xLower().getData();
        const RealT* upper = system.xUpper().getData();

        success *= isEqual(lower[0], PMIN) && isEqual(upper[0], PMAX);
        success *= isEqual(lower[1], QMIN) && isEqual(upper[1], QMAX);

        setVariables(system, GENERATOR_X);
        const RealT p = GENERATOR_X[0];

        success *= generator.evaluateObjective() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(generator.f(), C0 + C1 * p + C2 * p * p);

        success *= system.evaluateGradient() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(system.gradient().getData()[0], C1 + TWO<RealT> * C2 * p);

        success *= generator.evaluateJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= matches(sum(generator.jacobianEntries()), EntryMap{{{0, 0}, ONE<RealT>}, {{1, 1}, ONE<RealT>}}, kTol);

        const std::vector<RealT> lambda(system.sizeConstraints(), ONE<RealT>);
        setMultipliers(system, lambda);
        success *= generator.evaluateHessian(SIGMA) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= matches(sum(generator.hessianEntries()), EntryMap{{{0, 0}, TWO<RealT> * C2 * SIGMA}}, kTol);

        return success.report(__func__);
      }

      /// An out-of-service generator is fixed at zero with no cost and the same entries
      TestOutcome outOfServiceGenerator()
      {
        TestStatus success = true;

        auto data                                             = generatorData();
        data.parameters[OPF::GeneratorParameters::in_service] = false;

        auto  bus_owner        = std::make_unique<OPF::Bus<RealT, IdxT>>(busData(1));
        auto& bus              = *bus_owner;
        auto  generator_owner  = std::make_unique<OPF::Generator<RealT, IdxT>>(data);
        auto& generator        = *generator_owner;
        success               *= generator.connect(0, bus) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        auto system  = makeSystem<RealT>(std::move(generator_owner), std::move(bus_owner));
        success     *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT* lower = system.xLower().getData();
        const RealT* upper = system.xUpper().getData();

        success *= isEqual(lower[0], ZERO<RealT>) && isEqual(upper[0], ZERO<RealT>);
        success *= isEqual(lower[1], ZERO<RealT>) && isEqual(upper[1], ZERO<RealT>);

        setVariables(system, GENERATOR_X);
        success *= generator.evaluateObjective() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(generator.f(), ZERO<RealT>);

        success *= generator.evaluateJacobian() == 0;
        success *= generator.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= generator.jacobianEntries().nnz == 2;
        success *= generator.hessianEntries().nnz == 1;

        return success.report(__func__);
      }

      /// The bus shunt adds its power to the balance, with derivatives in closed form
      TestOutcome busShunt()
      {
        TestStatus success = true;

        auto data                              = busData(1);
        data.parameters[OPF::BusParameters::G] = SHUNT_G;
        data.parameters[OPF::BusParameters::B] = SHUNT_B;

        auto  bus_owner  = std::make_unique<OPF::Bus<RealT, IdxT>>(data);
        auto& bus        = *bus_owner;
        auto  system     = makeSystem<RealT>(std::move(bus_owner));
        success         *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* g = system.g().getData();

        setVariables(system, SHUNT_X);
        const RealT vr = SHUNT_X[0];
        const RealT vi = SHUNT_X[1];

        success *= system.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(g[0], -SHUNT_G * (vr * vr + vi * vi));
        success *= isEqual(g[1], SHUNT_B * (vr * vr + vi * vi));

        success *= bus.evaluateJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= matches(sum(bus.jacobianEntries()),
                           EntryMap{{{0, 0}, -TWO<RealT> * SHUNT_G * vr},
                                    {{0, 1}, -TWO<RealT> * SHUNT_G * vi},
                                    {{1, 0}, TWO<RealT> * SHUNT_B * vr},
                                    {{1, 1}, TWO<RealT> * SHUNT_B * vi},
                                    {{2, 0}, TWO<RealT> * vr},
                                    {{2, 1}, TWO<RealT> * vi}},
                           kTol);

        setMultipliers(system, SHUNT_LAMBDA);
        success *= bus.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT hessian = TWO<RealT> * (-SHUNT_G * SHUNT_LAMBDA[0] + SHUNT_B * SHUNT_LAMBDA[1]);

        success *= matches(sum(bus.hessianEntries()), EntryMap{{{0, 0}, hessian}, {{1, 1}, hessian}}, kTol);

        return success.report(__func__);
      }

      /// Devices supply the bus balance, the squared magnitude has the limits,
      /// and a reference bus fixes its imaginary voltage at zero
      TestOutcome bus()
      {
        TestStatus success = true;

        auto data                                 = busData(3);
        data.parameters[OPF::BusParameters::Vmin] = VMIN;
        data.parameters[OPF::BusParameters::Vmax] = VMAX;

        auto  bus_owner  = std::make_unique<OPF::Bus<RealT, IdxT>>(data);
        auto& bus        = *bus_owner;
        auto  system     = makeSystem<RealT>(std::move(bus_owner));
        success         *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* x = system.x().getData();

        success *= system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(x[0], ONE<RealT>) && isEqual(x[1], ZERO<RealT>);
        for (IdxT j = 0; j < system.size(); ++j)
        {
          success *= std::isinf(system.xLower().getData()[j]) && std::isinf(system.xUpper().getData()[j]);
        }

        // The balance is an equality at zero, and the squared magnitude an inequality within the limits
        const RealT* g_lower  = system.gLower().getData();
        const RealT* g_upper  = system.gUpper().getData();
        success              *= system.sizeConstraints() == 3;
        success              *= system.tag() == std::vector<bool>{false, false, true};
        success              *= isEqual(g_lower[0], ZERO<RealT>) && isEqual(g_upper[0], ZERO<RealT>);
        success              *= isEqual(g_lower[1], ZERO<RealT>) && isEqual(g_upper[1], ZERO<RealT>);
        success              *= isEqual(g_lower[2], VMIN * VMIN) && isEqual(g_upper[2], VMAX * VMAX);

        setVariables(system, BUS_X);
        const RealT vr = BUS_X[0];
        const RealT vi = BUS_X[1];

        success *= system.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* g  = system.g().getData();
        success        *= isEqual(g[0], ZERO<RealT>) && isEqual(g[1], ZERO<RealT>);
        success        *= isEqual(g[2], vr * vr + vi * vi);

        success *= bus.evaluateJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        // A bus without a shunt contributes nothing to its balance
        success *= matches(sum(bus.jacobianEntries()),
                           EntryMap{{{0, 0}, ZERO<RealT>},
                                    {{0, 1}, ZERO<RealT>},
                                    {{1, 0}, ZERO<RealT>},
                                    {{1, 1}, ZERO<RealT>},
                                    {{2, 0}, TWO<RealT> * vr},
                                    {{2, 1}, TWO<RealT> * vi}},
                           kTol);

        setMultipliers(system, BUS_LAMBDA);
        success *= bus.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT hessian = TWO<RealT> * BUS_LAMBDA[2];

        success *= matches(sum(bus.hessianEntries()), EntryMap{{{0, 0}, hessian}, {{1, 1}, hessian}}, kTol);

        data.parameters[OPF::BusParameters::reference]  = true;
        auto reference_owner                            = std::make_unique<OPF::Bus<RealT, IdxT>>(data);
        auto reference_system                           = makeSystem<RealT>(std::move(reference_owner));
        success                                        *= reference_system.allocate() == 0 && reference_system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= std::isinf(reference_system.xLower().getData()[0]) && std::isinf(reference_system.xUpper().getData()[0]);
        success *= isEqual(reference_system.xLower().getData()[1], ZERO<RealT>) && isEqual(reference_system.xUpper().getData()[1], ZERO<RealT>);

        return success.report(__func__);
      }

      /// Constant-power demand has zero derivatives and is zero out of service
      TestOutcome loadConstantPower()
      {
        TestStatus success = true;

        OPF::LoadData<RealT, IdxT> data;
        data.id                                 = "load_2";
        data.ports[OPF::LoadPorts::bus]         = 2;
        data.parameters[OPF::LoadParameters::P] = LOAD_P;
        data.parameters[OPF::LoadParameters::Q] = LOAD_Q;

        auto  bus_owner                                   = std::make_unique<OPF::Bus<RealT, IdxT>>(busData(2));
        auto& bus                                         = *bus_owner;
        auto  load_owner                                  = std::make_unique<OPF::Load<RealT, IdxT>>(data);
        auto& load                                        = *load_owner;
        data.id                                           = "load_2_out_of_service";
        data.parameters[OPF::LoadParameters::in_service]  = false;
        auto  out_of_service_owner                        = std::make_unique<OPF::Load<RealT, IdxT>>(data);
        auto& out_of_service                              = *out_of_service_owner;
        success                                          *= load.connect(0, bus) == 0;
        success                                          *= out_of_service.connect(0, bus) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        auto system  = makeSystem<RealT>(std::move(load_owner), std::move(out_of_service_owner), std::move(bus_owner));
        success     *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* g = system.g().getData();

        setVariables(system, LOAD_X);

        // The out-of-service load adds no demand
        success *= system.evaluateConstraints() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= isEqual(g[0], -LOAD_P) && isEqual(g[1], -LOAD_Q);

        const std::vector<RealT> lambda(system.sizeConstraints(), ONE<RealT>);
        setMultipliers(system, lambda);
        success *= load.evaluateJacobian() == 0;
        success *= load.evaluateHessian(ONE<RealT>) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        for (const RealT value : values(load.jacobianEntries()))
        {
          success *= isEqual(value, ZERO<RealT>);
        }
        for (const RealT value : values(load.hessianEntries()))
        {
          success *= isEqual(value, ZERO<RealT>);
        }

        return success.report(__func__);
      }

      /// ZIP demand at each voltage, in closed form, and none out of service
      TestOutcome loadEqualities()
      {
        TestStatus success = true;
        for (const bool in_service : {true, false})
        {
          for (const ZipParts& zip : ZIP_PARTS)
          {
            auto system  = loadSystem<RealT>(loadData(in_service, zip));
            success     *= system.allocate() == 0 && system.initialize() == 0;
            if (!success)
            {
              return success.report(__func__);
            }
            for (const auto& voltage : ZIP_VOLTAGES)
            {
              setVariables(system, voltage);
              success *= system.evaluateConstraints() == 0;
              if (!success)
              {
                return success.report(__func__);
              }
              const RealT v2 = voltage[0] * voltage[0] + voltage[1] * voltage[1];
              const RealT v  = std::sqrt(v2);
              const RealT u  = status(in_service);

              success *= isEqual(system.g().getData()[0], -u * (LOAD_P + zip.Pi * v + zip.Pz * v2), kTol);
              success *= isEqual(system.g().getData()[1], -u * (LOAD_Q + zip.Qi * v + zip.Qz * v2), kTol);
            }
          }
        }
        return success.report(__func__);
      }

      /// Enzyme Jacobian entries equal the dependency-tracking and closed-form derivatives
      TestOutcome loadJacobian()
      {
        TestStatus success = true;
        for (const bool in_service : {true, false})
        {
          for (const ZipParts& zip : ZIP_PARTS)
          {
            const auto data     = loadData(in_service, zip);
            auto       system   = loadSystem<RealT>(data);
            auto       tracked  = loadSystem<VariableT>(data);
            success            *= system.allocate() == 0 && system.initialize() == 0;
            success            *= tracked.allocate() == 0 && tracked.initialize() == 0;
            if (!success)
            {
              return success.report(__func__);
            }
            for (const auto& voltage : ZIP_VOLTAGES)
            {
              setVariables(system, voltage);
              setVariables(tracked, voltage);
              success *= system.evaluateJacobian() == 0;
              success *= tracked.evaluateJacobian() == 0;
              if (!success)
              {
                return success.report(__func__);
              }
              const auto  actual   = MapFromCsr(system.getCsrJacobian());
              const auto  expected = MapFromCsr(tracked.getCsrJacobian());
              const RealT v        = std::hypot(voltage[0], voltage[1]);
              const RealT u        = status(in_service);

              // Derivative of the demand over each voltage component, divided by that component
              const std::array<RealT, 2> slope{u * (zip.Pi / v + TWO<RealT> * zip.Pz),
                                               u * (zip.Qi / v + TWO<RealT> * zip.Qz)};
              for (IdxT row = 0; row < slope.size(); ++row)
              {
                success *= isEqual(actual.at(row), expected.at(row), kTol);
                success *= isEqual(actual.at(row), DependencyMap{{0, -slope[row] * voltage[0]}, {1, -slope[row] * voltage[1]}}, kTol);
              }
            }
          }
        }
        return success.report(__func__);
      }

      /// Enzyme Hessian entries equal the closed-form second derivatives
      TestOutcome loadHessian()
      {
        TestStatus success = true;
        for (const bool in_service : {true, false})
        {
          for (const ZipParts& zip : ZIP_PARTS)
          {
            auto system  = loadSystem<RealT>(loadData(in_service, zip));
            success     *= system.allocate() == 0 && system.initialize() == 0;
            if (!success)
            {
              return success.report(__func__);
            }
            const auto [current, impedance] = loadWeights(in_service, zip);
            for (const auto& voltage : ZIP_VOLTAGES)
            {
              setVariables(system, voltage);
              setMultipliers(system, LOAD_LAMBDA);
              success *= system.evaluateHessian(ONE<RealT>) == 0;
              if (!success)
              {
                return success.report(__func__);
              }
              const RealT vr = voltage[0];
              const RealT vi = voltage[1];
              const RealT v2 = vr * vr + vi * vi;
              const RealT v3 = v2 * std::sqrt(v2);

              const DependencyMap row0{{0, current * vi * vi / v3 + TWO<RealT> * impedance}};
              const DependencyMap row1{{0, -current * vr * vi / v3}, {1, current * vr * vr / v3 + TWO<RealT> * impedance}};

              const auto actual  = MapFromCsr(system.getCsrHessian());
              success           *= isEqual(actual.at(0), row0, kTol);
              success           *= isEqual(actual.at(1), row1, kTol);
            }
          }
        }
        return success.report(__func__);
      }

      /// Demands are finite reals, integers count as reals, and the service status is a boolean
      TestOutcome loadConfiguration()
      {
        TestStatus                 success = true;
        OPF::LoadData<RealT, IdxT> data;
        success *= OPF::Load<RealT, IdxT>(data).verify().passed();

        data.parameters  = {{OPF::LoadParameters::P, IdxT{1}}};
        success         *= OPF::Load<RealT, IdxT>(data).verify().passed();

        for (const auto parameter : {OPF::LoadParameters::P,
                                     OPF::LoadParameters::Q,
                                     OPF::LoadParameters::Pi,
                                     OPF::LoadParameters::Qi,
                                     OPF::LoadParameters::Pz,
                                     OPF::LoadParameters::Qz})
        {
          data.parameters  = {{parameter, std::numeric_limits<RealT>::infinity()}};
          success         *= configurationRejected([&]()
                                           { return OPF::Load<RealT, IdxT>(data).verify(); });
          data.parameters  = {{parameter, true}};
          success         *= configurationRejected([&]()
                                           { return OPF::Load<RealT, IdxT>(data).verify(); });
        }

        data.parameters  = {{OPF::LoadParameters::in_service, IdxT{1}}};
        success         *= configurationRejected([&]()
                                         { return OPF::Load<RealT, IdxT>(data).verify(); });

        return success.report(__func__);
      }

    private:
      /// Constant-current and constant-impedance demands of a ZIP load
      struct ZipParts
      {
        RealT Pi;
        RealT Qi;
        RealT Pz;
        RealT Qz;
      };

      /// Constant power only, constant current, constant impedance, and all three parts
      static constexpr std::array<ZipParts, 4> ZIP_PARTS{{{0.0, 0.0, 0.0, 0.0},
                                                          {0.5, -0.2, 0.0, 0.0},
                                                          {0.0, 0.0, 0.3, 0.1},
                                                          {0.4, 0.15, 0.25, -0.05}}};

      static constexpr std::array<std::array<RealT, 2>, 3> ZIP_VOLTAGES{{{1.1, 0.0}, {0.8, -0.6}, {1.15, 0.2}}};
      static constexpr std::array<RealT, 3>                LOAD_LAMBDA{0.7, -0.4, 0.0};

      static OPF::LoadData<RealT, IdxT> loadData(bool in_service, const ZipParts& zip)
      {
        OPF::LoadData<RealT, IdxT> data;
        data.id         = "zip";
        data.parameters = {{OPF::LoadParameters::P, LOAD_P},
                           {OPF::LoadParameters::Q, LOAD_Q},
                           {OPF::LoadParameters::Pi, zip.Pi},
                           {OPF::LoadParameters::Qi, zip.Qi},
                           {OPF::LoadParameters::Pz, zip.Pz},
                           {OPF::LoadParameters::Qz, zip.Qz},
                           {OPF::LoadParameters::in_service, in_service}};
        return data;
      }

      /// Service status as the factor the models apply
      static RealT status(bool in_service)
      {
        if (in_service)
        {
          return ONE<RealT>;
        }
        return ZERO<RealT>;
      }

      /// Multiplier-weighted constant-current and constant-impedance demands, signed as power into the bus
      static std::pair<RealT, RealT> loadWeights(bool in_service, const ZipParts& zip)
      {
        const RealT u = status(in_service);
        return {-u * (LOAD_LAMBDA[0] * zip.Pi + LOAD_LAMBDA[1] * zip.Qi),
                -u * (LOAD_LAMBDA[0] * zip.Pz + LOAD_LAMBDA[1] * zip.Qz)};
      }

      template <typename ScalarT>
      static OPF::SystemModel<ScalarT, IdxT> loadSystem(const OPF::LoadData<RealT, IdxT>& data)
      {
        auto bus  = std::make_unique<OPF::Bus<ScalarT, IdxT>>(busData(1));
        auto load = std::make_unique<OPF::Load<ScalarT, IdxT>>(data);
        if (load->connect(0, *bus) != 0)
        {
          throw std::invalid_argument("Invalid load test connection");
        }
        return makeSystem<ScalarT>(std::move(load), std::move(bus));
      }

      template <typename ScalarT>
      static OPF::SystemModel<ScalarT, IdxT>
      branchSystem(std::unique_ptr<OPF::Branch<ScalarT, IdxT>> branch)
      {
        auto bus1 = std::make_unique<OPF::Bus<ScalarT, IdxT>>(busData(1));
        auto bus2 = std::make_unique<OPF::Bus<ScalarT, IdxT>>(busData(2));
        if (branch->connect(0, *bus1) != 0 || branch->connect(1, *bus2) != 0)
        {
          throw std::invalid_argument("Invalid branch test connections");
        }
        return makeSystem<ScalarT>(std::move(branch), std::move(bus1), std::move(bus2));
      }

      static constexpr const char* BRANCH_ID = "branch_1_2";

      static constexpr RealT R     = 0.1;
      static constexpr RealT X     = 0.1;
      static constexpr RealT G     = 0.2;
      static constexpr RealT B     = 0.4;
      static constexpr RealT GMAG  = 0.05;
      static constexpr RealT BMAG  = 0.1;
      static constexpr RealT TAP   = 2.0;
      static constexpr RealT PHASE = std::numbers::pi / 2.0;
      static constexpr RealT SMAX  = 1.2;

      /// Admittances g_k1, b_k1, g_k2, b_k2 of the current into bus k: series
      /// admittance 5 - j5, tap 2 at 90 degrees, and the shunts
      static constexpr std::array<RealT, 4> BUS1_ADMITTANCE = {-1.4, 0.95, 2.5, 2.5};
      static constexpr std::array<RealT, 4> BUS2_ADMITTANCE = {-2.5, -2.5, -5.1, 4.8};

      /// Constraints of a branch system: S1 and S2, then P, Q, and VM2 of each bus
      static constexpr std::array<IdxT, 2> BRANCH_FLOW_ROWS{0, 1};
      static constexpr std::array<IdxT, 4> BRANCH_BALANCE_ROWS{2, 3, 5, 6};

      static constexpr IdxT BRANCH_JACOBIAN_NNZ = 24;
      static constexpr IdxT BRANCH_HESSIAN_NNZ  = 10;

      static constexpr RealT PMIN  = 0.1;
      static constexpr RealT PMAX  = 0.9;
      static constexpr RealT QMIN  = -0.4;
      static constexpr RealT QMAX  = 0.4;
      static constexpr RealT C0    = 3.0;
      static constexpr RealT C1    = 11.0;
      static constexpr RealT C2    = 0.7;
      static constexpr RealT SIGMA = 1.7;

      static constexpr RealT SHUNT_G = 0.05;
      static constexpr RealT SHUNT_B = 0.3;

      static constexpr RealT VMIN = 0.9;
      static constexpr RealT VMAX = 1.1;

      static constexpr RealT LOAD_P = 0.6;
      static constexpr RealT LOAD_Q = 0.2;

      static constexpr std::array<RealT, 4> POWER_X{1.0, 0.0, 0.0, 1.0};
      static constexpr std::array<RealT, 4> BRANCH_X{1.02, 0.12, 0.96, -0.08};
      static constexpr std::array<RealT, 8> BRANCH_LAMBDA{0.3, -0.2, 1.5, -0.7, 0.0, 2.1, 0.4, 0.0};
      static constexpr std::array<RealT, 4> FLAT_START{1.0, 0.0, 1.0, 0.0};
      static constexpr std::array<RealT, 4> GENERATOR_X{0.3, 0.1, 1.01, 0.05};
      static constexpr std::array<RealT, 2> SHUNT_X{0.98, 0.1};
      static constexpr std::array<RealT, 3> SHUNT_LAMBDA{1.3, -0.6, 0.0};
      static constexpr std::array<RealT, 2> BUS_X{1.02, -0.05};
      static constexpr std::array<RealT, 3> BUS_LAMBDA{0.0, 0.0, 1.3};
      static constexpr std::array<RealT, 2> LOAD_X{1.01, -0.05};

      static OPF::BusData<RealT, IdxT> busData(IdxT number)
      {
        OPF::BusData<RealT, IdxT> data;
        data.number = number;
        return data;
      }

      static OPF::BranchData<RealT, IdxT> branchData()
      {
        OPF::BranchData<RealT, IdxT> data;
        data.id                                       = BRANCH_ID;
        data.ports[OPF::BranchPorts::bus1]            = 1;
        data.ports[OPF::BranchPorts::bus2]            = 2;
        data.parameters[OPF::BranchParameters::R]     = R;
        data.parameters[OPF::BranchParameters::X]     = X;
        data.parameters[OPF::BranchParameters::G]     = G;
        data.parameters[OPF::BranchParameters::B]     = B;
        data.parameters[OPF::BranchParameters::Gmag]  = GMAG;
        data.parameters[OPF::BranchParameters::Bmag]  = BMAG;
        data.parameters[OPF::BranchParameters::tap]   = TAP;
        data.parameters[OPF::BranchParameters::phase] = PHASE;
        data.parameters[OPF::BranchParameters::Smax]  = SMAX;
        return data;
      }

      static OPF::GeneratorData<RealT, IdxT> generatorData()
      {
        OPF::GeneratorData<RealT, IdxT> data;
        data.id                                         = "gen_1";
        data.ports[OPF::GeneratorPorts::bus]            = 1;
        data.parameters[OPF::GeneratorParameters::Pmin] = PMIN;
        data.parameters[OPF::GeneratorParameters::Pmax] = PMAX;
        data.parameters[OPF::GeneratorParameters::Qmin] = QMIN;
        data.parameters[OPF::GeneratorParameters::Qmax] = QMAX;
        data.parameters[OPF::GeneratorParameters::c0]   = C0;
        data.parameters[OPF::GeneratorParameters::c1]   = C1;
        data.parameters[OPF::GeneratorParameters::c2]   = C2;
        return data;
      }

      /// Stored values
      static std::vector<RealT> values(const CooBufferT& entries)
      {
        return {entries.values.getData(), entries.values.getData() + entries.nnz};
      }

      /// Entries of each row
      static std::map<IdxT, DependencyMap> rowsOf(const EntryMap& entries)
      {
        std::map<IdxT, DependencyMap> rows;
        for (const auto& [coordinate, value] : entries)
        {
          rows[coordinate.first][coordinate.second] = value;
        }
        return rows;
      }

      /// Add `value` at the lower-triangle coordinate of (i, j)
      static void addEntry(EntryMap& hessian, IdxT i, IdxT j, RealT value)
      {
        hessian[{std::max(i, j), std::min(i, j)}] += value;
      }

      /// Add `weight` times the lower triangle of the outer product of `gradient`
      static void addOuterProduct(EntryMap& hessian, const DependencyMap& gradient, RealT weight)
      {
        for (const auto& [i, di] : gradient)
        {
          for (const auto& [j, dj] : gradient)
          {
            if (i >= j)
            {
              hessian[{i, j}] += weight * di * dj;
            }
          }
        }
      }

      /**
       * @brief Add the Hessian of \f$w_P P_k + w_Q Q_k\f$ for the current
       * \f$I_k = (g_{k1} + j b_{k1}) V_1 + (g_{k2} + j b_{k2}) V_2\f$
       *
       * \f$P_k + j Q_k = V_k I_k^*\f$ is bilinear in the voltages, so each
       * admittance adds a constant block. Local indices \f$2m\f$ and
       * \f$2m + 1\f$ hold the real and imaginary voltage of bus \f$m + 1\f$.
       *
       * @param[in] admittance - \f$g_{k1}, b_{k1}, g_{k2}, b_{k2}\f$
       */
      static void addPowerHessian(EntryMap&                   hessian,
                                  IdxT                        k,
                                  const std::array<RealT, 4>& admittance,
                                  RealT                       weight_p,
                                  RealT                       weight_q)
      {
        for (IdxT m = 0; m < 2; ++m)
        {
          const RealT g = admittance[2 * m];
          const RealT b = admittance[2 * m + 1];

          // Real and imaginary parts of the weighted admittance (w_P + j w_Q)(g + j b)
          const RealT re = weight_p * g - weight_q * b;
          const RealT im = weight_p * b + weight_q * g;

          addProduct(hessian, 2 * k, 2 * m, re);
          addProduct(hessian, 2 * k, 2 * m + 1, -im);
          addProduct(hessian, 2 * k + 1, 2 * m, im);
          addProduct(hessian, 2 * k + 1, 2 * m + 1, re);
        }
      }

      /// Add the Hessian of \f$c x_i x_j\f$
      static void addProduct(EntryMap& hessian, IdxT i, IdxT j, RealT c)
      {
        RealT scale = ONE<RealT>;
        if (i == j)
        {
          scale = TWO<RealT>;
        }
        addEntry(hessian, i, j, scale * c);
      }
    };
  } // namespace Testing
} // namespace GridKit
