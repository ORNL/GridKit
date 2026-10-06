#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakers.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class BranchBreakersTests
    {
    private:
      using RealT         = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DataT         = PhasorDynamics::BranchBreakersData<RealT, IdxT>;
      using Parameters    = PhasorDynamics::BranchBreakersParameters;
      using Inputs        = PhasorDynamics::BranchBreakersSignalInputs;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      static constexpr RealT kTbrk = 0.05;
      static constexpr RealT kTol  = 100.0 * std::numeric_limits<RealT>::epsilon();

      /// Off-nominal branch shared with BranchTests::offNominalResidual.
      static DataT makeData()
      {
        DataT data;
        data.parameters[Parameters::R]     = 2.0;
        data.parameters[Parameters::X]     = 4.0;
        data.parameters[Parameters::G]     = 0.4;
        data.parameters[Parameters::B]     = 0.8;
        data.parameters[Parameters::Gmag]  = 0.2;
        data.parameters[Parameters::Bmag]  = 1.2;
        data.parameters[Parameters::tap]   = 1.25;
        data.parameters[Parameters::phase] = 0.3;
        data.parameters[Parameters::Tbrk]  = kTbrk;
        return data;
      }

      /// Breaker-terminated branch with fixture-owned trip and reset commands.
      template <class T>
      struct Fixture
      {
        PhasorDynamics::Bus<T, IdxT>            bus1{T{10.0}, T{20.0}};
        PhasorDynamics::Bus<T, IdxT>            bus2{T{30.0}, T{40.0}};
        PhasorDynamics::BranchBreakers<T, IdxT> branch{&bus1, &bus2, makeData()};

        std::array<T, 4>                                   command{};
        std::array<IdxT, 4>                                command_index{6, 7, 8, 9};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, 4> command_node{};

        Fixture()
        {
          for (size_t i = 0; i < command.size(); ++i)
          {
            command[i] = T{0.0};
            command_node[i].link(&command[i], &command_index[i]);
          }
          branch.getPorts().in.template port<Inputs::trip1>().connect(&command_node[0]);
          branch.getPorts().in.template port<Inputs::reset1>().connect(&command_node[1]);
          branch.getPorts().in.template port<Inputs::trip2>().connect(&command_node[2]);
          branch.getPorts().in.template port<Inputs::reset2>().connect(&command_node[3]);

          bus1.allocate();
          bus2.allocate();
          branch.allocate();
          for (IdxT i = 0; i < 2; ++i)
          {
            bus1.setVariableIndex(i, 2 + i);
            bus1.setResidualIndex(i, 2 + i);
            bus2.setVariableIndex(i, 4 + i);
            bus2.setResidualIndex(i, 4 + i);
          }
          bus1.initialize();
          bus2.initialize();
          branch.initialize();
        }

        /// Residuals at latch states z1, z2 and commands trip1, reset1, trip2, reset2.
        const T* evaluate(RealT z1, RealT z2, std::array<RealT, 4> commands)
        {
          auto* y  = branch.y().getData();
          auto* yp = branch.yp().getData();
          y[0]     = z1;
          y[1]     = z2;
          yp[0]    = 0.0;
          yp[1]    = 0.0;
          for (size_t i = 0; i < command.size(); ++i)
          {
            command[i] = commands[i];
          }
          if constexpr (std::is_same_v<T, DependencyTracking::Variable>)
          {
            for (size_t i = 0; i < 2; ++i)
            {
              y[i].setVariableNumber(2 * i);
              yp[i].setVariableNumber(2 * i + 1);
            }
            for (size_t i = 0; i < command.size(); ++i)
            {
              command[i].setVariableNumber(2 * command_index[i]);
            }
          }
          branch.y().setDataUpdated();
          branch.yp().setDataUpdated();
          branch.updateTime(0.0, 1.0);
          bus1.evaluateResidual();
          bus2.evaluateResidual();
          branch.evaluateResidual();
          return branch.getResidual().getData();
        }
      };

    public:
      BranchBreakersTests()  = default;
      ~BranchBreakersTests() = default;

      TestOutcome residual()
      {
        // Open breakers reproduce the Kron-reduced sides of the pi model.
        TestStatus success = true;

        Fixture<ScalarT> fixture;

        auto currents = [&](RealT z1, RealT z2, RealT ir1, RealT ii1, RealT ir2, RealT ii2)
        {
          fixture.evaluate(z1, z2, {0.0, 0.0, 0.0, 0.0});
          success *= isEqual(fixture.bus1.Ir(), ir1, kTol);
          success *= isEqual(fixture.bus1.Ii(), ii1, kTol);
          success *= isEqual(fixture.bus2.Ir(), ir2, kTol);
          success *= isEqual(fixture.bus2.Ii(), ii2, kTol);
        };

        currents(1.0, 1.0, 33.679793434963472, -22.927960563981181, 2.821345956502423, -19.182080826645358);
        currents(1.0, 0.0, 24.553846153846159, -25.969230769230769, 0.0, 0.0);
        currents(0.0, 1.0, 0.0, 0.0, -1.8619022031166033, -18.576034390112845);
        currents(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        return success.report(__func__);
      }

      TestOutcome latch()
      {
        // Trip has priority over reset, a closed breaker holds, and reset closes.
        TestStatus success = true;

        Fixture<ScalarT> fixture;
        const RealT      rate = std::numbers::ln2_v<RealT> / kTbrk;

        const auto* f  = fixture.evaluate(1.0, 0.0, {1.0, 1.0, 0.0, 1.0});
        success       *= isEqual(f[0], -rate, kTol);
        success       *= isEqual(f[1], rate, kTol);

        f        = fixture.evaluate(1.0, 1.0, {0.0, 0.0, 0.0, 0.0});
        success *= isEqual(f[0], 0.0, kTol);
        success *= isEqual(f[1], 0.0, kTol);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        // Enzyme matches DependencyTracking for the latches and the terminal currents.
        TestStatus success = true;

        const auto reference = dependencyTrackingJacobian();
        const auto enzyme    = enzymeJacobian();

        success *= (reference.size() == enzyme.size());
        for (size_t i = 0; i < reference.size() && i < enzyme.size(); ++i)
        {
          success *= isEqual(reference[i], enzyme[i], kTol);
        }

        return success.report(__func__);
      }

      TestOutcome infiniteBusJacobian()
      {
        // Infinite-bus rows and columns are skipped, not stored.
        TestStatus success = true;

        PhasorDynamics::BusInfinite<ScalarT, IdxT>    bus1(1.0, 0.1);
        PhasorDynamics::Bus<ScalarT, IdxT>            bus2(0.9, -0.2);
        PhasorDynamics::BranchBreakers<ScalarT, IdxT> branch(&bus1, &bus2, makeData());

        bus1.allocate();
        bus2.allocate();
        branch.allocate();
        for (IdxT i = 0; i < 2; ++i)
        {
          bus2.setVariableIndex(i, 2 + i);
          bus2.setResidualIndex(i, 2 + i);
        }
        bus1.initialize();
        bus2.initialize();
        branch.initialize();
        branch.updateTime(0.0, 1.0);
        bus2.evaluateResidual();
        branch.evaluateResidual();
        branch.evaluateJacobian();

        // Latch diagonals in y and yp, and bus-2 currents in both latches and bus-2 voltages.
        const IdxT  expected_nnz = 12;
        auto*       jacobian     = branch.getCooJacobian();
        const auto  nnz          = jacobian->getNnz();
        const auto* rows         = jacobian->getRowData();
        const auto* cols         = jacobian->getColData();

        success *= (nnz == expected_nnz);
        for (IdxT i = 0; i < nnz; ++i)
        {
          success *= (rows[i] != INVALID_INDEX<IdxT>);
          success *= (cols[i] != INVALID_INDEX<IdxT>);
        }

        return success.report(__func__);
      }

    private:
      static constexpr RealT                kLatch    = 0.5;
      static constexpr std::array<RealT, 4> kCommands = {0.3, 0.2, 0.1, 0.4};

      /// Latch rows from the component CSR and bus-current rows from the bus residuals.
      std::vector<DependencyMap> dependencyTrackingJacobian() const
      {
        Fixture<DependencyTracking::Variable> fixture;
        fixture.evaluate(kLatch, kLatch, kCommands);
        fixture.branch.evaluateJacobian();

        std::vector<DependencyMap> rows = MapFromCsr(fixture.branch.getCsrJacobian());
        for (const auto* current : {&fixture.bus1.Ir(), &fixture.bus1.Ii(), &fixture.bus2.Ir(), &fixture.bus2.Ii()})
        {
          // Bus currents carry no derivative terms, so every variable number is even.
          DependencyMap columns;
          for (const auto& [number, value] : current->getDependencies())
          {
            columns[number / 2] += value;
          }
          rows.push_back(columns);
        }
        return rows;
      }

      std::vector<DependencyMap> enzymeJacobian() const
      {
        Fixture<ScalarT> fixture;
        fixture.evaluate(kLatch, kLatch, kCommands);
        fixture.branch.evaluateJacobian();
        fixture.branch.constructCsr();
        return MapFromCsr(fixture.branch.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
