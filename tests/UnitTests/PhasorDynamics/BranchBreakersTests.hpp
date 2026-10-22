#pragma once

#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakers.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class BranchBreakersTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      BranchBreakersTests()  = default;
      ~BranchBreakersTests() = default;

      // The worst entry, bus-2 ii2 against z2, sums block slopes that cancel
      // to -8.5 and matches at 6.6 eps.
      static constexpr RealT kTol = 8.0 * std::numeric_limits<RealT>::epsilon();

      /// Construction and every verify() error class: parameter types and
      /// values, the impedance and self-admittances, buses, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors from the invalid cases below.
        Log::setVerbosity(Log::Verbosity::NONE);

        // Without R and X the series impedance is zero.
        Fixture<ScalarT> empty(Data{});
        success *= (empty.breakers.size() == static_cast<IdxT>(Utilities::enum_size<Internal>()));
        success *= !empty.breakers.verify().passed();

        Fixture<ScalarT> configured(makeData());
        success *= configured.breakers.verify().passed();

        auto zero_impedance                   = makeData();
        zero_impedance.parameters[Params::R]  = 0.0;
        zero_impedance.parameters[Params::X]  = 0.0;
        success                              *= rejected(zero_impedance);

        // A shunt that cancels the series admittance leaves Y11 = Y22 = 0.
        auto cancelled                   = makeOpenEndData();
        cancelled.parameters[Params::B]  = 4.0;
        success                         *= rejected(cancelled);

        success *= invalidParameterCase(Params::tap, 0.0);
        success *= invalidParameterCase(Params::Tbrk, -0.1);
        success *= invalidParameterCase(Params::X, true);

        for (const Params parameter : {Params::R,
                                       Params::X,
                                       Params::G,
                                       Params::B,
                                       Params::Gmag,
                                       Params::Bmag,
                                       Params::tap,
                                       Params::phase,
                                       Params::Tbrk})
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // Integer JSON values are accepted for real parameters.
        auto integer_real                   = makeData();
        integer_real.parameters[Params::X]  = static_cast<IdxT>(1);
        success                            *= !rejected(integer_real);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        BranchBreakersT                    bus1less(nullptr, &bus, makeData());
        BranchBreakersT                    bus2less(&bus, nullptr, makeData());
        success *= !bus1less.verify().passed();
        success *= !bus2less.verify().passed();

        success *= unlinkedSignalRejected<Inputs::trip1>();
        success *= unlinkedSignalRejected<Inputs::close1>();
        success *= unlinkedSignalRejected<Inputs::trip2>();
        success *= unlinkedSignalRejected<Inputs::close2>();

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// With both breakers closed, the terminal currents are those of a
      /// plain Branch with the same data, and they enter the bus residuals.
      TestOutcome closedMatchesBranch()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus1(kVr1, kVi1);
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(kVr2, kVi2);
        success *= (bus1.allocate() == 0) && (bus1.initialize() == 0);
        success *= (bus2.allocate() == 0) && (bus2.initialize() == 0);

        PhasorDynamics::BranchData<RealT, IdxT> line;
        line.parameters[PhasorDynamics::BranchParameters::R]     = kR;
        line.parameters[PhasorDynamics::BranchParameters::X]     = kX;
        line.parameters[PhasorDynamics::BranchParameters::G]     = kG;
        line.parameters[PhasorDynamics::BranchParameters::B]     = kB;
        line.parameters[PhasorDynamics::BranchParameters::Gmag]  = kGmag;
        line.parameters[PhasorDynamics::BranchParameters::Bmag]  = kBmag;
        line.parameters[PhasorDynamics::BranchParameters::tap]   = kTap;
        line.parameters[PhasorDynamics::BranchParameters::phase] = kPhase;

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, line);
        success *= (bus1.evaluateResidual() == 0);
        success *= (bus2.evaluateResidual() == 0);
        success *= (branch.evaluateResidual() == 0);

        Fixture<ScalarT> closed(makeData());
        success *= closed.prepare();
        success *= (closed.evaluate() == 0);

        success *= currentsMatch(closed, {bus1.Ir(), bus1.Ii(), bus2.Ir(), bus2.Ii()}, "closed");
        success *= currentsMatch(closed,
                                 {closed.bus1.Ir(), closed.bus1.Ii(), closed.bus2.Ir(), closed.bus2.Ii()},
                                 "bus residual");

        return success.report(__func__);
      }

      /// With one breaker open, the other terminal sees the open-end
      /// admittance; with both open, the branch carries no current.
      TestOutcome openEnds()
      {
        TestStatus success = true;

        // R = 0, X = 1/2, and B = 2 give Y11 = Y22 = j and Y12 = Y21 = -2j,
        // so each open-end admittance is j - (-4)/j = -3j.
        Fixture<ScalarT> fixture(makeOpenEndData(), 1.0, 0.0, 1.0, 0.0);
        success *= fixture.prepare();

        fixture.setLatches(1.0, 0.0);
        success *= (fixture.evaluate() == 0);
        success *= currentsMatch(fixture, {0.0, -3.0, 0.0, 0.0}, "bus-2 breaker open");

        fixture.setLatches(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= currentsMatch(fixture, {0.0, 0.0, 0.0, -3.0}, "bus-1 breaker open");

        fixture.setLatches(0.0, 0.0);
        success *= (fixture.evaluate() == 0);
        success *= currentsMatch(fixture, {0.0, 0.0, 0.0, 0.0}, "both breakers open");

        return success.report(__func__);
      }

      /// A close command closes, a trip command opens and overrides a close,
      /// and without commands each latch holds its side of 1/2. The latch
      /// reaches 1/2 after Tbrk, so the latch time is Tbrk / ln 2.
      TestOutcome latchResidual()
      {
        TestStatus success = true;

        const RealT latch_time = kTbrk / std::numbers::ln2_v<RealT>;

        Fixture<ScalarT> fixture(makeData());
        fixture.attachCommands();
        success *= fixture.prepare();

        fixture.setLatches(0.25, 0.75);
        fixture.command(Inputs::close1)  = 1.0;
        fixture.command(Inputs::trip2)   = 1.0;
        success                         *= (fixture.evaluate() == 0);
        success                         *= latchRatesMatch(fixture, latch_time, {0.75, -0.75}, "close and trip");

        fixture.setLatches(0.75, 0.25);
        fixture.command(Inputs::trip1)  = 1.0;
        fixture.command(Inputs::trip2)  = 0.0;
        success                        *= (fixture.evaluate() == 0);
        success                        *= latchRatesMatch(fixture, latch_time, {-0.75, -0.25}, "trip priority and hold open");

        fixture.setLatches(0.75, 0.75);
        fixture.command(Inputs::close1)  = 0.0;
        fixture.command(Inputs::trip1)   = 0.0;
        success                         *= (fixture.evaluate() == 0);
        success                         *= latchRatesMatch(fixture, latch_time, {0.25, 0.25}, "hold closed");

        // A zero breaker time uses the 1 ms floor.
        auto instant                     = makeData();
        instant.parameters[Params::Tbrk] = 0.0;

        Fixture<ScalarT> floored(instant);
        floored.attachCommands();
        success                        *= floored.prepare();
        floored.command(Inputs::trip1)  = 1.0;
        success                        *= (floored.evaluate() == 0);
        success                        *= latchRatesMatch(floored,
                                   1.0e-3 / std::numbers::ln2_v<RealT>,
                                                          {-1.0, 0.0},
                                   "floored breaker time");

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Every Enzyme block matches dependency tracking at a partly open
      /// state: the latch rows, both bus-current rows, and the
      /// terminal-current output rows that readers chain through.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> tracked(makeData());
        tracked.attachCommands(kJacobianCommand);
        success *= tracked.prepare();
        tracked.setLatches(kJacobianZ1, kJacobianZ2);
        numberVariables(tracked);
        tracked.breakers.updateTime(0.0, 1.0);
        success *= (tracked.evaluate() == 0);
        success *= (tracked.breakers.evaluateJacobian() == 0);

        const auto latch_rows = MapFromCsr(tracked.breakers.getCsrJacobian());

        Fixture<ScalarT> enzyme(makeData());
        enzyme.attachCommands(kJacobianCommand);
        success *= enzyme.prepare();
        enzyme.setLatches(kJacobianZ1, kJacobianZ2);
        enzyme.breakers.updateTime(0.0, 1.0);
        success *= (enzyme.evaluate() == 0);
        success *= (enzyme.breakers.evaluateJacobian() == 0);
        success *= (enzyme.breakers.constructCsr() == 0);

        const auto rows = MapFromCsr(enzyme.breakers.getCsrJacobian());

        success *= (rows.size() == static_cast<size_t>(kOutputIndex) + Utilities::enum_size<Outputs>());

        for (size_t row = 0; row < latch_rows.size(); ++row)
        {
          success *= rowMatches(rows[row], latch_rows[row], "latch", row);
        }

        // The terminal currents are both the bus-residual contributions and
        // the computed outputs.
        const std::array<IdxT, Utilities::enum_size<Outputs>()> bus_rows{kBus1Index,
                                                                         kBus1Index + 1,
                                                                         kBus2Index,
                                                                         kBus2Index + 1};
        for (const auto output : Utilities::enum_values<Outputs>())
        {
          const auto k          = static_cast<size_t>(output);
          const auto bus_row    = static_cast<size_t>(bus_rows[k]);
          const auto output_row = static_cast<size_t>(kOutputIndex) + k;
          const auto expected   = variableColumns(tracked.current(output));

          success *= rowMatches(rows[bus_row], expected, "bus current", bus_row);
          success *= rowMatches(rows[output_row], expected, "output", output_row);
        }

        return success.report(__func__);
      }
#endif

    private:
      using BranchBreakersT = PhasorDynamics::BranchBreakers<ScalarT, IdxT>;
      using Data            = typename BranchBreakersT::ModelDataT;
      using Params          = typename Data::Parameters;
      using Inputs          = typename Data::SignalInputs;
      using Outputs         = typename Data::SignalOutputs;
      using Internal        = typename BranchBreakersT::InternalVariablesT;
      using DependencyMap   = DependencyTracking::Variable::DependencyMap;

      static constexpr RealT kR     = 0.125;
      static constexpr RealT kX     = 0.5;
      static constexpr RealT kG     = 0.0625;
      static constexpr RealT kB     = 0.25;
      static constexpr RealT kGmag  = 0.0625;
      static constexpr RealT kBmag  = -0.125;
      static constexpr RealT kTap   = 1.25;
      static constexpr RealT kPhase = 0.5;
      static constexpr RealT kTbrk  = 0.05;

      static constexpr RealT kVr1 = 1.0;
      static constexpr RealT kVi1 = 0.25;
      static constexpr RealT kVr2 = 0.75;
      static constexpr RealT kVi2 = -0.125;

      // Both latches at their switching point, where every latch column is large
      static constexpr RealT kJacobianZ1      = 0.5;
      static constexpr RealT kJacobianZ2      = 0.5;
      static constexpr RealT kJacobianCommand = 0.5;

      // System indices: latches first, then both buses, the commands, and the outputs
      static constexpr IdxT kBus1Index   = 2;
      static constexpr IdxT kBus2Index   = 4;
      static constexpr IdxT kInputIndex  = 6;
      static constexpr IdxT kOutputIndex = 10;

      /// Owns both terminal buses, the BranchBreakers, the command signals,
      /// and the nodes that read the terminal-current outputs. Signal storage
      /// is declared before the model so every referenced node outlives it.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, Utilities::enum_size<Inputs>()>    command_values_{};
        std::array<IdxT, Utilities::enum_size<Inputs>()> command_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>,
                   Utilities::enum_size<Inputs>()>
            command_nodes_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>,
                   Utilities::enum_size<Outputs>()>
            output_nodes_{};

      public:
        explicit Fixture(const Data& data,
                         RealT       vr1 = kVr1,
                         RealT       vi1 = kVi1,
                         RealT       vr2 = kVr2,
                         RealT       vi2 = kVi2)
          : bus1(static_cast<T>(vr1), static_cast<T>(vi1)),
            bus2(static_cast<T>(vr2), static_cast<T>(vi2)),
            breakers(&bus1, &bus2, data)
        {
          for (const auto output : Utilities::enum_values<Outputs>())
          {
            breakers.getPorts().out[output].connect(&output_nodes_[static_cast<size_t>(output)]);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every breaker command.
        void attachCommands(RealT value = 0.0)
        {
          for (const auto input : Utilities::enum_values<Inputs>())
          {
            const auto k        = static_cast<size_t>(input);
            command_values_[k]  = static_cast<T>(value);
            command_indices_[k] = kInputIndex + static_cast<IdxT>(k);
            command_nodes_[k].link(&command_values_[k], &command_indices_[k]);
            breakers.getPorts().in[input].connect(&command_nodes_[k]);
          }
        }

        /// Allocate and initialize the buses and the breakers, with the
        /// system indices a SystemModel would assign.
        bool prepare()
        {
          bool success = (bus1.allocate() == 0) && (bus2.allocate() == 0);
          for (IdxT i = 0; i < bus1.size(); ++i)
          {
            bus1.setVariableIndex(i, kBus1Index + i);
            bus1.setResidualIndex(i, kBus1Index + i);
            bus2.setVariableIndex(i, kBus2Index + i);
            bus2.setResidualIndex(i, kBus2Index + i);
          }
          success = success && (bus1.initialize() == 0) && (bus2.initialize() == 0)
                    && (breakers.allocate() == 0);
          for (IdxT k = 0; k < breakers.sizeOutputs(); ++k)
          {
            breakers.setOutputIndex(k, kOutputIndex + k);
          }
          success = success && (breakers.initialize() == 0);
          if (!success)
          {
            std::cout << "BranchBreakers fixture preparation failed\n";
          }
          return success;
        }

        /// Clear the bus residuals and evaluate the breakers.
        int evaluate()
        {
          bus1.evaluateResidual();
          bus2.evaluateResidual();
          return breakers.evaluateResidual();
        }

        void setLatches(RealT z1, RealT z2)
        {
          auto* y                              = breakers.y().getData();
          y[static_cast<size_t>(Internal::Z1)] = static_cast<T>(z1);
          y[static_cast<size_t>(Internal::Z2)] = static_cast<T>(z2);
          breakers.y().setDataUpdated();
        }

        T& command(Inputs input)
        {
          return command_values_[static_cast<size_t>(input)];
        }

        const T& current(Outputs output) const
        {
          return output_nodes_[static_cast<size_t>(output)].read();
        }

        PhasorDynamics::Bus<T, IdxT>            bus1;
        PhasorDynamics::Bus<T, IdxT>            bus2;
        PhasorDynamics::BranchBreakers<T, IdxT> breakers;
      };

      Data makeData() const
      {
        Data data;
        data.device_class              = "BranchBreakers";
        data.disambiguation_string     = "branch_breakers_test";
        data.parameters[Params::R]     = kR;
        data.parameters[Params::X]     = kX;
        data.parameters[Params::G]     = kG;
        data.parameters[Params::B]     = kB;
        data.parameters[Params::Gmag]  = kGmag;
        data.parameters[Params::Bmag]  = kBmag;
        data.parameters[Params::tap]   = kTap;
        data.parameters[Params::phase] = kPhase;
        data.parameters[Params::Tbrk]  = kTbrk;
        return data;
      }

      Data makeOpenEndData() const
      {
        Data data;
        data.parameters[Params::R] = 0.0;
        data.parameters[Params::X] = 0.5;
        data.parameters[Params::B] = 2.0;
        return data;
      }

      bool rejected(const Data& data) const
      {
        return configurationRejected(
            [&]()
            {
              Fixture<ScalarT> fixture(data);
              return fixture.breakers.verify();
            });
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        return rejected(data);
      }

      template <Inputs variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData());
        fixture.breakers.getPorts().in.template port<variable>().connect(&unlinked_node);
        return !fixture.breakers.verify().passed();
      }

      bool currentsMatch(const Fixture<ScalarT>&                                   fixture,
                         const std::array<RealT, Utilities::enum_size<Outputs>()>& expected,
                         const char*                                               context) const
      {
        bool success = true;
        for (const auto output : Utilities::enum_values<Outputs>())
        {
          const auto  k      = static_cast<size_t>(output);
          const RealT actual = fixture.current(output);
          if (!isEqual(actual, expected[k], kTol))
          {
            std::cout << "BranchBreakers " << context << " current " << k << " mismatch: "
                      << std::setprecision(16) << actual << " != " << expected[k] << '\n';
            success = false;
          }
        }
        return success;
      }

      /// Compare each latch residual, scaled by the latch time, to the
      /// expected latch rate.
      bool latchRatesMatch(const Fixture<ScalarT>&                                    fixture,
                           RealT                                                      latch_time,
                           const std::array<RealT, Utilities::enum_size<Internal>()>& expected,
                           const char*                                                context) const
      {
        bool        success = true;
        const auto* f       = fixture.breakers.getResidual().getData();
        for (size_t k = 0; k < expected.size(); ++k)
        {
          const RealT actual = f[k] * latch_time;
          if (!isEqual(actual, expected[k], kTol))
          {
            std::cout << "BranchBreakers " << context << " latch " << k << " mismatch: "
                      << std::setprecision(16) << actual << " != " << expected[k] << '\n';
            success = false;
          }
        }
        return success;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// @todo Remove and setup the test to not rely on explicit variable numbering
      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y  = fixture.breakers.y().getData();
        auto* yp = fixture.breakers.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.breakers.size()); ++i)
        {
          y[i].setVariableNumber(2 * i);
          yp[i].setVariableNumber(2 * i + 1);
        }
        for (const auto input : Utilities::enum_values<Inputs>())
        {
          const auto k = static_cast<size_t>(kInputIndex) + static_cast<size_t>(input);
          fixture.command(input).setVariableNumber(2 * k);
        }
        fixture.breakers.y().setDataUpdated();
        fixture.breakers.yp().setDataUpdated();
      }

      /// Dependencies of a tracked value, keyed by variable column.
      static DependencyMap variableColumns(const DependencyTracking::Variable& value)
      {
        DependencyMap columns;
        for (const auto& [number, derivative] : value.getDependencies())
        {
          columns[number / 2] += derivative;
        }
        return columns;
      }

      static bool rowMatches(const DependencyMap& actual,
                             const DependencyMap& expected,
                             const char*          what,
                             size_t               row)
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "BranchBreakers " << what << " row " << row
                  << " mismatch between dependency tracking and Enzyme\n";
        return false;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
