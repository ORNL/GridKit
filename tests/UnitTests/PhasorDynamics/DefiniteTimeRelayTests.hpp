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
#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakers.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/DefiniteTimeRelay/DefiniteTimeRelay.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
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
    class DefiniteTimeRelayTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      DefiniteTimeRelayTests()  = default;
      ~DefiniteTimeRelayTests() = default;

      static constexpr RealT kTol = std::numeric_limits<RealT>::epsilon();

      /// Construction and every verify() error class: required and invalid
      /// parameters, parameter types, required ports, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors from the invalid cases below.
        Log::setVerbosity(Log::Verbosity::NONE);

        Fixture<ScalarT> empty(Data{});
        success *= (empty.relay.size() == static_cast<IdxT>(Utilities::enum_size<Internal>()));
        success *= !empty.relay.verify().passed();

        Fixture<ScalarT> configured(makeOvercurrentData());
        success *= configured.relay.verify().passed();

        for (const Params parameter : {Params::Mpickup, Params::Mref, Params::Ttrip})
        {
          auto missing = makeOvercurrentData();
          missing.parameters.erase(parameter);
          success *= rejected(missing);
        }

        success *= invalidParameterCase(Params::Mpickup, 0.0);
        success *= invalidParameterCase(Params::Mref, -0.1);

        // A reference level equal to the pickup leaves no operating direction.
        success *= invalidParameterCase(Params::Mref, 1.0);

        success *= invalidParameterCase(Params::Ttrip, -0.1);
        success *= invalidParameterCase(Params::Mpickup, true);
        success *= invalidParameterCase(Params::Lockout, static_cast<IdxT>(1));

        for (const Params parameter : {Params::Mpickup, Params::Mref, Params::Ttrip})
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // Integer JSON values are accepted for real parameters.
        auto integer_real                       = makeOvercurrentData();
        integer_real.parameters[Params::Ttrip]  = static_cast<IdxT>(1);
        success                                *= !rejected(integer_real);

        PhasorDynamics::SignalNode<ScalarT, IdxT> trip_node;
        RelayT                                    no_input(makeOvercurrentData());
        no_input.getPorts().out.template port<Outputs::trip>().connect(&trip_node);
        success *= !no_input.verify().passed();

        ScalarT                                   value{0.0};
        IdxT                                      index{0};
        PhasorDynamics::SignalNode<ScalarT, IdxT> input_node;
        input_node.link(&value, &index);
        RelayT no_output(makeOvercurrentData());
        no_output.getPorts().in.template port<Inputs::mr>().connect(&input_node);
        success *= !no_output.verify().passed();

        success *= unlinkedSignalRejected<Inputs::mr>();
        success *= unlinkedSignalRejected<Inputs::mi>();

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// Overcurrent and undervoltage elements pick up on the correct side of
      /// the pickup level, lockout commits past 1/2, and the trip command
      /// rises at 3/4.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        const auto overcurrent  = makeOvercurrentData();
        const auto undervoltage = makeUndervoltageData();

        auto self_reset                        = makeOvercurrentData();
        self_reset.parameters[Params::Lockout] = false;

        success *= residualsMatch(overcurrent, 2.0, 0.0, 0.25, 0.75, 0.0, "overcurrent picked up");
        success *= residualsMatch(overcurrent, 2.0, 0.0, 1.0, 0.0, 1.0, "overcurrent tripped");
        success *= residualsMatch(overcurrent, 0.0, 0.0, 0.25, -0.25, 0.0, "overcurrent reset");
        success *= residualsMatch(overcurrent, 0.0, 0.0, 0.75, 0.25, 0.5, "lockout committed");
        success *= residualsMatch(self_reset, 0.0, 0.0, 0.75, -0.75, 0.5, "self-reset");
        success *= residualsMatch(undervoltage, 0.15, 0.2, 0.25, 0.75, 0.0, "undervoltage picked up");
        success *= residualsMatch(undervoltage, 0.6, 0.8, 0.25, -0.25, 0.0, "undervoltage healthy");

        // An unattached mi is zero, whatever its node holds.
        success *= residualsMatch(undervoltage, 0.25, 1.0, 0.25, 0.75, 0.0, "unattached mi", false);

        return success.report(__func__);
      }

      /// The relay starts reset with zero derivatives and zero residuals
      /// below pickup, and the trip output aliases the trip command.
      TestOutcome initialization()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeOvercurrentData());
        fixture.attachImaginary();
        fixture.input(Inputs::mr)  = 0.5;
        success                   *= fixture.prepare();
        success                   *= (fixture.evaluate() == 0);

        const auto* y  = fixture.relay.y().getData();
        const auto* yp = fixture.relay.yp().getData();
        const auto* f  = fixture.relay.getResidual().getData();
        for (const auto variable : Utilities::enum_values<Internal>())
        {
          const auto k  = static_cast<size_t>(variable);
          success      *= (y[k] == 0.0) && (yp[k] == 0.0) && (f[k] == 0.0);
        }

        fixture.setState(0.5, 0.25);
        success *= (fixture.trip() == 0.25);

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        Fixture<ScalarT> empty(Data{});
        success *= (empty.relay.allocate() == 0);
        success *= (empty.relay.initialize() != 0);
        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// At the pickup knee, with the reset drive depending on the pickup,
      /// every Enzyme row matches dependency tracking.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        for (const RealT x : {0.5, 0.75})
        {
          const auto tracked = dependencyTrackingJacobian(x, success);
          const auto enzyme  = enzymeJacobian(x, success);

          success *= (tracked.size() == enzyme.size());
          for (size_t row = 0; row < std::min(tracked.size(), enzyme.size()); ++row)
          {
            if (!isEqual(enzyme[row], tracked[row], kTol))
            {
              std::cout << "DefiniteTimeRelay Jacobian row " << row << " at x = " << x
                        << " mismatch between dependency tracking and Enzyme\n";
              success = false;
            }
          }
        }

        return success.report(__func__);
      }

      /// An overcurrent relay that reads BranchBreakers terminal-current
      /// outputs reaches the bus voltages through the chain rule, and the
      /// assembled system Jacobian matches dependency tracking.
      TestOutcome breakerCurrentJacobian()
      {
        TestStatus success = true;

        const auto tracked = relaySystemJacobian<DependencyTracking::Variable>(success);
        const auto enzyme  = relaySystemJacobian<ScalarT>(success);

        success *= (tracked.size() == enzyme.size());
        for (size_t row = 0; row < std::min(tracked.size(), enzyme.size()); ++row)
        {
          if (!isEqual(enzyme[row], tracked[row], kTol))
          {
            std::cout << "DefiniteTimeRelay system Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using RelayT        = PhasorDynamics::Relay::DefiniteTimeRelay<ScalarT, IdxT>;
      using Data          = typename RelayT::ModelDataT;
      using Params        = typename Data::Parameters;
      using Inputs        = typename Data::SignalInputs;
      using Outputs       = typename Data::SignalOutputs;
      using Internal      = typename RelayT::InternalVariablesT;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;
      using JacobianT     = std::vector<DependencyMap>;

      static constexpr RealT kTtrip = 0.1;

      // System indices: the relay variables, then the measured quantity
      static constexpr IdxT kInputIndex = 2;

      /// Owns the relay, the measured-quantity signals, and the node that
      /// reads the trip output. Signal storage is declared before the model
      /// so every referenced node outlives it.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, Utilities::enum_size<Inputs>()>    input_values_{};
        std::array<IdxT, Utilities::enum_size<Inputs>()> input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>,
                   Utilities::enum_size<Inputs>()>
                                            input_nodes_{};
        PhasorDynamics::SignalNode<T, IdxT> trip_node_;

      public:
        explicit Fixture(const Data& data)
          : relay(data)
        {
          for (const auto input : Utilities::enum_values<Inputs>())
          {
            const auto k      = static_cast<size_t>(input);
            input_indices_[k] = kInputIndex + static_cast<IdxT>(k);
            input_nodes_[k].link(&input_values_[k], &input_indices_[k]);
          }
          relay.getPorts().in.template port<Inputs::mr>().connect(
              &input_nodes_[static_cast<size_t>(Inputs::mr)]);
          relay.getPorts().out.template port<Outputs::trip>().connect(&trip_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachImaginary()
        {
          relay.getPorts().in.template port<Inputs::mi>().connect(
              &input_nodes_[static_cast<size_t>(Inputs::mi)]);
        }

        bool prepare()
        {
          const bool success = (relay.allocate() == 0) && (relay.initialize() == 0);
          if (!success)
          {
            std::cout << "DefiniteTimeRelay fixture preparation failed\n";
          }
          return success;
        }

        int evaluate()
        {
          return relay.evaluateResidual();
        }

        void setState(RealT x, RealT s)
        {
          auto* y                             = relay.y().getData();
          y[static_cast<size_t>(Internal::X)] = static_cast<T>(x);
          y[static_cast<size_t>(Internal::S)] = static_cast<T>(s);
          relay.y().setDataUpdated();
        }

        T& input(Inputs input)
        {
          return input_values_[static_cast<size_t>(input)];
        }

        const T& trip() const
        {
          return trip_node_.read();
        }

        PhasorDynamics::Relay::DefiniteTimeRelay<T, IdxT> relay;
      };

      /// Instantaneous-style overcurrent element with pickup at 1 p.u.
      Data makeOvercurrentData() const
      {
        Data data;
        data.device_class                = "DefiniteTimeRelay";
        data.disambiguation_string       = "definite_time_relay_test";
        data.parameters[Params::Mpickup] = 1.0;
        data.parameters[Params::Mref]    = 0.0;
        data.parameters[Params::Ttrip]   = kTtrip;
        return data;
      }

      /// Undervoltage element that picks up below half the reference level.
      Data makeUndervoltageData() const
      {
        auto data                        = makeOvercurrentData();
        data.parameters[Params::Mpickup] = 0.5;
        data.parameters[Params::Mref]    = 1.0;
        return data;
      }

      bool rejected(const Data& data) const
      {
        return configurationRejected(
            [&]()
            {
              Fixture<ScalarT> fixture(data);
              return fixture.relay.verify();
            });
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeOvercurrentData();
        data.parameters[parameter] = value;
        return rejected(data);
      }

      template <Inputs variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeOvercurrentData());
        fixture.relay.getPorts().in.template port<variable>().connect(&unlinked_node);
        return !fixture.relay.verify().passed();
      }

      /// Compare the latch residual, scaled by the latch time, and the trip
      /// residual with s = 0 to their expected values.
      bool residualsMatch(const Data& data,
                          RealT       mr,
                          RealT       mi,
                          RealT       x,
                          RealT       expected_latch,
                          RealT       expected_trip,
                          const char* label,
                          bool        attach_mi = true) const
      {
        Fixture<ScalarT> fixture(data);
        if (attach_mi)
        {
          fixture.attachImaginary();
        }
        if (!fixture.prepare())
        {
          return false;
        }

        fixture.input(Inputs::mr) = mr;
        fixture.input(Inputs::mi) = mi;
        fixture.setState(x, 0.0);
        if (fixture.evaluate() != 0)
        {
          return false;
        }

        const RealT latch_time = kTtrip / (2.0 * std::numbers::ln2_v<RealT>);
        const auto* f          = fixture.relay.getResidual().getData();
        const RealT latch      = f[static_cast<size_t>(Internal::X)] * latch_time;
        const RealT trip       = f[static_cast<size_t>(Internal::S)];

        const bool success = isEqual(latch, expected_latch, kTol) && isEqual(trip, expected_trip, kTol);
        if (!success)
        {
          std::cout << "DefiniteTimeRelay " << label << " mismatch: latch "
                    << std::setprecision(16) << latch << " != " << expected_latch
                    << ", trip " << trip << " != " << expected_trip << '\n';
        }
        return success;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Self-reset overcurrent element measuring |m| = 1, at its pickup knee.
      template <typename T>
      void setKneeState(Fixture<T>& fixture, RealT x) const
      {
        fixture.input(Inputs::mr) = 0.6;
        fixture.input(Inputs::mi) = 0.8;
        fixture.setState(x, 0.25);
        fixture.relay.yp().getData()[static_cast<size_t>(Internal::X)] = 0.125;
        fixture.relay.yp().setDataUpdated();
        fixture.relay.updateTime(0.0, 1.0);
      }

      Data makeKneeData() const
      {
        auto data                        = makeOvercurrentData();
        data.parameters[Params::Lockout] = false;
        return data;
      }

      JacobianT dependencyTrackingJacobian(RealT x, TestStatus& success) const
      {
        Fixture<DependencyTracking::Variable> fixture(makeKneeData());
        fixture.attachImaginary();
        success *= fixture.prepare();
        setKneeState(fixture, x);

        // @todo Remove and setup the test to not rely on explicit variable numbering
        auto* y  = fixture.relay.y().getData();
        auto* yp = fixture.relay.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.relay.size()); ++i)
        {
          y[i].setVariableNumber(2 * i);
          yp[i].setVariableNumber(2 * i + 1);
        }
        for (const auto input : Utilities::enum_values<Inputs>())
        {
          const auto k = static_cast<size_t>(kInputIndex) + static_cast<size_t>(input);
          fixture.input(input).setVariableNumber(2 * k);
        }

        success *= (fixture.evaluate() == 0);
        success *= (fixture.relay.evaluateJacobian() == 0);
        return MapFromCsr(fixture.relay.getCsrJacobian());
      }

      JacobianT enzymeJacobian(RealT x, TestStatus& success) const
      {
        Fixture<ScalarT> fixture(makeKneeData());
        fixture.attachImaginary();
        success *= fixture.prepare();
        setKneeState(fixture, x);

        success *= (fixture.evaluate() == 0);
        success *= (fixture.relay.evaluateJacobian() == 0);
        success *= (fixture.relay.constructCsr() == 0);
        return MapFromCsr(fixture.relay.getCsrJacobian());
      }

      /// Bus(1, 0) - BranchBreakers(X = 1/2) - Bus(1/2, 0) carries I1 = j,
      /// so the relay, with pickup at 1 p.u., sits at its pickup knee. The
      /// relay latch row must reach both bus voltages.
      template <typename T>
      JacobianT relaySystemJacobian(TestStatus& success) const
      {
        using BreakerParams  = PhasorDynamics::BranchBreakersParameters;
        using BreakerInputs  = PhasorDynamics::BranchBreakersSignalInputs;
        using BreakerOutputs = PhasorDynamics::BranchBreakersSignalOutputs;

        PhasorDynamics::SignalNode<T, IdxT> ir1_node;
        PhasorDynamics::SignalNode<T, IdxT> ii1_node;
        PhasorDynamics::SignalNode<T, IdxT> trip_node;

        PhasorDynamics::BranchBreakersData<RealT, IdxT> breaker_data;
        breaker_data.parameters[BreakerParams::X] = 0.5;

        PhasorDynamics::Bus<T, IdxT>                      bus1(static_cast<T>(1.0), static_cast<T>(0.0));
        PhasorDynamics::Bus<T, IdxT>                      bus2(static_cast<T>(0.5), static_cast<T>(0.0));
        PhasorDynamics::BranchBreakers<T, IdxT>           breakers(&bus1, &bus2, breaker_data);
        PhasorDynamics::Relay::DefiniteTimeRelay<T, IdxT> relay(makeKneeData());

        breakers.getPorts().out.template port<BreakerOutputs::ir1>().connect(&ir1_node);
        breakers.getPorts().out.template port<BreakerOutputs::ii1>().connect(&ii1_node);
        breakers.getPorts().in.template port<BreakerInputs::trip1>().connect(&trip_node);
        relay.getPorts().in.template port<Inputs::mr>().connect(&ir1_node);
        relay.getPorts().in.template port<Inputs::mi>().connect(&ii1_node);
        relay.getPorts().out.template port<Outputs::trip>().connect(&trip_node);

        PhasorDynamics::SystemModel<T, IdxT> system;
        system.addBus(&bus1);
        system.addBus(&bus2);
        system.addComponent(&breakers);
        system.addComponent(&relay);
        system.updateTime(0.0, 1.0);

        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.evaluateJacobian() == 0);

        const auto rows      = MapFromCsr(system.getCsrJacobian());
        const auto latch_row = static_cast<size_t>(relay.getResidualIndices()[static_cast<size_t>(Internal::X)]);
        for (const auto* bus : {&bus1, &bus2})
        {
          for (const auto column : bus->getVariableIndices())
          {
            success *= rows[latch_row].contains(static_cast<size_t>(column));
          }
        }
        return rows;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
