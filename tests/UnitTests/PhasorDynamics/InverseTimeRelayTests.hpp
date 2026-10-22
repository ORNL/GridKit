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
#include <GridKit/Model/PhasorDynamics/Relay/InverseTimeRelay/InverseTimeRelay.hpp>
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
    class InverseTimeRelayTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      InverseTimeRelayTests()  = default;
      ~InverseTimeRelayTests() = default;

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

        Fixture<ScalarT> configured(makeData());
        success *= configured.relay.verify().passed();

        for (const Params parameter : {Params::Ipickup, Params::TD, Params::A, Params::B, Params::alpha})
        {
          auto missing = makeData();
          missing.parameters.erase(parameter);
          success *= rejected(missing);
        }

        success *= invalidParameterCase(Params::Ipickup, 0.0);
        success *= invalidParameterCase(Params::TD, 0.0);
        success *= invalidParameterCase(Params::A, 0.0);
        success *= invalidParameterCase(Params::B, -0.1);
        success *= invalidParameterCase(Params::alpha, 0.0);
        success *= invalidParameterCase(Params::tr, -0.1);
        success *= invalidParameterCase(Params::Ipickup, true);
        success *= invalidParameterCase(Params::Lockout, static_cast<IdxT>(1));

        for (const Params parameter :
             {Params::Ipickup, Params::TD, Params::A, Params::B, Params::alpha, Params::tr})
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // Integer JSON values are accepted for real parameters.
        auto integer_real                    = makeData();
        integer_real.parameters[Params::TD]  = static_cast<IdxT>(1);
        success                             *= !rejected(integer_real);

        success *= missingPortRejected<Inputs::ir>();
        success *= missingPortRejected<Inputs::ii>();

        ScalarT                                   value{0.0};
        IdxT                                      index{0};
        PhasorDynamics::SignalNode<ScalarT, IdxT> input_node;
        input_node.link(&value, &index);
        RelayT no_output(makeData());
        no_output.getPorts().in.template port<Inputs::ir>().connect(&input_node);
        no_output.getPorts().in.template port<Inputs::ii>().connect(&input_node);
        success *= !no_output.verify().passed();

        success *= unlinkedSignalRejected<Inputs::ir>();
        success *= unlinkedSignalRejected<Inputs::ii>();

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// Above pickup the latch rate is ln 4 over the IEEE C37.112 operate
      /// time; below pickup it is ln 4 over the reset time, and lockout
      /// commits past 1/2 at that rate.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        // A = 1, B = 0, alpha = 2, TD = 1: the operate time at M = 2 is 1/3 s.
        const auto very_inverse = makeData();

        // B = 1 adds TD B to the operate time, giving 4/3 s at M = 2.
        auto offset                  = makeData();
        offset.parameters[Params::B] = 1.0;

        // TD tr = 2 gives a reset time of 2 s at M = 0.
        auto self_reset                        = makeData();
        self_reset.parameters[Params::Lockout] = false;

        success *= residualsMatch(very_inverse, 2.0, 0.0, 0.25, 3.0 * 0.75, 0.0, "operate");
        success *= residualsMatch(very_inverse, 0.0, 2.0, 0.25, 3.0 * 0.75, 0.0, "operate on ii");
        success *= residualsMatch(offset, 2.0, 0.0, 0.25, 0.75 * 0.75, 0.0, "offset curve");
        success *= residualsMatch(very_inverse, 0.0, 0.0, 0.25, 0.5 * -0.25, 0.0, "reset");
        success *= residualsMatch(very_inverse, 0.0, 0.0, 0.75, 0.5 * 0.25, 0.5, "lockout committed");
        success *= residualsMatch(self_reset, 0.0, 0.0, 0.75, 0.5 * -0.75, 0.5, "self-reset");

        return success.report(__func__);
      }

      /// The relay starts reset with zero derivatives and zero residuals
      /// below pickup, and the trip output aliases the trip command.
      TestOutcome initialization()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeData());
        fixture.input(Inputs::ir)  = 0.5;
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
      /// At the pickup knee, where the operate and reset rates and the
      /// reset drive all depend on the current, every Enzyme row matches
      /// dependency tracking.
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
              std::cout << "InverseTimeRelay Jacobian row " << row << " at x = " << x
                        << " mismatch between dependency tracking and Enzyme\n";
              success = false;
            }
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using RelayT        = PhasorDynamics::Relay::InverseTimeRelay<ScalarT, IdxT>;
      using Data          = typename RelayT::ModelDataT;
      using Params        = typename Data::Parameters;
      using Inputs        = typename Data::SignalInputs;
      using Outputs       = typename Data::SignalOutputs;
      using Internal      = typename RelayT::InternalVariablesT;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;
      using JacobianT     = std::vector<DependencyMap>;

      // System indices: the relay variables, then the measured current
      static constexpr IdxT kInputIndex = 2;

      /// Owns the relay, the measured-current signals, and the node that
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
            relay.getPorts().in[input].connect(&input_nodes_[k]);
          }
          relay.getPorts().out.template port<Outputs::trip>().connect(&trip_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        bool prepare()
        {
          const bool success = (relay.allocate() == 0) && (relay.initialize() == 0);
          if (!success)
          {
            std::cout << "InverseTimeRelay fixture preparation failed\n";
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

        PhasorDynamics::Relay::InverseTimeRelay<T, IdxT> relay;
      };

      /// IEEE C37.112-shaped curve with round coefficients and pickup at 1 p.u.
      Data makeData() const
      {
        Data data;
        data.device_class                = "InverseTimeRelay";
        data.disambiguation_string       = "inverse_time_relay_test";
        data.parameters[Params::Ipickup] = 1.0;
        data.parameters[Params::TD]      = 1.0;
        data.parameters[Params::A]       = 1.0;
        data.parameters[Params::B]       = 0.0;
        data.parameters[Params::alpha]   = 2.0;
        data.parameters[Params::tr]      = 2.0;
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
        auto data                  = makeData();
        data.parameters[parameter] = value;
        return rejected(data);
      }

      /// A relay with the trip output and only the other current input.
      template <Inputs missing>
      bool missingPortRejected() const
      {
        ScalarT                                   value{0.0};
        IdxT                                      index{0};
        PhasorDynamics::SignalNode<ScalarT, IdxT> input_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> trip_node;
        input_node.link(&value, &index);

        RelayT relay(makeData());
        for (const auto input : Utilities::enum_values<Inputs>())
        {
          if (input != missing)
          {
            relay.getPorts().in[input].connect(&input_node);
          }
        }
        relay.getPorts().out.template port<Outputs::trip>().connect(&trip_node);
        return !relay.verify().passed();
      }

      template <Inputs variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData());
        fixture.relay.getPorts().in.template port<variable>().connect(&unlinked_node);
        return !fixture.relay.verify().passed();
      }

      /// Compare the latch residual, scaled by 1 / ln 4, and the trip
      /// residual with s = 0 to their expected values.
      bool residualsMatch(const Data& data,
                          RealT       ir,
                          RealT       ii,
                          RealT       x,
                          RealT       expected_latch,
                          RealT       expected_trip,
                          const char* label) const
      {
        Fixture<ScalarT> fixture(data);
        if (!fixture.prepare())
        {
          return false;
        }

        fixture.input(Inputs::ir) = ir;
        fixture.input(Inputs::ii) = ii;
        fixture.setState(x, 0.0);
        if (fixture.evaluate() != 0)
        {
          return false;
        }

        const RealT ln4   = 2.0 * std::numbers::ln2_v<RealT>;
        const auto* f     = fixture.relay.getResidual().getData();
        const RealT latch = f[static_cast<size_t>(Internal::X)] / ln4;
        const RealT trip  = f[static_cast<size_t>(Internal::S)];

        const bool success = isEqual(latch, expected_latch, kTol) && isEqual(trip, expected_trip, kTol);
        if (!success)
        {
          std::cout << "InverseTimeRelay " << label << " mismatch: latch "
                    << std::setprecision(16) << latch << " != " << expected_latch
                    << ", trip " << trip << " != " << expected_trip << '\n';
        }
        return success;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Self-reset relay with an offset curve measuring |I| = 1, at its
      /// pickup knee.
      template <typename T>
      void setKneeState(Fixture<T>& fixture, RealT x) const
      {
        fixture.input(Inputs::ir) = 0.6;
        fixture.input(Inputs::ii) = 0.8;
        fixture.setState(x, 0.25);
        fixture.relay.yp().getData()[static_cast<size_t>(Internal::X)] = 0.125;
        fixture.relay.yp().setDataUpdated();
        fixture.relay.updateTime(0.0, 1.0);
      }

      Data makeKneeData() const
      {
        auto data                        = makeData();
        data.parameters[Params::B]       = 0.5;
        data.parameters[Params::Lockout] = false;
        return data;
      }

      JacobianT dependencyTrackingJacobian(RealT x, TestStatus& success) const
      {
        Fixture<DependencyTracking::Variable> fixture(makeKneeData());
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
        success *= fixture.prepare();
        setKneeState(fixture, x);

        success *= (fixture.evaluate() == 0);
        success *= (fixture.relay.evaluateJacobian() == 0);
        success *= (fixture.relay.constructCsr() == 0);
        return MapFromCsr(fixture.relay.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
