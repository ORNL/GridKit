/**
 * @file ForcedOscillationTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Unit tests for the forced-oscillation signal source.
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ForcedOscillation.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ForcedOscillationData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ForcedOscillationTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using SourceT = PhasorDynamics::ForcedOscillation<ScalarT, IdxT>;
      using DataT   = PhasorDynamics::ForcedOscillationData<typename SourceT::RealT, IdxT>;
      using RealT   = typename SourceT::RealT;

      ForcedOscillationTests()  = default;
      ~ForcedOscillationTests() = default;

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      /// Parameter types and domains, required output linkage, and empty DAE lifecycle.
      TestOutcome validation()
      {
        TestStatus success = true;

        SourceT empty;
        success *= (empty.size() == static_cast<IdxT>(0));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0); // required output assignment is absent

        Fixture configured(makeData());
        success *= configured.initialize();
        success *= configured.output.linked();
        success *= (configured.output.getVariableIndex() == INVALID_INDEX<IdxT>);
        success *= (configured.source.tagDifferentiable() == 0);
        success *= (configured.source.setAbsoluteTolerance(static_cast<RealT>(1.0e-8)) == 0);
        success *= (configured.source.evaluateResidual() == 0);
        success *= (configured.source.evaluateJacobian() == 0);
        success *= (configured.source.nnz() == static_cast<IdxT>(0));

        const std::array<Params, 9> real_parameters{{
            Params::A,
            Params::f,
            Params::Kf,
            Params::Phi,
            Params::Ton,
            Params::Toff,
            Params::Tr,
            Params::Tf,
            Params::Kd,
        }};
        for (const Params parameter : real_parameters)
        {
          success *= rejects(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= rejects(parameter, std::numeric_limits<RealT>::infinity());
        }

        const std::array<Params, 6> nonnegative_parameters{{
            Params::A,
            Params::f,
            Params::Kf,
            Params::Tr,
            Params::Tf,
            Params::Kd,
        }};
        for (const Params parameter : nonnegative_parameters)
        {
          success *= rejects(parameter, static_cast<RealT>(-0.1));
        }

        auto reversed                      = makeData();
        reversed.parameters[Params::Ton]   = static_cast<RealT>(2.0);
        reversed.parameters[Params::Toff]  = static_cast<RealT>(1.0);
        success                           *= (verifyData(reversed) > 0);

        auto no_stop                      = makeData();
        no_stop.parameters[Params::Ton]   = static_cast<RealT>(2.0);
        no_stop.parameters[Params::Toff]  = static_cast<RealT>(-0.5);
        success                          *= (verifyData(no_stop) == 0);

        auto integer_amplitude                   = makeData();
        integer_amplitude.parameters[Params::A]  = static_cast<IdxT>(2);
        success                                 *= (verifyData(integer_amplitude) == 0);

        auto boolean_amplitude                   = makeData();
        boolean_amplitude.parameters[Params::A]  = true;
        success                                 *= (verifyData(boolean_amplitude) > 0);

        return success.report(__func__);
      }

      /// Initialization restores the source-owned output to the documented zero-time value.
      TestOutcome initialization()
      {
        TestStatus success = true;

        auto data                    = makeData();
        data.parameters[Params::A]   = static_cast<RealT>(2.0);
        data.parameters[Params::Phi] = pi() / static_cast<RealT>(2.0);
        data.parameters[Params::Ton] = static_cast<RealT>(-1.0);
        data.parameters[Params::Kd]  = std::log(static_cast<RealT>(2.0));

        Fixture fixture(data);
        success *= (fixture.source.allocate() == 0);
        success *= (fixture.source.verify() == 0);

        fixture.source.updateTime(static_cast<RealT>(1.0), static_cast<RealT>(0.0));
        fixture.output.init(static_cast<ScalarT>(42.0));
        success *= isEqual(fixture.output.read(), static_cast<ScalarT>(42.0), kTol);
        success *= (fixture.source.initialize() == 0);
        success *= isEqual(fixture.output.read(), static_cast<ScalarT>(0.5), kTol);
        success *= (fixture.output.getVariableIndex() == INVALID_INDEX<IdxT>);

        return success.report(__func__);
      }

      /// Time updates reproduce an ordinary sinusoid and a chirped decaying waveform.
      TestOutcome waveform()
      {
        TestStatus success = true;

        auto sine_data                  = makeData();
        sine_data.parameters[Params::A] = static_cast<RealT>(2.0);
        sine_data.parameters[Params::f] = static_cast<RealT>(0.25);

        Fixture sine(sine_data);
        success *= sine.initialize();
        sine.source.updateTime(static_cast<RealT>(1.0), static_cast<RealT>(7.0));
        success *= isEqual(sine.output.read(), static_cast<ScalarT>(2.0), kTol);

        auto chirp_data                    = makeData();
        chirp_data.parameters[Params::A]   = static_cast<RealT>(2.0);
        chirp_data.parameters[Params::f]   = static_cast<RealT>(0.1);
        chirp_data.parameters[Params::Kf]  = static_cast<RealT>(0.2);
        chirp_data.parameters[Params::Ton] = static_cast<RealT>(1.0);
        chirp_data.parameters[Params::Kd]  = static_cast<RealT>(0.5);

        Fixture chirp(chirp_data);
        success *= chirp.initialize();

        const RealT tau      = static_cast<RealT>(2.0);
        const RealT phase    = static_cast<RealT>(1.2) * pi();
        const RealT expected = static_cast<RealT>(2.0)
                               * std::exp(static_cast<RealT>(-0.5) * tau)
                               * std::sin(phase);
        chirp.source.updateTime(static_cast<RealT>(3.0), static_cast<RealT>(0.0));
        success *= isEqual(chirp.output.read(), static_cast<ScalarT>(expected), kTol);

        return success.report(__func__);
      }

      /// Raised-cosine boundaries publish the expected output, envelope, and active flag.
      TestOutcome activationWindowAndMonitors()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(2.0);
        data.parameters[Params::Phi]  = pi() / static_cast<RealT>(2.0);
        data.parameters[Params::Ton]  = static_cast<RealT>(2.0);
        data.parameters[Params::Toff] = static_cast<RealT>(6.0);
        data.parameters[Params::Tr]   = static_cast<RealT>(2.0);
        data.parameters[Params::Tf]   = static_cast<RealT>(2.0);
        data.monitored_variables.insert(Monitor::output);
        data.monitored_variables.insert(Monitor::envelope);
        data.monitored_variables.insert(Monitor::active);

        Fixture fixture(data);
        success *= fixture.initialize();
        success *= (fixture.source.getMonitor() != nullptr);
        success *= (!fixture.source.getMonitor()->empty());

        const std::array<WindowCase, 6> cases{{
            {1.0, 0.0, 0.0, 0.0},
            {2.0, 0.0, 0.0, 1.0},
            {3.0, 1.0, 0.5, 1.0},
            {4.0, 2.0, 1.0, 1.0},
            {5.0, 1.0, 0.5, 1.0},
            {6.0, 0.0, 0.0, 0.0},
        }};

        for (const auto& test_case : cases)
        {
          const RealT time = static_cast<RealT>(test_case.time);
          fixture.source.updateTime(time, static_cast<RealT>(0.0));
          success *= isEqual(fixture.output.read(),
                             static_cast<ScalarT>(test_case.output),
                             kTol);

          const auto values  = monitorValues(fixture.source, time);
          success           *= (values.size() == 4);
          if (values.size() == 4)
          {
            success *= isEqual(values[0], time, kTol);
            success *= isEqual(values[1], static_cast<RealT>(test_case.output), kTol);
            success *= isEqual(values[2], static_cast<RealT>(test_case.envelope), kTol);
            success *= isEqual(values[3], static_cast<RealT>(test_case.active), kTol);
          }
        }

        return success.report(__func__);
      }

      /// DependencyTracking publishes the same exogenous value and invalid DAE index.
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        using Variable = DependencyTracking::Variable;
        using DtSource = PhasorDynamics::ForcedOscillation<Variable, IdxT>;
        using DtSignal = PhasorDynamics::SignalNode<Variable, IdxT>;

        auto data                  = makeData();
        data.parameters[Params::A] = static_cast<RealT>(3.0);
        data.parameters[Params::f] = static_cast<RealT>(0.25);

        DtSource source(data);
        DtSignal output;
        source.getSignals().template assignSignalNode<Internal::OUTPUT>(&output);

        success *= (source.allocate() == 0);
        success *= (source.verify() == 0);
        success *= (source.initialize() == 0);
        source.updateTime(static_cast<RealT>(1.0), static_cast<RealT>(0.0));
        success *= isEqual(output.read().getValue(), static_cast<RealT>(3.0), kTol);
        success *= (output.getVariableIndex() == INVALID_INDEX<IdxT>);
        success *= (source.evaluateResidual() == 0);
        success *= (source.evaluateJacobian() == 0);
        success *= (source.nnz() == static_cast<IdxT>(0));

        return success.report(__func__);
      }

    private:
      using Params   = PhasorDynamics::ForcedOscillationParameters;
      using Internal = PhasorDynamics::ForcedOscillationInternalVariables;
      using Monitor  = PhasorDynamics::ForcedOscillationMonitorableVariables;
      using SignalT  = PhasorDynamics::SignalNode<ScalarT, IdxT>;

      struct Fixture
      {
        explicit Fixture(const DataT& data)
          : source(data)
        {
          source.getSignals().template assignSignalNode<Internal::OUTPUT>(&output);
        }

        bool initialize()
        {
          return source.allocate() == 0
                 && source.verify() == 0
                 && source.initialize() == 0;
        }

        SourceT source;
        SignalT output;
      };

      struct WindowCase
      {
        RealT time;
        RealT output;
        RealT envelope;
        RealT active;
      };

      static RealT pi()
      {
        return std::acos(static_cast<RealT>(-1.0));
      }

      static DataT makeData()
      {
        DataT data;
        data.device_class          = "ForcedOscillation";
        data.disambiguation_string = "fo_test";
        return data;
      }

      static int verifyData(const DataT& data)
      {
        Fixture fixture(data);
        fixture.source.allocate();
        return fixture.source.verify();
      }

      static bool rejects(Params parameter, RealT value)
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        return verifyData(data) > 0;
      }

      static std::vector<RealT> monitorValues(const SourceT& source, const RealT& time)
      {
        Model::VariableMonitorController<ScalarT> controller(time);
        controller.addMonitor(source.getMonitor());

        std::stringstream output;
        controller.addSink({Model::VariableMonitorFormat::CSV}, output);
        controller.print();
        return Tokenizer<RealT>(output.str(), ',')();
      }
    };
  } // namespace Testing
} // namespace GridKit
