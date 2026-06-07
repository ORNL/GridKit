#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/Delay.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/DelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/SignalHistory.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class DelayTests
    {
    private:
      using RealT       = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DelayT      = PhasorDynamics::SignalBlock::Delay<ScalarT, IdxT>;
      using DelayDataT  = typename DelayT::ModelDataT;
      using DelayPorts  = typename DelayDataT::Ports;
      using DelayParams = typename DelayDataT::Parameters;
      using SignalT     = PhasorDynamics::SignalNode<ScalarT, IdxT>;

      using InternalVars = PhasorDynamics::SignalBlock::DelayInternalVariables;
      using ExternalVars = PhasorDynamics::SignalBlock::DelayExternalVariables;

    public:
      /// SignalHistory: linear interpolation, endpoint exactness, prehistory
      /// before the first sample, and hold after the last sample.
      TestOutcome historyInterpolates()
      {
        TestStatus success = true;

        PhasorDynamics::SignalBlock::SignalHistory<RealT> h;
        h.requireWindow(2.0);
        h.record(0.0, 0.0);
        h.record(1.0, 2.0);

        success *= isEqual(h.read(0.5, static_cast<RealT>(99.0)), static_cast<RealT>(1.0));
        success *= isEqual(h.read(0.0, static_cast<RealT>(99.0)), static_cast<RealT>(0.0));
        success *= isEqual(h.read(1.0, static_cast<RealT>(99.0)), static_cast<RealT>(2.0));
        success *= isEqual(h.read(-0.5, static_cast<RealT>(99.0)), static_cast<RealT>(99.0));
        success *= isEqual(h.read(2.0, static_cast<RealT>(99.0)), static_cast<RealT>(2.0));

        return success.report(__func__);
      }

      /// SignalHistory: samples older than the retention window are pruned, but
      /// interpolation across the window edge stays correct.
      TestOutcome historyPrunesToWindow()
      {
        TestStatus success = true;

        PhasorDynamics::SignalBlock::SignalHistory<RealT> h;
        h.requireWindow(1.0);
        h.record(0.0, 0.0);
        h.record(1.0, 1.0);
        h.record(2.0, 2.0); // earliest retained time becomes 1.0; sample at 0.0 pruned

        success *= isEqual(h.read(1.5, static_cast<RealT>(99.0)), static_cast<RealT>(1.5));
        success *= isEqual(h.read(0.5, static_cast<RealT>(99.0)), static_cast<RealT>(99.0));

        return success.report(__func__);
      }

      /// Delay output equals the input shifted by tau, interpolated between
      /// accepted steps.
      TestOutcome delayShiftsInput()
      {
        TestStatus success = true;

        DelayDataT data;
        data.device_class                   = "Delay";
        data.disambiguation_string          = "D";
        data.parameters[DelayParams::delay] = static_cast<RealT>(1.0);
        data.ports[DelayPorts::input]       = 1;
        data.ports[DelayPorts::output]      = 2;

        DelayT  delay(data);
        SignalT in;
        SignalT out;
        ScalarT in_value{5.0};
        IdxT    in_index{3};

        in.set(&in_value, &in_index);
        delay.getSignals().template attachSignalNode<ExternalVars::IN>(&in);
        delay.getSignals().template assignSignalNode<InternalVars::OUT>(&out);
        delay.allocate();
        success *= (delay.verify() == 0);

        delay.resetHistory(0.0); // record (0.0, 5.0); prehistory = 5.0
        in_value = 7.0;
        delay.stepAccepted(1.0); // record (1.0, 7.0); history {(0,5),(1,7)}

        delay.updateTime(1.5, 0.0); // output = history(0.5) interp(0->5,1->7) = 6.0
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(6.0));

        in_value = 9.0;
        delay.stepAccepted(2.0); // record (2.0, 9.0); (0,5) pruned past the window

        delay.updateTime(2.0, 0.0); // output = history(1.0) exact = 7.0
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(7.0));

        delay.updateTime(2.5, 0.0); // output = history(1.5) interp(1->7,2->9) = 8.0
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(8.0));

        return success.report(__func__);
      }

      /// Before t0 + tau the output holds the prehistory value (default = input
      /// at the initial condition).
      TestOutcome delayPrehistoryBeforeStart()
      {
        TestStatus success = true;

        DelayDataT data;
        data.device_class                   = "Delay";
        data.disambiguation_string          = "D";
        data.parameters[DelayParams::delay] = static_cast<RealT>(1.0);
        data.ports[DelayPorts::input]       = 1;
        data.ports[DelayPorts::output]      = 2;

        DelayT  delay(data);
        SignalT in;
        SignalT out;
        ScalarT in_value{5.0};
        IdxT    in_index{3};

        in.set(&in_value, &in_index);
        delay.getSignals().template attachSignalNode<ExternalVars::IN>(&in);
        delay.getSignals().template assignSignalNode<InternalVars::OUT>(&out);
        delay.allocate();
        delay.resetHistory(0.0);

        delay.updateTime(0.5, 0.0); // lookup -0.5 < first sample -> prehistory = 5.0
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(5.0));

        return success.report(__func__);
      }

      /// An explicit prehistory parameter overrides the input-at-t0 default.
      TestOutcome delayConstantPrehistory()
      {
        TestStatus success = true;

        DelayDataT data;
        data.device_class                        = "Delay";
        data.disambiguation_string               = "D";
        data.parameters[DelayParams::delay]      = static_cast<RealT>(1.0);
        data.parameters[DelayParams::prehistory] = static_cast<RealT>(-3.0);
        data.ports[DelayPorts::input]            = 1;
        data.ports[DelayPorts::output]           = 2;

        DelayT  delay(data);
        SignalT in;
        SignalT out;
        ScalarT in_value{5.0};
        IdxT    in_index{3};

        in.set(&in_value, &in_index);
        delay.getSignals().template attachSignalNode<ExternalVars::IN>(&in);
        delay.getSignals().template assignSignalNode<InternalVars::OUT>(&out);
        delay.allocate();
        delay.resetHistory(0.0);

        delay.updateTime(0.5, 0.0); // before start -> constant prehistory -3.0
        success *= isEqual(static_cast<RealT>(out.read()), static_cast<RealT>(-3.0));

        return success.report(__func__);
      }

      /// maxStepSize reports the delay so the solver caps its step there.
      TestOutcome delayMaxStepSize()
      {
        TestStatus success = true;

        DelayDataT data;
        data.parameters[DelayParams::delay] = static_cast<RealT>(0.25);
        data.ports[DelayPorts::input]       = 1;
        data.ports[DelayPorts::output]      = 2;

        DelayT delay(data);
        success *= isEqual(delay.maxStepSize(), static_cast<RealT>(0.25));

        return success.report(__func__);
      }

      /// A non-positive delay is a configuration error.
      TestOutcome delayRejectsNonpositiveDelay()
      {
        TestStatus success = true;

        DelayDataT data;
        data.parameters[DelayParams::delay] = static_cast<RealT>(0.0);
        data.ports[DelayPorts::input]       = 1;
        data.ports[DelayPorts::output]      = 2;

        DelayT  delay(data);
        SignalT in;
        SignalT out;
        ScalarT in_value{1.0};
        IdxT    in_index{3};

        in.set(&in_value, &in_index);
        delay.getSignals().template attachSignalNode<ExternalVars::IN>(&in);
        delay.getSignals().template assignSignalNode<InternalVars::OUT>(&out);
        delay.allocate();

        success *= (delay.verify() != 0);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
