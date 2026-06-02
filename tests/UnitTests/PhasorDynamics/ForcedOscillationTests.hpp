/**
 * @file ForcedOscillationTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Unit tests for the forced oscillation signal operator.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <variant>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOperator/ForcedOscillation/ForcedOscillation.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOperator/ForcedOscillation/ForcedOscillationData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ForcedOscillationTests
    {
    public:
      using RealT  = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DataT  = PhasorDynamics::SignalOperator::ForcedOscillationData<RealT, IdxT>;
      using ModelT = PhasorDynamics::SignalOperator::ForcedOscillation<ScalarT, IdxT>;

      ForcedOscillationTests()  = default;
      ~ForcedOscillationTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto data = makeData();
        data.monitored_variables.insert(Mon::out);

        auto* forced  = new ModelT(data);
        success      *= (forced != nullptr);
        success      *= (forced->getMonitor() != nullptr);

        delete forced;

        return success.report(__func__);
      }

      TestOutcome source_mode()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(2.0);
        data.parameters[Params::f]    = static_cast<RealT>(0.25);
        data.parameters[Params::Bias] = static_cast<RealT>(0.1);
        data.parameters[Params::u0]   = static_cast<RealT>(0.3);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        auto&            forced = harness.forced;

        forced.allocate();
        success *= (forced.verify() == 0);
        forced.initialize();

        const RealT t        = 1.0;
        const RealT expected = expectedOutput(data, t, 0.3);
        forced.updateTime(t, 0.0);
        forced.y()[0] = expected;
        forced.evaluateResidual();

        success *= isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);
        success *= isEqual(harness.output_node.read(), static_cast<ScalarT>(expected), kTol);

        return success.report(__func__);
      }

      TestOutcome additive_mode()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(0.5);
        data.parameters[Params::f]    = static_cast<RealT>(0.25);
        data.parameters[Params::Bias] = static_cast<RealT>(-0.1);
        data.parameters[Params::Kin]  = static_cast<RealT>(2.0);
        data.parameters[Params::u0]   = static_cast<RealT>(99.0);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        harness.attachInput(0.7);
        auto& forced = harness.forced;

        forced.allocate();
        success *= (forced.verify() == 0);
        forced.initialize();

        const RealT t        = 1.0;
        const RealT expected = expectedOutput(data, t, 0.7);
        forced.updateTime(t, 0.0);
        forced.y()[0] = expected;
        forced.evaluateResidual();

        success *= isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);
        success *= isEqual(harness.output_node.read(), static_cast<ScalarT>(expected), kTol);

        return success.report(__func__);
      }

      TestOutcome delayed_additive_identity()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(0.1);
        data.parameters[Params::f]    = static_cast<RealT>(0.7);
        data.parameters[Params::Ton]  = static_cast<RealT>(10.0);
        data.parameters[Params::Tr]   = static_cast<RealT>(1.0);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);
        data.monitored_variables.insert(Mon::in);
        data.monitored_variables.insert(Mon::force);
        data.monitored_variables.insert(Mon::out);
        data.monitored_variables.insert(Mon::active);

        ComponentHarness harness(data);
        harness.attachInput(0.0);
        auto& forced = harness.forced;

        forced.allocate();
        success *= (forced.verify() == 0);
        forced.initialize();

        forced.y()[0] = static_cast<ScalarT>(0.7);
        forced.initializeInputFromOutput();
        success *= isEqual(harness.input_value, static_cast<ScalarT>(0.7), kTol);

        const RealT t = 5.0;
        forced.updateTime(t, 0.0);
        forced.y()[0] = static_cast<ScalarT>(0.7);
        forced.evaluateResidual();

        success *= isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);
        success *= isEqual(harness.output_node.read(), static_cast<ScalarT>(0.7), kTol);

        Model::VariableMonitorController<ScalarT> monitor(t);
        monitor.addMonitor(forced.getMonitor());

        std::stringstream os;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, os);
        monitor.print();

        auto values = Tokenizer<RealT>(os.str(), ',')();
        if (values.size() == 5)
        {
          success *= isEqual(values[1], static_cast<RealT>(0.7), kRealTol);
          success *= isEqual(values[2], static_cast<RealT>(0.0), kRealTol);
          success *= isEqual(values[3], static_cast<RealT>(0.7), kRealTol);
          success *= isEqual(values[4], static_cast<RealT>(0.0), kRealTol);
        }
        else
        {
          std::cout << "Unexpected monitor value count: " << values.size() << "\n";
          success = false;
        }

        return success.report(__func__);
      }

      TestOutcome consistent_initial_condition()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(1.0);
        data.parameters[Params::f]    = static_cast<RealT>(0.25);
        data.parameters[Params::Bias] = static_cast<RealT>(0.2);
        data.parameters[Params::Kin]  = static_cast<RealT>(3.0);
        data.parameters[Params::u0]   = static_cast<RealT>(-1.0);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        harness.attachInput(0.4);
        auto& forced = harness.forced;

        forced.allocate();
        success *= (forced.verify() == 0);
        forced.updateTime(1.0, 0.0);
        forced.initialize();
        forced.evaluateResidual();
        forced.tagDifferentiable();

        const RealT expected  = expectedOutput(data, 1.0, 0.4);
        success              *= isEqual(forced.y()[0], static_cast<ScalarT>(expected), kTol);
        success              *= isEqual(forced.yp()[0], static_cast<ScalarT>(0.0), kTol);
        success              *= isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);
        success              *= (!forced.tag()[0]);

        return success.report(__func__);
      }

      TestOutcome sine_chirp_decay()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(2.0);
        data.parameters[Params::f]    = static_cast<RealT>(0.1);
        data.parameters[Params::Kf]   = static_cast<RealT>(0.2);
        data.parameters[Params::Kd]   = static_cast<RealT>(0.5);
        data.parameters[Params::Ton]  = static_cast<RealT>(1.0);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        auto&            forced = harness.forced;

        forced.allocate();
        forced.initialize();

        const RealT t        = 3.0;
        const RealT expected = expectedOutput(data, t, 0.0);
        forced.updateTime(t, 0.0);
        forced.y()[0] = expected;
        forced.evaluateResidual();

        success *= isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);

        return success.report(__func__);
      }

      TestOutcome raised_cosine_window()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(1.0);
        data.parameters[Params::Phi]  = static_cast<RealT>(0.5 * kPi);
        data.parameters[Params::Ton]  = static_cast<RealT>(2.0);
        data.parameters[Params::Toff] = static_cast<RealT>(6.0);
        data.parameters[Params::Tr]   = static_cast<RealT>(2.0);
        data.parameters[Params::Tf]   = static_cast<RealT>(2.0);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        auto&            forced = harness.forced;

        forced.allocate();
        forced.initialize();

        success *= checkOutput(forced, data, 1.0, 0.0);
        success *= checkOutput(forced, data, 3.0, 0.0);
        success *= checkOutput(forced, data, 5.0, 0.0);
        success *= checkOutput(forced, data, 6.0, 0.0);

        success *= isEqual(expectedEnvelope(data, 1.0), static_cast<RealT>(0.0), kRealTol);
        success *= isEqual(expectedEnvelope(data, 3.0), static_cast<RealT>(0.5), kRealTol);
        success *= isEqual(expectedEnvelope(data, 5.0), static_cast<RealT>(0.5), kRealTol);
        success *= isEqual(expectedEnvelope(data, 6.0), static_cast<RealT>(0.0), kRealTol);

        return success.report(__func__);
      }

      TestOutcome clamp_defaults_and_limits()
      {
        TestStatus success = true;

        auto wide_data                     = makeData();
        wide_data.parameters[Params::Bias] = static_cast<RealT>(2.0);

        ComponentHarness wide_harness(wide_data);
        auto&            wide = wide_harness.forced;
        wide.allocate();
        wide.initialize();
        success *= isEqual(wide.y()[0], static_cast<ScalarT>(2.0), static_cast<RealT>(1.0e-8));

        auto limited_data                     = makeData();
        limited_data.parameters[Params::Bias] = static_cast<RealT>(2.0);
        limited_data.parameters[Params::Lmin] = static_cast<RealT>(-0.5);
        limited_data.parameters[Params::Lmax] = static_cast<RealT>(0.75);

        ComponentHarness limited_harness(limited_data);
        auto&            limited = limited_harness.forced;
        limited.allocate();
        limited.initialize();
        success *= isEqual(limited.y()[0], static_cast<ScalarT>(0.75), static_cast<RealT>(1.0e-8));

        return success.report(__func__);
      }

      TestOutcome monitor_values()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::A]    = static_cast<RealT>(1.0);
        data.parameters[Params::f]    = static_cast<RealT>(0.25);
        data.parameters[Params::Bias] = static_cast<RealT>(0.1);
        data.parameters[Params::u0]   = static_cast<RealT>(0.2);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);
        data.monitored_variables.insert(Mon::in);
        data.monitored_variables.insert(Mon::env);
        data.monitored_variables.insert(Mon::force);
        data.monitored_variables.insert(Mon::vraw);
        data.monitored_variables.insert(Mon::out);
        data.monitored_variables.insert(Mon::active);

        ComponentHarness harness(data);
        auto&            forced = harness.forced;
        forced.allocate();
        forced.initialize();

        const RealT t        = 1.0;
        const RealT expected = expectedOutput(data, t, 0.2);
        forced.updateTime(t, 0.0);
        forced.y()[0] = expected;
        forced.evaluateResidual();

        Model::VariableMonitorController<ScalarT> monitor(t);
        monitor.addMonitor(forced.getMonitor());

        std::stringstream os;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, os);
        monitor.print();

        auto values = Tokenizer<RealT>(os.str(), ',')();
        if (values.size() == 7)
        {
          success *= isEqual(values[1], static_cast<RealT>(0.2), kRealTol);
          success *= isEqual(values[2], static_cast<RealT>(1.0), kRealTol);
          success *= isEqual(values[3], static_cast<RealT>(1.0), kRealTol);
          success *= isEqual(values[4], static_cast<RealT>(1.3), kRealTol);
          success *= isEqual(values[5], static_cast<RealT>(1.3), kRealTol);
          success *= isEqual(values[6], static_cast<RealT>(1.0), kRealTol);
        }
        else
        {
          std::cout << "Unexpected monitor value count: " << values.size() << "\n";
          success = false;
        }

        return success.report(__func__);
      }

      TestOutcome validation()
      {
        TestStatus success = true;

        auto   missing_output = makeData();
        ModelT no_output(missing_output);
        success *= (no_output.verify() > 0);

        auto   data = makeData();
        ModelT forced(data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> input_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> output_node;

        forced.getSignals().template attachSignalNode<Ext::IN>(&input_node);
        forced.getSignals().template assignSignalNode<Int::OUT>(&output_node);
        success *= (forced.verify() > 0);

        auto invalid                     = makeData();
        invalid.parameters[Params::A]    = static_cast<RealT>(-1.0);
        invalid.parameters[Params::f]    = static_cast<RealT>(-1.0);
        invalid.parameters[Params::Kf]   = static_cast<RealT>(-1.0);
        invalid.parameters[Params::Tr]   = static_cast<RealT>(-1.0);
        invalid.parameters[Params::Tf]   = static_cast<RealT>(-1.0);
        invalid.parameters[Params::Kd]   = static_cast<RealT>(-1.0);
        invalid.parameters[Params::Ton]  = static_cast<RealT>(5.0);
        invalid.parameters[Params::Toff] = static_cast<RealT>(4.0);
        invalid.parameters[Params::Lmin] = static_cast<RealT>(2.0);
        invalid.parameters[Params::Lmax] = static_cast<RealT>(1.0);

        ComponentHarness invalid_harness(invalid);
        success *= (invalid_harness.forced.verify() >= 8);

        return success.report(__func__);
      }

      TestOutcome json_and_system_wiring()
      {
        TestStatus success = true;

        std::istringstream input(R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "forced oscillation source",
    "case_description": "ForcedOscillation parser test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "infinite_bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
  ],
  "signals": [
    { "signal_id": 9001, "name": "Forced Oscillation Output" }
  ],
  "devices": [
    {
      "class": "ForcedOscillation",
      "id": "FO1",
      "ports": { "output": 9001 },
      "params": {
        "A": 0.5,
        "f": 0.25,
        "Kf": 0.1,
        "Phi": 0.0,
        "Bias": 0.1,
        "Kin": 1.0,
        "u0": 0.2,
        "Ton": 0.0,
        "Toff": -1.0,
        "Tr": 0.0,
        "Tf": 0.0,
        "Kd": 0.0,
        "Lmin": -10.0,
        "Lmax": 10.0
      },
      "mon": ["force", "out", "active"]
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.forced_oscillation.size() == 1);
        success   *= (data.forced_oscillation[0].ports.at(Ports::output) == 9001);
        success   *= (std::get_if<RealT>(&data.forced_oscillation[0].parameters.at(Params::A)) != nullptr);
        success   *= data.forced_oscillation[0].monitored_variables.contains(Mon::force);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.tagDifferentiable() == 0);
        system.updateTime(1.0, 0.0);
        success *= (system.evaluateResidual() == 0);
#ifdef GRIDKIT_ENABLE_ENZYME
        success *= (system.evaluateJacobian() == 0);
#endif
        success *= (system.size() == 1);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        success *= checkJacobian(false);
        success *= checkJacobian(true);

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::SignalOperator::ForcedOscillationParameters;
      using Ports  = PhasorDynamics::SignalOperator::ForcedOscillationPorts;
      using Mon    = PhasorDynamics::SignalOperator::ForcedOscillationMonitorableVariables;
      using Int    = PhasorDynamics::SignalOperator::ForcedOscillationInternalVariables;
      using Ext    = PhasorDynamics::SignalOperator::ForcedOscillationExternalVariables;

      static constexpr RealT kPi      = static_cast<RealT>(3.141592653589793238462643383279502884L);
      static constexpr RealT kRealTol = static_cast<RealT>(1.0e-9);
      static constexpr RealT kTol     = static_cast<RealT>(1.0e-8);

      struct ComponentHarness
      {
        explicit ComponentHarness(const DataT& data)
          : forced(data)
        {
          output_node.set(&output_value, &output_index);
          forced.getSignals().template assignSignalNode<Int::OUT>(&output_node);
        }

        void attachInput(RealT value, IdxT index = 7)
        {
          input_value = static_cast<ScalarT>(value);
          input_index = index;
          input_node.set(&input_value, &input_index);
          forced.getSignals().template attachSignalNode<Ext::IN>(&input_node);
        }

        ModelT forced;

        ScalarT                                   output_value{0.0};
        IdxT                                      output_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<ScalarT, IdxT> output_node;

        ScalarT                                   input_value{0.0};
        IdxT                                      input_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<ScalarT, IdxT> input_node;
      };

      static auto makeData() -> DataT
      {
        DataT data;
        data.device_class          = "ForcedOscillation";
        data.disambiguation_string = "fo_test";
        return data;
      }

      static RealT parameter(const DataT& data, Params key, RealT fallback)
      {
        if (!data.parameters.contains(key))
        {
          return fallback;
        }

        const auto& value = data.parameters.at(key);
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          return *real_value;
        }
        if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          return static_cast<RealT>(*index_value);
        }

        return fallback;
      }

      static RealT expectedEnvelope(const DataT& data, RealT time)
      {
        const RealT Ton  = parameter(data, Params::Ton, 0.0);
        const RealT Toff = parameter(data, Params::Toff, -1.0);
        const RealT Tr   = parameter(data, Params::Tr, 0.0);
        const RealT Tf   = parameter(data, Params::Tf, 0.0);

        bool active = (time >= Ton);
        if (Toff >= 0.0 && time >= Toff)
        {
          active = false;
        }

        RealT envelope = 0.0;
        if (active)
        {
          envelope = 1.0;

          if (Tr > 0.0 && time < Ton + Tr)
          {
            const RealT x  = (time - Ton) / Tr;
            envelope      *= 0.5 * (1.0 - std::cos(kPi * x));
          }

          if (Toff >= 0.0 && Tf > 0.0 && time > Toff - Tf)
          {
            const RealT x  = (Toff - time) / Tf;
            envelope      *= 0.5 * (1.0 - std::cos(kPi * x));
          }
        }

        return envelope;
      }

      static RealT expectedForce(const DataT& data, RealT time)
      {
        const RealT A    = parameter(data, Params::A, 0.0);
        const RealT f    = parameter(data, Params::f, 0.0);
        const RealT Kf   = parameter(data, Params::Kf, 0.0);
        const RealT Phi  = parameter(data, Params::Phi, 0.0);
        const RealT Ton  = parameter(data, Params::Ton, 0.0);
        const RealT Kd   = parameter(data, Params::Kd, 0.0);
        const RealT tau  = std::max(time - Ton, static_cast<RealT>(0.0));
        const RealT env  = expectedEnvelope(data, time);
        const RealT ph   = Phi + 2.0 * kPi * (f * tau + 0.5 * Kf * tau * tau);
        const RealT damp = std::exp(-Kd * tau);
        return A * env * damp * std::sin(ph);
      }

      static RealT expectedOutput(const DataT& data, RealT time, RealT input)
      {
        const RealT Bias = parameter(data, Params::Bias, 0.0);
        const RealT Kin  = parameter(data, Params::Kin, 1.0);
        const RealT Lmin = parameter(data, Params::Lmin, -1.0e6);
        const RealT Lmax = parameter(data, Params::Lmax, 1.0e6);
        const RealT raw  = Kin * input + Bias + expectedForce(data, time);
        return Math::clamp(raw, Lmin, Lmax);
      }

      static bool checkOutput(ModelT& forced, const DataT& data, RealT time, RealT input)
      {
        const RealT expected = expectedOutput(data, time, input);
        forced.updateTime(time, 0.0);
        forced.y()[0] = expected;
        forced.evaluateResidual();
        return isEqual(forced.getResidual()[0], static_cast<ScalarT>(0.0), kTol);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      bool checkJacobian(bool attach_input)
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::Bias] = static_cast<RealT>(0.1);
        data.parameters[Params::Kin]  = static_cast<RealT>(1.5);
        data.parameters[Params::u0]   = static_cast<RealT>(0.3);
        data.parameters[Params::Lmin] = static_cast<RealT>(-10.0);
        data.parameters[Params::Lmax] = static_cast<RealT>(10.0);

        ComponentHarness harness(data);
        if (attach_input)
        {
          harness.attachInput(0.4, 7);
        }

        auto& forced = harness.forced;
        forced.allocate();
        forced.initialize();
        forced.y()[0] = 0.8;
        forced.evaluateResidual();

        const RealT base_y     = forced.y()[0];
        const RealT base_input = static_cast<RealT>(harness.input_value);
        const RealT eps        = static_cast<RealT>(1.0e-6);

        forced.y()[0] = base_y + eps;
        forced.evaluateResidual();
        const RealT f_plus_y = static_cast<RealT>(forced.getResidual()[0]);

        forced.y()[0] = base_y - eps;
        forced.evaluateResidual();
        const RealT f_minus_y = static_cast<RealT>(forced.getResidual()[0]);

        const RealT fd_y = (f_plus_y - f_minus_y) / (2.0 * eps);

        RealT fd_input = 0.0;
        if (attach_input)
        {
          forced.y()[0]       = base_y;
          harness.input_value = base_input + eps;
          forced.evaluateResidual();
          const RealT f_plus_u = static_cast<RealT>(forced.getResidual()[0]);

          harness.input_value = base_input - eps;
          forced.evaluateResidual();
          const RealT f_minus_u = static_cast<RealT>(forced.getResidual()[0]);
          fd_input              = (f_plus_u - f_minus_u) / (2.0 * eps);

          harness.input_value = base_input;
        }

        forced.y()[0] = base_y;
        forced.evaluateResidual();
        forced.evaluateJacobian();

        auto jacobian = forced.getJacobian();
        jacobian.deduplicate();
        auto map = MapFromCOO(jacobian);

        success *= (!map.empty());
        if (!map.empty())
        {
          success *= isEqual(map[0][0], fd_y, static_cast<RealT>(1.0e-6));

          if (attach_input)
          {
            success *= map[0].contains(7);
            success *= isEqual(map[0][7], fd_input, static_cast<RealT>(1.0e-6));
          }
          else
          {
            success *= (!map[0].contains(INVALID_INDEX<IdxT>));
            success *= (map[0].size() == 1);
          }
        }

        return success;
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
