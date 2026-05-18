#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/Reeca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/ReecaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ConverterReecaTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterReecaTests()  = default;
      ~ConverterReecaTests() = default;

      static constexpr ScalarT kTol         = static_cast<ScalarT>(1.0e-10);
      static constexpr ScalarT kResidualTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);
        success *= (minimal.verify() > 0);

        auto                                            data = makeTestData();
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> from_data(&bus, data);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(from_data);

        success *= (from_data.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (from_data.getMonitor() != nullptr);
        success *= (from_data.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        auto missing = makeTestData();
        missing.parameters.erase(Params::Tiq);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> missing_model(&bus, missing);
        success *= verifyWithRequiredSignals(missing_model);

        auto bad_switch                    = makeTestData();
        bad_switch.parameters[Params::spf] = static_cast<IdxT>(2);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_switch_model(&bus, bad_switch);
        success *= verifyWithRequiredSignals(bad_switch_model);

        success *= invalidParameterCase(bus, Params::Sbase, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Trv, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Tp, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vdip, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vup, static_cast<RealT>(0.8));
        success *= invalidParameterCase(bus, Params::Dbd1, static_cast<RealT>(0.1));
        success *= invalidParameterCase(bus, Params::Dbd2, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Iql1, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Thld, static_cast<RealT>(1.0));
        success *= invalidParameterCase(bus, Params::Thld2, static_cast<RealT>(1.0));
        success *= invalidParameterCase(bus, Params::Qmin, static_cast<RealT>(1.0));
        success *= invalidParameterCase(bus, Params::Vmin, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Tiq, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Tpord, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::dPmin, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::dPmax, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Pmin, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Imax, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::vq2, static_cast<RealT>(0.5));
        success *= invalidParameterCase(bus, Params::lq1, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::vp2, static_cast<RealT>(0.5));
        success *= invalidParameterCase(bus, Params::lp1, static_cast<RealT>(-0.1));

        auto                                            valid = makeTestData();
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> valid_model(&bus, valid);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(valid_model);
        success *= (valid_model.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome initialization()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeTestData());
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(reeca);
        signals.assignOutputs(reeca);

        success *= (reeca.allocate() == 0);
        success *= (signals.iqcmd_node.linked());
        success *= (signals.ipcmd_node.linked());
        success *= (reeca.initialize() == 0);

        success *= isEqual(reeca.y()[index(Vars::VMEAS)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reeca.y()[index(Vars::PMEAS)], signals.pe, kTol);
        success *= isEqual(reeca.y()[index(Vars::VT)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reeca.y()[index(Vars::QREF)], signals.qgen, kTol);
        success *= isEqual(reeca.y()[index(Vars::PORD)], signals.pe, kTol);
        success *= isEqual(reeca.y()[index(Vars::IQCMD)], signals.qgen, kResidualTol);
        success *= isEqual(reeca.y()[index(Vars::IPCMD)], signals.pe, kResidualTol);
        success *= isEqual(signals.iqcmd_node.read(), signals.qgen, kResidualTol);
        success *= isEqual(signals.ipcmd_node.read(), signals.pe, kResidualTol);

        for (size_t i = 0; i < static_cast<size_t>(reeca.size()); ++i)
        {
          success *= isEqual(reeca.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeTestData());
        SignalFixture                                   signals;
        signals.qext = signals.qgen;
        signals.linkAllInputs();
        signals.attachAllInputs(reeca);
        signals.assignOutputs(reeca);

        reeca.allocate();
        success *= (reeca.initialize() == 0);
        success *= (reeca.evaluateResidual() == 0);

        success *= (reeca.getResidual().size() == static_cast<size_t>(Vars::MAXIMUM));
        for (const auto& residual : reeca.getResidual())
        {
          success *= std::isfinite(static_cast<RealT>(residual));
          success *= isEqual(residual, static_cast<ScalarT>(0.0), kResidualTol);
        }

        return success.report(__func__);
      }

      TestOutcome zeroTimeConstantTags()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeTestData();
        data.parameters[Params::Trv] = static_cast<RealT>(0.0);
        data.parameters[Params::Tp]  = static_cast<RealT>(0.0);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(reeca);

        reeca.allocate();
        success *= (reeca.tagDifferentiable() == 0);
        success *= (!reeca.tag()[index(Vars::VMEAS)]);
        success *= (!reeca.tag()[index(Vars::PMEAS)]);
        success *= (!reeca.tag()[index(Vars::VT)]);
        success *= (reeca.tag()[index(Vars::XPIQ)]);
        success *= (reeca.tag()[index(Vars::XPIV)]);
        success *= (reeca.tag()[index(Vars::QV)]);
        success *= (reeca.tag()[index(Vars::PORD)]);
        success *= (reeca.initialize() == 0);
        success *= (reeca.evaluateResidual() == 0);
        success *= residualsAreNearZero(reeca);

        return success.report(__func__);
      }

      TestOutcome outputAvailability()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        auto                               data = makeTestData();
        addAllMonitors(data);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(reeca);

        success *= (reeca.getMonitor() != nullptr);
        success *= (!reeca.getMonitor()->empty());

        return success.report(__func__);
      }

      TestOutcome signalVerification()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeTestData());
        SignalFixture                                   signals;

        success *= (reeca.verify() > 0);

        reeca.getSignals().template attachSignalNode<Ext::PE>(&signals.pe_node);
        success *= (reeca.verify() > 0);

        signals.linkPe();
        success *= (reeca.verify() > 0);

        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&signals.qgen_node);
        success *= (reeca.verify() > 0);

        signals.linkQgen();
        success *= (reeca.verify() == 0);

        reeca.getSignals().template attachSignalNode<Ext::OMEGA>(&signals.omega_node);
        success *= (reeca.verify() > 0);
        signals.linkOmega();
        success *= (reeca.verify() == 0);

        auto                                            spf_data = makeTestData(true);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> spf_model(&bus, spf_data);
        SignalFixture                                   spf_signals;
        spf_signals.linkRequired();
        spf_signals.attachRequired(spf_model);
        success *= (spf_model.verify() > 0);

        spf_model.getSignals().template attachSignalNode<Ext::PFAREF>(&spf_signals.pfaref_node);
        success *= (spf_model.verify() > 0);
        spf_signals.linkPfaref();
        success *= (spf_model.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome deviceInitFallbacks()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto init_data                      = makeTestData();
        init_data.parameters[Params::Sbase] = static_cast<RealT>(50.0);
        relaxCurrentLimits(init_data, static_cast<RealT>(2.0));
        init_data.external_initial_values[Ext::PE]   = static_cast<RealT>(0.8);
        init_data.external_initial_values[Ext::QGEN] = static_cast<RealT>(0.2);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> init_model(&bus, init_data);
        success *= (init_model.verify() == 0);
        success *= (init_model.allocate() == 0);
        success *= (init_model.initialize() == 0);
        success *= isEqual(init_model.y()[index(Vars::VT)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(init_model.y()[index(Vars::PMEAS)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(init_model.y()[index(Vars::QREF)], static_cast<ScalarT>(0.4), kTol);
        success *= isEqual(init_model.y()[index(Vars::PORD)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(init_model.y()[index(Vars::IQCMD)], static_cast<ScalarT>(0.4), kResidualTol);
        success *= isEqual(init_model.y()[index(Vars::IPCMD)], static_cast<ScalarT>(1.6), kResidualTol);
        success *= (init_model.evaluateResidual() == 0);
        success *= residualsAreNearZero(init_model);

        auto linked_data                             = init_data;
        linked_data.external_initial_values[Ext::PE] = static_cast<RealT>(0.1);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> linked_model(&bus, linked_data);
        SignalFixture                                   signals;
        signals.pe   = static_cast<ScalarT>(0.8);
        signals.qgen = static_cast<ScalarT>(0.2);
        signals.linkRequired();
        signals.attachRequired(linked_model);

        success *= (linked_model.verify() == 0);
        success *= (linked_model.allocate() == 0);
        success *= (linked_model.initialize() == 0);
        success *= isEqual(linked_model.y()[index(Vars::VT)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(linked_model.y()[index(Vars::PMEAS)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(linked_model.y()[index(Vars::QREF)], static_cast<ScalarT>(0.4), kTol);

        auto missing_qgen                             = makeTestData();
        missing_qgen.external_initial_values[Ext::PE] = static_cast<RealT>(0.8);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> missing_qgen_model(&bus, missing_qgen);
        success *= (missing_qgen_model.verify() > 0);

        auto bad_init_key                                  = init_data;
        bad_init_key.external_initial_values[Ext::MAXIMUM] = static_cast<RealT>(0.8);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_init_key_model(&bus, bad_init_key);
        success *= (bad_init_key_model.verify() > 0);

        auto bad_init_type                             = init_data;
        bad_init_type.external_initial_values[Ext::PE] = true;
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_init_type_model(&bus, bad_init_type);
        success *= (bad_init_type_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome priorityInitialization()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeTestData();
        data.parameters[Params::sPQ] = true;
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(reeca);
        signals.assignOutputs(reeca);

        success *= (reeca.allocate() == 0);
        success *= (reeca.initialize() == 0);
        success *= isEqual(reeca.y()[index(Vars::VT)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reeca.y()[index(Vars::IPCMD)], signals.pe, kResidualTol);
        success *= isEqual(reeca.y()[index(Vars::IQCMD)], signals.qgen, kResidualTol);
        success *= isEqual(reeca.y()[index(Vars::IPCIRC)], static_cast<ScalarT>(1.2), kTol);

        const auto circle =
            reeca.y()[index(Vars::IQCIRC)] * reeca.y()[index(Vars::IQCIRC)]
            + reeca.y()[index(Vars::IPCMD)] * reeca.y()[index(Vars::IPCMD)];
        success *= isEqual(circle, static_cast<ScalarT>(1.2 * 1.2), kResidualTol);

        success *= (reeca.evaluateResidual() == 0);
        success *= residualsAreNearZero(reeca);

        return success.report(__func__);
      }

      TestOutcome unsupportedInitializationRejects()
      {
        TestStatus success = true;

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeTestData());
          SignalFixture                                   signals;
          signals.linkRequired();
          signals.attachRequired(reeca);
          reeca.allocate();
          success *= (reeca.initialize() != 0);
        }

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
          bus.allocate();
          bus.initialize();

          auto data                     = makeTestData();
          data.parameters[Params::Pmax] = static_cast<RealT>(2.0);
          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
          SignalFixture                                   signals;
          signals.pe = static_cast<ScalarT>(1.1);
          signals.linkRequired();
          signals.attachRequired(reeca);
          reeca.allocate();
          success *= (reeca.initialize() != 0);
        }

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeTestData());
          SignalFixture                                   signals;
          signals.qext = signals.qgen + static_cast<ScalarT>(0.02);
          signals.linkAllInputs();
          signals.attachAllInputs(reeca);
          reeca.allocate();
          success *= (reeca.initialize() != 0);
        }

        return success.report(__func__);
      }

      TestOutcome jsonParseAndSystemAssembly()
      {
        TestStatus success = true;

        std::istringstream input(validSystemJson());
        auto               data = PhasorDynamics::parseSystemModelData(input);

        success *= (data.reeca.size() == 1);
        success *= (data.reeca[0].device_class == "Reeca");
        success *= (data.reeca[0].ports.at(Ports::pe) == 10);
        success *= (std::get_if<size_t>(&data.reeca[0].parameters.at(Params::Sbase)) != nullptr);
        success *= (std::get_if<bool>(&data.reeca[0].parameters.at(Params::spf)) != nullptr);
        success *= data.reeca[0].external_initial_values.empty();

        PhasorDynamics::SystemModel<double, size_t> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.tagDifferentiable() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.evaluateJacobian() == 0);
        success *= (system.size() == 2 + 2 * 3 + static_cast<size_t>(Vars::MAXIMUM));

        for (const auto& residual : system.getResidual())
        {
          success *= std::isfinite(residual);
        }

        std::istringstream invalid_key_input(
            initFallbackSystemJsonWithInit(R"json("init": { "pmeas": 0.8, "qgen": 0.2 })json"));
        auto invalid_key_data  = PhasorDynamics::parseSystemModelData(invalid_key_input);
        success               *= (invalid_key_data.reeca[0].input_error_count > 0);

        PhasorDynamics::Bus<double, size_t>              bad_key_bus(1.0, 0.0);
        PhasorDynamics::Converter::Reeca<double, size_t> bad_key_model(&bad_key_bus, invalid_key_data.reeca[0]);
        success *= (bad_key_model.verify() > 0);

        std::istringstream invalid_type_input(
            initFallbackSystemJsonWithInit(R"json("init": { "pe": "invalid", "qgen": 0.2 })json"));
        auto invalid_type_data  = PhasorDynamics::parseSystemModelData(invalid_type_input);
        success                *= (invalid_type_data.reeca[0].input_error_count > 0);

        PhasorDynamics::Bus<double, size_t>              bad_type_bus(1.0, 0.0);
        PhasorDynamics::Converter::Reeca<double, size_t> bad_type_model(&bad_type_bus, invalid_type_data.reeca[0]);
        success *= (bad_type_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome initFallbackSystemAssembly()
      {
        TestStatus success = true;

        std::istringstream input(initFallbackSystemJson());
        auto               data = PhasorDynamics::parseSystemModelData(input);

        success *= (data.reeca.size() == 1);
        success *= (data.reeca[0].ports.contains(Ports::pe) == false);
        success *= (data.reeca[0].ports.contains(Ports::qgen) == false);
        success *= (data.reeca[0].input_error_count == 0);
        success *= (std::get_if<double>(&data.reeca[0].external_initial_values.at(Ext::PE)) != nullptr);
        success *= (std::get_if<double>(&data.reeca[0].external_initial_values.at(Ext::QGEN)) != nullptr);

        PhasorDynamics::SystemModel<double, size_t> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.tagDifferentiable() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.evaluateJacobian() == 0);
        success *= (system.size() == 2 + static_cast<size_t>(Vars::MAXIMUM));

        for (const auto& residual : system.getResidual())
        {
          success *= std::isfinite(residual);
        }

        return success.report(__func__);
      }

      TestOutcome systemRejectsUnlinkedRequiredSignals()
      {
        TestStatus success = true;

        std::istringstream input(unlinkedSystemJson());
        auto               data = PhasorDynamics::parseSystemModelData(input);

        bool threw = false;
        try
        {
          PhasorDynamics::SystemModel<double, size_t> system(data);
          system.allocate();
        }
        catch (const std::runtime_error&)
        {
          threw = true;
        }
        success *= threw;

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto dependency_tracking_jacobian = DependencyTrackingJacobian();
        auto enzyme_jacobian              = EnzymeJacobian();

        success          *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        const auto nrows  = std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());
        for (size_t i = 0; i < nrows; ++i)
        {
          const auto dependency_tracking_row =
              pruneJacobianTails(dependency_tracking_jacobian[i], static_cast<RealT>(1.0e-8));
          const auto enzyme_row =
              pruneJacobianTails(enzyme_jacobian[i], static_cast<RealT>(1.0e-8));
          success *= isEqual(dependency_tracking_row, enzyme_row, static_cast<RealT>(1.0e-8));
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Converter::ReecaParameters;
      using Ports  = PhasorDynamics::Converter::ReecaPorts;
      using Vars   = PhasorDynamics::Converter::ReecaInternalVariables;
      using Ext    = PhasorDynamics::Converter::ReecaExternalVariables;
      using Mon    = PhasorDynamics::Converter::ReecaMonitorableVariables;

      static size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      struct SignalFixture
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qext_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pfaref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;

        ScalarT pe{0.75};
        ScalarT qgen{0.20};
        ScalarT omega{0.01};
        ScalarT qext{0.22};
        ScalarT pfaref{0.0};
        ScalarT pref{0.76};

        IdxT pe_index{101};
        IdxT qgen_index{102};
        IdxT omega_index{103};
        IdxT qext_index{104};
        IdxT pfaref_index{105};
        IdxT pref_index{106};

        void linkPe()
        {
          pe_node.set(&pe, &pe_index);
        }

        void linkQgen()
        {
          qgen_node.set(&qgen, &qgen_index);
        }

        void linkOmega()
        {
          omega_node.set(&omega, &omega_index);
        }

        void linkPfaref()
        {
          pfaref_node.set(&pfaref, &pfaref_index);
        }

        void linkRequired()
        {
          linkPe();
          linkQgen();
        }

        void linkAllInputs()
        {
          linkRequired();
          linkOmega();
          qext_node.set(&qext, &qext_index);
          linkPfaref();
          pref_node.set(&pref, &pref_index);
        }

        template <typename ModelT>
        void attachRequired(ModelT& model)
        {
          model.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
          model.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        }

        template <typename ModelT>
        void attachAllInputs(ModelT& model)
        {
          attachRequired(model);
          model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
          model.getSignals().template attachSignalNode<Ext::QEXT>(&qext_node);
          model.getSignals().template attachSignalNode<Ext::PFAREF>(&pfaref_node);
          model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        }

        template <typename ModelT>
        void assignOutputs(ModelT& model)
        {
          model.getSignals().template assignSignalNode<Vars::IQCMD>(&iqcmd_node);
          model.getSignals().template assignSignalNode<Vars::IPCMD>(&ipcmd_node);
        }
      };

      auto makeTestData(bool spf = false) -> PhasorDynamics::Converter::ReecaData<RealT, IdxT>
      {
        PhasorDynamics::Converter::ReecaData<RealT, IdxT> data;
        data.device_class          = "Reeca";
        data.disambiguation_string = "reeca_test";
        data.ports[Ports::bus]     = 1;
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);

        data.parameters[Params::Sbase] = static_cast<IdxT>(100);
        data.parameters[Params::spf]   = spf;
        data.parameters[Params::sV]    = true;
        data.parameters[Params::sQ]    = true;
        data.parameters[Params::sP]    = false;
        data.parameters[Params::sPQ]   = false;
        data.parameters[Params::Trv]   = static_cast<RealT>(0.02);
        data.parameters[Params::Tp]    = static_cast<RealT>(0.04);
        data.parameters[Params::Vdip]  = static_cast<RealT>(0.90);
        data.parameters[Params::Vup]   = static_cast<RealT>(1.10);
        data.parameters[Params::Dbd1]  = static_cast<RealT>(-0.02);
        data.parameters[Params::Dbd2]  = static_cast<RealT>(0.02);
        data.parameters[Params::Kqv]   = static_cast<RealT>(2.0);
        data.parameters[Params::Iql1]  = static_cast<RealT>(-1.0);
        data.parameters[Params::Iqh1]  = static_cast<RealT>(1.0);
        data.parameters[Params::Iqfrz] = static_cast<RealT>(0.0);
        data.parameters[Params::Thld]  = static_cast<RealT>(0.0);
        data.parameters[Params::Qmax]  = static_cast<RealT>(0.5);
        data.parameters[Params::Qmin]  = static_cast<RealT>(-0.5);
        data.parameters[Params::Kqp]   = static_cast<RealT>(0.1);
        data.parameters[Params::Kqi]   = static_cast<RealT>(0.2);
        data.parameters[Params::Vmax]  = static_cast<RealT>(1.1);
        data.parameters[Params::Vmin]  = static_cast<RealT>(0.9);
        data.parameters[Params::Vref1] = static_cast<RealT>(0.0);
        data.parameters[Params::Kvp]   = static_cast<RealT>(0.1);
        data.parameters[Params::Kvi]   = static_cast<RealT>(0.2);
        data.parameters[Params::Tiq]   = static_cast<RealT>(0.02);
        data.parameters[Params::Tpord] = static_cast<RealT>(0.02);
        data.parameters[Params::dPmax] = static_cast<RealT>(999.0);
        data.parameters[Params::dPmin] = static_cast<RealT>(-999.0);
        data.parameters[Params::Pmax]  = static_cast<RealT>(1.0);
        data.parameters[Params::Pmin]  = static_cast<RealT>(0.0);
        data.parameters[Params::Imax]  = static_cast<RealT>(1.2);
        data.parameters[Params::vq1]   = static_cast<RealT>(0.5);
        data.parameters[Params::lq1]   = static_cast<RealT>(1.0);
        data.parameters[Params::vq2]   = static_cast<RealT>(0.7);
        data.parameters[Params::lq2]   = static_cast<RealT>(1.0);
        data.parameters[Params::vq3]   = static_cast<RealT>(0.9);
        data.parameters[Params::lq3]   = static_cast<RealT>(1.0);
        data.parameters[Params::vq4]   = static_cast<RealT>(1.1);
        data.parameters[Params::lq4]   = static_cast<RealT>(1.0);
        data.parameters[Params::vp1]   = static_cast<RealT>(0.5);
        data.parameters[Params::lp1]   = static_cast<RealT>(1.0);
        data.parameters[Params::vp2]   = static_cast<RealT>(0.7);
        data.parameters[Params::lp2]   = static_cast<RealT>(1.0);
        data.parameters[Params::vp3]   = static_cast<RealT>(0.9);
        data.parameters[Params::lp3]   = static_cast<RealT>(1.0);
        data.parameters[Params::vp4]   = static_cast<RealT>(1.1);
        data.parameters[Params::lp4]   = static_cast<RealT>(1.0);
        data.parameters[Params::Thld2] = static_cast<RealT>(0.0);

        return data;
      }

      void addAllMonitors(PhasorDynamics::Converter::ReecaData<RealT, IdxT>& data)
      {
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);
        data.monitored_variables.insert(Mon::piq);
        data.monitored_variables.insert(Mon::piv);
        data.monitored_variables.insert(Mon::qv);
        data.monitored_variables.insert(Mon::pord);
        data.monitored_variables.insert(Mon::qref);
        data.monitored_variables.insert(Mon::sdip);
        data.monitored_variables.insert(Mon::iqmax);
        data.monitored_variables.insert(Mon::ipmax);
        data.monitored_variables.insert(Mon::iqv);
        data.monitored_variables.insert(Mon::vqctrl);
        data.monitored_variables.insert(Mon::iqbase);
      }

      void relaxCurrentLimits(PhasorDynamics::Converter::ReecaData<RealT, IdxT>& data, RealT limit)
      {
        data.parameters[Params::Pmax] = limit;
        data.parameters[Params::Imax] = static_cast<RealT>(limit + 0.2);
        data.parameters[Params::lq1]  = limit;
        data.parameters[Params::lq2]  = limit;
        data.parameters[Params::lq3]  = limit;
        data.parameters[Params::lq4]  = limit;
        data.parameters[Params::lp1]  = limit;
        data.parameters[Params::lp2]  = limit;
        data.parameters[Params::lp3]  = limit;
        data.parameters[Params::lp4]  = limit;
      }

      template <typename ModelT>
      bool residualsAreNearZero(const ModelT& model)
      {
        bool success = true;
        for (const auto& residual : model.getResidual())
        {
          success *= std::isfinite(static_cast<RealT>(residual));
          success *= isEqual(residual, static_cast<ScalarT>(0.0), kResidualTol);
        }
        return success;
      }

      bool verifyWithRequiredSignals(PhasorDynamics::Converter::Reeca<ScalarT, IdxT>& model)
      {
        SignalFixture signals;
        signals.linkRequired();
        signals.attachRequired(model);
        return model.verify() > 0;
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params param, RealT value)
      {
        auto data              = makeTestData();
        data.parameters[param] = value;
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> model(&bus, data);
        SignalFixture                                   signals;
        signals.linkRequired();
        signals.attachRequired(model);
        return model.verify() > 0;
      }

      std::string validSystemJson()
      {
        return R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "REECA physics",
    "case_description": "REECA parser and system physics smoke test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
  ],
  "signals": [
    { "signal_id": 10, "name": "pe fixture" },
    { "signal_id": 11, "name": "qgen fixture" },
    { "signal_id": 12, "name": "iqcmd" },
    { "signal_id": 13, "name": "ipcmd" }
  ],
  "devices": [
    {
      "class": "Tgov1",
      "ports": { "pmech": 10 },
      "id": "GP",
      "params": { "R": 0.05, "T1": 0.5, "T2": 2.5, "T3": 7.5, "Pvmax": 1.0, "Pvmin": 0.0, "Dt": 0.0 }
    },
    {
      "class": "Tgov1",
      "ports": { "pmech": 11 },
      "id": "GQ",
      "params": { "R": 0.05, "T1": 0.5, "T2": 2.5, "T3": 7.5, "Pvmax": 1.0, "Pvmin": 0.0, "Dt": 0.0 }
    },
    {
      "class": "Reeca",
      "ports": { "bus": 1, "pe": 10, "qgen": 11, "iqcmd": 12, "ipcmd": 13 },
      "id": "RC",
      "params": {
        "Sbase": 100,
        "spf": false,
        "sV": true,
        "sQ": true,
        "sP": false,
        "sPQ": false,
        "Trv": 0.02,
        "Tp": 0.04,
        "Vdip": 0.9,
        "Vup": 1.1,
        "Dbd1": -0.02,
        "Dbd2": 0.02,
        "Kqv": 2.0,
        "Iql1": -1.0,
        "Iqh1": 1.0,
        "Iqfrz": 0.0,
        "Thld": 0.0,
        "Qmax": 0.5,
        "Qmin": -0.5,
        "Kqp": 0.1,
        "Kqi": 0.2,
        "Vmax": 1.1,
        "Vmin": 0.9,
        "Vref1": 0.0,
        "Kvp": 0.1,
        "Kvi": 0.2,
        "Tiq": 0.02,
        "Tpord": 0.02,
        "dPmax": 999.0,
        "dPmin": -999.0,
        "Pmax": 1.0,
        "Pmin": 0.0,
        "Imax": 1.2,
        "vq1": 0.5,
        "lq1": 1.0,
        "vq2": 0.7,
        "lq2": 1.0,
        "vq3": 0.9,
        "lq3": 1.0,
        "vq4": 1.1,
        "lq4": 1.0,
        "vp1": 0.5,
        "lp1": 1.0,
        "vp2": 0.7,
        "lp2": 1.0,
        "vp3": 0.9,
        "lp3": 1.0,
        "vp4": 1.1,
        "lp4": 1.0,
        "Thld2": 0.0
      },
      "mon": ["vmeas", "pmeas", "iqcmd", "ipcmd"]
    }
  ]
}
)json";
      }

      std::string initFallbackSystemJson()
      {
        return R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "REECA init fallback",
    "case_description": "REECA accepts init feedback fallbacks",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
  ],
  "devices": [
    {
      "class": "Reeca",
      "ports": { "bus": 1 },
      "id": "RC",
      "params": {
        "Sbase": 50,
        "spf": false,
        "sV": true,
        "sQ": true,
        "sP": false,
        "sPQ": false,
        "Trv": 0.02,
        "Tp": 0.04,
        "Vdip": 0.9,
        "Vup": 1.1,
        "Dbd1": -0.02,
        "Dbd2": 0.02,
        "Kqv": 2.0,
        "Iql1": -1.0,
        "Iqh1": 1.0,
        "Iqfrz": 0.0,
        "Thld": 0.0,
        "Qmax": 0.5,
        "Qmin": -0.5,
        "Kqp": 0.1,
        "Kqi": 0.2,
        "Vmax": 1.1,
        "Vmin": 0.9,
        "Vref1": 0.0,
        "Kvp": 0.1,
        "Kvi": 0.2,
        "Tiq": 0.02,
        "Tpord": 0.02,
        "dPmax": 999.0,
        "dPmin": -999.0,
        "Pmax": 2.0,
        "Pmin": 0.0,
        "Imax": 2.2,
        "vq1": 0.5,
        "lq1": 1.0,
        "vq2": 0.7,
        "lq2": 1.0,
        "vq3": 0.9,
        "lq3": 1.0,
        "vq4": 1.1,
        "lq4": 1.0,
        "vp1": 0.5,
        "lp1": 2.0,
        "vp2": 0.7,
        "lp2": 2.0,
        "vp3": 0.9,
        "lp3": 2.0,
        "vp4": 1.1,
        "lp4": 2.0
      },
      "init": { "pe": 0.8, "qgen": 0.2 },
      "mon": ["vmeas", "pmeas"]
    }
  ]
}
)json";
      }

      std::string initFallbackSystemJsonWithInit(const std::string& replacement_init)
      {
        auto              text = initFallbackSystemJson();
        const std::string init = R"json("init": { "pe": 0.8, "qgen": 0.2 })json";
        text.replace(text.find(init), init.size(), replacement_init);
        return text;
      }

      std::string unlinkedSystemJson()
      {
        return R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "REECA missing source",
    "case_description": "REECA must reject unlinked required feedback signals",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
  ],
  "signals": [
    { "signal_id": 10, "name": "unlinked pe" },
    { "signal_id": 11, "name": "unlinked qgen" }
  ],
  "devices": [
    {
      "class": "Reeca",
      "ports": { "bus": 1, "pe": 10, "qgen": 11 },
      "id": "RC",
      "params": {
        "Sbase": 100,
        "spf": false,
        "sV": true,
        "sQ": true,
        "sP": false,
        "sPQ": false,
        "Trv": 0.02,
        "Tp": 0.04,
        "Vdip": 0.9,
        "Vup": 1.1,
        "Dbd1": -0.02,
        "Dbd2": 0.02,
        "Kqv": 2.0,
        "Iql1": -1.0,
        "Iqh1": 1.0,
        "Iqfrz": 0.0,
        "Thld": 0.0,
        "Qmax": 0.5,
        "Qmin": -0.5,
        "Kqp": 0.1,
        "Kqi": 0.2,
        "Vmax": 1.1,
        "Vmin": 0.9,
        "Vref1": 0.0,
        "Kvp": 0.1,
        "Kvi": 0.2,
        "Tiq": 0.02,
        "Tpord": 0.02,
        "dPmax": 999.0,
        "dPmin": -999.0,
        "Pmax": 1.0,
        "Pmin": 0.0,
        "Imax": 1.2,
        "vq1": 0.5,
        "lq1": 1.0,
        "vq2": 0.7,
        "lq2": 1.0,
        "vq3": 0.9,
        "lq3": 1.0,
        "vq4": 1.1,
        "lq4": 1.0,
        "vp1": 0.5,
        "lp1": 1.0,
        "vp2": 0.7,
        "lp2": 1.0,
        "vp3": 0.9,
        "lp3": 1.0,
        "vp4": 1.1,
        "lp4": 1.0,
        "Thld2": 0.0
      }
    }
  ]
}
)json";
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      DependencyTracking::Variable::DependencyMap pruneJacobianTails(
          DependencyTracking::Variable::DependencyMap dependencies,
          RealT                                       tolerance)
      {
        for (auto it = dependencies.begin(); it != dependencies.end();)
        {
          if (std::abs(it->second) <= tolerance)
          {
            it = dependencies.erase(it);
          }
          else
          {
            ++it;
          }
        }
        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data = makeTestData();

        PhasorDynamics::Bus<DepVar, IdxT>              bus(DepVar{0.96}, DepVar{0.28});
        PhasorDynamics::Converter::Reeca<DepVar, IdxT> reeca(&bus, data);

        PhasorDynamics::SignalNode<DepVar, IdxT> pe_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> qgen_node;
        DepVar                                   pe_value{0.8};
        DepVar                                   qgen_value{0.2};
        IdxT                                     pe_index   = static_cast<IdxT>(reeca.size() + bus.size());
        IdxT                                     qgen_index = static_cast<IdxT>(reeca.size() + bus.size() + 1);

        bus.allocate();
        reeca.allocate();
        bus.initialize();

        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);

        reeca.initialize();

        for (IdxT i = 0; i < reeca.size(); ++i)
        {
          reeca.y()[static_cast<size_t>(i)].setVariableNumber(i);
          reeca.yp()[static_cast<size_t>(i)].setVariableNumber(i);
        }
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.y()[static_cast<size_t>(i)].setVariableNumber(i + reeca.size());
        }
        pe_value.setVariableNumber(pe_index);
        qgen_value.setVariableNumber(qgen_index);

        reeca.evaluateResidual();

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(
            static_cast<size_t>(reeca.size()));
        for (IdxT i = 0; i < reeca.size(); ++i)
        {
          dependencies[static_cast<size_t>(i)] = reeca.getResidual()[static_cast<size_t>(i)].getDependencies();
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(0.96, 0.28);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        ScalarT                                   pe_value{0.8};
        ScalarT                                   qgen_value{0.2};
        IdxT                                      pe_index   = static_cast<IdxT>(reeca.size() + bus.size());
        IdxT                                      qgen_index = static_cast<IdxT>(reeca.size() + bus.size() + 1);

        bus.allocate();
        reeca.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + reeca.size());
          bus.setResidualIndex(i, i + reeca.size());
        }

        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);

        bus.initialize();
        reeca.initialize();
        reeca.updateTime(0.0, 1.0);

        reeca.evaluateResidual();
        reeca.evaluateJacobian();

        auto model_jacobian = reeca.getJacobian();
        model_jacobian.deduplicate();
        return MapFromCOO(model_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
