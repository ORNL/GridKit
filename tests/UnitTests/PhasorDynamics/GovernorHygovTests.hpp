#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class GovernorHygovTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Gov     = PhasorDynamics::Governor::Hygov<ScalarT, IdxT>;
      using Data    = PhasorDynamics::Governor::HygovData<RealT, IdxT>;
      using Var     = PhasorDynamics::Governor::HygovInternalVariables;
      using Ext     = PhasorDynamics::Governor::HygovExternalVariables;
      using Params  = PhasorDynamics::Governor::HygovParameters;
      using Mon     = PhasorDynamics::Governor::HygovMonitorableVariables;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        Gov hygov(makeHygovData());
        success *= (hygov.size() == static_cast<IdxT>(Var::MAXIMUM));
        success *= (hygov.getMonitor() != nullptr);
        success *= (hygov.verify() == 0);

        auto source_default_hygov = makeHygovSourceDefaultData();
        Gov  source_default_hygov_model(source_default_hygov);
        success *= (source_default_hygov_model.verify() == 0);

        auto bad_hygov_curve                    = makeHygovData();
        bad_hygov_curve.parameters[Params::Gv2] = static_cast<RealT>(0.2);
        Gov bad_hygov_curve_model(bad_hygov_curve);
        success *= (bad_hygov_curve_model.verify() > 0);

        auto bad_hygov_trate                      = makeHygovData();
        bad_hygov_trate.parameters[Params::Trate] = static_cast<RealT>(0.0);
        Gov bad_hygov_trate_model(bad_hygov_trate);
        success *= (bad_hygov_trate_model.verify() > 0);

        auto missing_hygov_trate = makeHygovData();
        missing_hygov_trate.parameters.erase(Params::Trate);
        Gov missing_hygov_trate_model(missing_hygov_trate);
        success *= (missing_hygov_trate_model.verify() > 0);

        auto bad_hygov_hdam                     = makeHygovData();
        bad_hygov_hdam.parameters[Params::Hdam] = static_cast<RealT>(0.0);
        Gov bad_hygov_hdam_model(bad_hygov_hdam);
        success *= (bad_hygov_hdam_model.verify() > 0);

        auto bad_hygov_rtemp                      = makeHygovData();
        bad_hygov_rtemp.parameters[Params::Rtemp] = static_cast<RealT>(0.0);
        Gov bad_hygov_rtemp_model(bad_hygov_rtemp);
        success *= (bad_hygov_rtemp_model.verify() > 0);

        auto bad_hygov_limits                     = makeHygovData();
        bad_hygov_limits.parameters[Params::Velm] = static_cast<RealT>(-0.1);
        bad_hygov_limits.parameters[Params::db1]  = static_cast<RealT>(-0.1);
        Gov bad_hygov_limits_model(bad_hygov_limits);
        success *= (bad_hygov_limits_model.verify() > 0);

        auto bad_hygov_gate_limits                     = makeHygovData();
        bad_hygov_gate_limits.parameters[Params::Gmin] = static_cast<RealT>(1.1);
        Gov bad_hygov_gate_limits_model(bad_hygov_gate_limits);
        success *= (bad_hygov_gate_limits_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signals()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        ScalarT                                   pmech_value{0.40};
        ScalarT                                   omega_value{0.0};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = 5;
        pmech_node.set(&pmech_value, &pmech_index);
        omega_node.set(&omega_value, &omega_index);

        Gov   hygov(makeHygovData());
        auto& hygov_signals = hygov.getSignals();
        hygov_signals.template assignSignalNode<Var::PMECH>(&pmech_node);
        hygov_signals.template attachSignalNode<Ext::OMEGA>(&omega_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech_value);
        success *= pmech_node.linked();
        success *= (pmech_node.getVariableIndex() == static_cast<IdxT>(Var::PMECH));
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);
        success *= (hygov.tagDifferentiable() == 0);
        success *= (hygov.evaluateResidual() == 0);

        const auto* y         = hygov.y().getData();
        const auto* yp        = hygov.yp().getData();
        const auto& residual  = hygov.getResidual();
        const auto* f         = residual.getData();
        success              *= isEqual(y[index(Var::Q)], static_cast<ScalarT>(0.50), kTol);
        success              *= isEqual(y[index(Var::G)], static_cast<ScalarT>(0.50), kTol);
        success              *= isEqual(y[index(Var::C)], static_cast<ScalarT>(0.50), kTol);
        success              *= (hygov.tag()[index(Var::G)] == true);

        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome sourceDefault()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.40};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        pmech_node.set(&pmech_value, &pmech_index);

        Gov hygov(makeHygovSourceDefaultData());
        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech_value);
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);
        success *= (hygov.tagDifferentiable() == 0);
        success *= (hygov.evaluateResidual() == 0);

        const auto* y         = hygov.y().getData();
        const auto* yp        = hygov.yp().getData();
        const auto& residual  = hygov.getResidual();
        const auto* f         = residual.getData();
        success              *= isEqual(y[index(Var::Q)], static_cast<ScalarT>(0.50), kTol);
        success              *= isEqual(y[index(Var::PGV)], static_cast<ScalarT>(0.50), kTol);
        success              *= isEqual(y[index(Var::G)], static_cast<ScalarT>(0.50), kTol);
        success              *= (hygov.tag()[index(Var::XN)] == true);

        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome zeroTimeConstants()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.40};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        pmech_node.set(&pmech_value, &pmech_index);

        auto data                    = makeHygovData();
        data.parameters[Params::Tr]  = static_cast<RealT>(0.0);
        data.parameters[Params::Tf]  = static_cast<RealT>(0.0);
        data.parameters[Params::Tg]  = static_cast<RealT>(0.0);
        data.parameters[Params::Tw]  = static_cast<RealT>(0.0);
        data.parameters[Params::Tnp] = static_cast<RealT>(0.0);

        Gov hygov(data);
        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech_value);
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);
        success *= (hygov.tagDifferentiable() == 0);
        success *= (hygov.evaluateResidual() == 0);

        success *= (hygov.tag()[index(Var::XN)] == true);
        success *= (hygov.tag()[index(Var::XF)] == true);
        success *= (hygov.tag()[index(Var::C)] == true);
        success *= (hygov.tag()[index(Var::G)] == true);
        success *= (hygov.tag()[index(Var::Q)] == true);

        const auto* yp       = hygov.yp().getData();
        const auto& residual = hygov.getResidual();
        const auto* f        = residual.getData();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome baseConversion()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.40};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        pmech_node.set(&pmech_value, &pmech_index);

        auto data                      = makeHygovData();
        data.parameters[Params::Trate] = static_cast<RealT>(50.0);

        Gov hygov(data);
        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech_value);
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);
        success *= (hygov.evaluateResidual() == 0);

        const auto* y         = hygov.y().getData();
        const auto* yp        = hygov.yp().getData();
        const auto& residual  = hygov.getResidual();
        const auto* f         = residual.getData();
        success              *= isEqual(y[index(Var::Q)], static_cast<ScalarT>(0.90), kTol);
        success              *= isEqual(y[index(Var::G)], static_cast<ScalarT>(0.90), kTol);
        success              *= isEqual(y[index(Var::PMECH)], static_cast<ScalarT>(0.40), kTol);
        success              *= isEqual(pmech_node.read(), static_cast<ScalarT>(0.40), kTol);

        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome absoluteTolerance()
      {
        TestStatus success = true;

        Gov hygov(makeHygovData());

        success             *= (hygov.allocate() == 0);
        success             *= (hygov.setAbsoluteTolerance(static_cast<RealT>(1.0e-7)) == 0);
        const auto& abs_tol  = hygov.absoluteTolerance();
        success             *= (abs_tol.getSize() == static_cast<size_t>(Var::MAXIMUM));

        const auto* tolerances = abs_tol.getData();
        for (size_t i = 0; i < abs_tol.getSize(); ++i)
        {
          success *= isEqual(tolerances[i], scalar(1.0e-7), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome prefSignal()
      {
        TestStatus success = true;

        Gov hygov(makeHygovData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> paux_node;
        const ScalarT                             pmech0 = scalar(kInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   pref_value  = scalar(99.0);
        ScalarT                                   paux_value  = scalar(kInitialPaux);
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      pref_index  = 7;
        IdxT                                      paux_index  = 8;
        pmech_node.set(&pmech_value, &pmech_index);
        pref_node.set(&pref_value, &pref_index);
        paux_node.set(&paux_value, &paux_index);

        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        hygov.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        hygov.getSignals().template attachSignalNode<Ext::PAUX>(&paux_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech0);
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);
        success *= isEqual(pref_node.read(),
                           prefForInitialPoint(pmech0,
                                               paux_value,
                                               scalar(kDefaultTrate),
                                               scalar(kDefaultSystemBase)),
                           kTol);
        success *= (hygov.evaluateResidual() == 0);

        const auto* yp       = hygov.yp().getData();
        const auto& residual = hygov.getResidual();
        const auto* f        = residual.getData();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        pref_value += scalar(kPrefStep);
        success    *= (hygov.evaluateResidual() == 0);
        success    *= isEqual(f[index(Var::EF)], scalar(kPrefStep), kTol);

        return success.report(__func__);
      }

      TestOutcome prefSignalBaseConversion()
      {
        TestStatus success = true;

        auto data                      = makeHygovData();
        data.parameters[Params::Trate] = static_cast<RealT>(kConversionTrate);

        Gov hygov(data);
        hygov.setSystemBase(static_cast<RealT>(kSystemFrequency),
                            static_cast<RealT>(kConversionSystemBase * 1.0e6));

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> paux_node;
        const ScalarT                             pmech0 = scalar(kInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   pref_value  = scalar(99.0);
        ScalarT                                   paux_value  = scalar(kInitialPaux);
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      pref_index  = 7;
        IdxT                                      paux_index  = 8;
        pmech_node.set(&pmech_value, &pmech_index);
        pref_node.set(&pref_value, &pref_index);
        paux_node.set(&paux_value, &paux_index);

        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        hygov.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        hygov.getSignals().template attachSignalNode<Ext::PAUX>(&paux_node);

        success *= (hygov.allocate() == 0);
        pmech_node.init(pmech0);
        success *= (hygov.verify() == 0);
        success *= (hygov.initialize() == 0);

        const ScalarT expected_pref =
            prefForInitialPoint(pmech0,
                                paux_value,
                                scalar(kConversionTrate),
                                scalar(kConversionSystemBase));

        success       *= isEqual(pref_node.read(), expected_pref, kTol);
        const auto* y  = hygov.y().getData();
        success       *= isEqual(y[index(Var::Q)], scalar(0.90), kTol);
        success       *= isEqual(y[index(Var::G)], scalar(0.90), kTol);
        success       *= (hygov.evaluateResidual() == 0);

        const auto* yp       = hygov.yp().getData();
        const auto& residual = hygov.getResidual();
        const auto* f        = residual.getData();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
        }

        pref_value += scalar(kPrefStep);
        success    *= (hygov.evaluateResidual() == 0);
        success    *= isEqual(f[index(Var::EF)],
                           scalar(kConversionSystemBase / kConversionTrate * kPrefStep),
                           kTol);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        auto invalid_trate                      = makeHygovData();
        invalid_trate.parameters[Params::Trate] = true;
        Gov invalid_trate_model(invalid_trate);
        success *= (invalid_trate_model.verify() > 0);

        auto negative_time                   = makeHygovData();
        negative_time.parameters[Params::Tf] = static_cast<RealT>(-0.1);
        negative_time.parameters[Params::Tn] = static_cast<RealT>(-0.1);
        Gov negative_time_model(negative_time);
        success *= (negative_time_model.verify() > 0);

        auto invalid_at                   = makeHygovData();
        invalid_at.parameters[Params::At] = static_cast<RealT>(0.0);
        Gov invalid_at_model(invalid_at);
        success *= (invalid_at_model.verify() > 0);

        auto invalid_damping                      = makeHygovData();
        invalid_damping.parameters[Params::Dturb] = static_cast<RealT>(-0.1);
        Gov invalid_damping_model(invalid_damping);
        success *= (invalid_damping_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signalValidation()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        Gov                                       omega_model(makeHygovData());
        omega_model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        success *= (omega_model.verify() > 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        Gov                                       pref_model(makeHygovData());
        pref_model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        success *= (pref_model.verify() > 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT> paux_node;
        Gov                                       paux_model(makeHygovData());
        paux_model.getSignals().template attachSignalNode<Ext::PAUX>(&paux_node);
        success *= (paux_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome jsonParseAndSystemAssembly()
      {
        TestStatus success = true;

        std::istringstream input(R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "hydro governor",
    "case_description": "HYGOV parser test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    {
      "number": 1,
      "class": "bus",
      "name": "Bus 1",
      "init": { "Vr": 1.0, "Vi": 0.0 },
      "params": { "kv": 1.0 }
    }
  ],
  "signals": [
    { "signal_id": 10, "name": "Pmech" }
  ],
  "devices": [
    {
      "class": "Genrou",
      "ports": { "bus": 1, "pmech": 10 },
      "id": "GEN1",
      "params": {
        "p0": 0.3, "q0": 0.0, "H": 3.0, "D": 0.0, "Ra": 0.0,
        "Tdop": 7.0, "Tdopp": 0.04, "Tqop": 0.75, "Tqopp": 0.05,
        "Xd": 2.1, "Xdp": 0.2, "Xdpp": 0.18, "Xq": 0.5, "Xqp": 0.5,
        "Xqpp": 0.18, "Xl": 0.15, "S10": 0.0, "S12": 0.0,
        "mva": 100.0
      }
    },
    {
      "class": "Hygov",
      "ports": { "pmech": 10 },
      "id": "HYG1",
      "params": {
        "Trate": 50.0, "Rperm": 0.05, "Rtemp": 0.4, "Tr": 5.0,
        "Tf": 0.2, "Tg": 0.0, "Velm": 0.5, "Gmax": 1.0, "Gmin": 0.0,
        "Tw": 1.0, "At": 1.0, "Dturb": 0.0, "Qnl": 0.1,
        "Tn": 0.0, "Tnp": 1.0, "db1": 0.0, "db2": 0.0, "Hdam": 1.0,
        "Gv0": 0.0, "Gv1": 0.2, "Gv2": 0.4, "Gv3": 0.6, "Gv4": 0.8, "Gv5": 1.0,
        "Pgv0": 0.0, "Pgv1": 0.2, "Pgv2": 0.4, "Pgv3": 0.6, "Pgv4": 0.8, "Pgv5": 1.0
      }
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.hygov.size() == 1);
        const auto trate_param =
            data.hygov[0].parameters.at(PhasorDynamics::Governor::HygovParameters::Trate);
        success            *= (std::get<RealT>(trate_param) == static_cast<RealT>(50.0));
        using SignalOutput  = typename Data::SignalOutputs;
        success            *= data.hygov[0].buses.empty();
        success            *= data.hygov[0].signal_inputs.empty();
        success            *= (data.hygov[0].signal_outputs.at(SignalOutput::pmech)
                    == static_cast<IdxT>(10));
        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        const auto hygov_size =
            static_cast<IdxT>(PhasorDynamics::Governor::HygovInternalVariables::MAXIMUM);
        const auto  hygov_offset  = static_cast<size_t>(system.size() - hygov_size);
        const auto* system_y      = system.y().getData();
        success                  *= isEqual(
            system_y[hygov_offset + index(PhasorDynamics::Governor::HygovInternalVariables::Q)],
            static_cast<ScalarT>(0.70),
            kTol);
        success *= isEqual(
            system_y[hygov_offset + index(PhasorDynamics::Governor::HygovInternalVariables::G)],
            static_cast<ScalarT>(0.70),
            kTol);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 33);

        return success.report(__func__);
      }

    private:
      static constexpr RealT kSystemFrequency      = 60.0;
      static constexpr RealT kDefaultSystemBase    = 100.0;
      static constexpr RealT kDefaultTrate         = 100.0;
      static constexpr RealT kConversionSystemBase = 100.0;
      static constexpr RealT kConversionTrate      = 50.0;
      static constexpr RealT kInitialPmech         = 0.40;
      static constexpr RealT kInitialPaux          = 0.02;
      static constexpr RealT kPrefStep             = 0.10;
      static constexpr RealT kRperm                = 0.05;
      static constexpr RealT kAt                   = 1.0;
      static constexpr RealT kQnl                  = 0.1;
      static constexpr RealT kHdam                 = 1.0;

      static size_t index(PhasorDynamics::Governor::HygovInternalVariables variable)
      {
        return static_cast<size_t>(variable);
      }

      static ScalarT scalar(RealT value)
      {
        return static_cast<ScalarT>(value);
      }

      static ScalarT prefForInitialPoint(ScalarT pmech0,
                                         ScalarT paux0,
                                         ScalarT trate,
                                         ScalarT system_base)
      {
        const ScalarT k_base = system_base / trate;
        const ScalarT q0     = scalar(kQnl) + k_base * pmech0 / (scalar(kAt) * scalar(kHdam));
        const ScalarT gate0  = q0 / std::sqrt(scalar(kHdam));
        return (scalar(kRperm) * gate0 - k_base * paux0) / k_base;
      }

      auto makeHygovData() -> Data
      {
        Data data;
        data.device_class          = "Hygov";
        data.disambiguation_string = "hygov_test";
        data.monitored_variables.insert(Mon::pmech);
        data.monitored_variables.insert(Mon::gate);

        data.parameters[Params::Trate] = static_cast<RealT>(100.0);
        data.parameters[Params::Rperm] = static_cast<RealT>(0.05);
        data.parameters[Params::Rtemp] = static_cast<RealT>(0.4);
        data.parameters[Params::Tr]    = static_cast<RealT>(5.0);
        data.parameters[Params::Tf]    = static_cast<RealT>(0.2);
        data.parameters[Params::Tg]    = static_cast<RealT>(0.0);
        data.parameters[Params::Velm]  = static_cast<RealT>(0.5);
        data.parameters[Params::Gmax]  = static_cast<RealT>(1.0);
        data.parameters[Params::Gmin]  = static_cast<RealT>(0.0);
        data.parameters[Params::Tw]    = static_cast<RealT>(1.0);
        data.parameters[Params::At]    = static_cast<RealT>(1.0);
        data.parameters[Params::Dturb] = static_cast<RealT>(0.0);
        data.parameters[Params::Qnl]   = static_cast<RealT>(0.1);
        data.parameters[Params::Tn]    = static_cast<RealT>(0.0);
        data.parameters[Params::Tnp]   = static_cast<RealT>(1.0);
        data.parameters[Params::db1]   = static_cast<RealT>(0.0);
        data.parameters[Params::db2]   = static_cast<RealT>(0.0);
        data.parameters[Params::Hdam]  = static_cast<RealT>(1.0);
        data.parameters[Params::Gv0]   = static_cast<RealT>(0.0);
        data.parameters[Params::Gv1]   = static_cast<RealT>(0.2);
        data.parameters[Params::Gv2]   = static_cast<RealT>(0.4);
        data.parameters[Params::Gv3]   = static_cast<RealT>(0.6);
        data.parameters[Params::Gv4]   = static_cast<RealT>(0.8);
        data.parameters[Params::Gv5]   = static_cast<RealT>(1.0);
        data.parameters[Params::Pgv0]  = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv1]  = static_cast<RealT>(0.2);
        data.parameters[Params::Pgv2]  = static_cast<RealT>(0.4);
        data.parameters[Params::Pgv3]  = static_cast<RealT>(0.6);
        data.parameters[Params::Pgv4]  = static_cast<RealT>(0.8);
        data.parameters[Params::Pgv5]  = static_cast<RealT>(1.0);

        return data;
      }

      auto makeHygovSourceDefaultData() -> Data
      {
        auto data                     = makeHygovData();
        data.parameters[Params::Tn]   = static_cast<RealT>(0.0);
        data.parameters[Params::Tnp]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv0]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv1]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv2]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv3]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv4]  = static_cast<RealT>(0.0);
        data.parameters[Params::Gv5]  = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv0] = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv1] = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv2] = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv3] = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv4] = static_cast<RealT>(0.0);
        data.parameters[Params::Pgv5] = static_cast<RealT>(0.0);

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
