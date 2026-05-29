#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>

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

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> hygov(makeHygovData());
        success *= (hygov.size()
                    == static_cast<IdxT>(
                        PhasorDynamics::Governor::HygovInternalVariables::MAXIMUM));
        success *= (hygov.getMonitor() != nullptr);
        success *= (hygov.verify() == 0);

        auto                                           source_default_hygov = makeHygovSourceDefaultData();
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> source_default_hygov_model(
            source_default_hygov);
        success *= (source_default_hygov_model.verify() == 0);

        auto bad_hygov = makeHygovData();
        bad_hygov.parameters[PhasorDynamics::Governor::HygovParameters::db2] =
            static_cast<RealT>(0.01);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_hygov_model(bad_hygov);
        success *= (bad_hygov_model.verify() > 0);

        auto bad_hygov_leadlag = makeHygovData();
        bad_hygov_leadlag.parameters[PhasorDynamics::Governor::HygovParameters::Tn] =
            static_cast<RealT>(0.1);
        bad_hygov_leadlag.parameters[PhasorDynamics::Governor::HygovParameters::Tnp] =
            static_cast<RealT>(0.0);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_hygov_leadlag_model(
            bad_hygov_leadlag);
        success *= (bad_hygov_leadlag_model.verify() > 0);

        auto bad_hygov_curve = makeHygovData();
        bad_hygov_curve.parameters[PhasorDynamics::Governor::HygovParameters::Gv2] =
            static_cast<RealT>(0.2);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_hygov_curve_model(bad_hygov_curve);
        success *= (bad_hygov_curve_model.verify() > 0);

        auto bad_hygov_trate = makeHygovData();
        bad_hygov_trate.parameters[PhasorDynamics::Governor::HygovParameters::Trate] =
            static_cast<RealT>(0.0);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_hygov_trate_model(bad_hygov_trate);
        success *= (bad_hygov_trate_model.verify() > 0);

        auto bad_hygov_base = makeHygovData();
        bad_hygov_base.parameters[PhasorDynamics::Governor::HygovParameters::mva] =
            static_cast<RealT>(0.0);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_hygov_base_model(bad_hygov_base);
        success *= (bad_hygov_base_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signals()
      {
        using Var = PhasorDynamics::Governor::HygovInternalVariables;

        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        ScalarT                                   pmech_value{0.40};
        ScalarT                                   omega_value{0.0};
        IdxT                                      pmech_index = 4;
        IdxT                                      omega_index = 5;
        pmech_node.set(&pmech_value, &pmech_index);
        omega_node.set(&omega_value, &omega_index);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> hygov(makeHygovData());
        auto&                                          hygov_signals = hygov.getSignals();
        hygov_signals
            .template assignSignalNode<PhasorDynamics::Governor::HygovInternalVariables::PMECH>(
                &pmech_node);
        hygov_signals
            .template attachSignalNode<PhasorDynamics::Governor::HygovExternalVariables::OMEGA>(
                &omega_node);

        success                      *= (hygov.allocate() == 0);
        hygov.y()[index(Var::PMECH)]  = pmech_value;
        success                      *= (hygov.verify() == 0);
        success                      *= (hygov.initialize() == 0);
        success                      *= (hygov.tagDifferentiable() == 0);
        success                      *= (hygov.evaluateResidual() == 0);

        success *= isEqual(hygov.y()[index(Var::Q)], static_cast<ScalarT>(0.50), kTol);
        success *= isEqual(hygov.y()[index(Var::G)], static_cast<ScalarT>(0.50), kTol);
        success *= isEqual(hygov.y()[index(Var::C)], static_cast<ScalarT>(0.50), kTol);
        success *= (hygov.tag()[index(Var::G)] == false);

        for (size_t i = 0; i < hygov.getResidual().size(); ++i)
        {
          success *= isEqual(hygov.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(hygov.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome sourceDefault()
      {
        using Var = PhasorDynamics::Governor::HygovInternalVariables;

        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.40};
        IdxT                                      pmech_index = 4;
        pmech_node.set(&pmech_value, &pmech_index);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> hygov(makeHygovSourceDefaultData());
        hygov.getSignals()
            .template assignSignalNode<PhasorDynamics::Governor::HygovInternalVariables::PMECH>(
                &pmech_node);

        success                      *= (hygov.allocate() == 0);
        hygov.y()[index(Var::PMECH)]  = pmech_value;
        success                      *= (hygov.verify() == 0);
        success                      *= (hygov.initialize() == 0);
        success                      *= (hygov.tagDifferentiable() == 0);
        success                      *= (hygov.evaluateResidual() == 0);

        success *= isEqual(hygov.y()[index(Var::Q)], static_cast<ScalarT>(0.50), kTol);
        success *= isEqual(hygov.y()[index(Var::PGV)], static_cast<ScalarT>(0.50), kTol);
        success *= isEqual(hygov.y()[index(Var::G)], static_cast<ScalarT>(0.50), kTol);
        success *= (hygov.tag()[index(Var::XN)] == false);

        for (size_t i = 0; i < hygov.getResidual().size(); ++i)
        {
          success *= isEqual(hygov.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(hygov.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome baseConversion()
      {
        using Var    = PhasorDynamics::Governor::HygovInternalVariables;
        using Params = PhasorDynamics::Governor::HygovParameters;

        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.40};
        IdxT                                      pmech_index = 6;
        pmech_node.set(&pmech_value, &pmech_index);

        auto data                      = makeHygovData();
        data.parameters[Params::Trate] = static_cast<RealT>(50.0);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> hygov(data);
        hygov.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);

        success                      *= (hygov.allocate() == 0);
        hygov.y()[index(Var::PMECH)]  = pmech_value;
        success                      *= (hygov.verify() == 0);
        success                      *= (hygov.initialize() == 0);
        success                      *= (hygov.evaluateResidual() == 0);

        success *= isEqual(hygov.y()[index(Var::Q)], static_cast<ScalarT>(0.90), kTol);
        success *= isEqual(hygov.y()[index(Var::G)], static_cast<ScalarT>(0.90), kTol);
        success *= isEqual(hygov.y()[index(Var::PMECH)], static_cast<ScalarT>(0.40), kTol);
        success *= isEqual(pmech_node.read(), static_cast<ScalarT>(0.40), kTol);

        for (size_t i = 0; i < hygov.getResidual().size(); ++i)
        {
          success *= isEqual(hygov.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(hygov.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

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
      "v_base": 1.0
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
        "Trate": 50.0, "mva": 100.0, "Rperm": 0.05, "Rtemp": 0.4, "Tr": 5.0,
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
        const auto mva_param =
            data.hygov[0].parameters.at(PhasorDynamics::Governor::HygovParameters::mva);
        success *= (std::get<RealT>(trate_param) == static_cast<RealT>(50.0));
        success *= (std::get<RealT>(mva_param) == static_cast<RealT>(100.0));
        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        const auto hygov_size =
            static_cast<IdxT>(PhasorDynamics::Governor::HygovInternalVariables::MAXIMUM);
        const auto hygov_offset  = static_cast<size_t>(system.size() - hygov_size);
        success                 *= isEqual(
            system.y()[hygov_offset + index(PhasorDynamics::Governor::HygovInternalVariables::Q)],
            static_cast<ScalarT>(0.70),
            kTol);
        success *= isEqual(
            system.y()[hygov_offset + index(PhasorDynamics::Governor::HygovInternalVariables::G)],
            static_cast<ScalarT>(0.70),
            kTol);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 35);

        return success.report(__func__);
      }

    private:
      static size_t index(PhasorDynamics::Governor::HygovInternalVariables variable)
      {
        return static_cast<size_t>(variable);
      }

      auto makeHygovData() -> PhasorDynamics::Governor::HygovData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Governor::HygovParameters;
        using Mon    = PhasorDynamics::Governor::HygovMonitorableVariables;

        PhasorDynamics::Governor::HygovData<RealT, IdxT> data;
        data.device_class          = "Hygov";
        data.disambiguation_string = "hygov_test";
        data.monitored_variables.insert(Mon::pmech);
        data.monitored_variables.insert(Mon::gate);

        data.parameters[Params::Trate] = static_cast<RealT>(100.0);
        data.parameters[Params::mva]   = static_cast<RealT>(100.0);
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

      auto makeHygovSourceDefaultData() -> PhasorDynamics::Governor::HygovData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Governor::HygovParameters;

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
