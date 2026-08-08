#pragma once

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/GensalData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    using namespace GridKit::PhasorDynamics;
    using namespace GridKit::Model;

    template <typename real_type, typename index_type>
    class PDIntegrationTests
    {
    public:
      using RealT            = real_type;
      using IdxT             = index_type;
      using BusDataT         = BusData<RealT, IdxT>;
      using SystemModelDataT = SystemModelData<RealT, IdxT>;

      TestStatus compare(const SystemModelDataT& a, const SystemModelDataT& b)
      {
        using namespace GridKit::PhasorDynamics;
        using namespace GridKit::PhasorDynamics::Governor;
        using namespace GridKit::PhasorDynamics::Exciter;

        TestStatus success = true;

        success *= a.bus.size() == b.bus.size();
        for (std::size_t i = 0; i < a.bus.size(); ++i)
        {
          const auto& abus = a.bus[i];
          const auto& bbus = b.bus[i];

          success *= abus.bus_id == bbus.bus_id;
          success *= abus.bus_type == bbus.bus_type;
          success *= abus.Vr0 == bbus.Vr0;
          success *= abus.Vi0 == bbus.Vi0;
        }

        success *= a.signal.size() == b.signal.size();
        for (std::size_t i = 0; i < a.signal.size(); ++i)
        {
          const auto& as = a.signal[i];
          const auto& bs = b.signal[i];

          success *= as.name == bs.name;
          success *= as.signal_id == bs.signal_id;
        }

        success *= a.branch.size() == b.branch.size();
        for (std::size_t i = 0; i < a.branch.size(); ++i)
        {
          using Param = BranchParameters;

          const auto& abr = a.branch[i];
          const auto& bbr = b.branch[i];

          success *= abr.buses.at(BranchBuses::bus1) == bbr.buses.at(BranchBuses::bus1);
          success *= abr.buses.at(BranchBuses::bus2) == bbr.buses.at(BranchBuses::bus2);
          success *= abr.parameters.at(Param::R) == bbr.parameters.at(Param::R);
          success *= abr.parameters.at(Param::X) == bbr.parameters.at(Param::X);
          success *= abr.parameters.at(Param::G) == bbr.parameters.at(Param::G);
          success *= abr.parameters.at(Param::B) == bbr.parameters.at(Param::B);
        }

        success *= a.genrou.size() == b.genrou.size();
        for (std::size_t i = 0; i < a.genrou.size(); ++i)
        {
          constexpr auto BUS = GenrouBuses::bus;
          using Param        = GenrouParameters;

          const auto& ag = a.genrou[i];
          const auto& bg = b.genrou[i];

          success *= ag.buses.at(BUS) == bg.buses.at(BUS);
          success *= ag.parameters.at(Param::p0) == bg.parameters.at(Param::p0);
          success *= ag.parameters.at(Param::q0) == bg.parameters.at(Param::q0);
          success *= ag.parameters.at(Param::H) == bg.parameters.at(Param::H);
          success *= ag.parameters.at(Param::D) == bg.parameters.at(Param::D);
          success *= ag.parameters.at(Param::Ra) == bg.parameters.at(Param::Ra);
          success *= ag.parameters.at(Param::Tdop) == bg.parameters.at(Param::Tdop);
          success *= ag.parameters.at(Param::Tdopp) == bg.parameters.at(Param::Tdopp);
          success *= ag.parameters.at(Param::Tqopp) == bg.parameters.at(Param::Tqopp);
          success *= ag.parameters.at(Param::Tqop) == bg.parameters.at(Param::Tqop);
          success *= ag.parameters.at(Param::Xd) == bg.parameters.at(Param::Xd);
          success *= ag.parameters.at(Param::Xdp) == bg.parameters.at(Param::Xdp);
          success *= ag.parameters.at(Param::Xdpp) == bg.parameters.at(Param::Xdpp);
          success *= ag.parameters.at(Param::Xq) == bg.parameters.at(Param::Xq);
          success *= ag.parameters.at(Param::Xqp) == bg.parameters.at(Param::Xqp);
          success *= ag.parameters.at(Param::Xqpp) == bg.parameters.at(Param::Xqpp);
          success *= ag.parameters.at(Param::Xl) == bg.parameters.at(Param::Xl);
          success *= ag.parameters.at(Param::S10) == bg.parameters.at(Param::S10);
          success *= ag.parameters.at(Param::S12) == bg.parameters.at(Param::S12);
        }

        success *= a.genclassical.size() == b.genclassical.size();
        for (std::size_t i = 0; i < a.genclassical.size(); ++i)
        {
          constexpr auto BUS = GenClassicalBuses::bus;
          using Param        = GenClassicalParameters;

          const auto& ag = a.genclassical[i];
          const auto& bg = b.genclassical[i];

          success *= ag.buses.at(BUS) == bg.buses.at(BUS);
          success *= ag.parameters.at(Param::p0) == bg.parameters.at(Param::p0);
          success *= ag.parameters.at(Param::q0) == bg.parameters.at(Param::q0);
          success *= ag.parameters.at(Param::H) == bg.parameters.at(Param::H);
          success *= ag.parameters.at(Param::D) == bg.parameters.at(Param::D);
          success *= ag.parameters.at(Param::Ra) == bg.parameters.at(Param::Ra);
          success *= ag.parameters.at(Param::Xdp) == bg.parameters.at(Param::Xdp);
        }

        success *= a.gov.size() == b.gov.size();
        for (std::size_t i = 0; i < a.gov.size(); ++i)
        {
          using SignalIn  = Tgov1SignalInputs;
          using SignalOut = Tgov1SignalOutputs;
          using Param     = Tgov1Parameters;

          const auto& ag = a.gov[i];
          const auto& bg = b.gov[i];

          success *= ag.signal_inputs.count(SignalIn::speed) == bg.signal_inputs.count(SignalIn::speed);
          if (ag.signal_inputs.count(SignalIn::speed))
          {
            success *= ag.signal_inputs.at(SignalIn::speed) == bg.signal_inputs.at(SignalIn::speed);
          }
          success *= ag.signal_outputs.count(SignalOut::pmech) == bg.signal_outputs.count(SignalOut::pmech);
          if (ag.signal_outputs.count(SignalOut::pmech))
          {
            success *= ag.signal_outputs.at(SignalOut::pmech) == bg.signal_outputs.at(SignalOut::pmech);
          }
          success *= ag.parameters.at(Param::Trate) == bg.parameters.at(Param::Trate);
          success *= ag.parameters.at(Param::R) == bg.parameters.at(Param::R);
          success *= ag.parameters.at(Param::Pvmin) == bg.parameters.at(Param::Pvmin);
          success *= ag.parameters.at(Param::Pvmax) == bg.parameters.at(Param::Pvmax);
          success *= ag.parameters.at(Param::T1) == bg.parameters.at(Param::T1);
          success *= ag.parameters.at(Param::T2) == bg.parameters.at(Param::T2);
          success *= ag.parameters.at(Param::T3) == bg.parameters.at(Param::T3);
          success *= ag.parameters.at(Param::Dt) == bg.parameters.at(Param::Dt);
        }

        success *= a.exciter.size() == b.exciter.size();
        for (std::size_t i = 0; i < a.exciter.size(); ++i)
        {
          using SignalIn  = Ieeet1SignalInputs;
          using SignalOut = Ieeet1SignalOutputs;
          using Param     = Ieeet1Parameters;

          const auto& ae = a.exciter[i];
          const auto& be = b.exciter[i];

          success *= ae.buses.at(Ieeet1Buses::bus) == be.buses.at(Ieeet1Buses::bus);
          success *= ae.signal_inputs.count(SignalIn::speed) == be.signal_inputs.count(SignalIn::speed);
          if (ae.signal_inputs.count(SignalIn::speed))
          {
            success *= ae.signal_inputs.at(SignalIn::speed) == be.signal_inputs.at(SignalIn::speed);
          }
          success *= ae.signal_outputs.count(SignalOut::efd) == be.signal_outputs.count(SignalOut::efd);
          if (ae.signal_outputs.count(SignalOut::efd))
          {
            success *= ae.signal_outputs.at(SignalOut::efd) == be.signal_outputs.at(SignalOut::efd);
          }
          success *= ae.parameters.at(Param::Tr) == be.parameters.at(Param::Tr);
          success *= ae.parameters.at(Param::Ka) == be.parameters.at(Param::Ka);
          success *= ae.parameters.at(Param::Ta) == be.parameters.at(Param::Ta);
          success *= ae.parameters.at(Param::Ke) == be.parameters.at(Param::Ke);
          success *= ae.parameters.at(Param::Te) == be.parameters.at(Param::Te);
          success *= ae.parameters.at(Param::Kf) == be.parameters.at(Param::Kf);
          success *= ae.parameters.at(Param::Tf) == be.parameters.at(Param::Tf);
          success *= ae.parameters.at(Param::Vrmin) == be.parameters.at(Param::Vrmin);
          success *= ae.parameters.at(Param::Vrmax) == be.parameters.at(Param::Vrmax);
          success *= ae.parameters.at(Param::E1) == be.parameters.at(Param::E1);
          success *= ae.parameters.at(Param::E2) == be.parameters.at(Param::E2);
          success *= ae.parameters.at(Param::Se1) == be.parameters.at(Param::Se1);
          success *= ae.parameters.at(Param::Se2) == be.parameters.at(Param::Se2);
          success *= ae.parameters.at(Param::Ispdlim) == be.parameters.at(Param::Ispdlim);
        }

        success *= a.loadz.size() == b.loadz.size();
        for (std::size_t i = 0; i < a.loadz.size(); ++i)
        {
          using Param = LoadZParameters;

          const auto& al = a.loadz[i];
          const auto& bl = b.loadz[i];

          success *= al.buses.at(LoadZBuses::bus) == bl.buses.at(LoadZBuses::bus);
          success *= al.parameters.at(Param::R) == bl.parameters.at(Param::R);
          success *= al.parameters.at(Param::X) == bl.parameters.at(Param::X);
        }

        success *= a.bus_fault.size() == b.bus_fault.size();
        for (std::size_t i = 0; i < a.bus_fault.size(); ++i)
        {
          using Param = BusFaultParameters;

          const auto& af = a.bus_fault[i];
          const auto& bf = b.bus_fault[i];

          success *= af.buses.at(BusFaultBuses::bus) == bf.buses.at(BusFaultBuses::bus);
          success *= af.parameters.at(Param::R) == bf.parameters.at(Param::R);
          success *= af.parameters.at(Param::X) == bf.parameters.at(Param::X);
          success *= af.parameters.at(Param::state0) == bf.parameters.at(Param::state0);
        }

        return success;
      }

      std::unique_ptr<ErrorSet> runSimulation(SystemModelDataT&  system_data,
                                              const std::string& ref_file)
      {
        system_data.monitor_sink.emplace_back(VariableMonitorFormat::CSV, "mon.csv");

        // Create system model
        SystemModel<RealT, IdxT> sys(system_data);
        sys.allocate();

        // Get access to the fault
        auto* fault = sys.getBusFault(0);

        // Set time step to 1/4 of a 60Hz cycle
        real_type dt = 1.0 / 4.0 / 60.0;

        // Set up simulation
        AnalysisManager::Sundials::Ida<RealT, size_t> ida(&sys);
        ida.configureSimulation();

        // Run for 1s
        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0, dt);

        // Introduce fault and run for the next 0.1s
        fault->setStatus(true);
        ida.initializeSimulation(1.0);
        ida.runSimulation(1.1, dt);

        // Clear the fault and run until t = 10s.
        fault->setStatus(false);
        ida.initializeSimulation(1.1);
        ida.runSimulation(10.0, dt);

        sys.stopMonitor();

        // Return error set
        return GridKit::Testing::compareCSV("mon.csv", ref_file);
      }

      TestOutcome twoBusBasic()
      {
        using namespace GridKit::PhasorDynamics;
        using BusType   = BusDataT::BusType;
        using BusVar    = BusMonitorableVariables;
        using GenrouVar = GenrouMonitorableVariables;

        SystemModelDataT set_data;

        set_data.bus.resize(2);
        set_data.bus[0].bus_id   = 0;
        set_data.bus[0].bus_type = BusType::DEFAULT;
        set_data.bus[0].Vr0      = 0.9949877346411762;
        set_data.bus[0].Vi0      = 0.09999703952427966;
        set_data.bus[0].monitored_variables.insert(BusVar::Vm);

        set_data.bus[1].bus_id   = 1;
        set_data.bus[1].bus_type = BusType::SLACK;
        set_data.bus[1].Vr0      = 1.0;
        set_data.bus[1].Vi0      = 0.0;

        set_data.branch.resize(1);
        set_data.branch[0].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[0].buses[BranchBuses::bus2]        = set_data.bus[1].bus_id;
        set_data.branch[0].parameters[BranchParameters::R] = 0.0;
        set_data.branch[0].parameters[BranchParameters::X] = 0.1;
        set_data.branch[0].parameters[BranchParameters::G] = 0.0;
        set_data.branch[0].parameters[BranchParameters::B] = 0.0;

        set_data.genrou.resize(1);
        set_data.genrou[0].buses[GenrouBuses::bus]             = 0;
        set_data.genrou[0].parameters[GenrouParameters::p0]    = 1.;
        set_data.genrou[0].parameters[GenrouParameters::q0]    = 0.05013;
        set_data.genrou[0].parameters[GenrouParameters::H]     = 3.;
        set_data.genrou[0].parameters[GenrouParameters::D]     = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Ra]    = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Tdop]  = 7.;
        set_data.genrou[0].parameters[GenrouParameters::Tdopp] = .04;
        set_data.genrou[0].parameters[GenrouParameters::Tqopp] = .05;
        set_data.genrou[0].parameters[GenrouParameters::Tqop]  = .75;
        set_data.genrou[0].parameters[GenrouParameters::Xd]    = 2.1;
        set_data.genrou[0].parameters[GenrouParameters::Xdp]   = 0.2;
        set_data.genrou[0].parameters[GenrouParameters::Xdpp]  = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xq]    = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqp]   = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqpp]  = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xl]    = 0.15;
        set_data.genrou[0].parameters[GenrouParameters::S10]   = 0.;
        set_data.genrou[0].parameters[GenrouParameters::S12]   = 0.;
        set_data.genrou[0].monitored_variables.insert(GenrouVar::speed);

        set_data.bus_fault.resize(1);
        set_data.bus_fault[0].buses[BusFaultBuses::bus]              = 0;
        set_data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
        set_data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-3;
        set_data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

        std::string      base_name = "TwoBus/Basic/TwoBusBasic";
        std::string      in_file   = base_name + ".case.json";
        SystemModelDataT file_data = parseSystemModelData(in_file);

        auto success = compare(set_data, file_data);

        auto error_set = runSimulation(set_data, base_name + ".ref.csv");

        RealT error_V_allowed = 2.01e-4;
        RealT error_w_allowed = 1e-4;

        success *= error_set->var_errors[0].max_value < error_V_allowed;
        success *= error_set->var_errors[1].max_value < error_w_allowed;

        return success.report(__func__);
      }

      TestOutcome twoBusIeeet1()
      {
        using namespace GridKit::PhasorDynamics;
        using namespace GridKit::PhasorDynamics::Governor;
        using namespace GridKit::PhasorDynamics::Exciter;
        using BusType   = BusDataT::BusType;
        using BusVar    = BusMonitorableVariables;
        using GenrouVar = GenrouMonitorableVariables;
        using Ieeet1Var = Ieeet1MonitorableVariables;

        SystemModelDataT set_data;

        set_data.bus.resize(2);
        set_data.bus[0].bus_id   = 0;
        set_data.bus[0].bus_type = BusType::DEFAULT;
        set_data.bus[0].Vr0      = 0.9949877346411762;
        set_data.bus[0].Vi0      = 0.09999703952427966;
        set_data.bus[0].monitored_variables.insert(BusVar::Vm);
        set_data.bus[1].bus_id   = 1;
        set_data.bus[1].bus_type = BusType::SLACK;
        set_data.bus[1].Vr0      = 1.0;
        set_data.bus[1].Vi0      = 0.0;

        set_data.signal.resize(3);
        set_data.signal[0].name      = "Machine Speed Deviation";
        set_data.signal[0].signal_id = 0;
        set_data.signal[1].name      = "Mechanical Power";
        set_data.signal[1].signal_id = 1;
        set_data.signal[2].name      = "Excitation Field";
        set_data.signal[2].signal_id = 2;

        set_data.branch.resize(1);
        set_data.branch[0].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[0].buses[BranchBuses::bus2]        = set_data.bus[1].bus_id;
        set_data.branch[0].parameters[BranchParameters::R] = 0.0;
        set_data.branch[0].parameters[BranchParameters::X] = 0.1;
        set_data.branch[0].parameters[BranchParameters::G] = 0.0;
        set_data.branch[0].parameters[BranchParameters::B] = 0.0;

        set_data.bus_fault.resize(1);
        set_data.bus_fault[0].buses[BusFaultBuses::bus]              = 0;
        set_data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
        set_data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-3;
        set_data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

        set_data.genrou.resize(1);
        set_data.genrou[0].buses[GenrouBuses::bus]                    = 0;
        set_data.genrou[0].signal_outputs[GenrouSignalOutputs::speed] = 0;
        set_data.genrou[0].signal_inputs[GenrouSignalInputs::pmech]   = 1;
        set_data.genrou[0].signal_inputs[GenrouSignalInputs::efd]     = 2;
        set_data.genrou[0].parameters[GenrouParameters::p0]           = 1.;
        set_data.genrou[0].parameters[GenrouParameters::q0]           = 0.05013;
        set_data.genrou[0].parameters[GenrouParameters::H]            = 3.;
        set_data.genrou[0].parameters[GenrouParameters::D]            = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Ra]           = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Tdop]         = 7.;
        set_data.genrou[0].parameters[GenrouParameters::Tdopp]        = .04;
        set_data.genrou[0].parameters[GenrouParameters::Tqopp]        = .05;
        set_data.genrou[0].parameters[GenrouParameters::Tqop]         = .75;
        set_data.genrou[0].parameters[GenrouParameters::Xd]           = 2.1;
        set_data.genrou[0].parameters[GenrouParameters::Xdp]          = 0.2;
        set_data.genrou[0].parameters[GenrouParameters::Xdpp]         = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xq]           = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqp]          = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqpp]         = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xl]           = 0.15;
        set_data.genrou[0].parameters[GenrouParameters::S10]          = 0.;
        set_data.genrou[0].parameters[GenrouParameters::S12]          = 0.;
        set_data.genrou[0].parameters[GenrouParameters::mva]          = 100.0;
        set_data.genrou[0].monitored_variables.insert(GenrouVar::speed);

        set_data.gov.resize(1);
        set_data.gov[0].signal_inputs[Tgov1SignalInputs::speed]   = 0;
        set_data.gov[0].signal_outputs[Tgov1SignalOutputs::pmech] = 1;
        set_data.gov[0].parameters[Tgov1Parameters::Trate]        = 100.0;
        set_data.gov[0].parameters[Tgov1Parameters::R]            = 0.05;
        set_data.gov[0].parameters[Tgov1Parameters::Pvmin]        = 0.0;
        set_data.gov[0].parameters[Tgov1Parameters::Pvmax]        = 1.0;
        set_data.gov[0].parameters[Tgov1Parameters::T1]           = 0.5;
        set_data.gov[0].parameters[Tgov1Parameters::T2]           = 2.5;
        set_data.gov[0].parameters[Tgov1Parameters::T3]           = 7.5;
        set_data.gov[0].parameters[Tgov1Parameters::Dt]           = 0.0;

        set_data.exciter.resize(1);
        set_data.exciter[0].buses[Ieeet1Buses::bus]                  = 0;
        set_data.exciter[0].signal_inputs[Ieeet1SignalInputs::speed] = 0;
        set_data.exciter[0].signal_outputs[Ieeet1SignalOutputs::efd] = 2;
        set_data.exciter[0].parameters[Ieeet1Parameters::Tr]         = 0.0;
        set_data.exciter[0].parameters[Ieeet1Parameters::Ka]         = 50.;
        set_data.exciter[0].parameters[Ieeet1Parameters::Ta]         = 0.04;
        set_data.exciter[0].parameters[Ieeet1Parameters::Ke]         = -0.06;
        set_data.exciter[0].parameters[Ieeet1Parameters::Te]         = 0.6;
        set_data.exciter[0].parameters[Ieeet1Parameters::Kf]         = 0.09;
        set_data.exciter[0].parameters[Ieeet1Parameters::Tf]         = 1.46;
        set_data.exciter[0].parameters[Ieeet1Parameters::Vrmin]      = -1.;
        set_data.exciter[0].parameters[Ieeet1Parameters::Vrmax]      = 1.;
        set_data.exciter[0].parameters[Ieeet1Parameters::E1]         = 2.8;
        set_data.exciter[0].parameters[Ieeet1Parameters::E2]         = 3.373;
        set_data.exciter[0].parameters[Ieeet1Parameters::Se1]        = 0.04;
        set_data.exciter[0].parameters[Ieeet1Parameters::Se2]        = 0.33;
        set_data.exciter[0].parameters[Ieeet1Parameters::Ispdlim]    = 0.;
        set_data.exciter[0].monitored_variables.insert(Ieeet1Var::efd);

        std::string      in_file   = "TwoBus/Ieeet1/TwoBusIeeet1.case.json";
        SystemModelDataT file_data = parseSystemModelData(in_file);

        auto success = compare(set_data, file_data);

        return success.report(__func__);
      }

      TestOutcome twoBusTgov1()
      {
        using namespace GridKit::PhasorDynamics;
        using namespace GridKit::PhasorDynamics::Governor;
        using BusType   = BusDataT::BusType;
        using BusVar    = BusMonitorableVariables;
        using GenrouVar = GenrouMonitorableVariables;

        SystemModelDataT set_data;

        set_data.bus.resize(2);
        set_data.bus[0].bus_id   = 0;
        set_data.bus[0].bus_type = BusType::DEFAULT;
        set_data.bus[0].Vr0      = 0.9949877346411762;
        set_data.bus[0].Vi0      = 0.09999703952427966;
        set_data.bus[0].monitored_variables.insert(BusVar::Vm);
        set_data.bus[1].bus_id   = 1;
        set_data.bus[1].bus_type = BusType::SLACK;
        set_data.bus[1].Vr0      = 1.0;
        set_data.bus[1].Vi0      = 0.0;

        set_data.signal.resize(2);
        set_data.signal[0].name      = "Machine Speed Deviation";
        set_data.signal[0].signal_id = 0;
        set_data.signal[1].name      = "Mechanical Power";
        set_data.signal[1].signal_id = 1;

        set_data.branch.resize(1);
        set_data.branch[0].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[0].buses[BranchBuses::bus2]        = set_data.bus[1].bus_id;
        set_data.branch[0].parameters[BranchParameters::R] = 0.0;
        set_data.branch[0].parameters[BranchParameters::X] = 0.1;
        set_data.branch[0].parameters[BranchParameters::G] = 0.0;
        set_data.branch[0].parameters[BranchParameters::B] = 0.0;

        set_data.bus_fault.resize(1);
        set_data.bus_fault[0].buses[BusFaultBuses::bus]              = 0;
        set_data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
        set_data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-3;
        set_data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

        set_data.genrou.resize(1);
        set_data.genrou[0].buses[GenrouBuses::bus]                    = 0;
        set_data.genrou[0].signal_outputs[GenrouSignalOutputs::speed] = 0;
        set_data.genrou[0].signal_inputs[GenrouSignalInputs::pmech]   = 1;
        set_data.genrou[0].parameters[GenrouParameters::p0]           = 1.;
        set_data.genrou[0].parameters[GenrouParameters::q0]           = 0.05013;
        set_data.genrou[0].parameters[GenrouParameters::H]            = 3.;
        set_data.genrou[0].parameters[GenrouParameters::D]            = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Ra]           = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Tdop]         = 7.;
        set_data.genrou[0].parameters[GenrouParameters::Tdopp]        = .04;
        set_data.genrou[0].parameters[GenrouParameters::Tqopp]        = .05;
        set_data.genrou[0].parameters[GenrouParameters::Tqop]         = .75;
        set_data.genrou[0].parameters[GenrouParameters::Xd]           = 2.1;
        set_data.genrou[0].parameters[GenrouParameters::Xdp]          = 0.2;
        set_data.genrou[0].parameters[GenrouParameters::Xdpp]         = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xq]           = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqp]          = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::Xqpp]         = 0.18;
        set_data.genrou[0].parameters[GenrouParameters::Xl]           = 0.15;
        set_data.genrou[0].parameters[GenrouParameters::S10]          = 0.;
        set_data.genrou[0].parameters[GenrouParameters::S12]          = 0.;
        set_data.genrou[0].parameters[GenrouParameters::mva]          = 100.0;
        set_data.genrou[0].monitored_variables.insert(GenrouVar::speed);

        set_data.gov.resize(1);
        set_data.gov[0].signal_inputs[Tgov1SignalInputs::speed]   = 0;
        set_data.gov[0].signal_outputs[Tgov1SignalOutputs::pmech] = 1;
        set_data.gov[0].parameters[Tgov1Parameters::Trate]        = 100.0;
        set_data.gov[0].parameters[Tgov1Parameters::R]            = 0.05;
        set_data.gov[0].parameters[Tgov1Parameters::Pvmin]        = 0.0;
        set_data.gov[0].parameters[Tgov1Parameters::Pvmax]        = 1.0;
        set_data.gov[0].parameters[Tgov1Parameters::T1]           = 0.5;
        set_data.gov[0].parameters[Tgov1Parameters::T2]           = 2.5;
        set_data.gov[0].parameters[Tgov1Parameters::T3]           = 7.5;
        set_data.gov[0].parameters[Tgov1Parameters::Dt]           = 0.0;

        std::string      in_file   = "TwoBus/Tgov1/TwoBusTgov1.case.json";
        SystemModelDataT file_data = parseSystemModelData(in_file);

        auto success = compare(set_data, file_data);
        return success.report(__func__);
      }

      TestOutcome threeBusBasic()
      {
        using namespace GridKit::PhasorDynamics;
        using BusType   = BusDataT::BusType;
        using BusVar    = BusMonitorableVariables;
        using GenrouVar = GenrouMonitorableVariables;

        SystemModelDataT set_data;

        set_data.bus.resize(3);
        set_data.bus[0].bus_id   = 0;
        set_data.bus[0].bus_type = BusType::SLACK;
        set_data.bus[0].Vr0      = 1.06;
        set_data.bus[0].Vi0      = 0.0;
        set_data.bus[0].monitored_variables.insert(BusVar::Vm);
        set_data.bus[1].bus_id   = 1;
        set_data.bus[1].bus_type = BusType::DEFAULT;
        set_data.bus[1].Vr0      = 1.0599558398065716;
        set_data.bus[1].Vi0      = -0.009675621941024773;
        set_data.bus[1].monitored_variables.insert(BusVar::Vm);
        set_data.bus[2].bus_id   = 2;
        set_data.bus[2].bus_type = BusType::DEFAULT;
        set_data.bus[2].Vr0      = 0.9610827543495831;
        set_data.bus[2].Vi0      = -0.13122476630506485;
        set_data.bus[2].monitored_variables.insert(BusVar::Vm);

        set_data.branch.resize(3);
        set_data.branch[0].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[0].buses[BranchBuses::bus2]        = set_data.bus[1].bus_id;
        set_data.branch[0].parameters[BranchParameters::R] = 0.05;
        set_data.branch[0].parameters[BranchParameters::X] = 0.21;
        set_data.branch[0].parameters[BranchParameters::G] = 0.0;
        set_data.branch[0].parameters[BranchParameters::B] = 0.1;
        set_data.branch[1].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[1].buses[BranchBuses::bus2]        = set_data.bus[2].bus_id;
        set_data.branch[1].parameters[BranchParameters::R] = 0.06;
        set_data.branch[1].parameters[BranchParameters::X] = 0.15;
        set_data.branch[1].parameters[BranchParameters::G] = 0.;
        set_data.branch[1].parameters[BranchParameters::B] = 0.12;
        set_data.branch[2].buses[BranchBuses::bus1]        = set_data.bus[1].bus_id;
        set_data.branch[2].buses[BranchBuses::bus2]        = set_data.bus[2].bus_id;
        set_data.branch[2].parameters[BranchParameters::R] = 0.08;
        set_data.branch[2].parameters[BranchParameters::X] = 0.27;
        set_data.branch[2].parameters[BranchParameters::G] = 0.;
        set_data.branch[2].parameters[BranchParameters::B] = 0.45;

        set_data.genrou.resize(2);
        set_data.genrou[0].buses[GenrouBuses::bus]             = 1;
        set_data.genrou[0].parameters[GenrouParameters::p0]    = 0.5;
        set_data.genrou[0].parameters[GenrouParameters::q0]    = -0.07588;
        set_data.genrou[0].parameters[GenrouParameters::H]     = 2.7;
        set_data.genrou[0].parameters[GenrouParameters::D]     = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Ra]    = 0.;
        set_data.genrou[0].parameters[GenrouParameters::Tdop]  = 7.;
        set_data.genrou[0].parameters[GenrouParameters::Tdopp] = .04;
        set_data.genrou[0].parameters[GenrouParameters::Tqopp] = .05;
        set_data.genrou[0].parameters[GenrouParameters::Tqop]  = .75;
        set_data.genrou[0].parameters[GenrouParameters::Xd]    = 1.9;
        set_data.genrou[0].parameters[GenrouParameters::Xdp]   = 0.17;
        set_data.genrou[0].parameters[GenrouParameters::Xdpp]  = 0.15;
        set_data.genrou[0].parameters[GenrouParameters::Xq]    = 0.4;
        set_data.genrou[0].parameters[GenrouParameters::Xqp]   = 0.35;
        set_data.genrou[0].parameters[GenrouParameters::Xqpp]  = 0.15;
        set_data.genrou[0].parameters[GenrouParameters::Xl]    = 0.14999;
        set_data.genrou[0].parameters[GenrouParameters::S10]   = 0.;
        set_data.genrou[0].parameters[GenrouParameters::S12]   = 0.;
        set_data.genrou[0].monitored_variables.insert(GenrouVar::speed);
        set_data.genrou[1].buses[GenrouBuses::bus]             = 2;
        set_data.genrou[1].parameters[GenrouParameters::p0]    = 0.25;
        set_data.genrou[1].parameters[GenrouParameters::q0]    = 0.26587;
        set_data.genrou[1].parameters[GenrouParameters::H]     = 1.6;
        set_data.genrou[1].parameters[GenrouParameters::D]     = 0.;
        set_data.genrou[1].parameters[GenrouParameters::Ra]    = 0.;
        set_data.genrou[1].parameters[GenrouParameters::Tdop]  = 7.5;
        set_data.genrou[1].parameters[GenrouParameters::Tdopp] = .04;
        set_data.genrou[1].parameters[GenrouParameters::Tqopp] = .05;
        set_data.genrou[1].parameters[GenrouParameters::Tqop]  = .75;
        set_data.genrou[1].parameters[GenrouParameters::Xd]    = 2.3;
        set_data.genrou[1].parameters[GenrouParameters::Xdp]   = 0.2;
        set_data.genrou[1].parameters[GenrouParameters::Xdpp]  = 0.18;
        set_data.genrou[1].parameters[GenrouParameters::Xq]    = 0.5;
        set_data.genrou[1].parameters[GenrouParameters::Xqp]   = 0.5;
        set_data.genrou[1].parameters[GenrouParameters::Xqpp]  = 0.18;
        set_data.genrou[1].parameters[GenrouParameters::Xl]    = 0.15;
        set_data.genrou[1].parameters[GenrouParameters::S10]   = 0.;
        set_data.genrou[1].parameters[GenrouParameters::S12]   = 0.;
        set_data.genrou[1].monitored_variables.insert(GenrouVar::speed);

        set_data.loadz.resize(1);
        set_data.loadz[0].buses[LoadZBuses::bus]         = 2;
        set_data.loadz[0].parameters[LoadZParameters::R] = 0.4447197839297772;
        set_data.loadz[0].parameters[LoadZParameters::X] = 0.20330047265361242;

        set_data.bus_fault.resize(1);
        set_data.bus_fault[0].buses[BusFaultBuses::bus]              = 2;
        set_data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
        set_data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-5;
        set_data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

        std::string      in_file   = "ThreeBus/Basic/ThreeBusBasic.case.json";
        SystemModelDataT file_data = parseSystemModelData(in_file);

        auto success = compare(set_data, file_data);
        return success.report(__func__);
      }

      TestOutcome threeBusClassical()
      {
        using namespace GridKit::PhasorDynamics;
        using BusType = BusDataT::BusType;
        using BusVar  = BusMonitorableVariables;
        using GenVar  = GenClassicalMonitorableVariables;

        SystemModelDataT set_data;

        set_data.bus.resize(3);
        set_data.bus[0].bus_id   = 0;
        set_data.bus[0].bus_type = BusType::SLACK;
        set_data.bus[0].Vr0      = 1.06;
        set_data.bus[0].Vi0      = 0.0;
        set_data.bus[1].bus_id   = 1;
        set_data.bus[1].bus_type = BusType::DEFAULT;
        set_data.bus[1].Vr0      = 1.0599558398065716;
        set_data.bus[1].Vi0      = -0.009675621941024773;
        set_data.bus[1].monitored_variables.insert(BusVar::Vm);
        set_data.bus[2].bus_id   = 2;
        set_data.bus[2].bus_type = BusType::DEFAULT;
        set_data.bus[2].Vr0      = 0.9610827543495831;
        set_data.bus[2].Vi0      = -0.13122476630506485;
        set_data.bus[2].monitored_variables.insert(BusVar::Vm);

        set_data.branch.resize(3);
        set_data.branch[0].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[0].buses[BranchBuses::bus2]        = set_data.bus[1].bus_id;
        set_data.branch[0].parameters[BranchParameters::R] = 0.05;
        set_data.branch[0].parameters[BranchParameters::X] = 0.21;
        set_data.branch[0].parameters[BranchParameters::G] = 0.0;
        set_data.branch[0].parameters[BranchParameters::B] = 0.1;
        set_data.branch[1].buses[BranchBuses::bus1]        = set_data.bus[0].bus_id;
        set_data.branch[1].buses[BranchBuses::bus2]        = set_data.bus[2].bus_id;
        set_data.branch[1].parameters[BranchParameters::R] = 0.06;
        set_data.branch[1].parameters[BranchParameters::X] = 0.15;
        set_data.branch[1].parameters[BranchParameters::G] = 0.;
        set_data.branch[1].parameters[BranchParameters::B] = 0.12;
        set_data.branch[2].buses[BranchBuses::bus1]        = set_data.bus[1].bus_id;
        set_data.branch[2].buses[BranchBuses::bus2]        = set_data.bus[2].bus_id;
        set_data.branch[2].parameters[BranchParameters::R] = 0.08;
        set_data.branch[2].parameters[BranchParameters::X] = 0.27;
        set_data.branch[2].parameters[BranchParameters::G] = 0.;
        set_data.branch[2].parameters[BranchParameters::B] = 0.45;

        set_data.genclassical.resize(2);
        set_data.genclassical[0].buses[GenClassicalBuses::bus]           = set_data.bus[1].bus_id;
        set_data.genclassical[0].parameters[GenClassicalParameters::p0]  = 0.5;
        set_data.genclassical[0].parameters[GenClassicalParameters::q0]  = -0.07588;
        set_data.genclassical[0].parameters[GenClassicalParameters::H]   = 2.7;
        set_data.genclassical[0].parameters[GenClassicalParameters::D]   = 0.;
        set_data.genclassical[0].parameters[GenClassicalParameters::Ra]  = 0.;
        set_data.genclassical[0].parameters[GenClassicalParameters::Xdp] = 0.17;
        set_data.genclassical[0].monitored_variables.insert(GenVar::speed);
        set_data.genclassical[1].buses[GenClassicalBuses::bus]           = set_data.bus[2].bus_id;
        set_data.genclassical[1].parameters[GenClassicalParameters::p0]  = 0.25;
        set_data.genclassical[1].parameters[GenClassicalParameters::q0]  = 0.26587;
        set_data.genclassical[1].parameters[GenClassicalParameters::H]   = 1.6;
        set_data.genclassical[1].parameters[GenClassicalParameters::D]   = 0.;
        set_data.genclassical[1].parameters[GenClassicalParameters::Ra]  = 0.;
        set_data.genclassical[1].parameters[GenClassicalParameters::Xdp] = 0.2;
        set_data.genclassical[1].monitored_variables.insert(GenVar::speed);

        set_data.loadz.resize(1);
        set_data.loadz[0].buses[LoadZBuses::bus]         = 2;
        set_data.loadz[0].parameters[LoadZParameters::R] = 0.4447197839297772;
        set_data.loadz[0].parameters[LoadZParameters::X] = 0.20330047265361242;

        set_data.bus_fault.resize(1);
        set_data.bus_fault[0].buses[BusFaultBuses::bus]              = 2;
        set_data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
        set_data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-5;
        set_data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

        std::string      in_file   = "ThreeBus/Classical/ThreeBusClassical.case.json";
        SystemModelDataT file_data = parseSystemModelData(in_file);

        auto success = compare(set_data, file_data);
        return success.report(__func__);
      }

      /// A finite plant-reference pulse propagates through the public REPCA,
      /// REECB, and REGCA signals before the closed loop returns to equilibrium.
      TestOutcome regcaReecbRepca()
      {
        using namespace GridKit::PhasorDynamics::Controller;
        using namespace GridKit::PhasorDynamics::Converter;

        // Signal IDs are local graph keys. Their enum values only provide a
        // compact set of unique identifiers; the numbers carry no model meaning.
        enum class LoopSignal : size_t
        {
          IPCMD,
          IQCMD,
          IBRANCHR,
          IBRANCHI,
          PBRANCH,
          QBRANCH,
          QEXT,
          PEXT,
          PLANT_PREF,
        };
        const auto signalId = [](LoopSignal signal)
        {
          return static_cast<IdxT>(signal);
        };

        constexpr IdxT  RENEWABLE_BUS_ID       = static_cast<IdxT>(0);
        constexpr RealT SYSTEM_MVA             = static_cast<RealT>(100.0);
        // The 2:1 base ratio exercises each controller power-base conversion.
        constexpr RealT COMPONENT_MVA          = SYSTEM_MVA / static_cast<RealT>(2.0);
        // A nonzero interior operating point stays clear of controller limits.
        constexpr RealT INITIAL_ACTIVE_POWER   = static_cast<RealT>(0.4);
        constexpr RealT INITIAL_REACTIVE_POWER = static_cast<RealT>(0.05);
        constexpr RealT FAST_TIME_CONSTANT     = static_cast<RealT>(0.02);
        constexpr RealT PLANT_TIME_CONSTANT    = static_cast<RealT>(0.5);
        constexpr RealT REFERENCE_PULSE        = static_cast<RealT>(0.05);
        // Five fast time constants let the perturbation propagate through the
        // converter and electrical controller. One fifth of the pulse proves a
        // material response without pinning an exact transient trajectory.
        constexpr RealT PULSE_DURATION         = static_cast<RealT>(5.0) * FAST_TIME_CONSTANT;
        constexpr RealT MINIMUM_RESPONSE       = REFERENCE_PULSE / static_cast<RealT>(5.0);
        // Preserve the conservative fifty-plant-time-constant recovery window.
        constexpr RealT RECOVERY_HORIZON       = static_cast<RealT>(50.0) * PLANT_TIME_CONSTANT;
        // This is an output cadence; IDA still selects its internal steps.
        constexpr RealT OUTPUT_INTERVAL        = ONE<RealT> / static_cast<RealT>(60.0);
        constexpr RealT RECOVERY_TOLERANCE     = static_cast<RealT>(1.0e-6);

        TestStatus success = true;

        SystemModelDataT data;
        data.va_base = SYSTEM_MVA * static_cast<RealT>(1.0e6);

        auto& bus    = data.bus.emplace_back();
        bus.bus_id   = RENEWABLE_BUS_ID;
        bus.bus_type = BusDataT::BusType::SLACK;
        bus.Vr0      = ONE<RealT>;
        bus.Vi0      = ZERO<RealT>;

        data.signal = {{"Active Current Command", signalId(LoopSignal::IPCMD)},
                       {"Reactive Current Command", signalId(LoopSignal::IQCMD)},
                       {"Branch Current Real", signalId(LoopSignal::IBRANCHR)},
                       {"Branch Current Imaginary", signalId(LoopSignal::IBRANCHI)},
                       {"Branch Active Power", signalId(LoopSignal::PBRANCH)},
                       {"Branch Reactive Power", signalId(LoopSignal::QBRANCH)},
                       {"Reactive Power Command", signalId(LoopSignal::QEXT)},
                       {"Active Power Command", signalId(LoopSignal::PEXT)},
                       {"Plant Active Power Reference", signalId(LoopSignal::PLANT_PREF)}};

        auto& converter                                        = data.regca.emplace_back();
        converter.buses[RegcaBuses::bus]                       = RENEWABLE_BUS_ID;
        converter.signal_inputs[RegcaSignalInputs::ipcmd]      = signalId(LoopSignal::IPCMD);
        converter.signal_inputs[RegcaSignalInputs::iqcmd]      = signalId(LoopSignal::IQCMD);
        converter.signal_outputs[RegcaSignalOutputs::ibranchr] = signalId(LoopSignal::IBRANCHR);
        converter.signal_outputs[RegcaSignalOutputs::ibranchi] = signalId(LoopSignal::IBRANCHI);
        converter.signal_outputs[RegcaSignalOutputs::pbranch]  = signalId(LoopSignal::PBRANCH);
        converter.signal_outputs[RegcaSignalOutputs::qbranch]  = signalId(LoopSignal::QBRANCH);
        converter.parameters[RegcaParameters::p0]              = INITIAL_ACTIVE_POWER;
        converter.parameters[RegcaParameters::q0]              = INITIAL_REACTIVE_POWER;
        converter.parameters[RegcaParameters::mva]             = COMPONENT_MVA;
        converter.parameters[RegcaParameters::Tg]              = FAST_TIME_CONSTANT;
        converter.parameters[RegcaParameters::TM]              = FAST_TIME_CONSTANT;
        converter.parameters[RegcaParameters::Rqmax]           = static_cast<RealT>(999.0);
        converter.parameters[RegcaParameters::Rqmin]           = static_cast<RealT>(-999.0);
        converter.parameters[RegcaParameters::Rpmax]           = static_cast<RealT>(999.0);
        converter.parameters[RegcaParameters::sL]              = true;
        converter.parameters[RegcaParameters::IL1]             = static_cast<RealT>(1.1);
        converter.parameters[RegcaParameters::VL0]             = static_cast<RealT>(0.4);
        converter.parameters[RegcaParameters::VL1]             = static_cast<RealT>(0.9);
        converter.parameters[RegcaParameters::VA0]             = static_cast<RealT>(0.4);
        converter.parameters[RegcaParameters::VA1]             = static_cast<RealT>(0.9);
        converter.parameters[RegcaParameters::Vhvmax]          = static_cast<RealT>(1.2);

        auto& controller                                     = data.reecb.emplace_back();
        controller.buses[ReecbBuses::bus]                    = RENEWABLE_BUS_ID;
        controller.signal_inputs[ReecbSignalInputs::pe]      = signalId(LoopSignal::PBRANCH);
        controller.signal_inputs[ReecbSignalInputs::qgen]    = signalId(LoopSignal::QBRANCH);
        controller.signal_inputs[ReecbSignalInputs::qext]    = signalId(LoopSignal::QEXT);
        controller.signal_inputs[ReecbSignalInputs::pref]    = signalId(LoopSignal::PEXT);
        controller.signal_outputs[ReecbSignalOutputs::ipcmd] = signalId(LoopSignal::IPCMD);
        controller.signal_outputs[ReecbSignalOutputs::iqcmd] = signalId(LoopSignal::IQCMD);
        controller.parameters[ReecbParameters::mva]          = COMPONENT_MVA;
        controller.parameters[ReecbParameters::Trv]          = FAST_TIME_CONSTANT;
        controller.parameters[ReecbParameters::Tp]           = FAST_TIME_CONSTANT;
        controller.parameters[ReecbParameters::Kvi]          = static_cast<RealT>(5.0);
        controller.parameters[ReecbParameters::QFlag]        = true;
        controller.parameters[ReecbParameters::VFlag]        = true;

        auto& plant                                    = data.repca.emplace_back();
        plant.buses[RepcaBuses::bus]                   = RENEWABLE_BUS_ID;
        plant.signal_inputs[RepcaSignalInputs::ir]     = signalId(LoopSignal::IBRANCHR);
        plant.signal_inputs[RepcaSignalInputs::ii]     = signalId(LoopSignal::IBRANCHI);
        plant.signal_inputs[RepcaSignalInputs::p]      = signalId(LoopSignal::PBRANCH);
        plant.signal_inputs[RepcaSignalInputs::q]      = signalId(LoopSignal::QBRANCH);
        plant.signal_inputs[RepcaSignalInputs::pref]   = signalId(LoopSignal::PLANT_PREF);
        plant.signal_outputs[RepcaSignalOutputs::qext] = signalId(LoopSignal::QEXT);
        plant.signal_outputs[RepcaSignalOutputs::pext] = signalId(LoopSignal::PEXT);
        plant.parameters[RepcaParameters::mva]         = COMPONENT_MVA;
        plant.parameters[RepcaParameters::Freqflag]    = true;
        plant.parameters[RepcaParameters::Ddn]         = ZERO<RealT>;
        plant.parameters[RepcaParameters::Dup]         = ZERO<RealT>;
        plant.parameters[RepcaParameters::Tp]          = FAST_TIME_CONSTANT;
        plant.parameters[RepcaParameters::Tlag]        = PLANT_TIME_CONSTANT;

        auto& reference                                          = data.constant_source.emplace_back();
        reference.parameters[ConstantSignalSourceParameters::Sr] = INITIAL_ACTIVE_POWER;
        reference.signal_outputs[ConstantSignalSourceSignalOutputs::sr] =
            signalId(LoopSignal::PLANT_PREF);

        SystemModel<RealT, IdxT> system(data);
        success *= system.allocate() == 0;

        auto* pref_signal    = system.getSignal(signalId(LoopSignal::PLANT_PREF));
        auto* pext_signal    = system.getSignal(signalId(LoopSignal::PEXT));
        auto* qext_signal    = system.getSignal(signalId(LoopSignal::QEXT));
        auto* ipcmd_signal   = system.getSignal(signalId(LoopSignal::IPCMD));
        auto* pbranch_signal = system.getSignal(signalId(LoopSignal::PBRANCH));
        if (pref_signal == nullptr || pext_signal == nullptr || qext_signal == nullptr
            || ipcmd_signal == nullptr || pbranch_signal == nullptr)
        {
          std::cout << "Renewable-control loop is missing a required signal\n";
          success = false;
          return success.report(__func__);
        }

        AnalysisManager::Sundials::Ida<RealT, IdxT> ida(&system);
        success *= ida.configureSimulation() == 0;
        success *= ida.initializeSimulation(ZERO<RealT>) == 0;

        const auto*              equilibrium_values = system.y().getData();
        const std::vector<RealT> equilibrium(
            equilibrium_values,
            equilibrium_values + static_cast<size_t>(system.y().getSize()));

        success *= isEqual(static_cast<RealT>(pref_signal->read()),
                           INITIAL_ACTIVE_POWER,
                           RECOVERY_TOLERANCE);
        success *= isEqual(static_cast<RealT>(pext_signal->read()),
                           INITIAL_ACTIVE_POWER,
                           RECOVERY_TOLERANCE);
        success *= isEqual(static_cast<RealT>(qext_signal->read()),
                           INITIAL_REACTIVE_POWER,
                           RECOVERY_TOLERANCE);

        const RealT pref0    = static_cast<RealT>(pref_signal->read());
        const RealT pext0    = static_cast<RealT>(pext_signal->read());
        const RealT ipcmd0   = static_cast<RealT>(ipcmd_signal->read());
        const RealT pbranch0 = static_cast<RealT>(pbranch_signal->read());

        pref_signal->init(pref0 + REFERENCE_PULSE);
        success *= ida.initializeSimulation(ZERO<RealT>) == 0;
        success *= ida.runSimulation(PULSE_DURATION, OUTPUT_INTERVAL) == 0;

        struct SignalResponse
        {
          const char* label;
          RealT       value;
        };

        const std::array<SignalResponse, 3> responses{{
            {"REPCA PEXT", static_cast<RealT>(pext_signal->read()) - pext0},
            {"REECB IPCMD", static_cast<RealT>(ipcmd_signal->read()) - ipcmd0},
            {"REGCA PBRANCH", static_cast<RealT>(pbranch_signal->read()) - pbranch0},
        }};
        for (const auto& response : responses)
        {
          if (response.value <= MINIMUM_RESPONSE)
          {
            std::cout << response.label << " responded by only " << response.value
                      << " during the reference pulse\n";
            success = false;
          }
        }

        pref_signal->init(pref0);
        success *= ida.initializeSimulation(PULSE_DURATION) == 0;
        success *= ida.runSimulation(PULSE_DURATION + RECOVERY_HORIZON,
                                     OUTPUT_INTERVAL)
                   == 0;

        const auto* final_values = system.y().getData();
        for (size_t entry = 0; entry < equilibrium.size(); ++entry)
        {
          const RealT deviation = final_values[entry] - equilibrium[entry];
          if (!isEqual(deviation, ZERO<RealT>, RECOVERY_TOLERANCE))
          {
            std::cout << "State " << entry << " remains " << deviation
                      << " from its equilibrium after recovery\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
