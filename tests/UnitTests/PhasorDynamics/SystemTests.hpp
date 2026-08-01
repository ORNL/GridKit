#pragma once

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/StateDataAdapter.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using GridKit::PhasorDynamics::BranchBuses;
    using GridKit::PhasorDynamics::BranchParameters;

    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class SystemTests
    {
    private:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      SystemTests()  = default;
      ~SystemTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = nullptr;

        // Create an empty system
        system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        if (system == nullptr)
        {
          std::cout << "Default constructor failed!\n";
          success = false;
          return success.report(__func__);
        }

        delete system;
        system = nullptr;

        PhasorDynamics::SystemModelData<ScalarT, IdxT> data;

        // Set bus data
        data.bus.resize(2);

        // Bus 0
        data.bus[0].bus_id            = 0;
        data.bus[0].bus_type          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].initial_state     = Model::BusState{};
        data.bus[0].initial_state->vr = 10.0;
        data.bus[0].initial_state->vi = 20.0;

        // Bus 1
        data.bus[1].bus_id            = 1;
        data.bus[1].bus_type          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].initial_state     = Model::BusState{};
        data.bus[1].initial_state->vr = 30.0;
        data.bus[1].initial_state->vi = 40.0;

        // Set branch data
        data.branch.resize(1);

        // Branch 0-1
        data.branch[0].buses[BranchBuses::bus1]        = data.bus[0].bus_id;
        data.branch[0].buses[BranchBuses::bus2]        = data.bus[1].bus_id;
        data.branch[0].parameters[BranchParameters::R] = 2.0;
        data.branch[0].parameters[BranchParameters::X] = 4.0;
        data.branch[0].parameters[BranchParameters::G] = 0.2;
        data.branch[0].parameters[BranchParameters::B] = 1.2;

        // Create an empty system model
        system = new PhasorDynamics::SystemModel<ScalarT, IdxT>(data);
        system->allocate();
        system->initialize();
        system->evaluateResidual();

        // Answer keys
        const ScalarT Ir0{17.0};  ///< Solution: real current entering bus-0
        const ScalarT Ii0{-10.0}; ///< Solution: imaginary current entering bus-0
        const ScalarT Ir1{15.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-20.0}; ///< Solution: imaginary current entering bus-1

        auto* bus0 = system->getBus(0);
        auto* bus1 = system->getBus(1);

        success *= isEqual(bus0->Ir(), Ir0);
        success *= isEqual(bus0->Ii(), Ii0);
        success *= isEqual(bus1->Ir(), Ir1);
        success *= isEqual(bus1->Ii(), Ii1);

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome stateDataAttachment()
      {
        using namespace PhasorDynamics;

        SystemModelData<> model_data;
        model_data.bus.resize(2);
        model_data.bus[0].bus_id = 7;
        model_data.bus[1].bus_id = 8;

        Model::StateData state_data;
        state_data.buses["bus_id_7"].vr = 1.025;
        state_data.buses.emplace("bus_id_999", Model::BusState{});

        double state_value = 1.0;
        auto   addDevice   = [&](auto& devices, const char* id)
        {
          devices.resize(1);
          devices[0].disambiguation_string = id;
          state_data.devices[id].p         = state_value++;
        };

        addDevice(model_data.branch, "branch_1");
        addDevice(model_data.genrou, "genrou_1");
        addDevice(model_data.gensal, "gensal_1");
        addDevice(model_data.genclassical, "genclassical_1");
        addDevice(model_data.loadz, "loadz_1");
        addDevice(model_data.loadzip, "loadzip_1");
        addDevice(model_data.bus_fault, "fault_1");
        state_data.devices["fault_1"].active                 = true;
        model_data.branch[0].parameters[BranchParameters::R] = 0.1;
        state_data.devices.emplace("unknown_1", Model::DeviceState{});

        applyState(model_data, state_data);

        TestStatus success  = true;
        success            *= model_data.bus[0].initial_state.has_value();
        if (model_data.bus[0].initial_state.has_value())
        {
          success *= model_data.bus[0].initial_state->vr == 1.025;
        }
        success *= !model_data.bus[1].initial_state.has_value();

        auto checkDevice = [&](const auto& devices, const char* id)
        {
          success *= devices[0].initial_state.has_value();
          if (devices[0].initial_state.has_value())
          {
            success *= devices[0].initial_state->p == state_data.devices.at(id).p;
          }
        };

        checkDevice(model_data.branch, "branch_1");
        checkDevice(model_data.genrou, "genrou_1");
        checkDevice(model_data.gensal, "gensal_1");
        checkDevice(model_data.genclassical, "genclassical_1");
        checkDevice(model_data.loadz, "loadz_1");
        checkDevice(model_data.loadzip, "loadzip_1");
        checkDevice(model_data.bus_fault, "fault_1");
        success *= model_data.bus_fault[0].initial_state->active == true;
        success *= std::get<RealT>(model_data.branch[0].parameters[BranchParameters::R]) == 0.1;

        return success.report(__func__);
      }

      TestOutcome busFaultStateInput()
      {
        using namespace PhasorDynamics;

        SystemModelData<ScalarT, IdxT> data;
        data.bus.resize(1);
        auto& bus             = data.bus[0];
        bus.bus_id            = 0;
        bus.bus_type          = BusData<ScalarT, IdxT>::BusType::SLACK;
        bus.initial_state     = Model::BusState{};
        bus.initial_state->vr = 1.0;
        bus.initial_state->vi = 0.0;

        data.bus_fault.resize(1);
        auto& fault                             = data.bus_fault[0];
        fault.disambiguation_string             = "fault_1";
        fault.buses[BusFaultBuses::bus]         = 0;
        fault.parameters[BusFaultParameters::R] = 0.0;
        fault.parameters[BusFaultParameters::X] = 1.0e-3;
        fault.initial_state                     = Model::DeviceState{};
        fault.initial_state->active             = false;

        SystemModel<ScalarT, IdxT> system(data);
        system.allocate();
        system.initialize();
        system.evaluateResidual();

        TestStatus success  = true;
        success            *= isEqual(system.getBus(0)->Ir(), ScalarT{0.0});
        success            *= isEqual(system.getBus(0)->Ii(), ScalarT{0.0});

        system.setInput("fault_1", "active", 1.0);
        system.initialize();
        system.evaluateResidual();
        success *= isEqual(system.getBus(0)->Ir(), ScalarT{0.0});
        success *= isEqual(system.getBus(0)->Ii(), ScalarT{1000.0});

        system.setInput("fault_1", "active", 0.0);
        system.initialize();
        system.evaluateResidual();
        success *= isEqual(system.getBus(0)->Ir(), ScalarT{0.0});
        success *= isEqual(system.getBus(0)->Ii(), ScalarT{0.0});

        return success.report(__func__);
      }

      TestOutcome busStateInputs()
      {
        using namespace PhasorDynamics;

        SystemModelData<ScalarT, IdxT> data;
        data.bus.resize(2);

        data.bus[0].bus_id            = 4;
        data.bus[0].bus_type          = BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].initial_state     = Model::BusState{};
        data.bus[0].initial_state->vr = 1.1;
        data.bus[0].initial_state->vi = -0.1;

        data.bus[1].bus_id            = 5;
        data.bus[1].bus_type          = BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].initial_state     = Model::BusState{};
        data.bus[1].initial_state->vr = 0.9;
        data.bus[1].initial_state->vi = 0.2;

        SystemModel<ScalarT, IdxT> system(data);
        system.allocate();
        system.initialize();

        TestStatus success  = true;
        success            *= isEqual(system.getBus(4)->Vr(), ScalarT{1.1});
        success            *= isEqual(system.getBus(4)->Vi(), ScalarT{-0.1});
        success            *= isEqual(system.getBus(5)->Vr(), ScalarT{0.9});
        success            *= isEqual(system.getBus(5)->Vi(), ScalarT{0.2});

        return success.report(__func__);
      }

      TestOutcome deviceInputSignals()
      {
        using namespace PhasorDynamics;

        SystemModelData<ScalarT, IdxT> data;
        data.bus.resize(2);
        data.bus[0].bus_id            = 0;
        data.bus[0].bus_type          = BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].initial_state     = Model::BusState{};
        data.bus[0].initial_state->vr = 10.0;
        data.bus[0].initial_state->vi = 20.0;
        data.bus[1].bus_id            = 1;
        data.bus[1].bus_type          = BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].initial_state     = Model::BusState{};
        data.bus[1].initial_state->vr = 30.0;
        data.bus[1].initial_state->vi = 40.0;

        data.branch.resize(1);
        auto& branch                           = data.branch[0];
        branch.disambiguation_string           = "branch_0_1";
        branch.buses[BranchBuses::bus1]        = 0;
        branch.buses[BranchBuses::bus2]        = 1;
        branch.parameters[BranchParameters::R] = 2.0;
        branch.parameters[BranchParameters::X] = 4.0;
        branch.parameters[BranchParameters::G] = 0.2;
        branch.parameters[BranchParameters::B] = 1.2;
        branch.initial_state                   = Model::DeviceState{};
        branch.initial_state->tap              = 1.0;
        branch.initial_state->phase            = 0.0;
        branch.initial_state->open             = true;

        SystemModel<ScalarT, IdxT> system(data);
        system.allocate();
        system.initialize();
        system.evaluateResidual();

        TestStatus success  = true;
        success            *= isEqual(system.getBus(0)->Ir(), ScalarT{0.0});
        success            *= isEqual(system.getBus(0)->Ii(), ScalarT{0.0});
        success            *= isEqual(system.getBus(1)->Ir(), ScalarT{0.0});
        success            *= isEqual(system.getBus(1)->Ii(), ScalarT{0.0});

        system.setInput("branch_0_1", "open", 0.0);
        system.evaluateResidual();

        // Reclosing uses the tap and phase supplied by State.
        success *= isEqual(system.getBus(0)->Ir(), ScalarT{17.0});
        success *= isEqual(system.getBus(0)->Ii(), ScalarT{-10.0});
        success *= isEqual(system.getBus(1)->Ir(), ScalarT{15.0});
        success *= isEqual(system.getBus(1)->Ii(), ScalarT{-20.0});

        system.setInput("branch_0_1", "open", 1.0);
        system.evaluateResidual();
        success *= isEqual(system.getBus(0)->Ir(), ScalarT{0.0});
        success *= isEqual(system.getBus(0)->Ii(), ScalarT{0.0});
        success *= isEqual(system.getBus(1)->Ir(), ScalarT{0.0});
        success *= isEqual(system.getBus(1)->Ii(), ScalarT{0.0});

        return success.report(__func__);
      }

      TestOutcome loadZIPStateInputs()
      {
        using namespace PhasorDynamics;

        auto makeData = []
        {
          SystemModelData<ScalarT, IdxT> data;
          data.bus.resize(1);
          data.bus[0].bus_id            = 9;
          data.bus[0].bus_type          = BusData<ScalarT, IdxT>::BusType::SLACK;
          data.bus[0].initial_state     = Model::BusState{};
          data.bus[0].initial_state->vr = 0.6;
          data.bus[0].initial_state->vi = 0.8;

          data.loadzip.resize(1);
          auto& load                                 = data.loadzip[0];
          load.disambiguation_string                 = "loadzip_9";
          load.buses[LoadZIPBuses::bus]              = 9;
          load.parameters[LoadZIPParameters::alphaI] = 0.2;
          load.parameters[LoadZIPParameters::alphaP] = 0.4;
          return data;
        };

        TestStatus success = true;

        // State p and q are terminal injections at the stored voltage.
        auto state_data                             = makeData();
        state_data.loadzip[0].initial_state         = Model::DeviceState{};
        state_data.loadzip[0].initial_state->p      = -1.0;
        state_data.loadzip[0].initial_state->q      = -0.25;
        state_data.loadzip[0].initial_state->online = true;
        state_data.loadzip[0].initial_state->open   = true;

        SystemModel<ScalarT, IdxT> state_system(state_data);
        state_system.allocate();
        state_system.initialize();
        state_system.evaluateResidual();
        success *= isEqual(state_system.getBus(9)->Ir(), ScalarT{-0.8});
        success *= isEqual(state_system.getBus(9)->Ii(), ScalarT{-0.65});

#ifdef GRIDKIT_ENABLE_ENZYME
        auto* jacobian  = state_system.getCsrJacobian();
        success        *= jacobian != nullptr;
        if (jacobian == nullptr)
        {
          return success.report(__func__);
        }

        const IdxT        nnz = jacobian->getNnz();
        std::vector<IdxT> rows(
            jacobian->getRowData(), jacobian->getRowData() + state_system.size() + 1);
        std::vector<IdxT> cols(
            jacobian->getColData(), jacobian->getColData() + nnz);

        state_system.setInput("loadzip_9", "p", 0.0);
        state_system.setInput("loadzip_9", "q", 0.0);
        state_system.initialize();
        state_system.evaluateResidual();
        state_system.evaluateJacobian();
        success *= jacobian->getNnz() == nnz;
        for (IdxT i = 0; i < state_system.size() + 1; ++i)
        {
          success *= jacobian->getRowData()[i] == rows[static_cast<size_t>(i)];
        }
        for (IdxT i = 0; i < nnz; ++i)
        {
          success *= jacobian->getColData()[i] == cols[static_cast<size_t>(i)];
        }

        state_system.setInput("loadzip_9", "p", -1.0);
        state_system.setInput("loadzip_9", "q", -0.25);
#endif

        state_system.setInput("loadzip_9", "online", 0.0);
        state_system.evaluateResidual();
        success *= isEqual(state_system.getBus(9)->Ir(), ScalarT{0.0});
        success *= isEqual(state_system.getBus(9)->Ii(), ScalarT{0.0});

        state_system.setInput("loadzip_9", "p", -2.0);
        state_system.setInput("loadzip_9", "q", -0.5);
        state_system.setInput("loadzip_9", "online", 1.0);
        state_system.initialize();
        state_system.evaluateResidual();
        success *= isEqual(state_system.getBus(9)->Ir(), ScalarT{-1.6});
        success *= isEqual(state_system.getBus(9)->Ii(), ScalarT{-1.3});

        return success.report(__func__);
      }

      TestOutcome generatorStateInputs()
      {
        using namespace PhasorDynamics;

        SystemModelData<ScalarT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id            = 3;
        data.bus[0].bus_type          = BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].initial_state     = Model::BusState{};
        data.bus[0].initial_state->vr = 1.0;
        data.bus[0].initial_state->vi = 0.0;

        data.genclassical.resize(1);
        auto& gen                                   = data.genclassical[0];
        gen.disambiguation_string                   = "genclassical_3";
        gen.buses[GenClassicalBuses::bus]           = 3;
        gen.parameters[GenClassicalParameters::H]   = 0.5;
        gen.parameters[GenClassicalParameters::D]   = 0.0;
        gen.parameters[GenClassicalParameters::Ra]  = 0.0;
        gen.parameters[GenClassicalParameters::Xdp] = 0.5;
        gen.initial_state                           = Model::DeviceState{};
        gen.initial_state->p                        = 0.4;
        gen.initial_state->q                        = -0.2;
        gen.initial_state->online                   = false;
        gen.initial_state->tap                      = 7.0; // Irrelevant state is harmless.

        SystemModel<ScalarT, IdxT> system(data);
        system.allocate();
        system.initialize();
        system.evaluateResidual();

        TestStatus success  = true;
        success            *= isEqual(system.getBus(3)->Ir(), ScalarT{0.0});
        success            *= isEqual(system.getBus(3)->Ii(), ScalarT{0.0});

        system.setInput("genclassical_3", "online", 1.0);
        system.evaluateResidual();
        success *= isEqual(system.getBus(3)->Ir(), ScalarT{0.4});
        success *= isEqual(system.getBus(3)->Ii(), ScalarT{0.2});

        return success.report(__func__);
      }

      TestOutcome composer()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        ScalarT Vr1{10.0}; ///< Bus-1 real voltage
        ScalarT Vi1{20.0}; ///< Bus-1 imaginary voltage
        ScalarT Vr2{30.0}; ///< Bus-2 real voltage
        ScalarT Vi2{40.0}; ///< Bus-2 imaginary voltage

        const ScalarT Ir1{17.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
        const ScalarT Ir2{15.0};  ///< Solution: real current entering bus-2
        const ScalarT Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

        // Create an empty system model
        PhasorDynamics::SystemModel<ScalarT, IdxT> system;

        // Add a bus
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(Vr1, Vi1);
        system.addBus(&bus1);

        // Add a bus
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(Vr2, Vi2);
        system.addBus(&bus2);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, R, X, G, B);
        system.addComponent(&branch);

        system.allocate();
        system.initialize();
        system.evaluateResidual();

        success *= isEqual(bus1.Ir(), Ir1);
        success *= isEqual(bus1.Ii(), Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        return success.report(__func__);
      }

      TestOutcome reallocateAfterTopologyChange()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus1(1.0, 0.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         bus2(1.0, 0.0);
        PhasorDynamics::BusFault<ScalarT, IdxT>    fault(&bus1);

        system.addBus(&bus1);
        system.addComponent(&fault);
        success                    *= system.allocate() == 0;
        const IdxT size_before_bus  = system.size();

        system.addBus(&bus2);
        success *= system.allocate() == 0;
        success *= system.size() == size_before_bus + bus2.size();

#ifdef GRIDKIT_ENABLE_ENZYME
        const auto* jacobian  = system.getCsrJacobian();
        success              *= jacobian != nullptr;

        IdxT nnz_without_branch = 0;
        if (jacobian != nullptr)
        {
          success            *= jacobian->getNumRows() == system.size();
          success            *= jacobian->getNumColumns() == system.size();
          nnz_without_branch  = jacobian->getNnz();
        }
#endif

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2);
        system.addComponent(&branch);
        success *= system.allocate() == 0;
        success *= system.evaluateJacobian() == 0;

#ifdef GRIDKIT_ENABLE_ENZYME
        jacobian  = system.getCsrJacobian();
        success  *= jacobian != nullptr;
        if (jacobian != nullptr)
        {
          success *= jacobian->getNumRows() == system.size();
          success *= jacobian->getNumColumns() == system.size();
          success *= jacobian->getNnz() > nnz_without_branch;
        }
#endif

        return success.report(__func__);
      }

      TestOutcome modelVectorsAliasSystemStorage()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus1(1.0, 2.0);
        PhasorDynamics::BusInfinite<ScalarT, IdxT> infinite_bus;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus2(3.0, 4.0);
        PhasorDynamics::Branch<ScalarT, IdxT>      branch(&bus1, &bus2);
        PhasorDynamics::LoadZ<ScalarT, IdxT>       load(&bus2, 1.0, 1.0);

        system.addBus(&bus1);
        system.addBus(&infinite_bus);
        system.addBus(&bus2);
        system.addComponent(&branch);
        system.addComponent(&load);

        if (system.allocate() != 0
            || system.setAbsoluteTolerance(1e-4) != 0)
        {
          success = false;
          return success.report(__func__);
        }

        auto checkAlias = [&](auto& system_vector, auto& model_vector, IdxT offset)
        {
          auto*      system_data = system_vector.getData();
          auto*      model_data  = model_vector.getData();
          const auto first       = static_cast<std::size_t>(offset);

          if (!system_data || model_data != system_data + first)
          {
            success = false;
            return;
          }

          success *= system_vector.setToConst(ScalarT{3.0}) == 0;
          success *= isEqual(model_data[0], ScalarT{3.0});

          success *= model_vector.setToConst(ScalarT{4.0}) == 0;
          success *= isEqual(system_data[first], ScalarT{4.0});
        };

        auto checkModel = [&](auto& model, IdxT offset)
        {
          success *= model.getVariableIndex(0) == offset;
          success *= model.getResidualIndex(0) == offset;

          checkAlias(system.y(), model.y(), offset);
          checkAlias(system.yp(), model.yp(), offset);
          checkAlias(system.getResidual(), model.getResidual(), offset);
          checkAlias(system.absoluteTolerance(), model.absoluteTolerance(), offset);
        };

        const IdxT bus2_offset = bus1.size();
        const IdxT load_offset = bus1.size() + bus2.size();
        const auto bus2_first  = static_cast<std::size_t>(bus2_offset);

        auto rebind = [&](auto& model, IdxT offset)
        {
          return model.bind(system.y(),
                            system.yp(),
                            system.getResidual(),
                            system.absoluteTolerance(),
                            offset);
        };

        // Rebinding the same slices is a no-op.
        success *= rebind(bus2, bus2_offset) == 0;
        success *= rebind(load, load_offset) == 0;

        checkModel(bus2, bus2_offset);
        checkModel(load, load_offset);

        // Tags remain model-owned and are collected separately.
        system.tag()[bus2_first]  = true;
        success                  *= system.tagDifferentiable() == 0;
        success                  *= !system.tag()[bus2_first];

        bus2.tag()[0]  = true;
        success       *= !system.tag()[bus2_first];

        return success.report(__func__);
      }

      /**
       * @brief Test for exception when signals are incorrectly configured
       */
      TestOutcome signalError()
      {
        using namespace std::filesystem;
        using namespace GridKit::PhasorDynamics;
        auto input_file = current_path() / "ThreeBusBasicBad.json";
        auto data       = parseSystemModelData(input_file);
        auto sys        = SystemModel<double, size_t>(data);

        TestStatus status{true};
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing for exceptions when signals are incorrectly configured. "
                    << "Logged errors are are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);
        status *= throws<std::runtime_error>(
            [&]()
            { sys.allocate(); });

        return status.report(__func__);
      }

      /**
       * @brief Test for exception when a child cannot bind to system storage
       */
      TestOutcome allocationError()
      {
        using namespace GridKit::PhasorDynamics;

        TestStatus                 status{true};
        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus(ScalarT{1.0}, ScalarT{0.0});

        status *= bus.allocate() == 0;
        system.addBus(&bus);
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing for exceptions when when a child cannot bind to system storage. "
                    << "Logged errors are are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);
        status *= throws<std::runtime_error>(
            [&]()
            { system.allocate(); });

        return status.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModelData<ScalarT, IdxT> data;

        // Set bus data
        data.bus.resize(2);

        // Bus 0
        data.bus[0].bus_id            = 0;
        data.bus[0].bus_type          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].initial_state     = Model::BusState{};
        data.bus[0].initial_state->vr = 10.0;
        data.bus[0].initial_state->vi = 20.0;

        // Bus 1
        data.bus[1].bus_id            = 1;
        data.bus[1].bus_type          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].initial_state     = Model::BusState{};
        data.bus[1].initial_state->vr = 30.0;
        data.bus[1].initial_state->vi = 40.0;

        // Set branch data
        data.branch.resize(1);

        // Branch 0-1
        data.branch[0].buses[BranchBuses::bus1]        = data.bus[0].bus_id;
        data.branch[0].buses[BranchBuses::bus2]        = data.bus[1].bus_id;
        data.branch[0].parameters[BranchParameters::R] = 2.0;
        data.branch[0].parameters[BranchParameters::X] = 4.0;
        data.branch[0].parameters[BranchParameters::G] = 0.2;
        data.branch[0].parameters[BranchParameters::B] = 1.2;

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian(data);

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian(data);

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          PhasorDynamics::SystemModelData<ScalarT, IdxT> data)
      {
        // Create an empty system model
        PhasorDynamics::SystemModel<DependencyTracking::Variable, IdxT> system(data);

        // Allocate and initialize the system
        system.allocate();
        system.initialize();

        // Set independent variables
        auto* y = system.y().getData();
        for (size_t i = 0; i < system.size(); ++i)
        {
          y[i].setVariableNumber(i);
        }
        system.y().setDataUpdated();

        // Evaluate and get the system residuals
        system.evaluateResidual();
        auto&       residual      = system.getResidual();
        const auto* residual_data = residual.getData();

        // Print the dependencies
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          std::cout << i << "th residual: ";
          residual_data[i].print(std::cout);
          std::cout << "\n";
        }

        // Extract the dependencies
        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual.getSize());
        for (IdxT i = 0; i < residual.getSize(); ++i)
        {
          dependencies[i] = residual_data[i].getDependencies();
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          PhasorDynamics::SystemModelData<ScalarT, IdxT> data)
      {
        // Create an empty system model
        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        // Allocate and initialize the system
        system.allocate();
        system.initialize();

        // Evaluate and get the system Jacobian
        system.evaluateResidual();
        system.evaluateJacobian();
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* system_jacobian = system.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: System Jacobian\n";
        system_jacobian->print();

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
