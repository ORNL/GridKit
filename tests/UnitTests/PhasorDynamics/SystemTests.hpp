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
        data.bus[0].bus_id   = 0;
        data.bus[0].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = 10.0;
        data.bus[0].Vi0      = 20.0;

        // Bus 1
        data.bus[1].bus_id   = 1;
        data.bus[1].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].Vr0      = 30.0;
        data.bus[1].Vi0      = 40.0;

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
        model_data.bus[0].Vr0    = 0.9;
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
        model_data.branch[0].parameters[BranchParameters::R] = 0.1;
        state_data.devices.emplace("unknown_1", Model::DeviceState{});

        applyState(model_data, state_data);

        TestStatus success  = true;
        success            *= model_data.bus[0].initial_state.has_value();
        if (model_data.bus[0].initial_state.has_value())
        {
          success *= model_data.bus[0].initial_state->vr == 1.025;
        }
        success *= model_data.bus[0].Vr0 == 0.9;
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
        success *= std::get<RealT>(model_data.branch[0].parameters[BranchParameters::R]) == 0.1;

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
        data.bus[0].bus_id   = 0;
        data.bus[0].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = 10.0;
        data.bus[0].Vi0      = 20.0;

        // Bus 1
        data.bus[1].bus_id   = 1;
        data.bus[1].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].Vr0      = 30.0;
        data.bus[1].Vi0      = 40.0;

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
