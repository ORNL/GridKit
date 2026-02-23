#pragma once

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using GridKit::PhasorDynamics::BranchParameters;
    using GridKit::PhasorDynamics::BranchPorts;

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
        data.branch[0].ports[BranchPorts::bus1]        = data.bus[0].bus_id;
        data.branch[0].ports[BranchPorts::bus2]        = data.bus[1].bus_id;
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
        status *= throws<std::runtime_error>(
            [&]()
            { sys.allocate(); });

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
        data.branch[0].ports[BranchPorts::bus1]        = data.bus[0].bus_id;
        data.branch[0].ports[BranchPorts::bus2]        = data.bus[1].bus_id;
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
        for (size_t i = 0; i < system.size(); ++i)
        {
          system.y()[i].setVariableNumber(i);
        }

        // Evaluate and get the system residuals
        system.evaluateResidual();
        std::vector<DependencyTracking::Variable> residual = system.getResidual();

        // Print the dependencies
        for (size_t i = 0; i < residual.size(); ++i)
        {
          std::cout << i << "th residual: ";
          (residual[i]).print(std::cout);
          std::cout << "\n";
        }

        // Extract the dependencies
        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual.size());
        for (IdxT i = 0; i < residual.size(); ++i)
        {
          dependencies[i] = (residual[i]).getDependencies();
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
