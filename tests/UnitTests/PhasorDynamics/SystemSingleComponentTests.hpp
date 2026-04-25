#include <iomanip>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Smoke test for components (single component connected to an infinite bus)
    /// through the system model with the minimal constructors
    template <class ScalarT, typename IdxT>
    class SystemSingleComponentTests
    {
    private:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      SystemSingleComponentTests()  = default;
      ~SystemSingleComponentTests() = default;

      TestOutcome branch()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1, bus2;
        system->addBus(&bus1);
        system->addBus(&bus2);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2);
        system->addComponent(&branch);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == 0;
        success *= system->size() == branch.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome branchTrip()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1, bus2;
        system->addBus(&bus1);
        system->addBus(&bus2);

        // Branch data with realistic R/X/G/B so setDerivedParams stays finite.
        PhasorDynamics::BranchData<RealT, IdxT> branch_data;
        branch_data.disambiguation_string                           = "br_0";
        branch_data.parameters[PhasorDynamics::BranchParameters::R] = 0.05;
        branch_data.parameters[PhasorDynamics::BranchParameters::X] = 0.21;
        branch_data.parameters[PhasorDynamics::BranchParameters::G] = 0.0;
        branch_data.parameters[PhasorDynamics::BranchParameters::B] = 0.1;

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, branch_data);
        system->addComponent(&branch, branch_data.disambiguation_string);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;

        // Default: branch is closed (in service)
        success *= branch.status() == true;

        // Open via cue
        try
        {
          system->cue("br_0", PhasorDynamics::Action::Off);
        }
        catch (...)
        {
          success *= false;
        }
        success *= branch.status() == false;

        // Close via cue
        try
        {
          system->cue("br_0", PhasorDynamics::Action::On);
        }
        catch (...)
        {
          success *= false;
        }
        success *= branch.status() == true;

        // state0 = false -> branch starts open without a cue.
        PhasorDynamics::BranchData<RealT, IdxT> open_data;
        open_data.parameters[PhasorDynamics::BranchParameters::R]      = 0.05;
        open_data.parameters[PhasorDynamics::BranchParameters::X]      = 0.21;
        open_data.parameters[PhasorDynamics::BranchParameters::state0] = false;
        PhasorDynamics::Branch<ScalarT, IdxT> open_branch(&bus1, &bus2, open_data);
        success *= open_branch.status() == false;

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome bus()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::Bus<ScalarT, IdxT> bus;
        system->addBus(&bus);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == bus.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome busFault()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::BusFaultData<RealT, IdxT> fault_data;
        fault_data.disambiguation_string = "fault_0";
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, fault_data);
        system->addComponent(&fault, fault_data.disambiguation_string);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == 0;
        success *= system->size() == fault.size();

        // Cue routes to the matching child via SystemModel::cue.
        try
        {
          system->cue("fault_0", PhasorDynamics::Action::On);
        }
        catch (...)
        {
          success *= false;
        }

        // Unknown target throws.
        bool threw = false;
        try
        {
          system->cue("no_such", PhasorDynamics::Action::On);
        }
        catch (const std::exception&)
        {
          threw = true;
        }
        success *= threw;

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome ieeet1()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus);
        system->addComponent(&exciter);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == exciter.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome load()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::Load<ScalarT, IdxT> load(&bus);
        system->addComponent(&load);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == load.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome genrou()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, 0); ///< is unit_id really necessary?
        system->addComponent(&gen);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == gen.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome genClassical()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 0); ///< is unit_id really necessary?
        system->addComponent(&gen);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == gen.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome tgov1()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> tgov1;
        system->addComponent(&tgov1);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == tgov1.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
