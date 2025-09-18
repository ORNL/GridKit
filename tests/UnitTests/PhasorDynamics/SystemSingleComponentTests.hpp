#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/ComponentLibrary.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Test the components (connected to an infinite bus) through the system model with the default
    /// constructors
    template <class ScalarT, typename IdxT>
    class SystemSingleComponentTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

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

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

        success *= system->size() == 0;
        success *= system->size() == branch.size();
        
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

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

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

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

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

        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, 0); // is unit_id really necessary?
        system->addComponent(&gen);

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

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

        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 0); // is unit_id really necessary?
        system->addComponent(&gen);

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

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

        system->allocate();
        system->initialize();
        system->evaluateResidual();
        system->evaluateJacobian();

        success *= system->size() == tgov1.size();
        
        delete system;
        system = nullptr;

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
