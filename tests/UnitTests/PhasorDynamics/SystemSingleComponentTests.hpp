#include <iomanip>
#include <iostream>

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

        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus);
        system->addFault(&fault);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == fault.size();

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

        PhasorDynamics::LoadZ<ScalarT, IdxT> load(&bus);
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

      TestOutcome loadZIP()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(1.0, 0.0);
        system->addBus(&bus);

        PhasorDynamics::LoadZIP<ScalarT, IdxT> load(&bus);
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

      TestOutcome regca()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(1.0, 0.0);
        system->addBus(&bus);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeRegcaData());
        system->addComponent(&regca);

        success *= system->allocate() == 0;
        success *= system->initialize() == 0;
        success *= system->evaluateResidual() == 0;
        success *= system->evaluateJacobian() == 0;
        success *= system->size() == regca.size();

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

    private:
      auto makeRegcaData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Converter::RegcaParameters;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class               = "Regca";
        data.disambiguation_string      = "regca_test";
        data.parameters[Params::P0]     = static_cast<RealT>(1.0);
        data.parameters[Params::Q0]     = static_cast<RealT>(0.0);
        data.parameters[Params::mva]    = static_cast<RealT>(100.0);
        data.parameters[Params::Tg]     = static_cast<RealT>(0.02);
        data.parameters[Params::TM]     = static_cast<RealT>(0.02);
        data.parameters[Params::Rqmax]  = static_cast<RealT>(999.0);
        data.parameters[Params::Rqmin]  = static_cast<RealT>(-999.0);
        data.parameters[Params::Rpmax]  = static_cast<RealT>(999.0);
        data.parameters[Params::sL]     = true;
        data.parameters[Params::IL1]    = static_cast<RealT>(1.1);
        data.parameters[Params::VL0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VL1]    = static_cast<RealT>(0.9);
        data.parameters[Params::VA0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VA1]    = static_cast<RealT>(0.9);
        data.parameters[Params::Vhvmax] = static_cast<RealT>(1.2);
        return data;
      }

    public:
      TestOutcome genrou()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus;
        system->addBus(&bus);

        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus);
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

        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus);
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
