#include <iomanip>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        bus1.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&bus1_vr);
        bus1.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&bus1_vi);
        bus2.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&bus2_vr);
        bus2.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&bus2_vi);

        PhasorDynamics::Branch<ScalarT, IdxT> branch;
        branch.getSignals().template attachSignalNode<PhasorDynamics::BranchExternalVariables::VR1>(&bus1_vr);
        branch.getSignals().template attachSignalNode<PhasorDynamics::BranchExternalVariables::VI1>(&bus1_vi);
        branch.getSignals().template attachSignalNode<PhasorDynamics::BranchExternalVariables::VR2>(&bus2_vr);
        branch.getSignals().template attachSignalNode<PhasorDynamics::BranchExternalVariables::VI2>(&bus2_vi);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::BusFault<ScalarT, IdxT> fault;
        fault.getSignals().template attachSignalNode<PhasorDynamics::BusFaultExternalVariables::VR>(&vr_signal);
        fault.getSignals().template attachSignalNode<PhasorDynamics::BusFaultExternalVariables::VI>(&vi_signal);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter;
        exciter.getSignals().template attachSignalNode<PhasorDynamics::Exciter::Ieeet1ExternalVariables::VREAL>(&vr_signal);
        exciter.getSignals().template attachSignalNode<PhasorDynamics::Exciter::Ieeet1ExternalVariables::VIMAG>(&vi_signal);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::LoadZ<ScalarT, IdxT> load;
        load.getSignals().template attachSignalNode<PhasorDynamics::LoadZExternalVariables::VR>(&vr_signal);
        load.getSignals().template attachSignalNode<PhasorDynamics::LoadZExternalVariables::VI>(&vi_signal);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::LoadZIP<ScalarT, IdxT> load;
        load.getSignals().template attachSignalNode<PhasorDynamics::LoadZIPExternalVariables::VR>(&vr_signal);
        load.getSignals().template attachSignalNode<PhasorDynamics::LoadZIPExternalVariables::VI>(&vi_signal);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::Genrou<ScalarT, IdxT> gen;
        gen.getSignals().template attachSignalNode<PhasorDynamics::GenrouExternalVariables::VR>(&vr_signal);
        gen.getSignals().template attachSignalNode<PhasorDynamics::GenrouExternalVariables::VI>(&vi_signal);
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

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_signal;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_signal;
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VR>(&vr_signal);
        bus.getSignals().template assignSignalNode<PhasorDynamics::BusInternalVariables::VI>(&vi_signal);

        PhasorDynamics::GenClassical<ScalarT, IdxT> gen;
        gen.getSignals().template attachSignalNode<PhasorDynamics::GenClassicalExternalVariables::VR>(&vr_signal);
        gen.getSignals().template attachSignalNode<PhasorDynamics::GenClassicalExternalVariables::VI>(&vi_signal);
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
