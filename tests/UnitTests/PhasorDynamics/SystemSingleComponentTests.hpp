#include <iomanip>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

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

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus);
        Log::setVerbosity(previous_verbosity);
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

      /// ESDC1A through the production path: model data to system
      /// construction to required bus and field-voltage signal wiring. UEL
      /// mode 2 makes zero field voltage an admissible standalone state.
      TestOutcome esdc1a()
      {
        using Data    = PhasorDynamics::Exciter::Esdc1aData<RealT, IdxT>;
        using Buses   = typename Data::Buses;
        using Outputs = typename Data::SignalOutputs;
        using Params  = typename Data::Parameters;
        using Vars    = PhasorDynamics::Exciter::Esdc1aInternalVariables;

        constexpr IdxT bus_id = static_cast<IdxT>(1);
        constexpr IdxT efd_id = static_cast<IdxT>(1);

        TestStatus success = true;

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = bus_id;
        data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = static_cast<RealT>(1.0);
        data.bus[0].Vi0      = static_cast<RealT>(0.0);

        data.signal.resize(1);
        data.signal[0].signal_id = efd_id;
        data.signal[0].name      = "Field Voltage";

        Data esdc1a_data;
        esdc1a_data.device_class                 = "Esdc1a";
        esdc1a_data.disambiguation_string        = "esdc1a_system";
        esdc1a_data.buses[Buses::bus]            = bus_id;
        esdc1a_data.parameters[Params::Tr]       = static_cast<RealT>(0.02);
        esdc1a_data.parameters[Params::Tb]       = static_cast<RealT>(0.5);
        esdc1a_data.parameters[Params::UEL]      = static_cast<IdxT>(2);
        esdc1a_data.signal_outputs[Outputs::efd] = efd_id;
        data.esdc1a.push_back(esdc1a_data);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        auto* efd  = system.getSignal(efd_id);
        success   *= efd->linked();
        success   *= efd->getVariableIndex() == static_cast<IdxT>(Vars::EFD);

        auto missing_bus_data          = data;
        missing_bus_data.bus[0].bus_id = static_cast<IdxT>(0);
        missing_bus_data.esdc1a[0].buses.clear();

        PhasorDynamics::SystemModel<ScalarT, IdxT> missing_bus_system(missing_bus_data);
        const auto                                 previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= missing_bus_system.verify() > 0;
        Log::setVerbosity(previous_verbosity);

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

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = static_cast<IdxT>(1);
        data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = static_cast<RealT>(1.0);
        data.bus[0].Vi0      = static_cast<RealT>(0.0);
        data.regca.push_back(makeRegcaData());

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size()
                   == static_cast<IdxT>(PhasorDynamics::Converter::RegcaInternalVariables::MAXIMUM);

        return success.report(__func__);
      }

      TestOutcome repca()
      {
        using Buses  = PhasorDynamics::Controller::RepcaBuses;
        using Inputs = PhasorDynamics::Controller::RepcaSignalInputs;
        using Params = PhasorDynamics::Controller::RepcaParameters;
        using Vars   = PhasorDynamics::Controller::RepcaInternalVariables;

        constexpr IdxT bus_id   = static_cast<IdxT>(1);
        constexpr IdxT input_id = static_cast<IdxT>(1);

        TestStatus success = true;

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = bus_id;
        data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = static_cast<RealT>(1.0);
        data.bus[0].Vi0      = static_cast<RealT>(0.0);

        data.signal.resize(1);
        data.signal[0].signal_id = input_id;

        typename PhasorDynamics::SystemModelData<RealT, IdxT>::RepcaDataT repca_data;
        repca_data.buses[Buses::bus]         = bus_id;
        repca_data.signal_inputs[Inputs::ir] = input_id;
        repca_data.signal_inputs[Inputs::ii] = input_id;
        repca_data.signal_inputs[Inputs::p]  = input_id;
        repca_data.signal_inputs[Inputs::q]  = input_id;
        repca_data.parameters[Params::Tp]    = static_cast<RealT>(0.02);
        data.repca.push_back(repca_data);

        ScalarT input_value{};
        IdxT    input_index = INVALID_INDEX<IdxT>;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        system.getSignal(input_id)->set(&input_value, &input_index);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        return success.report(__func__);
      }

      /// Construct GASTPTI through the production system-data path.
      TestOutcome gastpti()
      {
        using Outputs = PhasorDynamics::Governor::GastPtiSignalOutputs;
        using Vars    = PhasorDynamics::Governor::GastPtiInternalVariables;

        TestStatus success = true;

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.signal.resize(1);
        data.signal[0].signal_id = static_cast<IdxT>(1);
        data.signal[0].name      = "Mechanical Power";

        auto& gastpti                          = data.gastpti.emplace_back();
        gastpti.device_class                   = "GastPti";
        gastpti.disambiguation_string          = "gastpti_test";
        gastpti.signal_outputs[Outputs::pmech] = static_cast<IdxT>(1);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        return success.report(__func__);
      }

      /// REECB through the production data path.
      TestOutcome reecb()
      {
        using Data   = PhasorDynamics::Controller::ReecbData<RealT, IdxT>;
        using Buses  = typename Data::Buses;
        using Params = typename Data::Parameters;
        using Vars   = PhasorDynamics::Controller::ReecbInternalVariables;

        constexpr IdxT bus_id = static_cast<IdxT>(1);

        TestStatus success = true;

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = bus_id;
        data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = static_cast<RealT>(1.0);
        data.bus[0].Vi0      = static_cast<RealT>(0.0);

        Data reecb_data;
        reecb_data.buses[Buses::bus]        = bus_id;
        reecb_data.parameters[Params::Tp]   = static_cast<RealT>(0.02);
        reecb_data.parameters[Params::Pmin] = static_cast<RealT>(-1.0);
        data.reecb.push_back(reecb_data);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        return success.report(__func__);
      }

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

        PhasorDynamics::SignalNode<ScalarT, IdxT>      pmech;
        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> tgov1;
        tgov1.getSignals()
            .template assignSignalNode<PhasorDynamics::Governor::Tgov1InternalVariables::PM>(&pmech);
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

      /// HYGOV through the production path: model data to system
      /// construction to signal wiring. The declared pmech signal starts at
      /// zero mechanical power, an admissible operating point.
      TestOutcome hygov()
      {
        using Data    = PhasorDynamics::Governor::HygovData<RealT, IdxT>;
        using Outputs = typename Data::SignalOutputs;
        using Params  = typename Data::Parameters;
        using Vars    = PhasorDynamics::Governor::HygovInternalVariables;

        constexpr IdxT pmech_id = static_cast<IdxT>(2);

        TestStatus success = true;

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.freq_base = 60.0;
        data.va_base   = 100.0e6;
        data.signal.resize(1);
        data.signal[0].signal_id = pmech_id;
        data.signal[0].name      = "Mechanical Power";

        Data hygov_data;
        hygov_data.device_class                   = "Hygov";
        hygov_data.disambiguation_string          = "hygov_system";
        hygov_data.parameters[Params::Trate]      = static_cast<RealT>(100.0);
        hygov_data.parameters[Params::Tnp]        = static_cast<RealT>(1.0);
        hygov_data.signal_outputs[Outputs::pmech] = pmech_id;
        data.hygov.push_back(hygov_data);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        auto* pmech  = system.getSignal(pmech_id);
        success     *= pmech->linked();
        success     *= pmech->getVariableIndex() == static_cast<IdxT>(Vars::PMECH);

        auto missing_output_data = data;
        missing_output_data.hygov[0].signal_outputs.clear();

        PhasorDynamics::SystemModel<ScalarT, IdxT> missing_output_system(missing_output_data);
        const auto                                 previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= missing_output_system.verify() > 0;
        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

    private:
      auto makeRegcaData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Converter::RegcaParameters;
        using Buses  = PhasorDynamics::Converter::RegcaBuses;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class               = "Regca";
        data.disambiguation_string      = "regca_test";
        data.buses[Buses::bus]          = static_cast<IdxT>(1);
        data.parameters[Params::p0]     = static_cast<RealT>(1.0);
        data.parameters[Params::q0]     = static_cast<RealT>(0.0);
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
    };

  } // namespace Testing
} // namespace GridKit
