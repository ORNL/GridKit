#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Component smoke tests and focused integration tests through SystemModel.
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
        std::cout << "Testing expected ESDC1A missing-bus configuration error.\n";
        success *= missing_bus_system.verify() > 0;

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

      TestOutcome reecb()
      {
        using ConstantParams  = PhasorDynamics::ConstantSignalSourceParameters;
        using ConstantOutputs = PhasorDynamics::ConstantSignalSourceSignalOutputs;
        using Vars            = PhasorDynamics::Converter::ReecbInternalVariables;
        using Ext             = PhasorDynamics::Converter::ReecbExternalVariables;
        using ReecbT          = PhasorDynamics::Converter::Reecb<ScalarT, IdxT>;

        TestStatus success = true;

        // A missing required bus must remain missing so REECB verification can
        // reject it; SystemModel must not silently substitute bus zero.
        {
          std::cout << "Testing missing REECB bus mapping; logged errors are expected.\n";
          auto missing_bus_model = makeReecbData();
          missing_bus_model.buses.clear();
          missing_bus_model.signal_inputs.clear();
          missing_bus_model.signal_outputs.clear();

          PhasorDynamics::SystemModelData<RealT, IdxT> missing_bus_data;
          missing_bus_data.bus.resize(1);
          missing_bus_data.bus[0].bus_id = static_cast<IdxT>(0);
          missing_bus_data.bus[0].bus_type =
              PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
          missing_bus_data.reecb.push_back(missing_bus_model);

          bool rejected = false;
          try
          {
            PhasorDynamics::SystemModel<ScalarT, IdxT> missing_bus_system(missing_bus_data);
            missing_bus_system.allocate();
          }
          catch (const std::runtime_error&)
          {
            rejected = true;
          }
          success *= rejected;
        }

        PhasorDynamics::SystemModelData<RealT, IdxT> data;
        data.va_base = static_cast<RealT>(100.0e6);
        data.bus.resize(1);
        data.bus[0].bus_id   = static_cast<IdxT>(1);
        data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = static_cast<RealT>(1.0);
        data.bus[0].Vi0      = static_cast<RealT>(0.0);

        data.signal.resize(static_cast<size_t>(ipcmd_signal_id));
        for (size_t signal = 0; signal < data.signal.size(); ++signal)
        {
          data.signal[signal].signal_id = pe_signal_id + static_cast<IdxT>(signal);
        }

        data.reecb.push_back(makeReecbData());

        typename PhasorDynamics::SystemModelData<RealT, IdxT>::ConstantSourceT source;
        source.parameters[ConstantParams::Sr]      = static_cast<RealT>(0.25);
        source.parameters[ConstantParams::Si]      = static_cast<RealT>(0.05);
        source.signal_outputs[ConstantOutputs::sr] = pe_signal_id;
        source.signal_outputs[ConstantOutputs::si] = qgen_signal_id;
        data.constant_source.push_back(source);

        source.signal_outputs[ConstantOutputs::sr] = qext_signal_id;
        source.signal_outputs[ConstantOutputs::si] = pfaref_signal_id;
        data.constant_source.push_back(source);

        source.signal_outputs.erase(ConstantOutputs::si);
        source.signal_outputs[ConstantOutputs::sr] = pref_signal_id;
        data.constant_source.push_back(source);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        success *= system.allocate() == 0;
        system.getSignal(iqcmd_signal_id)->init(static_cast<ScalarT>(0.05));
        system.getSignal(ipcmd_signal_id)->init(static_cast<ScalarT>(0.25));
        success *= system.verify() == 0;
        for (IdxT signal_id = pe_signal_id; signal_id <= ipcmd_signal_id; ++signal_id)
        {
          success *= system.getSignal(signal_id)->linked();
        }
        success *= system.initialize() == 0;
        success *= system.tagDifferentiable() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;
        success *= system.size() == static_cast<IdxT>(Vars::MAXIMUM);

        auto* reecb  = dynamic_cast<ReecbT*>(system.getComponent(static_cast<IdxT>(0)));
        success     *= reecb != nullptr;
        if (reecb != nullptr)
        {
          auto& signals  = reecb->getSignals();
          success       *= signals.template readExternalVariableIndex<Ext::PE>()
                     == system.getSignal(pe_signal_id)->getVariableIndex();
          success *= signals.template readExternalVariableIndex<Ext::QGEN>()
                     == system.getSignal(qgen_signal_id)->getVariableIndex();
          success *= signals.template readExternalVariableIndex<Ext::QEXT>()
                     == system.getSignal(qext_signal_id)->getVariableIndex();
          success *= signals.template readExternalVariableIndex<Ext::PFAREF>()
                     == system.getSignal(pfaref_signal_id)->getVariableIndex();
          success *= signals.template readExternalVariableIndex<Ext::PREF>()
                     == system.getSignal(pref_signal_id)->getVariableIndex();
          success *= signals.template getSignalNode<Vars::IQCMD>()
                     == system.getSignal(iqcmd_signal_id);
          success *= signals.template getSignalNode<Vars::IPCMD>()
                     == system.getSignal(ipcmd_signal_id);

          success *= residualsAreZero(*reecb, "REECB SystemModel");

          // The constant sources drive qext and pfaref with 0.25/0.05, but
          // REECB overwrites attached unknown references during initialize(),
          // so the published values below are its resolved setpoints.
          const std::array<RealT, static_cast<size_t>(ipcmd_signal_id)> expected_signals{
              0.25, // pe from the constant source
              0.05, // qgen from the constant source
              0.05, // qext resolved by REECB
              0.0,  // pfaref resolved by REECB
              0.25, // pref resolved by REECB
              0.05, // iqcmd owned by REECB
              0.25, // ipcmd owned by REECB
          };
          for (size_t signal = 0; signal < expected_signals.size(); ++signal)
          {
            success *= isEqual(
                static_cast<RealT>(system.getSignal(pe_signal_id + static_cast<IdxT>(signal))->read()),
                expected_signals[signal],
                static_cast<RealT>(1.0e-9));
          }

          // The component/system base ratio is two. Perturbing only the
          // system-base command therefore changes its component-base residual
          // by twice the perturbation.
          const auto* residual = reecb->getResidual().getData();
          system.getSignal(iqcmd_signal_id)->init(static_cast<ScalarT>(0.06));
          success *= system.evaluateResidual() == 0;
          success *= isEqual(static_cast<RealT>(residual[static_cast<size_t>(Vars::IQCMD)]),
                             static_cast<RealT>(-0.02),
                             static_cast<RealT>(1.0e-9));
        }
        return success.report(__func__);
      }

      /// REGCA initializes the shared current commands and actual power
      /// feedback before REECB consumes them, with mixed component/system
      /// bases at off-nominal terminal voltages. Signal wiring identity for
      /// the pair is covered by ReecbIntegrationTests.
      TestOutcome regcaReecb()
      {
        using ConstantParams  = PhasorDynamics::ConstantSignalSourceParameters;
        using ConstantOutputs = PhasorDynamics::ConstantSignalSourceSignalOutputs;
        using RegcaParams     = PhasorDynamics::Converter::RegcaParameters;
        using RegcaInputs     = PhasorDynamics::Converter::RegcaSignalInputs;
        using RegcaOutputs    = PhasorDynamics::Converter::RegcaSignalOutputs;
        using RegcaVars       = PhasorDynamics::Converter::RegcaInternalVariables;
        using ReecbParams     = PhasorDynamics::Converter::ReecbParameters;
        using ReecbVars       = PhasorDynamics::Converter::ReecbInternalVariables;
        using RegcaT          = PhasorDynamics::Converter::Regca<ScalarT, IdxT>;
        using ReecbT          = PhasorDynamics::Converter::Reecb<ScalarT, IdxT>;

        TestStatus success = true;

        for (const RealT terminal_voltage :
             {static_cast<RealT>(0.9), static_cast<RealT>(1.19)})
        {
          PhasorDynamics::SystemModelData<RealT, IdxT> data;
          data.va_base = static_cast<RealT>(100.0e6);
          data.bus.resize(1);
          data.bus[0].bus_id   = static_cast<IdxT>(1);
          data.bus[0].bus_type = PhasorDynamics::BusData<RealT, IdxT>::BusType::SLACK;
          data.bus[0].Vr0      = terminal_voltage;
          data.bus[0].Vi0      = static_cast<RealT>(0.0);
          data.signal.resize(static_cast<size_t>(ipcmd_signal_id));
          for (size_t signal = 0; signal < data.signal.size(); ++signal)
          {
            data.signal[signal].signal_id = pe_signal_id + static_cast<IdxT>(signal);
          }

          auto regca_data                                  = makeRegcaData();
          regca_data.parameters[RegcaParams::p0]           = static_cast<RealT>(0.25);
          regca_data.parameters[RegcaParams::q0]           = static_cast<RealT>(0.05);
          regca_data.parameters[RegcaParams::mva]          = static_cast<RealT>(50.0);
          regca_data.signal_inputs[RegcaInputs::ipcmd]     = ipcmd_signal_id;
          regca_data.signal_inputs[RegcaInputs::iqcmd]     = iqcmd_signal_id;
          regca_data.signal_outputs[RegcaOutputs::pbranch] = pe_signal_id;
          regca_data.signal_outputs[RegcaOutputs::qbranch] = qgen_signal_id;
          data.regca.push_back(regca_data);

          auto reecb_data                           = makeReecbData();
          reecb_data.parameters[ReecbParams::QFlag] = true;
          reecb_data.parameters[ReecbParams::VFlag] = true;
          reecb_data.parameters[ReecbParams::Vmin]  = static_cast<RealT>(0.5);
          reecb_data.parameters[ReecbParams::Vmax]  = static_cast<RealT>(1.4);
          data.reecb.push_back(reecb_data);

          typename PhasorDynamics::SystemModelData<RealT, IdxT>::ConstantSourceT source;
          source.parameters[ConstantParams::Sr]      = static_cast<RealT>(0.0);
          source.parameters[ConstantParams::Si]      = static_cast<RealT>(0.0);
          source.signal_outputs[ConstantOutputs::sr] = qext_signal_id;
          source.signal_outputs[ConstantOutputs::si] = pfaref_signal_id;
          data.constant_source.push_back(source);
          source.signal_outputs.erase(ConstantOutputs::si);
          source.signal_outputs[ConstantOutputs::sr] = pref_signal_id;
          data.constant_source.push_back(source);

          PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
          success *= system.allocate() == 0;
          success *= system.verify() == 0;
          for (IdxT signal_id = pe_signal_id; signal_id <= ipcmd_signal_id; ++signal_id)
          {
            success *= system.getSignal(signal_id)->linked();
          }
          success *= system.initialize() == 0;
          success *= system.tagDifferentiable() == 0;
          success *= system.evaluateResidual() == 0;
          success *= system.evaluateJacobian() == 0;
          success *= system.size()
                     == static_cast<IdxT>(RegcaVars::MAXIMUM)
                            + static_cast<IdxT>(ReecbVars::MAXIMUM);

          auto* regca  = dynamic_cast<RegcaT*>(system.getComponent(static_cast<IdxT>(0)));
          auto* reecb  = dynamic_cast<ReecbT*>(system.getComponent(static_cast<IdxT>(1)));
          success     *= regca != nullptr;
          success     *= reecb != nullptr;
          if (regca == nullptr || reecb == nullptr)
          {
            continue;
          }

          const RealT pe     = static_cast<RealT>(system.getSignal(pe_signal_id)->read());
          const RealT qgen   = static_cast<RealT>(system.getSignal(qgen_signal_id)->read());
          const RealT ipcmd  = static_cast<RealT>(system.getSignal(ipcmd_signal_id)->read());
          const RealT iqcmd  = static_cast<RealT>(system.getSignal(iqcmd_signal_id)->read());
          success           *= isEqual(pe, static_cast<RealT>(0.25), static_cast<RealT>(1.0e-12));
          success           *= isEqual(qgen, static_cast<RealT>(0.05), static_cast<RealT>(1.0e-12));
          success           *= isEqual(
              static_cast<RealT>(reecb->y().getData()[static_cast<size_t>(ReecbVars::PMEAS)]),
              static_cast<RealT>(0.5),
              static_cast<RealT>(1.0e-12));

          if (terminal_voltage < static_cast<RealT>(1.0))
          {
            success *= std::abs(pe - ipcmd * terminal_voltage)
                       > static_cast<RealT>(1.0e-6);
          }
          else
          {
            success *= std::abs(qgen - iqcmd * terminal_voltage)
                       > static_cast<RealT>(1.0e-6);
          }

          success *= residualsAreZero(*regca, "REGCA integrated");
          success *= residualsAreZero(*reecb, "REECB integrated");
        }

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
        std::cout << "Testing expected HYGOV missing-output configuration error.\n";
        success *= missing_output_system.verify() > 0;

        return success.report(__func__);
      }

    private:
      template <typename ComponentT>
      bool residualsAreZero(ComponentT& component, const char* name) const
      {
        const auto* residual = component.getResidual().getData();
        bool        match    = true;
        for (size_t row = 0; row < static_cast<size_t>(component.size()); ++row)
        {
          if (!isEqual(static_cast<RealT>(residual[row]),
                       static_cast<RealT>(0.0),
                       static_cast<RealT>(1.0e-9)))
          {
            std::cout << name << " residual row " << row << " is not at steady state: "
                      << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                      << residual[row] << '\n';
            match = false;
          }
        }
        return match;
      }

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

      static constexpr IdxT pe_signal_id     = 1;
      static constexpr IdxT qgen_signal_id   = 2;
      static constexpr IdxT qext_signal_id   = 3;
      static constexpr IdxT pfaref_signal_id = 4;
      static constexpr IdxT pref_signal_id   = 5;
      static constexpr IdxT iqcmd_signal_id  = 6;
      static constexpr IdxT ipcmd_signal_id  = 7;

      auto makeReecbData() -> PhasorDynamics::Converter::ReecbData<RealT, IdxT>
      {
        using Params  = PhasorDynamics::Converter::ReecbParameters;
        using Buses   = PhasorDynamics::Converter::ReecbBuses;
        using Inputs  = PhasorDynamics::Converter::ReecbSignalInputs;
        using Outputs = PhasorDynamics::Converter::ReecbSignalOutputs;

        PhasorDynamics::Converter::ReecbData<RealT, IdxT> data;
        data.device_class                   = "Reecb";
        data.disambiguation_string          = "reecb_test";
        data.buses[Buses::bus]              = static_cast<IdxT>(1);
        data.parameters[Params::mva]        = static_cast<RealT>(50.0);
        data.signal_inputs[Inputs::pe]      = pe_signal_id;
        data.signal_inputs[Inputs::qgen]    = qgen_signal_id;
        data.signal_inputs[Inputs::qext]    = qext_signal_id;
        data.signal_inputs[Inputs::pfaref]  = pfaref_signal_id;
        data.signal_inputs[Inputs::pref]    = pref_signal_id;
        data.signal_outputs[Outputs::iqcmd] = iqcmd_signal_id;
        data.signal_outputs[Outputs::ipcmd] = ipcmd_signal_id;
        return data;
      }
    };

  } // namespace Testing
} // namespace GridKit
