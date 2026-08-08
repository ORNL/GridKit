#pragma once

#include <limits>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Connection tests for pairs of components that share a signal node.
    /// Each case checks that the node links both components and that they
    /// agree on its value after initialization. Solver-driven cases live
    /// in @ref PDIntegrationTests.
    template <typename scalar_type, typename index_type>
    class ComponentConnectionTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ComponentConnectionTests()  = default;
      ~ComponentConnectionTests() = default;

      // The tolerance only absorbs floating-point roundoff.
      static constexpr RealT kTol = std::numeric_limits<RealT>::epsilon();

      /// GENROU initializes first and writes the field voltage it needs to
      /// the shared node. ESDC1A then initializes around that value and
      /// must leave it unchanged at a steady state.
      TestOutcome genrouEsdc1a()
      {
        using MachineExternal = PhasorDynamics::GenrouExternalVariables;
        using ExciterInternal = PhasorDynamics::Exciter::Esdc1aInternalVariables;
        using ExciterParams   = PhasorDynamics::Exciter::Esdc1aParameters;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> efd;
        PhasorDynamics::Genrou<ScalarT, IdxT>     machine(&bus);

        PhasorDynamics::Exciter::Esdc1aData<RealT, IdxT> exciter_data;
        exciter_data.parameters[ExciterParams::Tr] = static_cast<RealT>(0.02);
        exciter_data.parameters[ExciterParams::Tb] = static_cast<RealT>(0.5);

        PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT> exciter(&bus, exciter_data);

        machine.getSignals().template attachSignalNode<MachineExternal::EFD>(&efd);
        exciter.getSignals().template assignSignalNode<ExciterInternal::EFD>(&efd);

        system.addBus(&bus);
        system.addComponent(&machine);
        system.addComponent(&exciter);

        success *= system.allocate() == 0;
        success *= efd.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        // At zero power the required field voltage equals the terminal voltage.
        success *= isEqual(efd.read(), static_cast<ScalarT>(1.0), kTol);

        const auto* residual = exciter.getResidual().getData();
        for (IdxT row = 0; row < exciter.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      /// GENROU initializes first and writes the mechanical power it needs to
      /// the shared node. HYGOV then initializes around that value and must
      /// leave it unchanged at a steady state. The speed port is left
      /// unattached, so HYGOV reads the exactly zero deviation its
      /// initialization requires.
      TestOutcome genrouHygov()
      {
        using MachineExternal  = PhasorDynamics::GenrouExternalVariables;
        using GovernorInternal = PhasorDynamics::Governor::HygovInternalVariables;
        using GovernorParams   = PhasorDynamics::Governor::HygovParameters;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech;
        PhasorDynamics::Genrou<ScalarT, IdxT>     machine(&bus);

        PhasorDynamics::Governor::HygovData<RealT, IdxT> governor_data;
        governor_data.parameters[GovernorParams::Tnp] = static_cast<RealT>(1.0);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> governor(governor_data);

        machine.getSignals().template attachSignalNode<MachineExternal::PM>(&pmech);
        governor.getSignals().template assignSignalNode<GovernorInternal::PMECH>(&pmech);

        system.addBus(&bus);
        system.addComponent(&machine);
        system.addComponent(&governor);

        success *= system.allocate() == 0;
        success *= pmech.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        // At zero machine power the required mechanical power is zero.
        success *= isEqual(pmech.read(), static_cast<ScalarT>(0.0), kTol);

        const auto* residual = governor.getResidual().getData();
        for (IdxT row = 0; row < governor.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      /// REGCA initializes first and publishes its branch current and power
      /// to the four shared nodes. REPCA then initializes around those
      /// measurements and must hold a steady state without a frequency input.
      TestOutcome regcaRepca()
      {
        using ConverterInternal = PhasorDynamics::Converter::RegcaInternalVariables;
        using ConverterParams   = PhasorDynamics::Converter::RegcaParameters;
        using PlantInternal     = PhasorDynamics::Controller::RepcaInternalVariables;
        using PlantExternal     = PhasorDynamics::Controller::RepcaExternalVariables;
        using PlantParams       = PhasorDynamics::Controller::RepcaParameters;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> ir;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> converter_data;
        converter_data.parameters[ConverterParams::p0]     = static_cast<RealT>(0.4);
        converter_data.parameters[ConverterParams::q0]     = static_cast<RealT>(0.1);
        converter_data.parameters[ConverterParams::mva]    = static_cast<RealT>(100.0);
        converter_data.parameters[ConverterParams::Tg]     = static_cast<RealT>(0.02);
        converter_data.parameters[ConverterParams::TM]     = static_cast<RealT>(0.02);
        converter_data.parameters[ConverterParams::Rqmax]  = static_cast<RealT>(999.0);
        converter_data.parameters[ConverterParams::Rqmin]  = static_cast<RealT>(-999.0);
        converter_data.parameters[ConverterParams::Rpmax]  = static_cast<RealT>(999.0);
        converter_data.parameters[ConverterParams::sL]     = true;
        converter_data.parameters[ConverterParams::IL1]    = static_cast<RealT>(1.1);
        converter_data.parameters[ConverterParams::VL0]    = static_cast<RealT>(0.4);
        converter_data.parameters[ConverterParams::VL1]    = static_cast<RealT>(0.9);
        converter_data.parameters[ConverterParams::VA0]    = static_cast<RealT>(0.4);
        converter_data.parameters[ConverterParams::VA1]    = static_cast<RealT>(0.9);
        converter_data.parameters[ConverterParams::Vhvmax] = static_cast<RealT>(1.2);

        PhasorDynamics::Controller::RepcaData<RealT, IdxT> plant_data;
        plant_data.parameters[PlantParams::mva] = static_cast<RealT>(50.0);
        plant_data.parameters[PlantParams::Tp]  = static_cast<RealT>(0.05);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT>  converter(&bus, converter_data);
        PhasorDynamics::Controller::Repca<ScalarT, IdxT> plant(&bus, plant_data);

        auto& converter_signals = converter.getSignals();
        converter_signals.template assignSignalNode<ConverterInternal::IR>(&ir);
        converter_signals.template assignSignalNode<ConverterInternal::II>(&ii);
        converter_signals.template assignSignalNode<ConverterInternal::PBR>(&p);
        converter_signals.template assignSignalNode<ConverterInternal::QBR>(&q);

        auto& plant_signals = plant.getSignals();
        plant_signals.template attachSignalNode<PlantExternal::IR>(&ir);
        plant_signals.template attachSignalNode<PlantExternal::II>(&ii);
        plant_signals.template attachSignalNode<PlantExternal::P>(&p);
        plant_signals.template attachSignalNode<PlantExternal::Q>(&q);

        system.addBus(&bus);
        system.addComponent(&converter);
        system.addComponent(&plant);

        success *= system.allocate() == 0;
        success *= ir.linked();
        success *= ii.linked();
        success *= p.linked();
        success *= q.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        success *= isEqual(ir.read(), static_cast<ScalarT>(0.4), kTol);
        success *= isEqual(ii.read(), static_cast<ScalarT>(-0.1), kTol);
        success *= isEqual(p.read(), static_cast<ScalarT>(0.4), kTol);
        success *= isEqual(q.read(), static_cast<ScalarT>(0.1), kTol);

        const auto* state = plant.y().getData();

        success *= isEqual(state[static_cast<size_t>(PlantInternal::PMEAS)],
                           static_cast<ScalarT>(0.8),
                           kTol);
        success *= isEqual(state[static_cast<size_t>(PlantInternal::QMEAS)],
                           static_cast<ScalarT>(0.2),
                           kTol);

        const auto* residual = plant.getResidual().getData();
        for (IdxT row = 0; row < plant.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      /// REGCA initializes first and publishes the current commands it
      /// resolves to the shared nodes, alongside its branch powers. REECB
      /// then initializes around all four published values and must leave
      /// them unchanged at a steady state.
      TestOutcome regcaReecb()
      {
        using ConverterExternal  = PhasorDynamics::Converter::RegcaExternalVariables;
        using ConverterInternal  = PhasorDynamics::Converter::RegcaInternalVariables;
        using ConverterParams    = PhasorDynamics::Converter::RegcaParameters;
        using ControllerExternal = PhasorDynamics::Controller::ReecbExternalVariables;
        using ControllerInternal = PhasorDynamics::Controller::ReecbInternalVariables;
        using ControllerParams   = PhasorDynamics::Controller::ReecbParameters;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pe;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen;

        // The operating point is exactly representable, so the pair rests at
        // an exact steady state.
        PhasorDynamics::Converter::RegcaData<RealT, IdxT> converter_data;
        converter_data.parameters[ConverterParams::p0]     = static_cast<RealT>(0.375);
        converter_data.parameters[ConverterParams::q0]     = static_cast<RealT>(0.0625);
        converter_data.parameters[ConverterParams::mva]    = static_cast<RealT>(100.0);
        converter_data.parameters[ConverterParams::Tg]     = static_cast<RealT>(0.02);
        converter_data.parameters[ConverterParams::TM]     = static_cast<RealT>(0.02);
        converter_data.parameters[ConverterParams::Rqmax]  = static_cast<RealT>(999.0);
        converter_data.parameters[ConverterParams::Rqmin]  = static_cast<RealT>(-999.0);
        converter_data.parameters[ConverterParams::Rpmax]  = static_cast<RealT>(999.0);
        converter_data.parameters[ConverterParams::sL]     = true;
        converter_data.parameters[ConverterParams::IL1]    = static_cast<RealT>(1.1);
        converter_data.parameters[ConverterParams::VL0]    = static_cast<RealT>(0.25);
        converter_data.parameters[ConverterParams::VL1]    = static_cast<RealT>(0.75);
        converter_data.parameters[ConverterParams::VA0]    = static_cast<RealT>(0.25);
        converter_data.parameters[ConverterParams::VA1]    = static_cast<RealT>(0.75);
        converter_data.parameters[ConverterParams::Vhvmax] = static_cast<RealT>(1.5);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> converter(&bus, converter_data);

        PhasorDynamics::Controller::ReecbData<RealT, IdxT> controller_data;
        controller_data.parameters[ControllerParams::mva]    = static_cast<RealT>(100.0);
        controller_data.parameters[ControllerParams::Tp]     = static_cast<RealT>(0.02);
        controller_data.parameters[ControllerParams::QFlag]  = true;
        controller_data.parameters[ControllerParams::VFlag]  = true;
        controller_data.parameters[ControllerParams::Pqflag] = true;
        controller_data.parameters[ControllerParams::Imax]   = static_cast<RealT>(0.625);
        controller_data.parameters[ControllerParams::Vmin]   = static_cast<RealT>(0.5);
        controller_data.parameters[ControllerParams::Vmax]   = static_cast<RealT>(1.5);

        PhasorDynamics::Controller::Reecb<ScalarT, IdxT> controller(&bus, controller_data);

        controller.getSignals().template assignSignalNode<ControllerInternal::IPCMD>(&ipcmd);
        controller.getSignals().template assignSignalNode<ControllerInternal::IQCMD>(&iqcmd);
        converter.getSignals().template attachSignalNode<ConverterExternal::IPCMD>(&ipcmd);
        converter.getSignals().template attachSignalNode<ConverterExternal::IQCMD>(&iqcmd);
        converter.getSignals().template assignSignalNode<ConverterInternal::PBR>(&pe);
        converter.getSignals().template assignSignalNode<ConverterInternal::QBR>(&qgen);
        controller.getSignals().template attachSignalNode<ControllerExternal::PE>(&pe);
        controller.getSignals().template attachSignalNode<ControllerExternal::QGEN>(&qgen);

        system.addBus(&bus);
        system.addComponent(&converter);
        system.addComponent(&controller);

        success *= system.allocate() == 0;
        success *= ipcmd.linked() && iqcmd.linked() && pe.linked() && qgen.linked();
        success *= ipcmd.getVariableIndex()
                   == controller.getVariableIndex(
                       static_cast<IdxT>(ControllerInternal::IPCMD));
        success *= iqcmd.getVariableIndex()
                   == controller.getVariableIndex(
                       static_cast<IdxT>(ControllerInternal::IQCMD));
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        // At unit terminal voltage the shared nodes carry the scheduled powers.
        success *= isEqual(ipcmd.read(), static_cast<ScalarT>(0.375), kTol);
        success *= isEqual(iqcmd.read(), static_cast<ScalarT>(0.0625), kTol);
        success *= isEqual(pe.read(), static_cast<ScalarT>(0.375), kTol);
        success *= isEqual(qgen.read(), static_cast<ScalarT>(0.0625), kTol);

        const auto* residual = controller.getResidual().getData();
        for (IdxT row = 0; row < controller.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      /// IEEEST reinitializes before IEEET1 so the exciter does not latch a
      /// stale stabilizer output.
      TestOutcome ieeestIeeet1()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Exciter;
        using namespace PhasorDynamics::Stabilizer;

        TestStatus success = true;

        constexpr IdxT               bus_id   = static_cast<IdxT>(1);
        constexpr IdxT               input_id = static_cast<IdxT>(20);
        constexpr IdxT               vs_id    = static_cast<IdxT>(21);
        SystemModelData<RealT, IdxT> data;

        auto& bus    = data.bus.emplace_back();
        bus.bus_id   = bus_id;
        bus.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
        bus.Vr0      = ONE<RealT>;
        bus.Vi0      = ZERO<RealT>;

        data.signal = {{"IEEEST input", input_id}, {"IEEEST output", vs_id}};

        auto& source                                                 = data.constant_source.emplace_back();
        source.signal_outputs[ConstantSignalSourceSignalOutputs::sr] = input_id;
        source.parameters[ConstantSignalSourceParameters::Sr]        = ZERO<RealT>;

        auto& stabilizer                                       = data.stabilizer.emplace_back();
        stabilizer.signal_inputs[IeeestSignalInputs::input]    = input_id;
        stabilizer.signal_outputs[IeeestSignalOutputs::output] = vs_id;

        auto& exciter                                 = data.exciter.emplace_back();
        exciter.buses[Ieeet1Buses::bus]               = bus_id;
        exciter.signal_inputs[Ieeet1SignalInputs::vs] = vs_id;
        exciter.parameters[Ieeet1Parameters::Tr]      = static_cast<RealT>(0.02);

        SystemModel<ScalarT, IdxT> system(data);
        success *= system.allocate() == 0;

        auto* vs = system.getSignal(vs_id);
        if (vs == nullptr || !vs->linked())
        {
          success = false;
          return success.report(__func__);
        }

        // Reinitialization must not let the exciter latch this stale value.
        vs->init(static_cast<ScalarT>(0.05));
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= isEqual(vs->read(), ZERO<ScalarT>, kTol);

        const auto& residual      = system.getResidual();
        const auto* residual_data = residual.getData();
        for (IdxT row = 0; row < residual.getSize(); ++row)
        {
          success *= isEqual(residual_data[row], ZERO<ScalarT>, kTol);
        }

        return success.report(__func__);
      }
    };

    /// Production SystemModel wiring tests for GASTPTI and synchronous machines.
    template <typename scalar_type, typename index_type>
    class GastPtiConnectionTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      GastPtiConnectionTests()  = default;
      ~GastPtiConnectionTests() = default;

      /// Check production GASTPTI signal wiring to GENROU.
      TestOutcome genrouGastPti()
      {
        using namespace PhasorDynamics;

        return checkConnection<Genrou<ScalarT, IdxT>,
                               GenrouInternalVariables::OMEGA,
                               GenrouExternalVariables::PM>(makeGenrouGastPtiCase(),
                                                            __func__);
      }

      /// Check production GASTPTI signal wiring to GENSAL.
      TestOutcome gensalGastPti()
      {
        using namespace PhasorDynamics;

        return checkConnection<Gensal<ScalarT, IdxT>,
                               GensalInternalVariables::OMEGA,
                               GensalExternalVariables::PM>(makeGensalGastPtiCase(),
                                                            __func__);
      }

    private:
      using SystemDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;
      using SystemT     = PhasorDynamics::SystemModel<ScalarT, IdxT>;
      using GastPtiT    = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;

      static constexpr IdxT kBusId               = static_cast<IdxT>(17);
      static constexpr IdxT kSpeedSignalId       = static_cast<IdxT>(101);
      static constexpr IdxT kPmechSignalId       = static_cast<IdxT>(102);
      static constexpr IdxT kPrefSignalId        = static_cast<IdxT>(103);
      static constexpr IdxT kMachineComponentId  = static_cast<IdxT>(0);
      static constexpr IdxT kGovernorComponentId = static_cast<IdxT>(1);

      static constexpr RealT kSystemBaseVa         = static_cast<RealT>(100.0e6);
      static constexpr RealT kMachineBaseMva       = static_cast<RealT>(100.0);
      static constexpr RealT kTurbineRatingMw      = static_cast<RealT>(100.0);
      static constexpr RealT kInitialActivePower   = static_cast<RealT>(0.4);
      static constexpr RealT kInitialReactivePower = static_cast<RealT>(0.05);

      template <typename MachineT, auto speed_variable, auto pmech_variable>
      TestOutcome checkConnection(const SystemDataT& data, const char* test_name)
      {
        using namespace PhasorDynamics::Governor;

        TestStatus success = true;
        SystemT    system(data);

        success *= system.allocate() == 0;

        auto* machine  = dynamic_cast<MachineT*>(system.getComponent(kMachineComponentId));
        auto* governor = dynamic_cast<GastPtiT*>(system.getComponent(kGovernorComponentId));
        auto* speed    = system.getSignal(kSpeedSignalId);
        auto* pmech    = system.getSignal(kPmechSignalId);
        auto* pref     = system.getSignal(kPrefSignalId);

        if (machine == nullptr || governor == nullptr || speed == nullptr
            || pmech == nullptr || pref == nullptr)
        {
          success = false;
          return success.report(test_name);
        }

        const bool signals_linked  = speed->linked() && pmech->linked() && pref->linked();
        success                   *= signals_linked;
        if (!signals_linked)
        {
          return success.report(test_name);
        }

        auto& machine_signals  = machine->getSignals();
        auto& governor_signals = governor->getSignals();

        const bool ports_connected =
            machine_signals.template isAssigned<speed_variable>()
            && machine_signals.template isAttached<pmech_variable>()
            && governor_signals.template isAttached<GastPtiExternalVariables::OMEGA>()
            && governor_signals.template isAttached<GastPtiExternalVariables::PREF>()
            && governor_signals.template isAssigned<GastPtiInternalVariables::PMECH>();
        success *= ports_connected;
        if (!ports_connected)
        {
          return success.report(test_name);
        }

        const IdxT speed_index = speed->getVariableIndex();
        const IdxT pmech_index = pmech->getVariableIndex();
        const IdxT pref_index  = pref->getVariableIndex();

        success *= speed_index != pmech_index;
        success *= speed_index
                   == machine->getVariableIndex(static_cast<IdxT>(speed_variable));
        success *= speed_index
                   == governor_signals.template readExternalVariableIndex<
                       GastPtiExternalVariables::OMEGA>();
        success *= pmech_index
                   == machine_signals.template readExternalVariableIndex<pmech_variable>();
        success *= pmech_index
                   == governor->getVariableIndex(
                       static_cast<IdxT>(GastPtiInternalVariables::PMECH));
        success *= pref_index
                   == governor_signals.template readExternalVariableIndex<
                       GastPtiExternalVariables::PREF>();

        return success.report(test_name);
      }

      static SystemDataT makeGastPtiCase()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Governor;

        SystemDataT data;
        data.va_base = kSystemBaseVa;

        auto& bus    = data.bus.emplace_back();
        bus.bus_id   = kBusId;
        bus.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
        bus.Vr0      = ONE<RealT>;
        bus.Vi0      = ZERO<RealT>;

        data.signal = {{"Machine Speed Deviation", kSpeedSignalId},
                       {"Mechanical Power", kPmechSignalId},
                       {"Active Power Reference", kPrefSignalId}};

        auto& governor                                       = data.gastpti.emplace_back();
        governor.signal_inputs[GastPtiSignalInputs::speed]   = kSpeedSignalId;
        governor.signal_inputs[GastPtiSignalInputs::pref]    = kPrefSignalId;
        governor.signal_outputs[GastPtiSignalOutputs::pmech] = kPmechSignalId;
        governor.parameters[GastPtiParameters::Trate]        = kTurbineRatingMw;

        auto& source                                                 = data.constant_source.emplace_back();
        source.signal_outputs[ConstantSignalSourceSignalOutputs::sr] = kPrefSignalId;
        source.parameters[ConstantSignalSourceParameters::Sr]        = ZERO<RealT>;

        return data;
      }

      static SystemDataT makeGenrouGastPtiCase()
      {
        using namespace PhasorDynamics;

        auto data = makeGastPtiCase();

        auto& machine                                      = data.genrou.emplace_back();
        machine.buses[GenrouBuses::bus]                    = kBusId;
        machine.signal_outputs[GenrouSignalOutputs::speed] = kSpeedSignalId;
        machine.signal_inputs[GenrouSignalInputs::pmech]   = kPmechSignalId;
        machine.parameters[GenrouParameters::p0]           = kInitialActivePower;
        machine.parameters[GenrouParameters::q0]           = kInitialReactivePower;
        machine.parameters[GenrouParameters::H]            = RealT{3.0};
        machine.parameters[GenrouParameters::D]            = ZERO<RealT>;
        machine.parameters[GenrouParameters::Ra]           = ZERO<RealT>;
        machine.parameters[GenrouParameters::Tdop]         = RealT{7.0};
        machine.parameters[GenrouParameters::Tdopp]        = RealT{0.04};
        machine.parameters[GenrouParameters::Tqop]         = RealT{0.75};
        machine.parameters[GenrouParameters::Tqopp]        = RealT{0.05};
        machine.parameters[GenrouParameters::Xd]           = RealT{2.1};
        machine.parameters[GenrouParameters::Xdp]          = RealT{0.2};
        machine.parameters[GenrouParameters::Xdpp]         = RealT{0.18};
        machine.parameters[GenrouParameters::Xq]           = RealT{0.5};
        machine.parameters[GenrouParameters::Xqp]          = RealT{0.5};
        machine.parameters[GenrouParameters::Xqpp]         = RealT{0.18};
        machine.parameters[GenrouParameters::Xl]           = RealT{0.15};
        machine.parameters[GenrouParameters::S10]          = ZERO<RealT>;
        machine.parameters[GenrouParameters::S12]          = ZERO<RealT>;
        machine.parameters[GenrouParameters::mva]          = kMachineBaseMva;

        return data;
      }

      static SystemDataT makeGensalGastPtiCase()
      {
        using namespace PhasorDynamics;

        auto data = makeGastPtiCase();

        auto& machine                                      = data.gensal.emplace_back();
        machine.buses[GensalBuses::bus]                    = kBusId;
        machine.signal_outputs[GensalSignalOutputs::speed] = kSpeedSignalId;
        machine.signal_inputs[GensalSignalInputs::pmech]   = kPmechSignalId;
        machine.parameters[GensalParameters::p0]           = kInitialActivePower;
        machine.parameters[GensalParameters::q0]           = kInitialReactivePower;
        machine.parameters[GensalParameters::mva]          = kMachineBaseMva;

        return data;
      }
    };

  } // namespace Testing
} // namespace GridKit
