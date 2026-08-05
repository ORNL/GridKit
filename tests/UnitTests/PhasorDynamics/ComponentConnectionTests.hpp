#pragma once

#include <limits>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
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

      /// GENROU initializes first and writes the mechanical power it needs to
      /// the shared node. GASTPTI then initializes around that value and
      /// must leave it unchanged at a steady state.
      TestOutcome genrouGastPti()
      {
        using MachineExternal  = PhasorDynamics::GenrouExternalVariables;
        using GovernorInternal = PhasorDynamics::Governor::GastPtiInternalVariables;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT>        pmech;
        PhasorDynamics::Genrou<ScalarT, IdxT>            machine(&bus);
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> governor;

        machine.getSignals().template attachSignalNode<MachineExternal::PM>(&pmech);
        governor.getSignals().template assignSignalNode<GovernorInternal::PMECH>(&pmech);

        system.addBus(&bus);
        system.addComponent(&machine);
        system.addComponent(&governor);

        success *= system.allocate() == 0;
        success *= pmech.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        // At zero machine power the governor holds no fuel flow.
        success *= isEqual(pmech.read(), static_cast<ScalarT>(0.0), kTol);

        const auto* residual = governor.getResidual().getData();
        for (IdxT row = 0; row < governor.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      /// The same connection with GENSAL as the machine.
      TestOutcome gensalGastPti()
      {
        using MachineExternal  = PhasorDynamics::GensalExternalVariables;
        using GovernorInternal = PhasorDynamics::Governor::GastPtiInternalVariables;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech;

        // GENSAL has no parameterless constructor; its documented defaults
        // are the zero-power machine this case needs.
        PhasorDynamics::GensalData<RealT, IdxT> machine_data;

        PhasorDynamics::Gensal<ScalarT, IdxT>            machine(&bus, machine_data);
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> governor;

        machine.getSignals().template attachSignalNode<MachineExternal::PM>(&pmech);
        governor.getSignals().template assignSignalNode<GovernorInternal::PMECH>(&pmech);

        system.addBus(&bus);
        system.addComponent(&machine);
        system.addComponent(&governor);

        success *= system.allocate() == 0;
        success *= pmech.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        // At zero machine power the governor holds no fuel flow.
        success *= isEqual(pmech.read(), static_cast<ScalarT>(0.0), kTol);

        const auto* residual = governor.getResidual().getData();
        for (IdxT row = 0; row < governor.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
