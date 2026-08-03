#pragma once

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Verify REECB signal attachment to REGCA, the coupled steady state, and
    /// time-domain recovery of the closed loop.
    template <typename real_type, typename index_type>
    class ReecbIntegrationTests
    {
    public:
      using RealT = real_type;
      using IdxT  = index_type;

      /// Full command and branch-power feedback loop between REECB and REGCA.
      TestOutcome regca()
      {
        return checkConnection(makeClosedLoopCase(), true, __func__);
      }

      /// Command-only wiring; REECB reconstructs feedback from its commands.
      TestOutcome regcaReconstructedFeedback()
      {
        return checkConnection(makeCommandOnlyCase(), false, __func__);
      }

      /// Perturbation probes confirm response direction, rate, gating, and
      /// priority coupling across the closed REGCA/REECB loop.
      TestOutcome regcaLoopResponse()
      {
        using I = PhasorDynamics::Converter::ReecbInternalVariables;
        using R = PhasorDynamics::Converter::RegcaInternalVariables;

        TestStatus success = true;
        SystemT    system(makeClosedLoopCase());

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;

        auto* converter  = dynamic_cast<RegcaT*>(system.getComponent(kConverterComponentId));
        auto* controller = dynamic_cast<ReecbT*>(system.getComponent(kControllerComponentId));

        if (converter == nullptr || controller == nullptr)
        {
          success = false;
          return success.report(__func__);
        }

        const auto reecb = [&](I variable)
        { return controller->getVariableIndex(static_cast<IdxT>(variable)); };
        const auto regca = [&](R variable)
        { return converter->getVariableIndex(static_cast<IdxT>(variable)); };

        const RealT d = kProbeDelta;

        success *= runProbes(
            system,
            {
                {"converter reactive current chases the command", {reecb(I::IQCMD), d}, {}, regca(R::IQ), d / kTg},
                {"converter active current chases the command", {reecb(I::IPCMD), d}, {}, regca(R::IP), d / kTg},
                {"power measurement tracks the branch feedback", {regca(R::PBR), d}, {}, reecb(I::PMEAS), d / kTp},
                {"voltage measurement tracks the terminal magnitude", {reecb(I::VT), d}, {}, reecb(I::VMEAS), d / kTrv},
                {"controller restores its published reactive command", {reecb(I::IQCMD), d}, {}, reecb(I::IQCMD), -d},
                {"reactive command consumes low-priority headroom", {reecb(I::IQCMD), d}, {}, reecb(I::ILMAX), -(TWO<RealT> * kInitialReactivePower + d) * d},
                {"active command follows the power order", {reecb(I::PORD), d}, {}, reecb(I::IPCMD), d},
                {"volt-var loop integrates a voltage rise downward", {reecb(I::VMEAS), d}, {}, reecb(I::XPIV), -kKvi * d},
                {"voltage dip freezes the volt-var integrator", {reecb(I::VT), -kDipDepth}, {reecb(I::VMEAS), d}, reecb(I::XPIV), ZERO<RealT>},
            });

        // The probes restore every state they touch: the system must still be
        // at the initialized equilibrium.
        success *= system.evaluateResidual() == 0;
        success *= allNearZero(system.getResidual());

        return success.report(__func__);
      }

      /// Displaced measurement and current lag states integrate back to the
      /// initialized equilibrium, closing the REGCA/REECB loop in time domain.
      TestOutcome regcaLoopRecovery()
      {
        using I = PhasorDynamics::Converter::ReecbInternalVariables;
        using R = PhasorDynamics::Converter::RegcaInternalVariables;

        TestStatus success = true;
        SystemT    system(makeClosedLoopCase());

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;

        auto* converter  = dynamic_cast<RegcaT*>(system.getComponent(kConverterComponentId));
        auto* controller = dynamic_cast<ReecbT*>(system.getComponent(kControllerComponentId));

        if (converter == nullptr || controller == nullptr)
        {
          success = false;
          return success.report(__func__);
        }

        const auto*              equilibrium_values = system.y().getData();
        const std::vector<RealT> equilibrium(
            equilibrium_values,
            equilibrium_values + static_cast<size_t>(system.y().getSize()));

        // The displacements stay strictly inside every limiter, deadband, and
        // voltage band, so the return path is a smooth interior trajectory.
        auto* y                                                       = system.y().getData();
        y[controller->getVariableIndex(static_cast<IdxT>(I::VMEAS))] += kRecoveryDelta;
        y[controller->getVariableIndex(static_cast<IdxT>(I::PMEAS))] -= kRecoveryDelta;
        y[converter->getVariableIndex(static_cast<IdxT>(R::IQ))]     += kRecoveryDelta;
        y[converter->getVariableIndex(static_cast<IdxT>(R::IP))]     -= kRecoveryDelta;
        system.y().setDataUpdated();

        AnalysisManager::Sundials::Ida<RealT, IdxT> ida(&system);
        success *= ida.configureSimulation() == 0;
        success *= ida.initializeSimulation(ZERO<RealT>) == 0;
        // The step callback keeps the model state current so the final point
        // can be compared against the stored equilibrium.
        success *= ida.runSimulation(kRecoveryHorizon, kRecoveryMonitorStep, [](RealT) {}) == 0;

        const auto* final_values = system.y().getData();
        for (size_t entry = 0; entry < equilibrium.size(); ++entry)
        {
          const RealT deviation = std::abs(static_cast<RealT>(final_values[entry]) - equilibrium[entry]);
          if (deviation > kRecoveryTolerance)
          {
            std::cout << "State " << entry << " remains " << deviation
                      << " from its equilibrium after recovery\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

    private:
      using SystemDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;
      using SystemT     = PhasorDynamics::SystemModel<RealT, IdxT>;
      using RegcaT      = PhasorDynamics::Converter::Regca<RealT, IdxT>;
      using ReecbT      = PhasorDynamics::Converter::Reecb<RealT, IdxT>;

      static constexpr IdxT kBusId                 = static_cast<IdxT>(23);
      static constexpr IdxT kIpcmdSignalId         = static_cast<IdxT>(201);
      static constexpr IdxT kIqcmdSignalId         = static_cast<IdxT>(202);
      static constexpr IdxT kPbranchSignalId       = static_cast<IdxT>(203);
      static constexpr IdxT kQbranchSignalId       = static_cast<IdxT>(204);
      static constexpr IdxT kConverterComponentId  = static_cast<IdxT>(0);
      static constexpr IdxT kControllerComponentId = static_cast<IdxT>(1);

      static constexpr RealT kSystemBaseVa         = static_cast<RealT>(100.0e6);
      static constexpr RealT kConverterBaseMva     = static_cast<RealT>(100.0);
      static constexpr RealT kInitialActivePower   = static_cast<RealT>(0.4);
      static constexpr RealT kInitialReactivePower = static_cast<RealT>(0.05);

      // Case parameters shared with the probe-response expectations.
      static constexpr RealT kTg  = static_cast<RealT>(0.02);
      static constexpr RealT kTp  = static_cast<RealT>(0.02);
      static constexpr RealT kTrv = static_cast<RealT>(0.02);
      static constexpr RealT kKvi = static_cast<RealT>(5.0);

      // Probes stay clear of every smoothing transition, so responses are
      // exact; the dip lands the terminal voltage well below Vdip.
      static constexpr RealT kProbeDelta = static_cast<RealT>(1.0e-3);
      static constexpr RealT kDipDepth   = static_cast<RealT>(0.5);

      // The slowest closed-loop mode of the case pairs the reactive and
      // voltage integrators with a time constant near three seconds; the
      // horizon leaves the displaced trajectory well inside the tolerance.
      static constexpr RealT kRecoveryDelta       = static_cast<RealT>(2.0e-3);
      static constexpr RealT kRecoveryHorizon     = static_cast<RealT>(25.0);
      static constexpr RealT kRecoveryMonitorStep = static_cast<RealT>(1.0 / 60.0);
      static constexpr RealT kRecoveryTolerance   = static_cast<RealT>(1.0e-6);

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      TestOutcome checkConnection(const SystemDataT& data,
                                  bool               feedback_attached,
                                  const char*        test_name)
      {
        using PhasorDynamics::Converter::ReecbExternalVariables;
        using PhasorDynamics::Converter::ReecbInternalVariables;
        using PhasorDynamics::Converter::RegcaExternalVariables;
        using PhasorDynamics::Converter::RegcaInternalVariables;

        TestStatus success = true;
        SystemT    system(data);

        success *= system.allocate() == 0;

        auto* converter  = dynamic_cast<RegcaT*>(system.getComponent(kConverterComponentId));
        auto* controller = dynamic_cast<ReecbT*>(system.getComponent(kControllerComponentId));
        auto* ipcmd      = system.getSignal(kIpcmdSignalId);
        auto* iqcmd      = system.getSignal(kIqcmdSignalId);
        auto* pbranch    = feedback_attached ? system.getSignal(kPbranchSignalId) : nullptr;
        auto* qbranch    = feedback_attached ? system.getSignal(kQbranchSignalId) : nullptr;

        if (converter == nullptr || controller == nullptr || ipcmd == nullptr || iqcmd == nullptr
            || (feedback_attached && (pbranch == nullptr || qbranch == nullptr)))
        {
          success = false;
          return success.report(test_name);
        }

        bool signals_linked = ipcmd->linked() && iqcmd->linked();
        if (feedback_attached)
        {
          signals_linked = signals_linked && pbranch->linked() && qbranch->linked();
        }
        success *= signals_linked;
        if (!signals_linked)
        {
          return success.report(test_name);
        }

        auto& converter_signals  = converter->getSignals();
        auto& controller_signals = controller->getSignals();

        bool ports_connected =
            controller_signals.template isAssigned<ReecbInternalVariables::IPCMD>()
            && controller_signals.template isAssigned<ReecbInternalVariables::IQCMD>()
            && converter_signals.template isAttached<RegcaExternalVariables::IPCMD>()
            && converter_signals.template isAttached<RegcaExternalVariables::IQCMD>();
        if (feedback_attached)
        {
          ports_connected = ports_connected
                            && converter_signals.template isAssigned<RegcaInternalVariables::PBR>()
                            && converter_signals.template isAssigned<RegcaInternalVariables::QBR>()
                            && controller_signals.template isAttached<ReecbExternalVariables::PE>()
                            && controller_signals.template isAttached<ReecbExternalVariables::QGEN>();
        }
        else
        {
          ports_connected = ports_connected
                            && !controller_signals.template isAttached<ReecbExternalVariables::PE>()
                            && !controller_signals.template isAttached<ReecbExternalVariables::QGEN>();
        }
        success *= ports_connected;
        if (!ports_connected)
        {
          return success.report(test_name);
        }

        // Optional references stay unattached and latch their initialized setpoints.
        success *= !controller_signals.template isAttached<ReecbExternalVariables::QEXT>();
        success *= !controller_signals.template isAttached<ReecbExternalVariables::PFAREF>();
        success *= !controller_signals.template isAttached<ReecbExternalVariables::PREF>();

        // Each linked signal is one shared global unknown: the node, the
        // publishing state, and the subscribing port agree on its index.
        const IdxT ipcmd_index = ipcmd->getVariableIndex();
        const IdxT iqcmd_index = iqcmd->getVariableIndex();

        success *= ipcmd_index != iqcmd_index;
        success *= ipcmd_index
                   == controller->getVariableIndex(
                       static_cast<IdxT>(ReecbInternalVariables::IPCMD));
        success *= ipcmd_index
                   == converter_signals.template readExternalVariableIndex<
                       RegcaExternalVariables::IPCMD>();
        success *= iqcmd_index
                   == controller->getVariableIndex(
                       static_cast<IdxT>(ReecbInternalVariables::IQCMD));
        success *= iqcmd_index
                   == converter_signals.template readExternalVariableIndex<
                       RegcaExternalVariables::IQCMD>();

        if (feedback_attached)
        {
          const IdxT pbranch_index = pbranch->getVariableIndex();
          const IdxT qbranch_index = qbranch->getVariableIndex();

          success *= pbranch_index != qbranch_index;
          success *= pbranch_index != ipcmd_index;
          success *= pbranch_index != iqcmd_index;
          success *= qbranch_index != ipcmd_index;
          success *= qbranch_index != iqcmd_index;
          success *= pbranch_index
                     == converter->getVariableIndex(
                         static_cast<IdxT>(RegcaInternalVariables::PBR));
          success *= pbranch_index
                     == controller_signals.template readExternalVariableIndex<
                         ReecbExternalVariables::PE>();
          success *= qbranch_index
                     == converter->getVariableIndex(
                         static_cast<IdxT>(RegcaInternalVariables::QBR));
          success *= qbranch_index
                     == controller_signals.template readExternalVariableIndex<
                         ReecbExternalVariables::QGEN>();
        }

        // The published REGCA operating point initializes the coupled pair to
        // an exact steady state.
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= allNearZero(system.yp());
        success *= allNearZero(system.getResidual());

        return success.report(test_name);
      }

      /// One state write applied before a probe evaluation.
      struct Write
      {
        IdxT  state{INVALID_INDEX<IdxT>};
        RealT delta{};
      };

      /// One perturbation and the exact residual response it must produce.
      struct Probe
      {
        const char* label;
        Write       first;
        Write       second{};
        IdxT        row{};
        RealT       expected{};
      };

      /// Apply each probe to the initialized system, check the responding
      /// residual row against its exact expectation, and restore the state.
      TestStatus runProbes(SystemT& system, std::initializer_list<Probe> probes)
      {
        TestStatus success = true;
        auto*      y       = system.y().getData();

        for (const auto& probe : probes)
        {
          const RealT first_base  = static_cast<RealT>(y[probe.first.state]);
          RealT       second_base = ZERO<RealT>;

          y[probe.first.state] = first_base + probe.first.delta;
          if (probe.second.state != INVALID_INDEX<IdxT>)
          {
            second_base           = static_cast<RealT>(y[probe.second.state]);
            y[probe.second.state] = second_base + probe.second.delta;
          }
          system.y().setDataUpdated();
          success *= system.evaluateResidual() == 0;

          const RealT response = static_cast<RealT>(system.getResidual().getData()[probe.row]);
          if (!isEqual(response, probe.expected, kTol))
          {
            std::cout << "Probe '" << probe.label << "' expected " << probe.expected
                      << " but produced " << response << '\n';
            success = false;
          }

          y[probe.first.state] = first_base;
          if (probe.second.state != INVALID_INDEX<IdxT>)
          {
            y[probe.second.state] = second_base;
          }
        }
        system.y().setDataUpdated();

        return success;
      }

      static SystemDataT makeCommandOnlyCase()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Converter;

        SystemDataT data;
        data.va_base = kSystemBaseVa;

        auto& bus    = data.bus.emplace_back();
        bus.bus_id   = kBusId;
        bus.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
        bus.Vr0      = ONE<RealT>;
        bus.Vi0      = ZERO<RealT>;

        data.signal = {{"Active Current Command", kIpcmdSignalId},
                       {"Reactive Current Command", kIqcmdSignalId}};

        auto& converter                                   = data.regca.emplace_back();
        converter.buses[RegcaBuses::bus]                  = kBusId;
        converter.signal_inputs[RegcaSignalInputs::ipcmd] = kIpcmdSignalId;
        converter.signal_inputs[RegcaSignalInputs::iqcmd] = kIqcmdSignalId;
        converter.parameters[RegcaParameters::p0]         = kInitialActivePower;
        converter.parameters[RegcaParameters::q0]         = kInitialReactivePower;
        converter.parameters[RegcaParameters::mva]        = kConverterBaseMva;
        converter.parameters[RegcaParameters::Tg]         = kTg;
        converter.parameters[RegcaParameters::TM]         = static_cast<RealT>(0.02);
        converter.parameters[RegcaParameters::Rqmax]      = static_cast<RealT>(999.0);
        converter.parameters[RegcaParameters::Rqmin]      = static_cast<RealT>(-999.0);
        converter.parameters[RegcaParameters::Rpmax]      = static_cast<RealT>(999.0);
        converter.parameters[RegcaParameters::sL]         = true;
        converter.parameters[RegcaParameters::IL1]        = static_cast<RealT>(1.1);
        converter.parameters[RegcaParameters::VL0]        = static_cast<RealT>(0.4);
        converter.parameters[RegcaParameters::VL1]        = static_cast<RealT>(0.9);
        converter.parameters[RegcaParameters::VA0]        = static_cast<RealT>(0.4);
        converter.parameters[RegcaParameters::VA1]        = static_cast<RealT>(0.9);
        converter.parameters[RegcaParameters::Vhvmax]     = static_cast<RealT>(1.2);

        auto& controller                                     = data.reecb.emplace_back();
        controller.buses[ReecbBuses::bus]                    = kBusId;
        controller.signal_outputs[ReecbSignalOutputs::ipcmd] = kIpcmdSignalId;
        controller.signal_outputs[ReecbSignalOutputs::iqcmd] = kIqcmdSignalId;
        controller.parameters[ReecbParameters::mva]          = kConverterBaseMva;
        controller.parameters[ReecbParameters::Trv]          = kTrv;
        controller.parameters[ReecbParameters::Tp]           = kTp;
        controller.parameters[ReecbParameters::Kvi]          = kKvi;
        controller.parameters[ReecbParameters::QFlag]        = true;
        controller.parameters[ReecbParameters::VFlag]        = true;

        return data;
      }

      static SystemDataT makeClosedLoopCase()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Converter;

        auto data = makeCommandOnlyCase();

        data.signal.push_back({"Branch Active Power", kPbranchSignalId});
        data.signal.push_back({"Branch Reactive Power", kQbranchSignalId});

        auto& converter                                       = data.regca.front();
        converter.signal_outputs[RegcaSignalOutputs::pbranch] = kPbranchSignalId;
        converter.signal_outputs[RegcaSignalOutputs::qbranch] = kQbranchSignalId;

        auto& controller                                  = data.reecb.front();
        controller.signal_inputs[ReecbSignalInputs::pe]   = kPbranchSignalId;
        controller.signal_inputs[ReecbSignalInputs::qgen] = kQbranchSignalId;

        return data;
      }

      template <typename VectorT>
      static bool allNearZero(const VectorT& vector)
      {
        const auto* values = vector.getData();
        for (IdxT entry = 0; entry < vector.getSize(); ++entry)
        {
          if (!isEqual(values[entry], ZERO<RealT>, kTol))
          {
            return false;
          }
        }
        return true;
      }
    };
  } // namespace Testing
} // namespace GridKit
