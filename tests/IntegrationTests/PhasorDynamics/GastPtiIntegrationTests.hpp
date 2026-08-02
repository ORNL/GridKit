#pragma once

#include <limits>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Verify GASTPTI signal attachment to GENROU and GENSAL.
    template <typename real_type, typename index_type>
    class GastPtiIntegrationTests
    {
    public:
      using RealT = real_type;
      using IdxT  = index_type;

      TestOutcome genrou()
      {
        using namespace PhasorDynamics;

        return checkConnection<Genrou<RealT, IdxT>,
                               GenrouInternalVariables::OMEGA,
                               GenrouExternalVariables::PM>(makeGenrouCase(), __func__);
      }

      TestOutcome gensal()
      {
        using namespace PhasorDynamics;

        return checkConnection<Gensal<RealT, IdxT>,
                               GensalInternalVariables::OMEGA,
                               GensalExternalVariables::PM>(makeGensalCase(), __func__);
      }

    private:
      using SystemDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;
      using SystemT     = PhasorDynamics::SystemModel<RealT, IdxT>;
      using GastPtiT    = PhasorDynamics::Governor::GastPti<RealT, IdxT>;

      static constexpr IdxT kBusId               = static_cast<IdxT>(17);
      static constexpr IdxT kSpeedSignalId       = static_cast<IdxT>(101);
      static constexpr IdxT kPmechSignalId       = static_cast<IdxT>(102);
      static constexpr IdxT kMachineComponentId  = static_cast<IdxT>(0);
      static constexpr IdxT kGovernorComponentId = static_cast<IdxT>(1);

      static constexpr RealT kSystemBaseVa         = static_cast<RealT>(100.0e6);
      static constexpr RealT kMachineBaseMva       = static_cast<RealT>(100.0);
      static constexpr RealT kTurbineRatingMw      = static_cast<RealT>(100.0);
      static constexpr RealT kInitialActivePower   = static_cast<RealT>(0.4);
      static constexpr RealT kInitialReactivePower = static_cast<RealT>(0.05);
      static constexpr RealT kBehaviorTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

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

        if (machine == nullptr || governor == nullptr || speed == nullptr || pmech == nullptr)
        {
          success = false;
          return success.report(test_name);
        }

        const bool signals_linked  = speed->linked() && pmech->linked();
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
            && governor_signals.template isAssigned<GastPtiInternalVariables::PMECH>();
        success *= ports_connected;
        if (!ports_connected)
        {
          return success.report(test_name);
        }

        const IdxT speed_index = speed->getVariableIndex();
        const IdxT pmech_index = pmech->getVariableIndex();

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

        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= allNearZero(system.yp());
        success *= allNearZero(system.getResidual());

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
                       {"Mechanical Power", kPmechSignalId}};

        auto& governor                                       = data.gastpti.emplace_back();
        governor.signal_inputs[GastPtiSignalInputs::speed]   = kSpeedSignalId;
        governor.signal_outputs[GastPtiSignalOutputs::pmech] = kPmechSignalId;
        governor.parameters[GastPtiParameters::Trate]        = kTurbineRatingMw;

        return data;
      }

      static SystemDataT makeGenrouCase()
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

      static SystemDataT makeGensalCase()
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

      template <typename VectorT>
      static bool allNearZero(const VectorT& vector)
      {
        const auto* values = vector.getData();
        for (IdxT entry = 0; entry < vector.getSize(); ++entry)
        {
          if (!isEqual(values[entry], ZERO<RealT>, kBehaviorTol))
          {
            return false;
          }
        }
        return true;
      }
    };
  } // namespace Testing
} // namespace GridKit
