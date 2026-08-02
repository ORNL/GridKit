#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename real_type, typename index_type>
    class PDComponentIntegrationTests
    {
    public:
      using RealT = real_type;
      using IdxT  = index_type;

      /// Verify the GENROU-to-GASTPTI signal mapping and initialization order.
      TestOutcome genrouGastPtiInitialization()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Governor;

        using GenrouT  = Genrou<RealT, IdxT>;
        using GastPtiT = GastPti<RealT, IdxT>;

        TestStatus               success = true;
        SystemModel<RealT, IdxT> system(makeGenrouGastPtiCase());

        success *= system.allocate() == 0;

        auto* generator = dynamic_cast<GenrouT*>(system.getComponent(kGenrouComponentId));
        auto* governor  = dynamic_cast<GastPtiT*>(system.getComponent(kGastPtiComponentId));
        auto* speed     = system.getSignal(kSpeedSignalId);
        auto* pmech     = system.getSignal(kPmechSignalId);
        auto* pref      = system.getSignal(kPrefSignalId);

        success *= generator != nullptr;
        success *= governor != nullptr;
        success *= speed != nullptr;
        success *= pmech != nullptr;
        success *= pref != nullptr;
        if (generator == nullptr || governor == nullptr || speed == nullptr
            || pmech == nullptr || pref == nullptr)
        {
          return success.report(__func__);
        }

        success *= speed->linked();
        success *= pmech->linked();
        success *= pref->linked();
        success *= governor->getSignals()
                       .template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>()
                   == generator->getVariableIndex(
                       static_cast<IdxT>(GenrouInternalVariables::OMEGA));
        success *= generator->getSignals()
                       .template readExternalVariableIndex<GenrouExternalVariables::PM>()
                   == governor->getVariableIndex(
                       static_cast<IdxT>(GastPtiInternalVariables::PMECH));

        success                        *= generator->initialize() == 0;
        const RealT machine_pmech_seed  = pmech->read();

        auto* governor_y  = governor->y().getData();
        auto* governor_yp = governor->yp().getData();
        for (std::size_t row = 0; row < governor->y().getSize(); ++row)
        {
          governor_y[row]  = static_cast<RealT>(-9.0);
          governor_yp[row] = static_cast<RealT>(7.0);
        }
        pref->init(static_cast<RealT>(-0.25));

        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        success *= pmech->read() == machine_pmech_seed;
        success *= isEqual(machine_pmech_seed, kSystemPmechSeed, kTol);
        success *= isEqual(pref->read(), machine_pmech_seed, kTol);

        const std::array<RealT, 7> expected_governor_state{{
            kComponentFuelFlow,
            kComponentFuelFlow,
            kComponentFuelFlow,
            kComponentFuelFlow,
            static_cast<RealT>(1.4),
            kComponentFuelFlow,
            machine_pmech_seed,
        }};
        for (std::size_t row = 0; row < expected_governor_state.size(); ++row)
        {
          success *= isEqual(governor_y[row], expected_governor_state[row], kTol);
          success *= governor_yp[row] == ZERO<RealT>;
        }

        const auto& residual = system.getResidual();
        for (std::size_t row = 0; row < residual.getSize(); ++row)
        {
          success *= isEqual(residual.getData()[row], ZERO<RealT>, kTol);
        }

        return success.report(__func__);
      }

    private:
      using SystemDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;

      static constexpr IdxT kBusId              = static_cast<IdxT>(1);
      static constexpr IdxT kGenrouComponentId  = static_cast<IdxT>(0);
      static constexpr IdxT kGastPtiComponentId = static_cast<IdxT>(1);
      static constexpr IdxT kSpeedSignalId      = static_cast<IdxT>(1);
      static constexpr IdxT kPmechSignalId      = static_cast<IdxT>(2);
      static constexpr IdxT kPrefSignalId       = static_cast<IdxT>(3);

      static constexpr RealT kSystemPmechSeed   = static_cast<RealT>(0.4);
      static constexpr RealT kComponentFuelFlow = static_cast<RealT>(0.8);
      static constexpr RealT kComponentBaseMva  = static_cast<RealT>(50.0);
      static constexpr RealT kTol               = static_cast<RealT>(1.0e-12);

      static SystemDataT makeGenrouGastPtiCase()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Governor;

        SystemDataT data;

        data.bus.resize(1);
        data.bus[0].bus_id   = kBusId;
        data.bus[0].bus_type = BusData<RealT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = ONE<RealT>;
        data.bus[0].Vi0      = ZERO<RealT>;

        data.signal = {{"Machine Speed Deviation", kSpeedSignalId},
                       {"Mechanical Power", kPmechSignalId},
                       {"Load Reference", kPrefSignalId}};

        auto& generator                                      = data.genrou.emplace_back();
        generator.buses[GenrouBuses::bus]                    = kBusId;
        generator.signal_outputs[GenrouSignalOutputs::speed] = kSpeedSignalId;
        generator.signal_inputs[GenrouSignalInputs::pmech]   = kPmechSignalId;
        generator.parameters[GenrouParameters::p0]           = kSystemPmechSeed;
        generator.parameters[GenrouParameters::q0]           = RealT{0.05};
        generator.parameters[GenrouParameters::H]            = RealT{3.0};
        generator.parameters[GenrouParameters::D]            = ZERO<RealT>;
        generator.parameters[GenrouParameters::Ra]           = ZERO<RealT>;
        generator.parameters[GenrouParameters::Tdop]         = RealT{7.0};
        generator.parameters[GenrouParameters::Tdopp]        = RealT{0.04};
        generator.parameters[GenrouParameters::Tqop]         = RealT{0.75};
        generator.parameters[GenrouParameters::Tqopp]        = RealT{0.05};
        generator.parameters[GenrouParameters::Xd]           = RealT{2.1};
        generator.parameters[GenrouParameters::Xdp]          = RealT{0.2};
        generator.parameters[GenrouParameters::Xdpp]         = RealT{0.18};
        generator.parameters[GenrouParameters::Xq]           = RealT{0.5};
        generator.parameters[GenrouParameters::Xqp]          = RealT{0.5};
        generator.parameters[GenrouParameters::Xqpp]         = RealT{0.18};
        generator.parameters[GenrouParameters::Xl]           = RealT{0.15};
        generator.parameters[GenrouParameters::S10]          = ZERO<RealT>;
        generator.parameters[GenrouParameters::S12]          = ZERO<RealT>;
        generator.parameters[GenrouParameters::mva]          = kComponentBaseMva;

        auto& governor                                       = data.gastpti.emplace_back();
        governor.signal_inputs[GastPtiSignalInputs::speed]   = kSpeedSignalId;
        governor.signal_inputs[GastPtiSignalInputs::pref]    = kPrefSignalId;
        governor.signal_outputs[GastPtiSignalOutputs::pmech] = kPmechSignalId;
        governor.parameters[GastPtiParameters::Trate]        = kComponentBaseMva;

        auto& reference                                                 = data.constant_source.emplace_back();
        reference.signal_outputs[ConstantSignalSourceSignalOutputs::sr] = kPrefSignalId;
        reference.parameters[ConstantSignalSourceParameters::Sr]        = RealT{-0.25};

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
