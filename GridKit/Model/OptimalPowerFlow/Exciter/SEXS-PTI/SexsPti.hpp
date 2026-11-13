/**
 * @file SexsPti.hpp
 * @brief Declaration of an optimal power flow SEXS-PTI exciter.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      /// External variables of a `SexsPti`, the field voltage of its signal
      enum class SexsPtiExternalVariables : size_t
      {
        EFD, ///< \f$E_\mathrm{fd}\f$ Field voltage [p.u.]
      };

      /// Internal constraints of a `SexsPti`
      enum class SexsPtiInternalConstraints : size_t
      {
        EFD, ///< \f$E_\mathrm{fd}\f$ Exciter field voltage output [p.u.]
      };

      /// Simplified excitation system at steady state: the limits on the field voltage, which are linear, so it has no Hessian
      template <typename scalar_type, typename index_type>
      class SexsPti final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = SexsPtiData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = SexsPtiExternalVariables;
        using InternalConstraintsT = SexsPtiInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<SexsPtiSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<SexsPtiSignals>()> NUM_PORT_CONSTRAINTS{0};

        explicit SexsPti(const ModelDataT& data);

        Model::ConfigurationChecks verify() const override;
        int                        tagInequality() override;
        int                        initialize() override;

        int evaluateGradient() override;
        int evaluateJacobian() override;

        [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

      private:
        void initializeParameters(const ModelDataT& data);

        RealT Efdmin_{ZERO<RealT>};
        RealT Efdmax_{ZERO<RealT>};
      };
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
