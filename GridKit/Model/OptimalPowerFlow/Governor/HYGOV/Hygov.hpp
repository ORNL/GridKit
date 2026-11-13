/**
 * @file Hygov.hpp
 * @brief Declaration of an optimal power flow HYGOV governor.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      /// External variables of a `Hygov`, the mechanical power of its signal
      enum class HygovExternalVariables : size_t
      {
        PM, ///< \f$P_\mathrm{m}\f$ Mechanical power [p.u.]
      };

      /// Internal constraints of a `Hygov`
      enum class HygovInternalConstraints : size_t
      {
        GATE, ///< \f$g\f$ Gate position [p.u.]
      };

      /// Hydro turbine governor at steady state: the gate limits on the mechanical power, which are linear, so it has no Hessian
      template <typename scalar_type, typename index_type>
      class Hygov final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = HygovData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = HygovExternalVariables;
        using InternalConstraintsT = HygovInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<HygovSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<HygovSignals>()> NUM_PORT_CONSTRAINTS{0};

        explicit Hygov(const ModelDataT& data);

        Model::ConfigurationChecks verify() const override;
        int                        tagInequality() override;
        int                        initialize() override;

        int evaluateGradient() override;
        int evaluateJacobian() override;

        [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

      private:
        void initializeParameters(const ModelDataT& data);
        void setDerivedParameters();

        RealT Trate_{100.0};
        RealT Gmax_{ONE<RealT>};
        RealT Gmin_{ZERO<RealT>};
        RealT At_{1.2};
        RealT Qnl_{0.05};
        RealT Hdam_{ONE<RealT>};
      };
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
