/**
 * @file GastPti.hpp
 * @brief Declaration of an optimal power flow GASTPTI governor.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      /// External variables of a `GastPti`, the mechanical power of its signal
      enum class GastPtiExternalVariables : size_t
      {
        PM, ///< \f$P_\mathrm{m}\f$ Mechanical power [p.u.]
      };

      /// Internal constraints of a `GastPti`
      enum class GastPtiInternalConstraints : size_t
      {
        XF, ///< \f$x_F\f$ Fuel flow [p.u.]
      };

      /// Gas turbine governor at steady state: the valve and temperature limits on the mechanical power, which are linear, so it has no Hessian
      template <typename scalar_type, typename index_type>
      class GastPti final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = GastPtiData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = GastPtiExternalVariables;
        using InternalConstraintsT = GastPtiInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<GastPtiSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<GastPtiSignals>()> NUM_PORT_CONSTRAINTS{0};

        explicit GastPti(const ModelDataT& data);

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
        RealT Vmax_{ONE<RealT>};
        RealT Vmin_{ZERO<RealT>};
        RealT At_{ONE<RealT>};
      };
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
