/**
 * @file Esdc1a.hpp
 * @brief Declaration of an optimal power flow ESDC1A exciter.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      /// External variables of an `Esdc1a`, the field voltage of its signal
      enum class Esdc1aExternalVariables : size_t
      {
        EFD, ///< \f$E_\mathrm{fd}\f$ Field voltage [p.u.]
      };

      /// Internal constraints of an `Esdc1a`
      enum class Esdc1aInternalConstraints : size_t
      {
        VR,   ///< \f$V_R\f$ Voltage-regulator output [p.u.]
        EFDP, ///< \f$E_\mathrm{fd}'\f$ Exciter field-voltage state [p.u.]
      };

      /// IEEE DC1A excitation system at steady state: the regulator and field-voltage limits
      template <typename scalar_type, typename index_type>
      class Esdc1a final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = Esdc1aData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = Esdc1aExternalVariables;
        using InternalConstraintsT = Esdc1aInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<Esdc1aSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<Esdc1aSignals>()> NUM_PORT_CONSTRAINTS{0};

        explicit Esdc1a(const ModelDataT& data);

        Model::ConfigurationChecks verify() const override;
        int                        tagInequality() override;
        int                        initialize() override;

        int evaluateGradient() override;
        int evaluateJacobian() override;
        int evaluateHessian(RealT sigma) override;

        [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

      private:
        void initializeParameters(const ModelDataT& data);
        void setDerivedParameters();

        RealT Ke_{0.1};
        RealT Vrmax_{1.0};
        RealT Vrmin_{-1.0};
        RealT E1_{2.8};
        RealT Se1_{0.08};
        RealT E2_{3.7};
        RealT Se2_{0.33};
        bool  exclim_{true};

        RealT SA_{ZERO<RealT>};
        RealT SB_{ZERO<RealT>};
        RealT ke_on_{ONE<RealT>};  ///< Configured, not automatic, Ke as a factor
        RealT lim_on_{ONE<RealT>}; ///< Field-voltage lower limit as a factor
      };
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
