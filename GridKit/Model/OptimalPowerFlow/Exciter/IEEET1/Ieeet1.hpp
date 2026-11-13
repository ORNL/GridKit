/**
 * @file Ieeet1.hpp
 * @brief Declaration of an optimal power flow IEEET1 exciter.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      /// External variables of an `Ieeet1`, the field voltage of its signal
      enum class Ieeet1ExternalVariables : size_t
      {
        EFD, ///< \f$E_\mathrm{fd}\f$ Field voltage [p.u.]
      };

      /// Internal constraints of an `Ieeet1`
      enum class Ieeet1InternalConstraints : size_t
      {
        VR, ///< \f$V_R\f$ Voltage regulator [p.u.]
      };

      /// IEEE Type 1 excitation system at steady state: the regulator limits on the field voltage
      template <typename scalar_type, typename index_type>
      class Ieeet1 final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = Ieeet1Data<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = Ieeet1ExternalVariables;
        using InternalConstraintsT = Ieeet1InternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<Ieeet1Signals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<Ieeet1Signals>()> NUM_PORT_CONSTRAINTS{0};

        explicit Ieeet1(const ModelDataT& data);

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

        RealT Ke_{-0.06};
        RealT Vrmin_{-1.0};
        RealT Vrmax_{1.0};
        RealT E1_{2.8};
        RealT E2_{3.73};
        RealT Se1_{0.04};
        RealT Se2_{0.33};

        RealT SA_{ZERO<RealT>};
        RealT SB_{ZERO<RealT>};
        RealT ke_on_{ONE<RealT>}; ///< Configured, not automatic, Ke as a factor
      };
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
