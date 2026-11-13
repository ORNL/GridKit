/**
 * @file Tgov1.hpp
 * @brief Declaration of an optimal power flow TGOV1 governor.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      /// External variables of a `Tgov1`, the mechanical power of its signal
      enum class Tgov1ExternalVariables : size_t
      {
        PM, ///< \f$P_\mathrm{m}\f$ Mechanical power [p.u.]
      };

      /// Internal constraints of a `Tgov1`
      enum class Tgov1InternalConstraints : size_t
      {
        PV, ///< \f$P_\mathrm{v}\f$ Valve position [p.u.]
      };

      /// Steam turbine governor at steady state: the valve limits on the mechanical power, which are linear, so it has no Hessian
      template <typename scalar_type, typename index_type>
      class Tgov1 final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = Tgov1Data<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = Tgov1ExternalVariables;
        using InternalConstraintsT = Tgov1InternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<Tgov1Signals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<Tgov1Signals>()> NUM_PORT_CONSTRAINTS{0};

        explicit Tgov1(const ModelDataT& data);

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
        RealT Pvmax_{ONE<RealT>};
        RealT Pvmin_{ZERO<RealT>};
      };
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
