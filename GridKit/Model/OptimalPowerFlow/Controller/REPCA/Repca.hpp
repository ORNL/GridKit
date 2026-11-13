/**
 * @file Repca.hpp
 * @brief Declaration of an optimal power flow REPCA controller.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      /// External variables of a `Repca`, the values of its signals
      enum class RepcaExternalVariables : size_t
      {
        PBR, ///< \f$P\f$ Branch active power [p.u.]
        QBR, ///< \f$Q\f$ Branch reactive power [p.u.]
      };

      /// Internal constraints of a `Repca`
      enum class RepcaInternalConstraints : size_t
      {
        Q, ///< Reactive-power command [p.u.]
        P, ///< Active-power command [p.u.]
      };

      /// Renewable plant controller at steady state: the plant limits on the branch power, which are linear, so it has no Hessian
      template <typename scalar_type, typename index_type>
      class Repca final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = RepcaData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = RepcaExternalVariables;
        using InternalConstraintsT = RepcaInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<RepcaSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<RepcaSignals>()> NUM_PORT_CONSTRAINTS{0, 0};

        explicit Repca(const ModelDataT& data);

        Model::ConfigurationChecks verify() const override;
        int                        tagInequality() override;
        int                        initialize() override;

        int evaluateGradient() override;
        int evaluateJacobian() override;

        [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

      private:
        void initializeParameters(const ModelDataT& data);
        void setDerivedParameters();

        RealT mva_{100.0};
        bool  Freqflag_{false};
        RealT Qmax_{ONE<RealT>};
        RealT Qmin_{-ONE<RealT>};
        RealT Pmax_{TWO<RealT>};
        RealT Pmin_{ZERO<RealT>};

        RealT freq_on_{ZERO<RealT>}; ///< Active-power output selector as a factor
      };
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
