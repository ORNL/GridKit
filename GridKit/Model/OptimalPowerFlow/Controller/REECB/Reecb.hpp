/**
 * @file Reecb.hpp
 * @brief Declaration of an optimal power flow REECB controller.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      /// External variables of a `Reecb`, the values of its signals
      enum class ReecbExternalVariables : size_t
      {
        IPCMD, ///< \f$I_p^\mathrm{cmd}\f$ Active-current command [p.u.]
        IQCMD, ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command [p.u.]
        PE,    ///< \f$P_\mathrm{e}\f$ Active-power feedback [p.u.]
        QGEN,  ///< \f$Q^\mathrm{gen}\f$ Reactive-power feedback [p.u.]
      };

      /// Internal constraints of a `Reecb`
      enum class ReecbInternalConstraints : size_t
      {
        P,  ///< Active-power order [p.u.]
        Q,  ///< Reactive power [p.u.]
        I2, ///< Squared converter current [p.u.]
      };

      /// Renewable electrical controller at steady state: the power and current limits on the converter commands
      template <typename scalar_type, typename index_type>
      class Reecb final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = ReecbData<RealT, IdxT>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;

      public:
        using InternalVariablesT   = Utilities::EmptyEnum;
        using ExternalVariablesT   = ReecbExternalVariables;
        using InternalConstraintsT = ReecbInternalConstraints;
        using ExternalConstraintsT = Utilities::EmptyEnum;

        static constexpr std::array<IdxT, Utilities::enum_size<ReecbSignals>()> NUM_PORT_VARIABLES{
            SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<ReecbSignals>()> NUM_PORT_CONSTRAINTS{0, 0, 0, 0};

        explicit Reecb(const ModelDataT& data);

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

        RealT mva_{ZERO<RealT>};
        RealT Pmax_{ONE<RealT>};
        RealT Pmin_{ZERO<RealT>};
        RealT Qmax_{0.436};
        RealT Qmin_{-0.436};
        RealT Imax_{1.3};
      };
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
