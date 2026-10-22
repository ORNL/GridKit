/**
 * @file InverseTimeRelay.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the InverseTimeRelay phasor-dynamics relay model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <numbers>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/InverseTimeRelay/InverseTimeRelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Internal variables of an `InverseTimeRelay`.
      enum class InverseTimeRelayInternalVariables : size_t
      {
        X, ///< \f$x\f$ Differential relay latch [-]
        S, ///< \f$s\f$ Algebraic trip command [-]
      };

      /// External signal variables read by an `InverseTimeRelay`.
      enum class InverseTimeRelayExternalVariables : size_t
      {
        IR, ///< \f$I_r\f$ Known measured current, real component [p.u.]
        II, ///< \f$I_i\f$ Known measured current, imaginary component [p.u.]
      };

      /**
       * @brief Inverse-time overcurrent relay with optional lockout.
       *
       * Above pickup, the trip rises once the integral of the inverse IEEE
       * C37.112 operate time reaches one.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class InverseTimeRelay : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT         = InverseTimeRelayData<RealT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<InverseTimeRelay, InverseTimeRelayData>;
        using InternalVariablesT = InverseTimeRelayInternalVariables;
        using ExternalVariablesT = InverseTimeRelayExternalVariables;

        InverseTimeRelay(const ModelDataT& data);
        ~InverseTimeRelay();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        Model::ConfigurationChecks verify() const override final;

        SignalPortsT& getPorts()
        {
          return ports_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static constexpr RealT TRIP_LEVEL            = static_cast<RealT>(0.75);
        static constexpr RealT LN4                   = TWO<RealT> * std::numbers::ln2_v<RealT>;

        RealT Ipickup_{ZERO<RealT>};
        RealT TD_{ZERO<RealT>};
        RealT A_{ZERO<RealT>};
        RealT B_{ZERO<RealT>};
        RealT exponent_{ZERO<RealT>}; ///< \f$\alpha\f$; Component::alpha_ is the solver coefficient
        RealT tr_{ZERO<RealT>};
        bool  lockout_{true};
        bool  has_required_parameters_{false};

        // Derived parameters
        RealT Treset_{ZERO<RealT>};
        RealT self_reset_{ZERO<RealT>}; ///< \f$1 - s_{\mathrm{lock}}\f$

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
