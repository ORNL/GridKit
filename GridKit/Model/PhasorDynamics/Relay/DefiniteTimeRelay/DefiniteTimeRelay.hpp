/**
 * @file DefiniteTimeRelay.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the DefiniteTimeRelay phasor-dynamics relay model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <numbers>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/DefiniteTimeRelay/DefiniteTimeRelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Internal variables of a `DefiniteTimeRelay`.
      enum class DefiniteTimeRelayInternalVariables : size_t
      {
        X, ///< \f$x\f$ Differential relay latch [-]
        S, ///< \f$s\f$ Algebraic trip command [-]
      };

      /// External signal variables read by a `DefiniteTimeRelay`.
      enum class DefiniteTimeRelayExternalVariables : size_t
      {
        MR, ///< \f$m_r\f$ Known measured quantity, real component [p.u.]
        MI, ///< \f$m_i\f$ Known measured quantity, imaginary component [p.u.]
      };

      /**
       * @brief Definite-time relay with optional lockout.
       *
       * One relay models one element: it operates above the pickup level
       * when the pickup exceeds the reference level, and below it otherwise.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class DefiniteTimeRelay : public Component<scalar_type, index_type>
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
        using ModelDataT         = DefiniteTimeRelayData<RealT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<DefiniteTimeRelay, DefiniteTimeRelayData>;
        using InternalVariablesT = DefiniteTimeRelayInternalVariables;
        using ExternalVariablesT = DefiniteTimeRelayExternalVariables;

        DefiniteTimeRelay(const ModelDataT& data);
        ~DefiniteTimeRelay();

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

        RealT Mpickup_{ZERO<RealT>};
        RealT Mref_{ZERO<RealT>};
        RealT Ttrip_{ZERO<RealT>};
        bool  lockout_{true};
        bool  has_required_parameters_{false};

        // Derived parameters
        RealT Tlatch_{ZERO<RealT>};
        RealT self_reset_{ZERO<RealT>}; ///< \f$1 - s_{\mathrm{lock}}\f$

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
