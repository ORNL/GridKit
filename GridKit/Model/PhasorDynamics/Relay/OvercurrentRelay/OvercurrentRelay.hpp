/**
 * @file OvercurrentRelay.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the overcurrent relay model.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/OvercurrentRelay/OvercurrentRelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Internal variables of an `OvercurrentRelay`.
      enum class OvercurrentRelayInternalVariables : size_t
      {
        X,    ///< \f$x\f$ Differential lockout latch state [-]
        TRIP, ///< \f$s\f$ Algebraic trip command [-]
      };

      /// External signal variables read by an `OvercurrentRelay`.
      enum class OvercurrentRelayExternalVariables : size_t
      {
        IR, ///< \f$I_{\mathrm{r}}\f$ Known measured current, real component [p.u.]
        II, ///< \f$I_{\mathrm{i}}\f$ Known measured current, imaginary component [p.u.]
      };

      /**
       * @brief Definite-time overcurrent relay with lockout.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class OvercurrentRelay : public Component<scalar_type, index_type>
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
        using ModelDataT         = OvercurrentRelayData<RealT, IdxT>;
        using SignalNodeSetT     = SignalNodeSet<ScalarT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<OvercurrentRelay, OvercurrentRelayData>;
        using InternalVariablesT = OvercurrentRelayInternalVariables;
        using ExternalVariablesT = OvercurrentRelayExternalVariables;

        OvercurrentRelay();
        explicit OvercurrentRelay(const ModelDataT& data);
        ~OvercurrentRelay();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

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

        /// Trip above the x = 1/2 latch commit, so the relay may read the current it interrupts
        static constexpr RealT TRIP_LEVEL = THREE<RealT> * QUARTER<RealT>;

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

        RealT Ipickup_{ONE<RealT>};
        RealT Ttrip_{ONE<RealT>};
        RealT inv_Ipickup2_{ONE<RealT>};
        RealT Tlatch_{ONE<RealT>};

        IdxT parameter_error_count_{0};

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
