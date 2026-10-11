#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Load/PERC1/Perc1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Differential variables precede algebraic variables, in README order.
    enum class Perc1InternalVariables : size_t
    {
      XPLAG,  ///< Active-power lead-lag state
      XQLAG,  ///< Reactive-power lead-lag state
      XPWO,   ///< Active-path voltage-washout state
      XQWO,   ///< Reactive-path voltage-washout state
      DFMEAS, ///< Filtered frequency deviation
      VMEAS,  ///< Filtered voltage magnitude
      IP,     ///< Active-current output on component base
      IQ,     ///< Reactive-current output on component base
      XD,     ///< Committed-delay memory
      XC,     ///< Ceased-load memory
      XR,     ///< Reconnection-release memory
      XH,     ///< First-cessation memory
      TC,     ///< Undervoltage qualification clock
      TD,     ///< Committed-delay clock
      TR,     ///< Recovery qualification clock
      TA,     ///< Reconnection ramp clock
      V,      ///< Regularized terminal-voltage magnitude
      VSAFE,  ///< Safe filtered voltage magnitude
      DP,     ///< Active-power washout output
      DQ,     ///< Reactive-power washout output
      UP,     ///< Active-power lead-lag input
      UQ,     ///< Reactive-power lead-lag input
      PCTRL,  ///< Controlled active power
      QCTRL,  ///< Controlled reactive power
      IPV,    ///< Limited active-current target
      IQV,    ///< Limited controlled reactive-current target
      FRACON, ///< Connected-load multiplier
      IR,     ///< Terminal real current on system base, positive for injection
      II,     ///< Terminal imaginary current on system base, positive for injection
    };

    enum class Perc1ExternalVariables : size_t
    {
      FREQ, ///< Absolute per-unit frequency
    };

    /**
     * @brief Power-electronic load with voltage/frequency response and reconnection.
     *
     * Internal currents are positive for consumption on component base.
     * Terminal currents are positive for injection on system base.
     */
    template <typename scalar_type, typename index_type>
    class Perc1 : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::va_system_base_;
      using Component<scalar_type, index_type>::va_component_base_;
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::ws_;
      using Component<scalar_type, index_type>::ws_indices_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;

    public:
      using ScalarT            = scalar_type;
      using IdxT               = index_type;
      using RealT              = typename Component<ScalarT, IdxT>::RealT;
      using BusT               = BusBase<ScalarT, IdxT>;
      using ModelDataT         = Perc1Data<RealT, IdxT>;
      using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
      using MonitorT           = Model::VariableMonitor<Perc1, Perc1Data>;
      using InternalVariablesT = Perc1InternalVariables;
      using ExternalVariablesT = Perc1ExternalVariables;

      explicit Perc1(BusT* bus);
      Perc1(BusT* bus, const ModelDataT& data);
      ~Perc1();

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
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* ws, ScalarT* f);
      __attribute__((always_inline)) inline int evaluateBusResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParameters();

      /// Reset-priority smooth memory rate; see README Appendix A.
      static __attribute__((always_inline)) inline ScalarT latch(
          ScalarT x, ScalarT set, ScalarT reset)
      {
        return (ONE<RealT> - reset)
                   * (set + (ONE<RealT> - set) * Math::sigmoid(x - HALF<RealT>))
               - x;
      }

      /// Inverse CommonMath clamp for initialization; NaN means infeasible.
      static RealT inverseClamp(RealT value, RealT lower, RealT upper);

      static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
      static constexpr RealT VOLTAGE_EPSILON       = static_cast<RealT>(1.0e-6);
      static constexpr RealT VOLTAGE_MINIMUM       = static_cast<RealT>(0.01);
      static constexpr RealT TAU                   = ONE<RealT> / (Math::MU<RealT> * Math::MU<RealT>);
      static_assert(Math::MU<RealT> > static_cast<RealT>(4), "PERC1 latch requires MU > 4");

      BusT* bus_{nullptr};

      // Input parameters, with the defaults documented in the README.
      RealT Pnom_{0};
      RealT Qnom_{0};
      RealT Lfm_{0.80};
      RealT QPratio_{0.66};
      RealT Dbfl_{0};
      RealT Dbfh_{0};
      RealT Kdroop_{0};
      RealT Kvp_{0};
      RealT Tvp_{0.10};
      RealT Kvq_{0};
      RealT Tvq_{0.10};
      RealT Tap_{0};
      RealT Tbp_{0};
      RealT Taq_{0};
      RealT Tbq_{0};
      RealT nP_{0};
      RealT nQ_{1};
      RealT Ipmax_{1};
      RealT Ipmin_{0};
      RealT Iqmax_{0.66};
      RealT Iqmin_{-0.66};
      RealT Fcease_{1};
      RealT Vcease_{0.50};
      RealT Tcease_{0.01};
      RealT Tdelay_{0};
      RealT Vrecon_{0.60};
      RealT Trecon_{0.05};
      RealT Tramp_{1};
      RealT Frecon_{1};
      RealT Tt_{0.02};
      RealT Tv_{0.02};
      RealT Tf_{0.02};

      int parameter_error_count_{0};

      // Fixed at initialization. Static loads retain the same vector layout.
      RealT dynamic_on_{1};
      RealT static_on_{0};
      RealT timed_cease_{0};
      RealT timed_delay_{0};
      RealT timed_recon_{0};
      RealT timed_ramp_{0};
      RealT inv_Tcease_{0};
      RealT inv_Tdelay_{0};
      RealT inv_Trecon_{0};
      RealT inv_Tramp_{0};
      RealT G_{0};
      RealT B_{0};
      RealT Pref_{0};
      RealT Qref_{0};
      RealT Iqextra_{0};
      RealT Vnom_safe_{1};
      RealT cease_on_{1};

      SignalPortsT              ports_;
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
