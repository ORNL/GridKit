/**
 * @file Regca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Internal variables of a `Regca`
      enum class RegcaInternalVariables : size_t
      {
        VM,      ///< \f$V_M\f$ Filtered terminal voltage
        IQ,      ///< \f$I_q\f$ Reactive-current state on component base
        IP,      ///< \f$I_p\f$ Active-current state on component base
        VT,      ///< \f$V_T\f$ Terminal voltage magnitude
        IR,      ///< \f$I_\mathrm{r}\f$ Branch-current real component on system base
        II,      ///< \f$I_\mathrm{i}\f$ Branch-current imaginary component on system base
        IQEXTRA, ///< \f$I_q^\mathrm{extra}\f$ Extra inductive current from HVRCM on component base
        IL,      ///< \f$I_L\f$ LVPL upper-limit current curve on component base
        PBR,     ///< \f$P^\mathrm{br}\f$ Branch active power on system base
        QBR,     ///< \f$Q^\mathrm{br}\f$ Branch reactive power on system base
        MAXIMUM,
      };

      /// External variables of a `Regca`
      enum class RegcaExternalVariables : size_t
      {
        IPCMD, ///< \f$I_p^\mathrm{cmd}\f$ Active-current command on system base
        IQCMD, ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command on system base
        MAXIMUM,
      };

      /**
       * @brief First-generation WECC renewable generator/converter model
       *        (REGCA) for inverter-coupled resources.
       *
       * Tracks the active- and reactive-current commands through the
       * converter current-control lag and the REGCA limit logic (LVPL,
       * LVACM, HVRCM, and the recovery rate limits) and injects the
       * resulting branch current into the terminal bus. The model README
       * derives the equations and initialization realized here.
       *
       * Internal current states and limiter quantities are on component
       * base. Signal ports, monitor outputs, branch currents, and branch
       * powers are on system base.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       *
       * @see RegcaData for the parameter keys, ports, and monitor selections.
       */
      template <typename scalar_type, typename index_type>
      class Regca : public Component<scalar_type, index_type>
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
        using Component<scalar_type, index_type>::va_system_base_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using BusT               = BusBase<ScalarT, IdxT>;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = RegcaData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Regca, RegcaData>;
        using InternalVariablesT = RegcaInternalVariables;
        using ExternalVariablesT = RegcaExternalVariables;

        /**
         * @brief Construct a REGCA converter without model data.
         *
         * @param[in] bus Terminal bus the converter injects into.
         */
        Regca(BusT* bus);

        /**
         * @brief Construct a REGCA converter from model data.
         *
         * @param[in] bus Terminal bus the converter injects into.
         * @param[in] data Parameters and monitored-variable selections.
         *
         * Model-data errors are recorded for verify() rather than thrown.
         */
        Regca(BusT* bus, const ModelDataT& data);
        ~Regca();

        int setGridKitComponentID(IdxT component_id) override final;

        /// Allocate model storage and connect assigned output signals.
        int allocate() override final;

        /**
         * @brief Validate the REGCA configuration.
         *
         * Applies the parameter conditions in the README, requires a
         * terminal bus, and requires attached command ports to be linked.
         * Operating-point admissibility is checked by initialize().
         *
         * @return Number of configuration errors, zero when valid.
         */
        int verify() const override final;

        /**
         * @brief Initialize REGCA from the power-flow operating point.
         *
         * Resolves the internal states and current-command setpoints from
         * the initialized terminal-bus voltage and the system-base P0/Q0
         * injections, and publishes the resolved commands to attached
         * command ports.
         *
         * @pre allocate() has completed, verify() reports no errors, and
         *      the terminal bus has been initialized.
         * @pre \f$V_{A1} \le V_{T,0} < V_\mathrm{hv}^{\max}\f$, and with
         *      LVPL enabled \f$I_{p,0} \le I_{L,0}\f$.
         * @post All internal derivatives are zero, and unattached command
         *       ports latch the resolved setpoints as constant commands.
         * @return 0 on success, nonzero when the operating point is
         *         rejected.
         */
        int initialize() override final;

        /// Tag VM, IQ, and IP as differential variables.
        int tagDifferentiable() override final;

        /**
         * @brief Set every REGCA absolute tolerance to @p rel_tol.
         *
         * @param[in] rel_tol Solver relative tolerance.
         */
        int setAbsoluteTolerance(RealT rel_tol) override final;

        /**
         * @brief Evaluate the residuals and accumulate the branch current
         *        into the terminal bus.
         *
         * @pre The terminal-bus residual has been zeroed this evaluation.
         */
        int evaluateResidual() override final;

        /**
         * @brief Assemble the sparse component Jacobian when built with Enzyme.
         *
         * Plain and dependency-tracking instantiations do not assemble a
         * separate Jacobian.
         *
         * @pre evaluateResidual() has run at the current state.
         */
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                RegcaInternalVariables,
                                RegcaExternalVariables>&
        {
          return signals_;
        }

        /**
         * @brief Access the monitor.
         * @return Monitor for this model, or nullptr when the model was
         *         constructed without data.
         */
        const Model::VariableMonitorBase* getMonitor() const override;

        /**
         * @brief Internal residual: the README differential and algebraic
         *        equations in RegcaInternalVariables order.
         *
         * @param[in] y Internal variables.
         * @param[in] yp Internal variable derivatives.
         * @param[in] wb Terminal-bus voltage components.
         * @param[in] ws Current-command signal values on system base.
         * @param[out] f Internal residuals.
         */
        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* ws, ScalarT* f);

        /**
         * @brief Bus residual: the system-base branch current injected
         *        into the terminal bus.
         *
         * @param[in] y Internal variables.
         * @param[in] yp Internal variable derivatives, unused.
         * @param[in] wb Terminal-bus voltage components, unused.
         * @param[out] h Current injected into the terminal bus.
         */
        __attribute__((always_inline)) inline int evaluateBusResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        /// Convert a system-base per-unit value to the component base.
        ScalarT toComponentBase(ScalarT value) const
        {
          return value * va_system_base_ / va_converter_base_;
        }

        /// Convert a component-base per-unit value to the system base.
        ScalarT toSystemBase(ScalarT value) const
        {
          return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
        }

        /**
         * @brief Smooth approximation of the REGCA `rrpwr` rate limiter.
         *
         * Limits motion that increases \f$|x|\f$ while smoothly releasing
         * restoring motion toward zero. At \f$x = 0\f$, the result is the
         * symmetric smooth slew limit applied to @p f.
         *
         * @param[in] x State whose sign determines the limited direction.
         * @param[in] f Unconstrained derivative.
         * @param[in] rate Nonnegative symmetric rate limit.
         * @return Rate-limited derivative.
         *
         * @todo Move this reusable limiter to CommonMath.
         */
        static __attribute__((always_inline)) inline ScalarT rrpwr(
            const ScalarT x,
            const ScalarT f,
            const RealT   rate)
        {
          assert(rate >= ZERO<RealT>);

          const ScalarT t     = std::tanh(HALF<RealT> * Math::MU<RealT> * x);
          const ScalarT abs_t = std::abs(t);
          const ScalarT w_pos = HALF<RealT> * t * (t + abs_t);
          const ScalarT w_neg = HALF<RealT> * t * (t - abs_t);

          return f + (ONE<RealT> - w_pos) * Math::ramp(-f - rate)
                 - (ONE<RealT> - w_neg) * Math::ramp(f - rate);
        }

        /**
         * @brief Smooth upper-limit anti-windup derivative.
         *
         * Passes dynamics below the algebraic upper bound. At or above the
         * bound, negative restoring motion passes while positive motion that
         * increases the violation is smoothly blocked.
         *
         * @param[in] x Limited state.
         * @param[in] f Unconstrained derivative.
         * @param[in] limit_max Algebraic upper bound.
         * @return Anti-windup-limited derivative.
         *
         * @todo Move this one-sided anti-windup helper to CommonMath.
         */
        static __attribute__((always_inline)) inline ScalarT awmax(
            const ScalarT x,
            const ScalarT f,
            const ScalarT limit_max)
        {
          const ScalarT below = Math::sigmoid(limit_max - x);

          return (below
                  + (ONE<RealT> - below) * Math::sigmoid(-f))
                 * f;
        }

        // Nonnegative root of q = ramp(q - margin). The root diverges as the
        // positive margin approaches zero.
        ScalarT smoothConstraintCorrection(ScalarT margin) const;

        // Terminal-bus accessors. Ir() and Ii() are accumulation targets,
        // not assignment targets.
        ScalarT& Vr()
        {
          return bus_->Vr();
        }

        ScalarT& Vi()
        {
          return bus_->Vi();
        }

        ScalarT& Ir()
        {
          return bus_->Ir();
        }

        ScalarT& Ii()
        {
          return bus_->Ii();
        }

        // Well-posedness floor for the current-control and voltage-sensor lags
        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        BusT* bus_{nullptr}; ///< Terminal bus the converter injects into

        // Input parameters
        RealT p0_{0};       ///< Initial active power injection on system base
        RealT q0_{0};       ///< Initial reactive power injection on system base
        RealT mva_base_{0}; ///< REGCA component power base
        RealT Tg_{0};       ///< Converter current-control lag time constant
        RealT TM_{0};       ///< Terminal voltage sensor time constant
        RealT Rqmax_{0};    ///< Reactive-current recovery positive rate limit
        RealT Rqmin_{0};    ///< Reactive-current recovery negative rate limit
        RealT Rpmax_{0};    ///< Active-current magnitude recovery rate limit
        bool  sL_{false};   ///< LVPL switch
        RealT IL1_{0};      ///< LVPL upper-current ceiling
        RealT VL0_{0};      ///< LVPL zero-crossing voltage
        RealT VL1_{0};      ///< LVPL upper breakpoint voltage
        RealT VA0_{0};      ///< LVACM lower breakpoint voltage
        RealT VA1_{0};      ///< LVACM upper breakpoint voltage
        RealT Vhvmax_{0};   ///< Terminal-voltage ceiling for HV reactive management

        IdxT parameter_error_count_{0}; ///< Data errors counted for verify()

        // Derived parameters. The complementary runtime LVPL masks keep a
        // configuration-independent Jacobian sparsity pattern: both AD paths
        // retain the IL expression even when its contribution is zero.
        RealT va_converter_base_{0}; ///< Component power base in VA
        RealT use_lvpl_{0};          ///< LVPL mask, complements bypass_lvpl_
        RealT bypass_lvpl_{1};       ///< LVPL bypass mask, complements use_lvpl_
        RealT iq_use_upper_{0};      ///< Rqmax-branch mask, complements iq_use_lower_
        RealT iq_use_lower_{1};      ///< Rqmin-branch mask, complements iq_use_upper_

        // Command setpoints latched by initialize(), used when the matching
        // signal port is unattached
        ScalarT ipcmd_set_{0}; ///< Active-current command setpoint on system base
        ScalarT iqcmd_set_{0}; ///< Reactive-current command setpoint on system base

        ComponentSignals<ScalarT, IdxT, RegcaInternalVariables, RegcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        // Local copies of signal variables
        std::vector<ScalarT> ws_;         ///< Command signal values on system base
        std::vector<IdxT>    ws_indices_; ///< Global indices of attached command signals
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
