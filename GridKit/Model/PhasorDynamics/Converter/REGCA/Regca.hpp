/**
 * @file Regca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

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
        LP,      ///< \f$\ell_p\f$ Active-current lower rate bound on component base
        UP,      ///< \f$u_p\f$ Active-current upper rate bound on component base
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
         * @brief Construct a sized but unconfigured REGCA converter.
         *
         * @param[in] bus Terminal bus the converter injects into.
         * @post Every parameter keeps its zero default and no monitor is
         *       created, so verify() reports configuration errors until
         *       the data constructor is used instead.
         */
        Regca(BusT* bus);

        /**
         * @brief Construct a REGCA converter from model data.
         *
         * @param[in] bus Terminal bus the converter injects into.
         * @param[in] data Parameters and monitored-variable selections.
         * @post Parameters are loaded and the monitor is created. Data
         *       errors are counted and reported by verify() rather than
         *       thrown.
         */
        Regca(BusT* bus, const ModelDataT& data);
        ~Regca();

        /**
         * @brief Set the component ID assigned by the system model.
         * @return 0 on success.
         */
        int setGridKitComponentID(IdxT component_id) override final;

        /**
         * @brief Allocate the model vectors and wire the output signals.
         *
         * @post State, residual, tolerance, and interface buffers are
         *       sized, the identity index maps are seeded, and each
         *       assigned output signal node aliases the internal state it
         *       publishes.
         * @note Repeated calls reuse the already-allocated vectors.
         * @return 0 on success.
         */
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
         *      LVPL enabled \f$I_{p,0} < I_{L,0}\f$.
         * @post All internal derivatives are zero, and unattached command
         *       ports latch the resolved setpoints as constant commands.
         * @return 0 on success, nonzero when the operating point is
         *         rejected.
         */
        int initialize() override final;

        /**
         * @brief Tag VM, IQ, and IP differential and all others algebraic.
         * @return 0 on success.
         */
        int tagDifferentiable() override final;

        /**
         * @brief Share the relative tolerance as every variable's absolute
         *        floor.
         *
         * @note All REGCA variables are per-unit currents, voltages, or
         *       powers of the same order.
         * @param[in] rel_tol Solver relative tolerance.
         * @return 0 on success.
         */
        int setAbsoluteTolerance(RealT rel_tol) override final;

        /**
         * @brief Evaluate the residuals and accumulate the branch current
         *        into the terminal bus.
         *
         * @pre The terminal-bus residual has been zeroed this evaluation.
         * @return 0 on success.
         */
        int evaluateResidual() override final;

        /**
         * @brief Assemble the sparse component Jacobian.
         *
         * Implemented through Enzyme. The plain and dependency-tracking
         * builds log that no Jacobian is assembled.
         *
         * @pre evaluateResidual() has run at the current state.
         * @return 0 on success.
         */
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Regca`
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
         * @warning The body must stay free of branches and loops so sparse
         *          automatic differentiation resolves one fixed structure.
         * @invariant The limiter selections enter as complementary runtime
         *            masks, so the sparsity pattern is independent of the
         *            sL configuration.
         *
         * @param[in] y Internal variables.
         * @param[in] yp Internal variable derivatives.
         * @param[in] wb Terminal-bus voltage components.
         * @param[in] ws Current-command signal values on system base.
         * @param[out] f Internal residuals.
         * @return 0 on success.
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
         * @return 0 on success.
         */
        __attribute__((always_inline)) inline int evaluateBusResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);

      private:
        /// Load the required parameters, counting errors for verify().
        void initializeParameters(const ModelDataT& data);
        /// Bind the monitorable variables to their internal states.
        void initializeMonitor();
        /// Resolve the derived constants and complementary limiter masks.
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

        /// Active-current lower rate bound \f$L_p(I_p)\f$, diagram `Rdown`.
        ScalarT lpTarget(ScalarT ip) const;
        /// Active-current upper rate bound \f$U_p(I_p)\f$, diagram `Rup`.
        ScalarT upTarget(ScalarT ip) const;
        /**
         * @brief Nonnegative root of \f$q = \mathrm{ramp}(q - m)\f$, the
         *        correction holding a smooth constraint at margin @p margin.
         * @pre margin > 0. The root diverges as the margin approaches zero,
         *      so callers must reject a non-positive margin first.
         */
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

        /// Well-posedness floor for the current-control and voltage-sensor lags
        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        /// Multiple of Rpmax standing in for an inactive active-current rate limit
        static constexpr RealT INACTIVE_RATE_LIMIT_FACTOR = static_cast<RealT>(100.0);

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
        RealT Mp_{0};                ///< \f$M_p\f$ Finite surrogate for inactive rate limits
        RealT va_converter_base_{0}; ///< Component power base in VA
        RealT use_lvpl_{0};          ///< LVPL mask, complements bypass_lvpl_
        RealT bypass_lvpl_{1};       ///< LVPL bypass mask, complements use_lvpl_
        RealT iq_use_upper_{0};      ///< Rqmax-branch mask, complements iq_use_lower_
        RealT iq_use_lower_{1};      ///< Rqmin-branch mask, complements iq_use_upper_

        // Command setpoints latched by initialize(), used when the matching
        // signal port is unattached
        ScalarT ipcmd_set_{0}; ///< Active-current command setpoint on system base
        ScalarT iqcmd_set_{0}; ///< Reactive-current command setpoint on system base

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, RegcaInternalVariables, RegcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        // Local copies of signal variables
        std::vector<ScalarT> ws_;         ///< Command signal values on system base
        std::vector<IdxT>    ws_indices_; ///< Global indices of attached command signals
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
