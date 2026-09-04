/**
 * @file Reecb.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REECB electrical-control model.
 */

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    namespace Controller
    {
      /// Internal variables and residual rows of a `Reecb`.
      enum class ReecbInternalVariables : size_t
      {
        VMEAS,  ///< \f$V^\mathrm{meas}\f$ Differential filtered terminal voltage [p.u.]
        PMEAS,  ///< \f$P^\mathrm{meas}\f$ Differential filtered electrical power on component base [p.u.]
        XPIQ,   ///< \f$x_Q^\mathrm{PI}\f$ Differential reactive-power PI state [p.u.]
        XPIV,   ///< \f$x_V^\mathrm{PI}\f$ Differential voltage-control PI state on component base [p.u.]
        QV,     ///< \f$Q_V\f$ Differential reactive-current command lag state on component base [p.u.]
        PORD,   ///< \f$P^\mathrm{ord}\f$ Differential filtered active-power order on component base [p.u.]
        VT,     ///< \f$V_T\f$ Algebraic terminal-voltage magnitude [p.u.]
        VSAFE,  ///< \f$V_\mathrm{safe}^\mathrm{meas}\f$ Algebraic safe measured voltage [p.u.]
        SDIP,   ///< \f$s_\mathrm{dip}\f$ Algebraic voltage-band gate [-]
        IQV,    ///< \f$I_q^\mathrm{inj}\f$ Algebraic reactive-current injection on component base [p.u.]
        QREF,   ///< \f$Q^\mathrm{ref}\f$ Algebraic reactive-power reference on component base [p.u.]
        EQ,     ///< \f$e_Q\f$ Algebraic reactive-power error on component base [p.u.]
        VPIQ,   ///< \f$V_Q^\mathrm{PI}\f$ Algebraic reactive-power PI output [p.u.]
        EPIV,   ///< \f$e_V^\mathrm{PI}\f$ Algebraic voltage-control error [p.u.]
        RPORD,  ///< \f$r_P^\mathrm{ord}\f$ Algebraic limited active-power order rate [p.u./s]
        ILMAX,  ///< \f$I_L^\max\f$ Algebraic current-circle continuation state on component base [p.u.]
        ILCAP,  ///< \f$I_L^\mathrm{cap}\f$ Algebraic off-axis current capacity on component base [p.u.]
        IQMAX,  ///< \f$I_q^\max\f$ Algebraic reactive-current limit on component base [p.u.]
        IPMAX,  ///< \f$I_p^\max\f$ Algebraic active-current limit on component base [p.u.]
        IQBASE, ///< \f$I_q^\mathrm{base}\f$ Algebraic voltage-controller current on component base [p.u.]
        IQRAW,  ///< \f$I_q^\mathrm{raw}\f$ Algebraic pre-limit reactive-current command on component base [p.u.]
        IQCMD,  ///< \f$I_q^\mathrm{cmd}\f$ Algebraic reactive-current command output on system base [p.u.]
        IPCMD,  ///< \f$I_p^\mathrm{cmd}\f$ Algebraic active-current command output on system base [p.u.]
        MAXIMUM ///< Number of REECB internal variables and residual rows
      };

      /// External variables read by a `Reecb`.
      enum class ReecbExternalVariables : size_t
      {
        VR,     ///< \f$V_\mathrm{r}\f$ Terminal-bus real voltage [p.u.]
        VI,     ///< \f$V_\mathrm{i}\f$ Terminal-bus imaginary voltage [p.u.]
        PE,     ///< \f$P_e\f$ Optional Known active-power feedback input on system base [p.u.]
        QGEN,   ///< \f$Q^\mathrm{gen}\f$ Optional Known reactive-power feedback input on system base [p.u.]
        QEXT,   ///< \f$Q^\mathrm{ext}\f$ Optional Unknown Volt/VAr reference input: system-base reactive power [p.u.], or the terminal-voltage reference [p.u.] when \f$s_Q=1\f$ and \f$s_V=0\f$
        PFAREF, ///< \f$\phi^\mathrm{ref}\f$ Optional Unknown power-factor angle-reference input [rad]
        PREF,   ///< \f$P^\mathrm{ref}\f$ Optional Unknown active-power reference input on system base [p.u.]
        MAXIMUM ///< Number of REECB external variables
      };

      /**
       * @brief WECC renewable electrical controller with reactive-power,
       *        voltage, and active-power command paths (REECB).
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class Reecb : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::va_system_base_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::variable_indices_ext_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::y_ext_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using BusT               = BusBase<ScalarT, IdxT>;
        using ModelDataT         = ReecbData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Reecb, ReecbData>;
        using InternalVariablesT = ReecbInternalVariables;
        using ExternalVariablesT = ReecbExternalVariables;

        /// Current-circle regularization and initialization reconstruction tolerance.
        static constexpr RealT INITIALIZATION_TOLERANCE =
            static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

        Reecb(BusT* bus);
        Reecb(BusT* bus, const ModelDataT& data);
        ~Reecb();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateInternalResidual() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                ReecbInternalVariables,
                                ReecbExternalVariables>&;

        const Model::VariableMonitorBase* getMonitor() const override;

        [[gnu::always_inline]] inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* y_ext,
            ScalarT*       f);

      private:
        void gatherExternalVariables();

        struct InitialPoint;

        struct InitialCurrentLimit
        {
          RealT total_limit;
          RealT continuation;
          RealT off_axis_capacity;
        };

        /// Smooth asymmetric slew-rate limiter.
        [[gnu::always_inline]] static inline ScalarT aslew(ScalarT rate, RealT lower, RealT upper);

        /// Smooth anti-windup derivative within a moving symmetric band.
        [[gnu::always_inline]] static inline ScalarT awband(ScalarT state, ScalarT rate, ScalarT band);

        /// Current-circle continuation state for an initial component-base limit.
        static RealT circleState(RealT imax, RealT high);

        /// Off-axis component-base capacity provided by a continuation state.
        static RealT capacity(RealT ilmax);

        /// Bisect an initial-limit bracket to its first upper-side point.
        template <typename FuncT>
        static RealT bisect(RealT a, RealT b, FuncT below);

        /// Solve the smallest feasible initial limit at or above `lower`.
        static std::optional<InitialCurrentLimit> solveInitialLimit(RealT lower, RealT high, RealT low);

        bool buildInitialPoint(InitialPoint& point);
        void commitInitialPoint(const InitialPoint& point);

        void loadRealParameter(const ModelDataT& data,
                               ReecbParameters   parameter,
                               RealT&            target,
                               const char*       name);
        void loadBooleanParameter(const ModelDataT& data,
                                  ReecbParameters   parameter,
                                  bool&             target,
                                  const char*       name);
        bool floorTimeConstant(RealT& value, const char* name);
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        static RealT logOneMinusExp(RealT x);
        bool         iclamp(RealT output, RealT lower, RealT upper, RealT& input) const;
        RealT        componentPowerBase() const;

        template <typename ValueT>
        [[gnu::always_inline]] inline ValueT toComponentBase(ValueT value) const;

        template <typename ValueT>
        ValueT toSystemBase(ValueT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

        static constexpr RealT VMEAS_MINIMUM = static_cast<RealT>(0.01);

        BusT* bus_{nullptr};

        // Input parameters
        RealT mva_base_{0};
        bool  PfFlag_{false};
        bool  VFlag_{false};
        bool  QFlag_{false};
        bool  Pqflag_{false};
        RealT Trv_{0.02};
        RealT Tp_{0};
        RealT Vref0_{0};
        RealT Vdip_{0.85};
        RealT Vup_{1.15};
        RealT dbd1_{0};
        RealT dbd2_{0};
        RealT kqv_{5.0};
        RealT Iql1_{-1.1};
        RealT Iqh1_{1.1};
        RealT Qmax_{0.436};
        RealT Qmin_{-0.436};
        RealT Kqp_{0};
        RealT Kqi_{0.1};
        RealT Vmax_{1.1};
        RealT Vmin_{0.9};
        RealT Kvp_{18.0};
        RealT Kvi_{5.0};
        RealT Tiq_{0.02};
        RealT Tpord_{0.02};
        RealT dPmax_{99.0};
        RealT dPmin_{-99.0};
        RealT Pmax_{1};
        RealT Pmin_{0};
        RealT Imax_{1.3};

        bool mva_given_{false};
        bool Vref0_given_{false};
        IdxT parameter_error_count_{0};

        // Derived parameters
        RealT va_component_base_{0};
        RealT pf_on_{0};
        RealT pf_off_{1};
        RealT q_on_{0};
        RealT q_off_{1};
        RealT q_pi_on_{0};
        RealT v_ref_on_{0};
        RealT q_ref_on_{1};
        RealT pq_on_{0};
        RealT pq_off_{1};

        // Unattached signal setpoints
        ScalarT pe_set_{0};
        ScalarT qgen_set_{0};
        ScalarT qext_set_{0};
        ScalarT pfaref_set_{0};
        ScalarT pref_set_{0};

        ComponentSignals<ScalarT, IdxT, ReecbInternalVariables, ReecbExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
