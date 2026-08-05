/**
 * @file Reecb.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REECB electrical-control model.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

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
        ILMAX,  ///< \f$I_L^\max\f$ Algebraic current limit available to the low-priority command on component base [p.u.]
        IQCMD,  ///< \f$I_q^\mathrm{cmd}\f$ Algebraic reactive-current command output on system base [p.u.]
        IPCMD,  ///< \f$I_p^\mathrm{cmd}\f$ Algebraic active-current command output on system base [p.u.]
        MAXIMUM ///< Number of REECB internal variables and residual rows
      };

      /// External signal variables read or initialized by a `Reecb`.
      enum class ReecbExternalVariables : size_t
      {
        PE,     ///< \f$P_e\f$ Optional Known active-power feedback input on system base [p.u.]
        QGEN,   ///< \f$Q^\mathrm{gen}\f$ Optional Known reactive-power feedback input on system base [p.u.]
        QEXT,   ///< \f$Q^\mathrm{ext}\f$ Optional Unknown Volt/VAr reference input: system-base reactive power [p.u.], or the terminal-voltage reference [p.u.] when \f$s_Q=1\f$ and \f$s_V=0\f$
        PFAREF, ///< \f$\phi^\mathrm{ref}\f$ Optional Unknown power-factor angle-reference input [rad]
        PREF,   ///< \f$P^\mathrm{ref}\f$ Optional Unknown active-power reference input on system base [p.u.]
        MAXIMUM ///< Number of REECB external signal variables
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
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
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

        Reecb(BusT* bus);
        Reecb(BusT* bus, const ModelDataT& data);
        ~Reecb();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                ReecbInternalVariables,
                                ReecbExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
        static constexpr size_t index(ReecbInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        }

        static constexpr size_t index(ReecbExternalVariables variable)
        {
          return static_cast<size_t>(variable);
        }

        static void checkConfiguration(bool condition, const char* message, int& errors);
        void        loadRealParameter(const ModelDataT& data,
                                      ReecbParameters   parameter,
                                      RealT&            target,
                                      const char*       name);
        void        loadBooleanParameter(const ModelDataT& data,
                                         ReecbParameters   parameter,
                                         bool&             target,
                                         const char*       name);
        bool        floorTimeConstant(RealT& value, const char* name);
        void        initializeParameters(const ModelDataT& data);
        void        initializeMonitor();
        void        setDerivedParameters();

        RealT logOneMinusExp(RealT x) const;
        RealT unclamp(RealT output, RealT lower, RealT upper) const;
        RealT componentPowerBase() const;

        template <typename ValueT>
        __attribute__((always_inline)) inline ValueT toComponentBase(ValueT value) const;

        template <typename ValueT>
        ValueT toSystemBase(ValueT value) const;

        /**
         * @brief Smooth asymmetric slew-rate limiter.
         *
         * @param[in] f Unconstrained rate.
         * @param[in] rate_min Negative rate limit.
         * @param[in] rate_max Positive rate limit.
         * @return Limited rate.
         */
        static __attribute__((always_inline)) inline ScalarT aslew(
            const ScalarT f,
            const RealT   rate_min,
            const RealT   rate_max)
        {
          assert(rate_min < ZERO<RealT> && ZERO<RealT> < rate_max);

          return f
                 / (ONE<RealT>
                    + Math::ramp(f / rate_max - ONE<RealT>)
                    + Math::ramp(f / rate_min - ONE<RealT>));
        }

        /**
         * @brief Smooth anti-windup derivative within a moving symmetric band.
         *
         * Math::antiwindup over [-band, band] with a band edge that is an
         * algebraic quantity, so differentiation carries the band's own
         * contributions through the gate.
         *
         * @param[in] x Limited PI state.
         * @param[in] f Pre-limit derivative of x.
         * @param[in] band Nonnegative symmetric band edge.
         * @return Anti-windup-limited derivative.
         *
         * @todo Fold moving-limit support into Math::antiwindup in CommonMath.
         */
        static __attribute__((always_inline)) inline ScalarT awband(
            const ScalarT x,
            const ScalarT f,
            const ScalarT band)
        {
          const ScalarT above_min = Math::sigmoid(x + band);
          const ScalarT below_max = Math::sigmoid(band - x);

          return (above_min * below_max                          //
                  + (ONE<RealT> - below_max) * Math::sigmoid(-f) //
                  + (ONE<RealT> - above_min) * Math::sigmoid(f))
                 * f;
        }

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static constexpr RealT VMEAS_MINIMUM         = static_cast<RealT>(0.01);

        /// Accepted reconstruction error where an initialization reference
        /// round-trips through a transcendental function.
        static constexpr RealT INITIALIZATION_TOLERANCE =
            static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

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

        // Local copies of signal variables
        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
