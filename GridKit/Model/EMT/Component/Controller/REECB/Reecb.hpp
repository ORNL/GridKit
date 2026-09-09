/**
 * @file Reecb.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Renewable electrical controller.
 */
#pragma once

#include <array>
#include <limits>
#include <memory>
#include <optional>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class ReecbInternalVariables : size_t
      {
        VMEAS,  ///< \f$V^\mathrm{meas}\f$ Differential filtered terminal voltage [p.u.]
        PMEAS,  ///< \f$P^\mathrm{meas}\f$ Differential filtered electrical power on component base [p.u.]
        XPIQ,   ///< \f$x_Q^\mathrm{PI}\f$ Differential reactive-power PI state [p.u.]
        XPIV,   ///< \f$x_V^\mathrm{PI}\f$ Differential voltage-control PI state on component base [p.u.]
        QV,     ///< \f$Q_V\f$ Differential reactive-current command lag state on component base [p.u.]
        PORD,   ///< \f$P^\mathrm{ord}\f$ Differential filtered active-power order on component base [p.u.]
        VT,     ///< \f$V_t\f$ Algebraic terminal-voltage magnitude [p.u.]
        VSAFE,  ///< \f$V_\mathrm{safe}^\mathrm{meas}\f$ Algebraic safe measured voltage [p.u.]
        SDIP,   ///< \f$s_\mathrm{dip}\f$ Algebraic voltage-band gate [-]
        IQV,    ///< \f$I_q^\mathrm{inj}\f$ Algebraic reactive-current injection on component base [p.u.]
        QREF,   ///< \f$Q_c^{\mathrm{ref}}\f$ Algebraic selected reactive-power reference on component base [p.u.]
        EQ,     ///< \f$e_Q\f$ Algebraic reactive-power error on component base [p.u.]
        VPIQ,   ///< \f$V_Q^\mathrm{PI}\f$ Algebraic reactive-power PI output [p.u.]
        EPIV,   ///< \f$e_V^\mathrm{PI}\f$ Algebraic voltage-control error [p.u.]
        RPORD,  ///< \f$r_P^\mathrm{ord}\f$ Algebraic limited active-power order rate [p.u./s]
        ILCAP,  ///< \f$I_L^\mathrm{cap}\f$ Algebraic off-axis current capacity on component base [p.u.]
        IQMAX,  ///< \f$I_q^\max\f$ Algebraic reactive-current limit on component base [p.u.]
        IPMAX,  ///< \f$I_p^\max\f$ Algebraic active-current limit on component base [p.u.]
        IQBASE, ///< \f$I_q^\mathrm{base}\f$ Algebraic voltage-controller current on component base [p.u.]
        IQRAW,  ///< \f$I_q^\mathrm{raw}\f$ Algebraic pre-limit reactive-current command on component base [p.u.]
        IQCMD,  ///< \f$I_q^\mathrm{cmd}\f$ Algebraic reactive-current command on component base [p.u.]
        IPCMD,  ///< \f$I_p^\mathrm{cmd}\f$ Algebraic active-current command on component base [p.u.]
        ICMDD,  ///< \f$i_d^{\mathrm{cmd}}\f$ Terminal d-axis current command [A]
        ICMDQ,  ///< \f$i_q^{\mathrm{cmd}}\f$ Terminal q-axis current command [A]
        MAXIMUM ///< Number of REECB internal variables and residual rows
      };

      enum class ReecbExternalVariables : size_t
      {
        VD,
        VQ,
        ID,
        IQ,
        PREF,
        QREF,
        VREF,
        PFAREF,
        MAXIMUM
      };

      template <typename scalar_type, typename index_type>
      class Reecb : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::y_ext_;
        using Component<scalar_type, index_type>::yp_ext_;
        using Component<scalar_type, index_type>::variable_indices_ext_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::allocated_;

      public:
        using ScalarT      = scalar_type;
        using IdxT         = index_type;
        using RealT        = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT   = ReecbData<RealT, IdxT>;
        using Outputs      = typename ModelDataT::Outputs;
        using SignalT      = Signal<ScalarT, IdxT>;
        using MonitorT     = Model::VariableMonitor<Reecb, ReecbData>;
        using InputSignals = std::array<SignalT*, 8>;

        Reecb();
        explicit Reecb(const ModelDataT& data);
        ~Reecb();
        void     attachInput(InputSignals inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

        SignalT& inputSignal(ReecbInputs input)
        {
          return *signals_.getAttachedSignal(static_cast<ReecbExternalVariables>(input));
        }

        int                                                     setGridKitComponentID(IdxT) override final;
        int                                                     allocate() override final;
        int                                                     verify() const override final;
        int                                                     initialize(const std::map<Outputs, RealT>& outputs = {});
        int                                                     initializeState(const std::map<std::string, RealT>& values) override;
        void                                                    validateInitialState(const std::map<std::string, RealT>& values) const override;
        typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;
        void                                                    prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial) override;
        int                                                     setAbsoluteTolerance(RealT) override final;
        int                                                     evaluateInternalResidual() override final;
        int                                                     evaluateResidual() override final;
        int                                                     assembleJacobian(RealT y_scale, RealT yp_scale) override final;

        auto& getSignals()
        {
          return signals_;
        }

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void                              initializeParameters(const ModelDataT& data);
        void                              initializeMonitor();
        const Model::VariableMonitorBase* getMonitor() const override;

        struct OperatingPoint
        {
          std::array<RealT, 24> state{};
          std::array<RealT, 4>  reference{};
          RealT                 qmin, qmax, vmin, vmax, pmin, pmax, imax, vref;
        };

        struct InitialCurrentLimit
        {
          RealT total_limit;
          RealT off_axis_capacity;
        };

        /// Smooth asymmetric slew-rate limiter.
        [[gnu::always_inline]] static inline ScalarT aslew(ScalarT rate, RealT lower, RealT upper);

        /// Smooth anti-windup derivative within a moving symmetric band.
        [[gnu::always_inline]] static inline ScalarT awband(ScalarT state, ScalarT rate, ScalarT band);

        /// Smooth nonnegative root of a current-circle squared radius.
        template <typename ValueT>
        [[gnu::always_inline]] static inline ValueT sqrtramp(ValueT x);

        /// Analytic nonnegative seed for inverting sqrtramp().
        static RealT isqrtramp(RealT y);

        /// Overflow-resistant difference of squares.
        template <typename ValueT>
        [[gnu::always_inline]] static inline ValueT circleSquare(RealT limit, ValueT high);

        /// Solve a feasible initial limit at or above `lower`.
        static std::optional<InitialCurrentLimit> solveInitialLimit(RealT lower, RealT high, RealT low);

        OperatingPoint operatingPoint(const std::map<Outputs, RealT>& outputs, const std::array<RealT, 8>& input) const;
        bool           iclamp(RealT output, RealT lower, RealT upper, RealT& input) const;
        static RealT   logOneMinusExp(RealT x);
        bool           referenceActive(size_t n) const;

        static constexpr RealT TIME_CONSTANT_MINIMUM    = RealT{1e-3};
        static constexpr RealT VMEAS_MINIMUM            = RealT{0.01};
        static constexpr RealT INITIALIZATION_TOLERANCE = RealT{100} * std::numeric_limits<RealT>::epsilon();
        RealT                  S_{0.0};
        RealT                  V_{0.0};
        bool                   PfFlag_{false};
        bool                   VFlag_{false};
        bool                   QFlag_{false};
        bool                   Pqflag_{false};
        RealT                  Trv_{0.02};
        RealT                  Tp_{0};
        RealT                  Vref0_{0};
        RealT                  Vdip_{0.85};
        RealT                  Vup_{1.15};
        RealT                  dbd1_{0};
        RealT                  dbd2_{0};
        RealT                  kqv_{5.0};
        RealT                  Iql1_{-1.1};
        RealT                  Iqh1_{1.1};
        RealT                  Qmax_{0.436};
        RealT                  Qmin_{-0.436};
        RealT                  Kqp_{0};
        RealT                  Kqi_{0.1};
        RealT                  Vmax_{1.1};
        RealT                  Vmin_{0.9};
        RealT                  Kvp_{18.0};
        RealT                  Kvi_{5.0};
        RealT                  Tiq_{0.02};
        RealT                  Tpord_{0.02};
        RealT                  dPmax_{99.0};
        RealT                  dPmin_{-99.0};
        RealT                  Pmax_{1};
        RealT                  Pmin_{0};
        RealT                  Imax_{1.3};

        bool Vref0_given_{false};

        RealT pf_on_{0};
        RealT pf_off_{1};
        RealT q_on_{0};
        RealT q_off_{1};
        RealT q_pi_on_{0};
        RealT v_ref_on_{0};
        RealT q_ref_on_{1};
        RealT pq_on_{0};
        RealT pq_off_{1};

        std::array<RealT, 4>                                                            reference_{};
        ComponentSignals<ScalarT, IdxT, ReecbInternalVariables, ReecbExternalVariables> signals_;
        std::array<SignalT, 2>                                                          output_;
        std::array<SignalT*, 2>                                                         alias_{};
        std::unique_ptr<MonitorT>                                                       monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
