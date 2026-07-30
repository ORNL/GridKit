/**
 * @file Reecb.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REECB electrical-control model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/ReecbData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace Converter
    {
      /// Internal variables of a `Reecb`
      enum class ReecbInternalVariables : size_t
      {
        VMEAS,     ///< \f$V^\mathrm{meas}\f$ Filtered terminal voltage
        PMEAS,     ///< \f$P^\mathrm{meas}\f$ Filtered active power on component base
        XPIQ,      ///< \f$x_Q^\mathrm{PI}\f$ Reactive-power PI state
        XPIV,      ///< \f$x_V^\mathrm{PI}\f$ Voltage-control PI state
        QV,        ///< \f$Q_V\f$ Reactive-current command lag state
        PORD,      ///< \f$P^\mathrm{ord}\f$ Filtered active-power order
        VT,        ///< \f$V_T\f$ Terminal voltage magnitude
        VMEASSAFE, ///< \f$V_\mathrm{safe}^\mathrm{meas}\f$ Safe filtered voltage
        SDIP,      ///< \f$s_\mathrm{dip}\f$ Voltage inside-band control gate
        VERR,      ///< \f$e_V^\mathrm{db}\f$ Deadbanded voltage error
        IQV,       ///< \f$I_q^\mathrm{inj}\f$ Reactive-current injection candidate
        QREF,      ///< \f$Q^\mathrm{ref}\f$ Selected reactive-power reference
        EQ,        ///< \f$e_Q\f$ Reactive-power control error
        VPIQ,      ///< \f$V_Q^\mathrm{PI}\f$ Reactive-power PI output
        EPIV,      ///< \f$e_V^\mathrm{PI}\f$ Voltage-control PI error
        FPORD,     ///< \f$f_P^\mathrm{ord}\f$ Active-power derivative target
        RPORD,     ///< \f$r_P^\mathrm{ord}\f$ Ramp-limited active-power derivative
        IQCIRC,    ///< \f$I_q^\mathrm{circ}\f$ Reactive-current circle limit
        IPCIRC,    ///< \f$I_p^\mathrm{circ}\f$ Active-current circle limit
        IQMAX,     ///< \f$I_q^\max\f$ Reactive-current upper limit
        IPMAX,     ///< \f$I_p^\max\f$ Active-current upper limit
        IQBASE,    ///< \f$I_q^\mathrm{base}\f$ Base reactive-current command
        IQRAW,     ///< \f$I_q^\mathrm{raw}\f$ Raw reactive-current command
        IQCMD,     ///< \f$I_q^\mathrm{cmd}\f$ Command output on system base
        IPCMD,     ///< \f$I_p^\mathrm{cmd}\f$ Command output on system base
        MAXIMUM,
      };

      /// External variables of a `Reecb`
      enum class ReecbExternalVariables : size_t
      {
        PE,     ///< \f$P_e\f$ Active-power feedback on system base
        QGEN,   ///< \f$Q^\mathrm{gen}\f$ Reactive-power feedback on system base
        QEXT,   ///< \f$Q^\mathrm{ext}\f$ Reactive-power command on system base
        PFAREF, ///< \f$\phi^\mathrm{ref}\f$ Power-factor angle reference in radians
        PREF,   ///< \f$P^\mathrm{ref}\f$ Active-power reference on system base
        MAXIMUM,
      };

      /**
       * @brief Second-generation WECC renewable electrical-control model (REECB).
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
        using SignalT            = SignalNode<ScalarT, IdxT>;
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
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        template <typename LowerT, typename UpperT>
        bool solveLimiterInput(ScalarT requested_output, LowerT lower_limit, UpperT upper_limit, ScalarT& limiter_input) const;

        template <typename LowerT, typename UpperT>
        ScalarT steadyAntiWindupInput(ScalarT nominal_input, ScalarT rate, LowerT lower_limit, UpperT upper_limit) const;

        RealT logOneMinusExp(RealT x) const;

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM       = static_cast<RealT>(1.0e-3);
        static constexpr RealT VMEAS_MINIMUM               = static_cast<RealT>(0.01);
        static constexpr RealT INITIALIZATION_TOLERANCE    = static_cast<RealT>(1.0e-10);
        static constexpr RealT INITIALIZATION_LIMIT_OFFSET = static_cast<RealT>(0.1);

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

        bool Vref0_given_{false};
        IdxT parameter_error_count_{0};

        // Derived parameters
        RealT va_converter_base_{0};
        RealT pf_on_{0};
        RealT pf_off_{1};
        RealT v_on_{0};
        RealT v_off_{1};
        RealT q_on_{0};
        RealT q_off_{1};
        RealT p_priority_{0};
        RealT q_priority_{1};

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
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
