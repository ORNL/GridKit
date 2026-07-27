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
      /// Internal variables of a `Reecb`.
      enum class ReecbInternalVariables : size_t
      {
        VMEAS,     ///< Filtered terminal voltage
        PMEAS,     ///< Filtered electrical power
        XPIQ,      ///< Reactive-power PI state
        XPIV,      ///< Voltage PI state
        QV,        ///< Reactive-current command lag state
        PORD,      ///< Filtered active-power order
        VT,        ///< Terminal voltage magnitude
        VMEASSAFE, ///< Safe filtered terminal voltage
        SDIP,      ///< Voltage inside-band control gate
        VERR,      ///< Deadbanded voltage error
        IQV,       ///< Reactive-current injection candidate
        QREF,      ///< Selected reactive-power reference
        EQ,        ///< Reactive-power control error
        VPIQ,      ///< Reactive-power PI output
        EPIV,      ///< Voltage-control PI error
        FPORD,     ///< Active-power order derivative before ramp limiting
        RPORD,     ///< Ramp-rate-limited active-power derivative
        IQCIRC,    ///< Reactive-current limit from current circle
        IPCIRC,    ///< Active-current limit from current circle
        IQMAX,     ///< Final reactive-current upper limit
        IPMAX,     ///< Final active-current upper limit
        IQBASE,    ///< Base reactive-current command
        IQRAW,     ///< Raw reactive-current command
        IQCMD,     ///< Reactive-current command output
        IPCMD,     ///< Active-current command output
        MAXIMUM,
      };

      /// External variables of a `Reecb`.
      enum class ReecbExternalVariables : size_t
      {
        PE,     ///< Electrical active-power signal
        QGEN,   ///< Reactive-power signal
        QEXT,   ///< External reactive-power command
        PFAREF, ///< Power-factor angle reference
        PREF,   ///< External active-power reference
        MAXIMUM,
      };

      /// Indices into the REECB state, derivative, and residual vectors.
      struct ReecbIdx
      {
        static constexpr size_t VMEAS     = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        static constexpr size_t PMEAS     = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        static constexpr size_t XPIQ      = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        static constexpr size_t XPIV      = static_cast<size_t>(ReecbInternalVariables::XPIV);
        static constexpr size_t QV        = static_cast<size_t>(ReecbInternalVariables::QV);
        static constexpr size_t PORD      = static_cast<size_t>(ReecbInternalVariables::PORD);
        static constexpr size_t VT        = static_cast<size_t>(ReecbInternalVariables::VT);
        static constexpr size_t VMEASSAFE = static_cast<size_t>(ReecbInternalVariables::VMEASSAFE);
        static constexpr size_t SDIP      = static_cast<size_t>(ReecbInternalVariables::SDIP);
        static constexpr size_t VERR      = static_cast<size_t>(ReecbInternalVariables::VERR);
        static constexpr size_t IQV       = static_cast<size_t>(ReecbInternalVariables::IQV);
        static constexpr size_t QREF      = static_cast<size_t>(ReecbInternalVariables::QREF);
        static constexpr size_t EQ        = static_cast<size_t>(ReecbInternalVariables::EQ);
        static constexpr size_t VPIQ      = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        static constexpr size_t EPIV      = static_cast<size_t>(ReecbInternalVariables::EPIV);
        static constexpr size_t FPORD     = static_cast<size_t>(ReecbInternalVariables::FPORD);
        static constexpr size_t RPORD     = static_cast<size_t>(ReecbInternalVariables::RPORD);
        static constexpr size_t IQCIRC    = static_cast<size_t>(ReecbInternalVariables::IQCIRC);
        static constexpr size_t IPCIRC    = static_cast<size_t>(ReecbInternalVariables::IPCIRC);
        static constexpr size_t IQMAX     = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        static constexpr size_t IPMAX     = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        static constexpr size_t IQBASE    = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        static constexpr size_t IQRAW     = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        static constexpr size_t IQCMD     = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        static constexpr size_t IPCMD     = static_cast<size_t>(ReecbInternalVariables::IPCMD);
        static constexpr size_t MAXIMUM   = static_cast<size_t>(ReecbInternalVariables::MAXIMUM);
      };

      /// Indices into the REECB external-signal buffers.
      struct ReecbExt
      {
        static constexpr size_t PE      = static_cast<size_t>(ReecbExternalVariables::PE);
        static constexpr size_t QGEN    = static_cast<size_t>(ReecbExternalVariables::QGEN);
        static constexpr size_t QEXT    = static_cast<size_t>(ReecbExternalVariables::QEXT);
        static constexpr size_t PFAREF  = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        static constexpr size_t PREF    = static_cast<size_t>(ReecbExternalVariables::PREF);
        static constexpr size_t MAXIMUM = static_cast<size_t>(ReecbExternalVariables::MAXIMUM);
      };

      template <typename scalar_type, typename index_type>
      class Reecb : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
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
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using BusT       = BusBase<ScalarT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using ModelDataT = ReecbData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Reecb, ReecbData>;

        Reecb(BusT* bus);
        Reecb(BusT* bus, const ModelDataT& data);
        ~Reecb();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT) override final;
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
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        /// Solve the input required to produce a requested smooth-limiter output.
        /// The limits may be constant Real parameters or algebraic variables.
        template <typename LowerT, typename UpperT>
        bool solveLimiterInput(ScalarT requested_output, LowerT lower_limit, UpperT upper_limit, ScalarT& limiter_input) const;

        /// Select a limiter input that zeros an anti-windup rate to initialization tolerance.
        /// The limits may be constant Real parameters or algebraic variables.
        template <typename LowerT, typename UpperT>
        ScalarT steadyAntiWindupInput(ScalarT nominal_input, ScalarT rate, LowerT lower_limit, UpperT upper_limit) const;

        /// Evaluate log(1 - exp(-x)) without cancellation for positive x.
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

        RealT mva_base_{ZERO<RealT>};
        bool  PfFlag_{false};
        bool  VFlag_{false};
        bool  QFlag_{false};
        bool  Pqflag_{false};
        RealT Trv_{static_cast<RealT>(0.02)};
        RealT Tp_{ZERO<RealT>};
        RealT Vref0_{ZERO<RealT>};
        RealT Vdip_{static_cast<RealT>(0.85)};
        RealT Vup_{static_cast<RealT>(1.15)};
        RealT dbd1_{ZERO<RealT>};
        RealT dbd2_{ZERO<RealT>};
        RealT kqv_{static_cast<RealT>(5.0)};
        RealT Iql1_{static_cast<RealT>(-1.1)};
        RealT Iqh1_{static_cast<RealT>(1.1)};
        RealT Qmax_{static_cast<RealT>(0.436)};
        RealT Qmin_{static_cast<RealT>(-0.436)};
        RealT Kqp_{ZERO<RealT>};
        RealT Kqi_{static_cast<RealT>(0.1)};
        RealT Vmax_{static_cast<RealT>(1.1)};
        RealT Vmin_{static_cast<RealT>(0.9)};
        RealT Kvp_{static_cast<RealT>(18.0)};
        RealT Kvi_{static_cast<RealT>(5.0)};
        RealT Tiq_{static_cast<RealT>(0.02)};
        RealT Tpord_{static_cast<RealT>(0.02)};
        RealT dPmax_{static_cast<RealT>(99.0)};
        RealT dPmin_{static_cast<RealT>(-99.0)};
        RealT Pmax_{ONE<RealT>};
        RealT Pmin_{ZERO<RealT>};
        RealT Imax_{static_cast<RealT>(1.3)};
        RealT va_converter_base_{0};
        RealT pf_on_{0};
        RealT pf_off_{1};
        RealT v_on_{0};
        RealT v_off_{1};
        RealT q_on_{0};
        RealT q_off_{1};
        RealT p_priority_{0};
        RealT q_priority_{1};

        bool Vref0_given_{false};
        IdxT parameter_error_count_{0};

        ScalarT qext_set_{0};
        ScalarT pe_set_{0};
        ScalarT qgen_set_{0};
        ScalarT pfaref_set_{0};
        ScalarT pref_set_{0};

        ComponentSignals<ScalarT, IdxT, ReecbInternalVariables, ReecbExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
