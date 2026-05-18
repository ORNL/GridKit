/**
 * @file Reeca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REECA phasor-dynamics converter-control model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/ReecaData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
    class SignalNode;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Internal variables of a `Reeca`.
      enum class ReecaInternalVariables : size_t
      {
        VMEAS,      ///< Filtered terminal voltage
        PMEAS,      ///< Filtered active-power feedback
        XPIQ,       ///< Reactive-power PI controller state
        XPIV,       ///< Voltage PI controller state
        QV,         ///< Reactive-current command lag state
        PORD,       ///< Filtered active-power order
        VT,         ///< Terminal voltage magnitude
        VMEAS_SAFE, ///< Voltage magnitude guarded for division
        SDIP,       ///< Voltage dip/overvoltage indicator
        VERR,       ///< Voltage-control deadband output
        IQV,        ///< Voltage-support reactive-current command
        QREF,       ///< Reactive-power reference
        EQ,         ///< Reactive-power error
        VPIQ,       ///< Reactive-power PI output
        EPIV,       ///< Voltage-control error
        FPORD,      ///< Active-power order pre-filter derivative
        RPORD,      ///< Active-power order rate-limiter output
        IQCIRC,     ///< Reactive-current circle limit
        IPCIRC,     ///< Active-current circle limit
        IQMAX,      ///< Reactive-current maximum
        IPMAX,      ///< Active-current maximum
        IQBASE,     ///< Base reactive-current command
        IQRAW,      ///< Raw reactive-current command
        IQCMD,      ///< Reactive-current command output
        IPCMD,      ///< Active-current command output
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Reeca : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::h_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::mva_system_base_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::residual_indices_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using bus_type        = BusBase<ScalarT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using model_data_type = ReecaData<RealT, IdxT>;
        using MonitorT        = Model::VariableMonitor<Reeca, ReecaData>;

        Reeca(bus_type* bus);
        Reeca(bus_type* bus, const model_data_type& data);
        ~Reeca();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                ReecaInternalVariables,
                                ReecaExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        void initializeParameters(const model_data_type& data);
        void initializeMonitor();
        void setDerivedParameters();

        RealT vaSystemBase() const
        {
          return mva_system_base_ * static_cast<RealT>(1.0e6);
        }

        RealT vaReecaBase() const
        {
          return Sbase_ * static_cast<RealT>(1.0e6);
        }

        template <typename ValueT>
        ValueT toReecaBase(ValueT value) const
        {
          return value * static_cast<ValueT>(vaSystemBase() / vaReecaBase());
        }

        template <typename ValueT>
        ValueT toSystemBase(ValueT value) const
        {
          return value / static_cast<ValueT>(vaSystemBase() / vaReecaBase());
        }

        template <typename ValueT>
        __attribute__((always_inline)) inline ValueT vdlQ(ValueT voltage) const
        {
          return static_cast<ValueT>(lq1_)
                 + Math::linseg(voltage, vq1_, vq2_, lq2_ - lq1_)
                 + Math::linseg(voltage, vq2_, vq3_, lq3_ - lq2_)
                 + Math::linseg(voltage, vq3_, vq4_, lq4_ - lq3_);
        }

        template <typename ValueT>
        __attribute__((always_inline)) inline ValueT vdlP(ValueT voltage) const
        {
          return static_cast<ValueT>(lp1_)
                 + Math::linseg(voltage, vp1_, vp2_, lp2_ - lp1_)
                 + Math::linseg(voltage, vp2_, vp3_, lp3_ - lp2_)
                 + Math::linseg(voltage, vp3_, vp4_, lp4_ - lp3_);
        }

        template <ReecaExternalVariables variable>
        bool isSignalLinked() const
        {
          return signals_.template isAttached<variable>() && signals_.template isLinked<variable>();
        }

        ScalarT& Vr()
        {
          return bus_->Vr();
        }

        ScalarT& Vi()
        {
          return bus_->Vi();
        }

        const ScalarT& Vr() const
        {
          return bus_->Vr();
        }

        const ScalarT& Vi() const
        {
          return bus_->Vi();
        }

        bus_type* bus_{nullptr};

        RealT Sbase_{0};
        bool  spf_{false};
        bool  sV_{false};
        bool  sQ_{false};
        bool  sP_{false};
        bool  sPQ_{false};
        RealT Trv_{0};
        RealT Tp_{0};
        RealT Vref0_{0};
        RealT Vdip_{0};
        RealT Vup_{0};
        RealT Dbd1_{0};
        RealT Dbd2_{0};
        RealT Kqv_{0};
        RealT Iql1_{0};
        RealT Iqh1_{0};
        RealT Iqfrz_{0};
        RealT Thld_{0};
        RealT Qmax_{0};
        RealT Qmin_{0};
        RealT Kqp_{0};
        RealT Kqi_{0};
        RealT Vmax_{0};
        RealT Vmin_{0};
        RealT Vref1_{0};
        RealT Kvp_{0};
        RealT Kvi_{0};
        RealT Tiq_{0};
        RealT Tpord_{0};
        RealT dPmax_{0};
        RealT dPmin_{0};
        RealT Pmax_{0};
        RealT Pmin_{0};
        RealT Imax_{0};
        RealT vq1_{0};
        RealT lq1_{0};
        RealT vq2_{0};
        RealT lq2_{0};
        RealT vq3_{0};
        RealT lq3_{0};
        RealT vq4_{0};
        RealT lq4_{0};
        RealT vp1_{0};
        RealT lp1_{0};
        RealT vp2_{0};
        RealT lp2_{0};
        RealT vp3_{0};
        RealT lp3_{0};
        RealT vp4_{0};
        RealT lp4_{0};
        RealT Thld2_{0};

        IdxT bus_id_{0};
        IdxT parameter_error_count_{0};

        RealT spf_on_{0};
        RealT spf_off_{1};
        RealT sV_on_{0};
        RealT sV_off_{1};
        RealT sQ_on_{0};
        RealT sQ_off_{1};
        RealT sP_on_{0};
        RealT sPQ_on_{0};
        RealT sPQ_off_{1};

        bool has_Vref0_{false};
        bool has_pe_init_{false};
        bool has_qgen_init_{false};
        bool has_omega_init_{false};
        bool has_qext_init_{false};
        bool has_pfaref_init_{false};
        bool has_pref_init_{false};

        RealT pe_init_{0};
        RealT qgen_init_{0};
        RealT omega_init_{0};
        RealT qext_init_{0};
        RealT pfaref_init_{0};
        RealT pref_init_{0};

        ScalarT pe_set_{0};
        ScalarT qgen_set_{0};
        ScalarT omega_set_{0};
        ScalarT qext_set_{0};
        ScalarT pfaref_set_{0};
        ScalarT pref_set_{0};

        ComponentSignals<ScalarT, IdxT, ReecaInternalVariables, ReecaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
