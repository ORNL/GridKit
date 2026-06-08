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
        SDIP,      ///< Voltage-dip/overvoltage freeze indicator
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

      template <typename scalar_type, typename index_type>
      class Reecb : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::J_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
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
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        void initModelParams(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        ScalarT toComponentBase(ScalarT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        BusT* bus_{nullptr};

        RealT mva_base_{0};
        RealT PfFlag_{0};
        RealT VFlag_{0};
        RealT QFlag_{0};
        RealT Pqflag_{0};
        RealT Trv_{0};
        RealT Tp_{0};
        RealT Vref0_{0};
        RealT Vdip_{0};
        RealT Vup_{0};
        RealT dbd1_{0};
        RealT dbd2_{0};
        RealT kqv_{0};
        RealT Iql1_{0};
        RealT Iqh1_{0};
        RealT Qmax_{0};
        RealT Qmin_{0};
        RealT Kqp_{0};
        RealT Kqi_{0};
        RealT Vmax_{0};
        RealT Vmin_{0};
        RealT Kvp_{0};
        RealT Kvi_{0};
        RealT Tiq_{0};
        RealT Tpord_{0};
        RealT dPmax_{0};
        RealT dPmin_{0};
        RealT Pmax_{0};
        RealT Pmin_{0};
        RealT Imax_{0};
        RealT va_converter_base_{0};
        RealT pf_off_{1};
        RealT v_off_{1};
        RealT q_off_{1};

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
