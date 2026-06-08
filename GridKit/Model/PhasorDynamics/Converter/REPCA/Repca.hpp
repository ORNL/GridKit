/**
 * @file Repca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REPCA plant-control model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
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
      /// Internal variables of a `Repca`.
      enum class RepcaInternalVariables : size_t
      {
        VMEAS,  ///< Filtered regulated voltage
        QMEAS,  ///< Filtered reactive-power signal
        XQ,     ///< Reactive PI state
        XQEXT,  ///< Reactive-command lead-lag state
        PMEAS,  ///< Filtered active-power signal
        XP,     ///< Active-power PI state
        PREF,   ///< Active-power command lag state
        VREG,   ///< Regulated-bus voltage magnitude
        VLDC,   ///< Line-drop compensated voltage magnitude
        VDROOP, ///< Reactive-droop-compensated voltage
        VCTRL,  ///< Selected voltage-measurement input
        SFRZ,   ///< Reactive PI voltage-enable indicator
        ERQ,    ///< Selected reactive-loop error
        ERQDB,  ///< Deadbanded reactive-loop error
        ERQLIM, ///< Limited reactive-loop error
        QPI,    ///< Reactive PI output
        QEXT,   ///< Reactive-power command output
        EF,     ///< Frequency error after deadband
        EP,     ///< Active-power control error
        EPLIM,  ///< Limited active-power control error
        PPI,    ///< Active-power PI output
        PEXT,   ///< Active-power command output
        MAXIMUM,
      };

      /// External variables of a `Repca`.
      enum class RepcaExternalVariables : size_t
      {
        IBRANCHR,  ///< Branch current real component on component base
        IBRANCHI,  ///< Branch current imaginary component on component base
        QBRANCH,   ///< Branch reactive-power signal on system base
        PBRANCH,   ///< Branch active-power signal on system base
        VREF,      ///< Voltage reference
        QREF,      ///< Reactive-power reference
        PPLANTREF, ///< Plant active-power reference
        FREQ,      ///< Optional frequency input
        FREQREF,   ///< Optional frequency reference
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class Repca : public Component<scalar_type, index_type>
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
        using ModelDataT = RepcaData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Repca, RepcaData>;

        Repca(BusT* bus);
        Repca(BusT* bus, const ModelDataT& data);
        ~Repca();

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
                                RepcaInternalVariables,
                                RepcaExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        void initModelParams(const ModelDataT& data);
        void initializeMonitor();

        int     verifyBranchSignalPorts() const;
        ScalarT readExternalOrDefault(RepcaExternalVariables variable, ScalarT default_value) const;
        void    setDerivedParameters();
        ScalarT toComponentBase(ScalarT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        BusT* bus_{nullptr};

        RealT mva_base_{0};
        RealT VcompFlag_{0};
        RealT RefFlag_{0};
        RealT Freqflag_{0};
        RealT Tfltr_{0};
        RealT Tft_{0};
        RealT Tfv_{0};
        RealT Tp_{0};
        RealT Tlag_{0};
        RealT Vfrz_{0};
        RealT Rc_{0};
        RealT Xc_{0};
        RealT Kc_{0};
        RealT dbdlow_{0};
        RealT dbdupper_{0};
        RealT emax_{0};
        RealT emin_{0};
        RealT Kp_{0};
        RealT Ki_{0};
        RealT Qmax_{0};
        RealT Qmin_{0};
        RealT fdbd1_{0};
        RealT fdbd2_{0};
        RealT Ddn_{0};
        RealT Dup_{0};
        RealT femax_{0};
        RealT femin_{0};
        RealT Kpg_{0};
        RealT Kig_{0};
        RealT Pmax_{0};
        RealT Pmin_{0};

        IdxT  parameter_error_count_{0};
        RealT system_to_component_base_{1};
        RealT vcomp_off_{1};
        RealT ref_off_{1};

        ScalarT vref_set_{0};
        ScalarT qref_set_{0};
        ScalarT pplantref_set_{0};
        ScalarT freq_set_{0};
        ScalarT freqref_set_{0};

        ComponentSignals<ScalarT, IdxT, RepcaInternalVariables, RepcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
