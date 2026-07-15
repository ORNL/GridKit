/**
 * @file Repca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REPCA phasor-dynamics plant-control model.
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
        XQPI,   ///< Reactive PI controller state
        XQLAG,  ///< Reactive-command lead-lag state
        PMEAS,  ///< Filtered active-power signal
        XPPI,   ///< Active-power PI controller state
        PREF,   ///< Active-power command lag state
        V,      ///< Regulated-bus voltage magnitude
        VLDC,   ///< Line-drop compensated voltage magnitude
        VDROOP, ///< Reactive-droop-compensated voltage
        VCTRL,  ///< Selected voltage-measurement input
        SFRZ,   ///< Reactive-PI voltage-enable indicator
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
        IBRANCHR,  ///< Branch current real component
        IBRANCHI,  ///< Branch current imaginary component
        PBRANCH,   ///< Branch active power
        QBRANCH,   ///< Branch reactive power
        FREQ,      ///< Frequency input
        FREQREF,   ///< Frequency reference
        VREF,      ///< Voltage-control reference
        QREF,      ///< Reactive-power reference
        PPLANTREF, ///< Plant active-power reference
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class Repca : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha_;
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
        int setAbsoluteTolerance(RealT) override final;
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
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        ScalarT& Vr()
        {
          return bus_->Vr();
        }

        ScalarT& Vi()
        {
          return bus_->Vi();
        }

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        BusT* bus_{nullptr};

        RealT mva_{static_cast<RealT>(100.0)};
        bool  VcompFlag_{true};
        bool  RefFlag_{true};
        bool  Freqflag_{false};
        RealT Tfltr_{static_cast<RealT>(0.05)};
        RealT Vfrz_{static_cast<RealT>(0.7)};
        RealT Rc_{ZERO<RealT>};
        RealT Xc_{ZERO<RealT>};
        RealT Kc_{ONE<RealT>};
        RealT dbdlow_{ZERO<RealT>};
        RealT dbdupper_{ZERO<RealT>};
        RealT emax_{ONE<RealT>};
        RealT emin_{-ONE<RealT>};
        RealT Kp_{static_cast<RealT>(10.0)};
        RealT Ki_{static_cast<RealT>(10.0)};
        RealT Qmax_{ONE<RealT>};
        RealT Qmin_{-ONE<RealT>};
        RealT Tft_{ZERO<RealT>};
        RealT Tfv_{static_cast<RealT>(3.0)};
        RealT Tp_{ZERO<RealT>};
        RealT fdbd1_{ZERO<RealT>};
        RealT fdbd2_{ZERO<RealT>};
        RealT Ddn_{static_cast<RealT>(20.0)};
        RealT Dup_{ZERO<RealT>};
        RealT femax_{ONE<RealT>};
        RealT femin_{-ONE<RealT>};
        RealT Kpg_{static_cast<RealT>(10.0)};
        RealT Kig_{static_cast<RealT>(10.0)};
        RealT Pmax_{static_cast<RealT>(2.0)};
        RealT Pmin_{ZERO<RealT>};
        RealT Tlag_{static_cast<RealT>(3.0)};

        IdxT  parameter_error_count_{0};
        RealT va_component_base_{ZERO<RealT>};
        RealT vcomp_on_{ONE<RealT>};
        RealT vcomp_off_{ZERO<RealT>};
        RealT ref_on_{ONE<RealT>};
        RealT ref_off_{ZERO<RealT>};
        RealT freq_on_{ZERO<RealT>};

        ScalarT freqref_set_{ONE<RealT>};
        ScalarT vref_set_{ONE<RealT>};
        ScalarT qref_set_{ZERO<RealT>};
        ScalarT pplantref_set_{ZERO<RealT>};

        ComponentSignals<ScalarT, IdxT, RepcaInternalVariables, RepcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
