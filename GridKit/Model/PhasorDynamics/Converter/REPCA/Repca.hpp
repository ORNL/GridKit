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
        XQPI,   ///< Reactive-power PI controller state
        XQLAG,  ///< Reactive-command lead-lag state
        PMEAS,  ///< Filtered active-power signal
        XPPI,   ///< Active-power PI controller state
        PREF,   ///< Active-power command lag state
        V,      ///< Regulated-bus voltage magnitude
        VLDC,   ///< Line-drop compensated voltage magnitude
        VDROOP, ///< Reactive-droop-compensated voltage
        VCTRL,  ///< Selected voltage-measurement input
        SFRZ,   ///< Reactive-power PI voltage-enable gate
        ERQ,    ///< Selected reactive-loop error
        ERQDB,  ///< Deadbanded reactive-loop error
        ERQLIM, ///< Limited reactive-loop error
        QPI,    ///< Reactive-power PI output
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
        IBRANCHR,  ///< Required branch-current real component
        IBRANCHI,  ///< Required branch-current imaginary component
        PBRANCH,   ///< Required branch active power
        QBRANCH,   ///< Required branch reactive power
        FREQ,      ///< Required frequency input
        FREQREF,   ///< Optional frequency reference
        VREF,      ///< Optional voltage-control reference
        QREF,      ///< Optional reactive-power reference
        PPLANTREF, ///< Optional plant active-power reference
        MAXIMUM,
      };

      /// Indices into the REPCA state, derivative, and residual vectors.
      struct RepcaIdx
      {
        static constexpr size_t VMEAS   = static_cast<size_t>(RepcaInternalVariables::VMEAS);
        static constexpr size_t QMEAS   = static_cast<size_t>(RepcaInternalVariables::QMEAS);
        static constexpr size_t XQPI    = static_cast<size_t>(RepcaInternalVariables::XQPI);
        static constexpr size_t XQLAG   = static_cast<size_t>(RepcaInternalVariables::XQLAG);
        static constexpr size_t PMEAS   = static_cast<size_t>(RepcaInternalVariables::PMEAS);
        static constexpr size_t XPPI    = static_cast<size_t>(RepcaInternalVariables::XPPI);
        static constexpr size_t PREF    = static_cast<size_t>(RepcaInternalVariables::PREF);
        static constexpr size_t V       = static_cast<size_t>(RepcaInternalVariables::V);
        static constexpr size_t VLDC    = static_cast<size_t>(RepcaInternalVariables::VLDC);
        static constexpr size_t VDROOP  = static_cast<size_t>(RepcaInternalVariables::VDROOP);
        static constexpr size_t VCTRL   = static_cast<size_t>(RepcaInternalVariables::VCTRL);
        static constexpr size_t SFRZ    = static_cast<size_t>(RepcaInternalVariables::SFRZ);
        static constexpr size_t ERQ     = static_cast<size_t>(RepcaInternalVariables::ERQ);
        static constexpr size_t ERQDB   = static_cast<size_t>(RepcaInternalVariables::ERQDB);
        static constexpr size_t ERQLIM  = static_cast<size_t>(RepcaInternalVariables::ERQLIM);
        static constexpr size_t QPI     = static_cast<size_t>(RepcaInternalVariables::QPI);
        static constexpr size_t QEXT    = static_cast<size_t>(RepcaInternalVariables::QEXT);
        static constexpr size_t EF      = static_cast<size_t>(RepcaInternalVariables::EF);
        static constexpr size_t EP      = static_cast<size_t>(RepcaInternalVariables::EP);
        static constexpr size_t EPLIM   = static_cast<size_t>(RepcaInternalVariables::EPLIM);
        static constexpr size_t PPI     = static_cast<size_t>(RepcaInternalVariables::PPI);
        static constexpr size_t PEXT    = static_cast<size_t>(RepcaInternalVariables::PEXT);
        static constexpr size_t MAXIMUM = static_cast<size_t>(RepcaInternalVariables::MAXIMUM);
      };

      /// Indices into the REPCA external-signal buffers.
      struct RepcaExt
      {
        static constexpr size_t IBRANCHR  = static_cast<size_t>(RepcaExternalVariables::IBRANCHR);
        static constexpr size_t IBRANCHI  = static_cast<size_t>(RepcaExternalVariables::IBRANCHI);
        static constexpr size_t PBRANCH   = static_cast<size_t>(RepcaExternalVariables::PBRANCH);
        static constexpr size_t QBRANCH   = static_cast<size_t>(RepcaExternalVariables::QBRANCH);
        static constexpr size_t FREQ      = static_cast<size_t>(RepcaExternalVariables::FREQ);
        static constexpr size_t FREQREF   = static_cast<size_t>(RepcaExternalVariables::FREQREF);
        static constexpr size_t VREF      = static_cast<size_t>(RepcaExternalVariables::VREF);
        static constexpr size_t QREF      = static_cast<size_t>(RepcaExternalVariables::QREF);
        static constexpr size_t PPLANTREF = static_cast<size_t>(RepcaExternalVariables::PPLANTREF);
        static constexpr size_t MAXIMUM   = static_cast<size_t>(RepcaExternalVariables::MAXIMUM);
      };

      template <typename scalar_type, typename index_type>
      class Repca : public Component<scalar_type, index_type>
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

        bool solveLimiterInput(ScalarT  requested_output,
                               RealT    lower_limit,
                               RealT    upper_limit,
                               ScalarT& limiter_input) const;

        /// Evaluate log(1 - exp(-x)) without cancellation for positive x.
        RealT logOneMinusExp(RealT x) const;

        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM    = static_cast<RealT>(1.0e-3);
        static constexpr RealT INITIALIZATION_TOLERANCE = static_cast<RealT>(1.0e-10);

        /// Distance past a smooth-clamp bound; with CommonMath MU = 240,
        /// this keeps exact-limit initialization residuals below 2e-13.
        static constexpr RealT INITIALIZATION_LIMIT_OFFSET = static_cast<RealT>(0.1);

        BusT* bus_{nullptr};

        RealT mva_base_{static_cast<RealT>(100.0)};
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
        RealT va_converter_base_{ZERO<RealT>};
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
