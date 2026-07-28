/**
 * @file Repca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REPCA phasor-dynamics plant-control model.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    namespace Controller
    {
      /// Internal variables of `Repca`.
      enum class RepcaInternalVariables : size_t
      {
        VMEAS,  ///< \f$V^\mathrm{meas}\f$ Differential filtered regulated voltage [p.u.]
        QMEAS,  ///< \f$Q^\mathrm{meas}\f$ Differential filtered reactive power on component base [p.u.]
        XQPI,   ///< \f$x_Q^\mathrm{PI}\f$ Differential reactive-power PI state on component base [p.u.]
        XQLAG,  ///< \f$x_Q^\mathrm{lag}\f$ Differential reactive-command lead-lag state on component base [p.u.]
        PMEAS,  ///< \f$P^\mathrm{meas}\f$ Differential filtered active power on component base [p.u.]
        XPPI,   ///< \f$x_P^\mathrm{PI}\f$ Differential active-power PI state on component base [p.u.]
        PREF,   ///< \f$P^\mathrm{ref}\f$ Differential active-power command lag state on component base [p.u.]
        V,      ///< \f$V\f$ Algebraic regulated-bus voltage magnitude [p.u.]
        VLDC,   ///< \f$V^\mathrm{ldc}\f$ Algebraic line-drop compensated voltage magnitude [p.u.]
        VDROOP, ///< \f$V^\mathrm{droop}\f$ Algebraic reactive-droop-compensated voltage [p.u.]
        VCTRL,  ///< \f$V^\mathrm{ctrl}\f$ Algebraic selected voltage-measurement input [p.u.]
        SFRZ,   ///< \f$s_\mathrm{frz}\f$ Algebraic reactive-power PI voltage-enable gate [-]
        ERQ,    ///< \f$e_\mathrm{RQ}\f$ Algebraic selected reactive-loop error [p.u.]
        ERQDB,  ///< \f$e_\mathrm{RQ}^\mathrm{db}\f$ Algebraic deadbanded reactive-loop error [p.u.]
        ERQLIM, ///< \f$e_\mathrm{RQ}^\mathrm{lim}\f$ Algebraic limited reactive-loop error [p.u.]
        QPI,    ///< \f$Q^\mathrm{PI}\f$ Algebraic reactive-power PI output on component base [p.u.]
        QEXT,   ///< \f$Q^\mathrm{ext}\f$ Algebraic reactive-power command on system base [p.u.]
        EF,     ///< \f$e_f\f$ Algebraic frequency error after deadband [p.u.]
        EP,     ///< \f$e_P\f$ Algebraic active-power control error on component base [p.u.]
        EPLIM,  ///< \f$e_P^\mathrm{lim}\f$ Algebraic limited active-power control error on component base [p.u.]
        PPI,    ///< \f$P^\mathrm{PI}\f$ Algebraic active-power PI output on component base [p.u.]
        PEXT,   ///< \f$P^\mathrm{ext}\f$ Algebraic active-power command on system base [p.u.]
        MAXIMUM ///< Number of internal variables
      };

      /// External variables of `Repca`.
      enum class RepcaExternalVariables : size_t
      {
        IR,      ///< \f$I_\mathrm{r}\f$ Required branch-current real component on system base [p.u.]
        II,      ///< \f$I_\mathrm{i}\f$ Required branch-current imaginary component on system base [p.u.]
        P,       ///< \f$P\f$ Required branch active power on system base [p.u.]
        Q,       ///< \f$Q\f$ Required branch reactive power on system base [p.u.]
        FREQ,    ///< \f$f\f$ Optional absolute frequency input [p.u.]
        VREF,    ///< \f$V^\mathrm{ref}\f$ Optional voltage-control reference [p.u.]
        PREF,    ///< \f$P_\mathrm{plant}^\mathrm{ref}\f$ Optional plant active-power reference on system base [p.u.]
        QREF,    ///< \f$Q^\mathrm{ref}\f$ Optional reactive-power reference on system base [p.u.]
        FREQREF, ///< \f$f^\mathrm{ref}\f$ Optional absolute frequency reference [p.u.]
        MAXIMUM  ///< Number of external variables
      };

      /**
       * @class Repca
       * @brief WECC renewable plant controller with reactive-power and
       *        active-power control paths.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class Repca : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha;
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
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using BusT               = BusBase<ScalarT, IdxT>;
        using SignalT            = SignalNode<ScalarT, IdxT>;
        using ModelDataT         = RepcaData<RealT, IdxT>;
        using MonitorT           = Model::VariableMonitor<Repca, RepcaData>;
        using InternalVariablesT = RepcaInternalVariables;
        using ExternalVariablesT = RepcaExternalVariables;

        static constexpr RealT INITIALIZATION_TOLERANCE = static_cast<RealT>(1.0e-12);

        Repca(BusT* bus);
        Repca(BusT* bus, const ModelDataT& data);
        ~Repca();

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
                                RepcaInternalVariables,
                                RepcaExternalVariables>&;

        const Model::VariableMonitorBase* getMonitor() const override;

        [[gnu::always_inline]] inline int evaluateInternalResidual(
            const ScalarT* y,
            const ScalarT* yp,
            const ScalarT* wb,
            const ScalarT* ws,
            ScalarT*       f);

      private:
        /// Smooth asymmetric frequency-droop response.
        static __attribute__((always_inline)) inline ScalarT droop(ScalarT error, RealT down, RealT up);

        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();

        bool invertClamp(ScalarT output, RealT lower, RealT upper, ScalarT& input) const;

        bool invertDeadband(ScalarT output, RealT lower, RealT upper, ScalarT& input) const;

        static RealT logOneMinusExp(RealT x);

        [[gnu::always_inline]] inline ScalarT toComponentBase(ScalarT value) const;
        ScalarT                               toSystemBase(ScalarT value) const;

        ScalarT& Vr();
        ScalarT& Vi();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

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
        RealT va_component_base_{ZERO<RealT>};
        RealT vcomp_on_{ONE<RealT>};
        RealT vcomp_off_{ZERO<RealT>};
        RealT ref_on_{ONE<RealT>};
        RealT ref_off_{ZERO<RealT>};
        RealT freq_on_{ZERO<RealT>};

        ScalarT freqref_set_{ONE<RealT>};
        ScalarT vref_set_{ONE<RealT>};
        ScalarT qref_set_{ZERO<RealT>};
        ScalarT pref_set_{ZERO<RealT>};

        ComponentSignals<ScalarT, IdxT, RepcaInternalVariables, RepcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
