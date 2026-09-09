/**
 * @file Repca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Renewable plant controller.
 */
#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class RepcaInternalVariables : size_t
      {
        VMEAS,  ///< \f$V^\mathrm{meas}\f$ Differential filtered regulated voltage [p.u.]
        QMEAS,  ///< \f$Q^\mathrm{meas}\f$ Differential filtered reactive power on component base [p.u.]
        XQPI,   ///< \f$x_Q^\mathrm{PI}\f$ Differential reactive-power PI state on component base [p.u.]
        XQLAG,  ///< \f$x_Q^\mathrm{lag}\f$ Differential reactive-command lead-lag state on component base [p.u.]
        PMEAS,  ///< \f$P^\mathrm{meas}\f$ Differential filtered active power on component base [p.u.]
        XPPI,   ///< \f$x_P^\mathrm{PI}\f$ Differential active-power PI state on component base [p.u.]
        PREF,   ///< \f$P^\mathrm{ref}\f$ Differential active-power command lag state on component base [p.u.]
        V,      ///< \f$V_t\f$ Algebraic regulated-bus voltage magnitude [p.u.]
        VLDC,   ///< \f$V^\mathrm{ldc}\f$ Algebraic line-drop compensated voltage magnitude [p.u.]
        VDROOP, ///< \f$V^\mathrm{droop}\f$ Algebraic reactive-droop-compensated voltage [p.u.]
        VCTRL,  ///< \f$V^\mathrm{ctrl}\f$ Algebraic selected voltage-measurement input [p.u.]
        SFRZ,   ///< \f$s_\mathrm{frz}\f$ Algebraic reactive-power PI voltage-enable gate [-]
        ERQ,    ///< \f$e_\mathrm{RQ}\f$ Algebraic selected reactive-loop error [p.u.]
        ERQDB,  ///< \f$e_\mathrm{RQ}^\mathrm{db}\f$ Algebraic deadbanded reactive-loop error [p.u.]
        ERQLIM, ///< \f$e_\mathrm{RQ}^\mathrm{lim}\f$ Algebraic limited reactive-loop error [p.u.]
        QPI,    ///< \f$Q^\mathrm{PI}\f$ Algebraic reactive-power PI output on component base [p.u.]
        QEXT,   ///< \f$Q^\mathrm{ext}\f$ Algebraic reactive-power command [var]
        EF,     ///< \f$e_f\f$ Algebraic frequency error after deadband [p.u.]
        EP,     ///< \f$e_P\f$ Algebraic active-power control error on component base [p.u.]
        EPLIM,  ///< \f$e_P^\mathrm{lim}\f$ Algebraic limited active-power control error on component base [p.u.]
        PPI,    ///< \f$P^\mathrm{PI}\f$ Algebraic active-power PI output on component base [p.u.]
        PEXT,   ///< \f$P^\mathrm{ext}\f$ Algebraic active-power command [W]
        MAXIMUM ///< Number of internal variables
      };

      enum class RepcaExternalVariables : size_t
      {
        VD,
        VQ,
        ID,
        IQ,
        FREQ,
        VREF,
        PREF,
        QREF,
        FREQREF,
        MAXIMUM
      };

      template <typename scalar_type, typename index_type>
      class Repca : public Component<scalar_type, index_type>
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
        using ModelDataT   = RepcaData<RealT, IdxT>;
        using Outputs      = typename ModelDataT::Outputs;
        using SignalT      = Signal<ScalarT, IdxT>;
        using MonitorT     = Model::VariableMonitor<Repca, RepcaData>;
        using InputSignals = std::array<SignalT*, 9>;

        Repca();
        explicit Repca(const ModelDataT& data);
        ~Repca();
        void     attachInput(InputSignals inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

        SignalT& inputSignal(RepcaInputs input)
        {
          return *signals_.getAttachedSignal(static_cast<RepcaExternalVariables>(input));
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
          std::array<RealT, 22> state;
          std::array<RealT, 4>  reference;
          RealT                 qmin, qmax, pmin, pmax;
        };

        OperatingPoint                                                                  operatingPoint(const std::map<Outputs, RealT>& outputs, const std::array<RealT, 9>& input) const;
        bool                                                                            invertClamp(RealT output, RealT lower, RealT upper, RealT& input) const;
        bool                                                                            invertDeadband(RealT output, RealT lower, RealT upper, RealT& input) const;
        static RealT                                                                    logOneMinusExp(RealT x);
        static constexpr RealT                                                          TIME_CONSTANT_MINIMUM       = RealT{1e-3};
        static constexpr RealT                                                          INITIALIZATION_TOLERANCE    = RealT{1e-12};
        static constexpr RealT                                                          INITIALIZATION_LIMIT_OFFSET = RealT{0.1};
        RealT                                                                           S_{0.0};
        RealT                                                                           V_{0.0};
        bool                                                                            VcompFlag_{true};
        bool                                                                            RefFlag_{true};
        bool                                                                            Freqflag_{false};
        RealT                                                                           Tfltr_{static_cast<RealT>(0.05)};
        RealT                                                                           Vfrz_{static_cast<RealT>(0.7)};
        RealT                                                                           Rc_{ZERO<RealT>};
        RealT                                                                           Xc_{ZERO<RealT>};
        RealT                                                                           Kc_{ONE<RealT>};
        RealT                                                                           dbdlow_{ZERO<RealT>};
        RealT                                                                           dbdupper_{ZERO<RealT>};
        RealT                                                                           emax_{ONE<RealT>};
        RealT                                                                           emin_{-ONE<RealT>};
        RealT                                                                           Kp_{static_cast<RealT>(10.0)};
        RealT                                                                           Ki_{static_cast<RealT>(10.0)};
        RealT                                                                           Qmax_{ONE<RealT>};
        RealT                                                                           Qmin_{-ONE<RealT>};
        RealT                                                                           Tft_{ZERO<RealT>};
        RealT                                                                           Tfv_{static_cast<RealT>(3.0)};
        RealT                                                                           Tp_{ZERO<RealT>};
        RealT                                                                           fdbd1_{ZERO<RealT>};
        RealT                                                                           fdbd2_{ZERO<RealT>};
        RealT                                                                           Ddn_{static_cast<RealT>(20.0)};
        RealT                                                                           Dup_{ZERO<RealT>};
        RealT                                                                           femax_{ONE<RealT>};
        RealT                                                                           femin_{-ONE<RealT>};
        RealT                                                                           Kpg_{static_cast<RealT>(10.0)};
        RealT                                                                           Kig_{static_cast<RealT>(10.0)};
        RealT                                                                           Pmax_{static_cast<RealT>(2.0)};
        RealT                                                                           Pmin_{ZERO<RealT>};
        RealT                                                                           Tlag_{static_cast<RealT>(3.0)};
        RealT                                                                           vcomp_on_{ONE<RealT>};
        RealT                                                                           vcomp_off_{ZERO<RealT>};
        RealT                                                                           ref_on_{ONE<RealT>};
        RealT                                                                           ref_off_{ZERO<RealT>};
        RealT                                                                           freq_on_{ZERO<RealT>};
        std::array<RealT, 4>                                                            reference_{1.0, 0.0, 0.0, 1.0};
        ComponentSignals<ScalarT, IdxT, RepcaInternalVariables, RepcaExternalVariables> signals_;
        std::array<SignalT, 2>                                                          output_;
        std::array<SignalT*, 2>                                                         alias_{};
        std::unique_ptr<MonitorT>                                                       monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
