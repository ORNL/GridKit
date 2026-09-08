#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Source/REGFMA/RegfmaData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class RegfmaInternalVariables : size_t
    {
      PF,    ///< \f$P_f\f$ Filtered active power [p.u.]
      QF,    ///< \f$Q_f\f$ Filtered reactive power [p.u.]
      VF,    ///< \f$V_f\f$ Filtered voltage magnitude [p.u.]
      XPMAX, ///< \f$x_P^{\max}\f$ Upper active-power-limit integral [p.u.]
      XPMIN, ///< \f$x_P^{\min}\f$ Lower active-power-limit integral [p.u.]
      XQMAX, ///< \f$x_Q^{\max}\f$ Upper reactive-power-limit integral [p.u.]
      XQMIN, ///< \f$x_Q^{\min}\f$ Lower reactive-power-limit integral [p.u.]
      XV,    ///< \f$x_V\f$ Voltage-control integral [p.u.]
      DELTA, ///< \f$\delta\f$ Internal angle deviation [rad]
      IA,    ///< \f$i_a\f$ Phase-a current injection [A]
      IB,    ///< \f$i_b\f$ Phase-b current injection [A]
      IC,    ///< \f$i_c\f$ Phase-c current injection [A]
      MAXIMUM,
    };

    enum class RegfmaExternalVariables : size_t
    {
      VA,   ///< \f$v_a\f$ Phase-a terminal voltage [V]
      VB,   ///< \f$v_b\f$ Phase-b terminal voltage [V]
      VC,   ///< \f$v_c\f$ Phase-c terminal voltage [V]
      PREF, ///< \f$P_\mathrm{ref}\f$ Active-power reference [p.u.]
      QREF, ///< \f$Q_\mathrm{ref}\f$ Reactive-power reference [p.u.]
      VREF, ///< \f$V_\mathrm{ref}\f$ Voltage reference [p.u.]
      MAXIMUM,
    };

    /// Averaged EMT realization of WECC REGFM_A1 with physical RL coupling.
    template <typename scalar_type, typename index_type>
    class Regfma : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using Base       = Component<ScalarT, IdxT>;
      using RealT      = typename Base::RealT;
      using SignalT    = typename Base::SignalT;
      using ModelDataT = RegfmaData<RealT, IdxT>;
      using Outputs    = typename ModelDataT::Outputs;
      using MonitorT   = Model::VariableMonitor<Regfma, RegfmaData>;

      explicit Regfma(const ModelDataT& data);
      ~Regfma() override;

      int setGridKitComponentID(IdxT id) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize(const std::map<Outputs, RealT>& outputs = {});

      int  initializeState(const std::map<std::string, RealT>& values) override;
      void validateInitialState(const std::map<std::string, RealT>& values) const override;

      typename Base::InitializationPortsT initializationPorts() override;
      int                                 setAbsoluteTolerance(RealT tolerance) override final;
      int                                 evaluateInternalResidual() override final;
      int                                 evaluateResidual() override final;
      int                                 assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      SignalT& currentSignal(size_t phase)
      {
        return current_.at(phase);
      }

      void assignOutput(Outputs output, SignalT* signal);

      auto getSignals() -> ComponentSignals<ScalarT, IdxT, RegfmaInternalVariables, RegfmaExternalVariables>&
      {
        return signals_;
      }

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* ye, const ScalarT* ype, ScalarT* f);

    protected:
      void gatherExternalVariables() override;

    private:
      using I = RegfmaInternalVariables;
      using E = RegfmaExternalVariables;

      static constexpr size_t PF    = static_cast<size_t>(I::PF);
      static constexpr size_t QF    = static_cast<size_t>(I::QF);
      static constexpr size_t VF    = static_cast<size_t>(I::VF);
      static constexpr size_t XPMAX = static_cast<size_t>(I::XPMAX);
      static constexpr size_t XPMIN = static_cast<size_t>(I::XPMIN);
      static constexpr size_t XQMAX = static_cast<size_t>(I::XQMAX);
      static constexpr size_t XQMIN = static_cast<size_t>(I::XQMIN);
      static constexpr size_t XV    = static_cast<size_t>(I::XV);
      static constexpr size_t DELTA = static_cast<size_t>(I::DELTA);
      static constexpr size_t IA    = static_cast<size_t>(I::IA);
      static constexpr size_t IB    = static_cast<size_t>(I::IB);
      static constexpr size_t IC    = static_cast<size_t>(I::IC);
      static constexpr size_t VA    = static_cast<size_t>(E::VA);
      static constexpr size_t VB    = static_cast<size_t>(E::VB);
      static constexpr size_t VC    = static_cast<size_t>(E::VC);
      static constexpr size_t PREF  = static_cast<size_t>(E::PREF);
      static constexpr size_t QREF  = static_cast<size_t>(E::QREF);
      static constexpr size_t VREF  = static_cast<size_t>(E::VREF);

      void                              initializeParameters(const ModelDataT& data);
      void                              initializeMonitor();
      const Model::VariableMonitorBase* getMonitor() const override;

      /// P-limit correction, Q-limit correction, voltage error, droop voltage, and frequency deviation.
      __attribute__((always_inline)) inline std::array<ScalarT, 5> controls(const ScalarT* y, const ScalarT* ye) const;
      /// Terminal alpha-beta voltage/current and terminal active/reactive power, in per unit.
      __attribute__((always_inline)) inline std::array<ScalarT, 6> measurements(const ScalarT* y, const ScalarT* ye) const;
      /// Limited alpha-beta current reference in per unit.
      __attribute__((always_inline)) inline std::array<ScalarT, 2> source(const ScalarT* y, const ScalarT* ye) const;
      ScalarT                                                      monitorValue(RegfmaMonitorableVariables variable);
      RealT                                                        initialVoltageCommand(RealT voltage) const;

      static constexpr RealT VOLTAGE_EPSILON = static_cast<RealT>(1.0e-8);

      RealT S_, V_, omega0_, XL_, RL_, mp_, mq_, kpv_, kiv_;
      RealT Emin_, Emax_, Pmin_, Pmax_, Qmin_, Qmax_;
      RealT kppmax_, kipmax_, kpqmax_, kiqmax_, TPf_, TQf_, TVf_, ImaxF_;
      RealT voltage_control_;
      bool  QVFlag_;
      RealT current_base_;
      RealT kpP_, kiP_;

      ScalarT                                                                           pref_set_{0}, qref_set_{0}, vref_set_{1};
      std::array<SignalT, 3>                                                            current_;
      std::array<SignalT*, 3>                                                           alias_{};
      ComponentSignals<ScalarT, IdxT, RegfmaInternalVariables, RegfmaExternalVariables> signals_;
      std::unique_ptr<MonitorT>                                                         monitor_;
    };
  } // namespace EMT
} // namespace GridKit
