/**
 * @file Repca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT renewable energy plant control.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Repca/RepcaData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Repca`, each of size 1 and in per unit on the device rating unless noted
    enum class RepcaInternalVariables : size_t
    {
      VMEAS,  ///< \f$V^\mathrm{meas}\f$ Filtered regulated voltage
      QMEAS,  ///< \f$Q^\mathrm{meas}\f$ Filtered reactive power
      XQPI,   ///< \f$x_Q^\mathrm{PI}\f$ Reactive-power PI state
      XQLAG,  ///< \f$x_Q^\mathrm{lag}\f$ Reactive-command lead-lag state
      PMEAS,  ///< \f$P^\mathrm{meas}\f$ Filtered active power
      XPPI,   ///< \f$x_P^\mathrm{PI}\f$ Active-power PI state
      PREF,   ///< \f$P^\mathrm{ref}\f$ Active-power command lag state
      VLDC,   ///< \f$V^\mathrm{ldc}\f$ Line-drop compensated voltage magnitude
      VDROOP, ///< \f$V^\mathrm{droop}\f$ Reactive-droop compensated voltage
      VCTRL,  ///< \f$V^\mathrm{ctrl}\f$ Selected voltage measurement
      SFRZ,   ///< \f$s_\mathrm{frz}\f$ Voltage-enable gate of the reactive-power PI
      ERQ,    ///< \f$e_\mathrm{RQ}\f$ Selected reactive-loop error
      ERQDB,  ///< \f$e_\mathrm{RQ}^\mathrm{db}\f$ Deadbanded reactive-loop error
      ERQLIM, ///< \f$e_\mathrm{RQ}^\mathrm{lim}\f$ Limited reactive-loop error
      QPI,    ///< \f$Q^\mathrm{PI}\f$ Reactive-power PI output
      QEXT,   ///< \f$Q^\mathrm{ext}\f$ Reactive-power command on the system power base
      EF,     ///< \f$e_f\f$ Frequency error after the deadband
      EP,     ///< \f$e_P\f$ Active-power control error
      EPLIM,  ///< \f$e_P^\mathrm{lim}\f$ Limited active-power control error
      PPI,    ///< \f$P^\mathrm{PI}\f$ Active-power PI output
      PEXT,   ///< \f$P^\mathrm{ext}\f$ Active-power command on the system power base
    };

    /// Renewable energy plant control
    template <typename scalar_type, typename index_type>
    class Repca final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = RepcaData<RealT, IdxT>;
      using InternalVariablesT = RepcaInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = RepcaInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Repca";

      static constexpr std::array<IdxT, Utilities::enum_size<RepcaInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>             NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<RepcaInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit Repca(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in RepcaInternalVariables order
      enum : size_t
      {
        VMEAS,
        QMEAS,
        XQPI,
        XQLAG,
        PMEAS,
        XPPI,
        PREF,
        VLDC,
        VDROOP,
        VCTRL,
        SFRZ,
        ERQ,
        ERQDB,
        ERQLIM,
        QPI,
        QEXT,
        EF,
        EP,
        EPLIM,
        PPI,
        PEXT,
      };

      // Offsets of the external values, in RepcaInputs order
      enum : size_t
      {
        V_IN = 0,
        P_IN,
        Q_IN,
        OMEGA,
        VREF,
        PREF_IN,
        QREF,
        FREQREF,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `Tfltr`, `Tfv`, `Tp` and `Tlag` [s]
      static constexpr RealT V_MINIMUM{0.01};               ///< Floor of the voltage that divides [p.u.]

      template <typename T>
      [[gnu::always_inline]] static T droop(T error, RealT down, RealT up);

      RealT mva_{ZERO<RealT>};
      RealT Tfltr_{ZERO<RealT>};
      RealT Vfrz_{ZERO<RealT>};
      RealT Rc_{ZERO<RealT>};
      RealT Xc_{ZERO<RealT>};
      RealT Kc_{ZERO<RealT>};
      RealT dbdlow_{ZERO<RealT>};
      RealT dbdupper_{ZERO<RealT>};
      RealT emax_{ZERO<RealT>};
      RealT emin_{ZERO<RealT>};
      RealT Kp_{ZERO<RealT>};
      RealT Ki_{ZERO<RealT>};
      RealT Qmax_{ZERO<RealT>};
      RealT Qmin_{ZERO<RealT>};
      RealT Tft_{ZERO<RealT>};
      RealT Tfv_{ZERO<RealT>};
      RealT Tp_{ZERO<RealT>};
      RealT fdbd1_{ZERO<RealT>};
      RealT fdbd2_{ZERO<RealT>};
      RealT Ddn_{ZERO<RealT>};
      RealT Dup_{ZERO<RealT>};
      RealT femax_{ZERO<RealT>};
      RealT femin_{ZERO<RealT>};
      RealT Kpg_{ZERO<RealT>};
      RealT Kig_{ZERO<RealT>};
      RealT Pmax_{ZERO<RealT>};
      RealT Pmin_{ZERO<RealT>};
      RealT Tlag_{ZERO<RealT>};
      RealT V_{ZERO<RealT>};

      RealT vcomp_on_{ZERO<RealT>};  ///< One with line-drop compensation
      RealT vcomp_off_{ZERO<RealT>}; ///< One with reactive droop
      RealT ref_on_{ZERO<RealT>};    ///< One with voltage control
      RealT ref_off_{ZERO<RealT>};   ///< One with reactive-power control
      RealT freq_on_{ZERO<RealT>};   ///< One if the active-power command is enabled

      RealT q_feedback_gain_{ZERO<RealT>}; ///< Gain that tracks the limited output of the reactive-power PI
      RealT q_integral_gain_{ZERO<RealT>}; ///< Integral gain of the reactive-power PI without a proportional gain
      RealT p_feedback_gain_{ZERO<RealT>}; ///< Gain that tracks the limited output of the active-power PI
      RealT p_integral_gain_{ZERO<RealT>}; ///< Integral gain of the active-power PI without a proportional gain
    };
  } // namespace EMT
} // namespace GridKit
