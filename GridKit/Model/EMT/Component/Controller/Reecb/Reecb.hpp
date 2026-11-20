/**
 * @file Reecb.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT renewable energy electrical control.
 */

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <string_view>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reecb/ReecbData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Reecb`, each of size 1 and in per unit on the device rating
    enum class ReecbInternalVariables : size_t
    {
      VMEAS,  ///< \f$V^\mathrm{meas}\f$ Filtered terminal voltage
      PMEAS,  ///< \f$P^\mathrm{meas}\f$ Filtered electrical power
      XPIQ,   ///< \f$x_Q^\mathrm{PI}\f$ Reactive-power PI state
      XPIV,   ///< \f$x_V^\mathrm{PI}\f$ Voltage-control PI state
      QV,     ///< \f$Q_V\f$ Reactive-current command lag state
      PORD,   ///< \f$P^\mathrm{ord}\f$ Filtered active-power order
      VSAFE,  ///< \f$V_\mathrm{safe}^\mathrm{meas}\f$ Safe measured voltage
      SDIP,   ///< \f$s_\mathrm{dip}\f$ Voltage-band gate
      IQV,    ///< \f$I_q^\mathrm{inj}\f$ Reactive-current injection
      QREF,   ///< \f$Q^\mathrm{ref}\f$ Reactive-power reference
      EQ,     ///< \f$e_Q\f$ Reactive-power error
      VPIQ,   ///< \f$V_Q^\mathrm{PI}\f$ Reactive-power PI output
      EPIV,   ///< \f$e_V^\mathrm{PI}\f$ Voltage-control error
      RPORD,  ///< \f$r_P^\mathrm{ord}\f$ Limited active-power order rate [p.u./s]
      ILCAP,  ///< \f$I_L^\mathrm{cap}\f$ Off-axis current capacity
      IQMAX,  ///< \f$I_q^{\max}\f$ Reactive-current limit
      IPMAX,  ///< \f$I_p^{\max}\f$ Active-current limit
      IQBASE, ///< \f$I_q^\mathrm{base}\f$ Voltage-controller current
      IQRAW,  ///< \f$I_q^\mathrm{raw}\f$ Reactive-current command before the limit
      IQCMD,  ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command
      IPCMD,  ///< \f$I_p^\mathrm{cmd}\f$ Active-current command
    };

    /// Renewable energy electrical control
    template <typename scalar_type, typename index_type>
    class Reecb final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = ReecbData<RealT, IdxT>;
      using InternalVariablesT = ReecbInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = ReecbInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Reecb";

      static constexpr std::array<IdxT, Utilities::enum_size<ReecbInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>             NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<ReecbInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit Reecb(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in ReecbInternalVariables order
      enum : size_t
      {
        VMEAS,
        PMEAS,
        XPIQ,
        XPIV,
        QV,
        PORD,
        VSAFE,
        SDIP,
        IQV,
        QREF,
        EQ,
        VPIQ,
        EPIV,
        RPORD,
        ILCAP,
        IQMAX,
        IPMAX,
        IQBASE,
        IQRAW,
        IQCMD,
        IPCMD,
      };

      // Offsets of the external values, in ReecbInputs order
      enum : size_t
      {
        VT_IN = 0,
        P_IN,
        Q_IN,
        QEXT,
        PFAREF,
        PREF,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `Trv`, `Tp`, `Tiq` and `Tpord` [s]
      static constexpr RealT VMEAS_MINIMUM{0.01};           ///< Floor of the measured voltage that divides [p.u.]

      /// Width of the knee of the current circle
      static constexpr RealT CURRENT_CIRCLE_KNEE = 100.0 * std::numeric_limits<RealT>::epsilon() / Math::MU<RealT>;

      template <typename T>
      [[gnu::always_inline]] static T aslew(T rate, RealT lower, RealT upper);

      template <typename T>
      [[gnu::always_inline]] static T awband(T state, T rate, T band);

      template <typename T>
      [[gnu::always_inline]] static T sqrtramp(T x);

      template <typename T>
      [[gnu::always_inline]] static T circleSquare(RealT limit, T high);

      RealT Trv_{ZERO<RealT>};
      RealT Tp_{ZERO<RealT>};
      RealT Vref0_{ZERO<RealT>};
      RealT Vdip_{ZERO<RealT>};
      RealT Vup_{ZERO<RealT>};
      RealT dbd1_{ZERO<RealT>};
      RealT dbd2_{ZERO<RealT>};
      RealT kqv_{ZERO<RealT>};
      RealT Iql1_{ZERO<RealT>};
      RealT Iqh1_{ZERO<RealT>};
      RealT Qmax_{ZERO<RealT>};
      RealT Qmin_{ZERO<RealT>};
      RealT Kqp_{ZERO<RealT>};
      RealT Kqi_{ZERO<RealT>};
      RealT Vmax_{ZERO<RealT>};
      RealT Vmin_{ZERO<RealT>};
      RealT Kvp_{ZERO<RealT>};
      RealT Kvi_{ZERO<RealT>};
      RealT Tiq_{ZERO<RealT>};
      RealT Tpord_{ZERO<RealT>};
      RealT dPmax_{ZERO<RealT>};
      RealT dPmin_{ZERO<RealT>};
      RealT Pmax_{ZERO<RealT>};
      RealT Pmin_{ZERO<RealT>};
      RealT Imax_{ZERO<RealT>};
      RealT V_{ZERO<RealT>};
      RealT mva_{ZERO<RealT>};

      RealT pf_on_{ZERO<RealT>};    ///< One with power-factor control
      RealT pf_off_{ZERO<RealT>};   ///< One without it
      RealT q_on_{ZERO<RealT>};     ///< One with Volt/VAr control
      RealT q_off_{ZERO<RealT>};    ///< One with the reactive-current lag
      RealT q_pi_on_{ZERO<RealT>};  ///< One with the cascaded reactive-power PI
      RealT v_ref_on_{ZERO<RealT>}; ///< One if input `qext` is a voltage command
      RealT q_ref_on_{ZERO<RealT>}; ///< One if it is a reactive-power command
      RealT pq_on_{ZERO<RealT>};    ///< One with active-current priority
      RealT pq_off_{ZERO<RealT>};   ///< One with reactive-current priority
    };
  } // namespace EMT
} // namespace GridKit
