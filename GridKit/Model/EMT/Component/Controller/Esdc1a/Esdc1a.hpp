/**
 * @file Esdc1a.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT IEEE type DC1A excitation system.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esdc1a/Esdc1aData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of an `Esdc1a`, each of size 1 and in per unit
    enum class Esdc1aInternalVariables : size_t
    {
      EFDP, ///< \f$E_{\mathrm{fd}}'\f$ Exciter field-voltage state
      VC,   ///< \f$V_C\f$ Filtered terminal-voltage magnitude
      VR,   ///< \f$V_R\f$ Voltage-regulator output
      VF,   ///< \f$V_F\f$ Stabilizing feedback state
      XLL,  ///< \f$x_{\mathrm{LL}}\f$ Input lead-lag denominator state
      EV,   ///< \f$e_V\f$ Voltage-error summing output
      VLL,  ///< \f$V_{\mathrm{LL}}\f$ Input lead-lag output
      VHV,  ///< \f$V_{\mathrm{HV}}\f$ High-value gate output
      SE,   ///< \f$s_e\f$ Saturation contribution
      VFE,  ///< \f$V_{\mathrm{FE}}\f$ Exciter feedback drive
      EFD,  ///< \f$E_{\mathrm{fd}}\f$ Field voltage
    };

    /// IEEE type DC1A excitation system
    template <typename scalar_type, typename index_type>
    class Esdc1a final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Esdc1aData<RealT, IdxT>;
      using InternalVariablesT = Esdc1aInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Esdc1aInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Esdc1a";

      static constexpr std::array<IdxT, Utilities::enum_size<Esdc1aInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>              NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<Esdc1aInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                NUM_OUTPUT_VARIABLES{};

      explicit Esdc1a(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in Esdc1aInternalVariables order
      enum : size_t
      {
        EFDP,
        VC,
        VR,
        VF,
        XLL,
        EV,
        VLL,
        VHV,
        SE,
        VFE,
        EFD,
      };

      // Offsets of the external values, in Esdc1aInputs order
      enum : size_t
      {
        VT,
        OMEGA,
        VREF,
        VS,
        VUEL,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `Tr`, `Ta`, `Tb`, `Te` and `Tf1` [s]
      static constexpr IdxT  SUMMING_UEL{2};                ///< Least `UEL` that adds the limiter signal to the voltage error
      static constexpr IdxT  MAX_UEL{3};                    ///< Greatest `UEL`

      template <typename T>
      [[gnu::always_inline]] static T awmin(T x, T f, RealT xmin);

      RealT V_{ZERO<RealT>};
      RealT Tr_{ZERO<RealT>};
      RealT Ka_{ZERO<RealT>};
      RealT Ta_{ZERO<RealT>};
      RealT Tb_{ZERO<RealT>};
      RealT Tc_{ZERO<RealT>};
      RealT Vrmax_{ZERO<RealT>};
      RealT Vrmin_{ZERO<RealT>};
      RealT Ke_{ZERO<RealT>};
      RealT Te_{ZERO<RealT>};
      RealT Kf_{ZERO<RealT>};
      RealT Tf1_{ZERO<RealT>};
      RealT E1_{ZERO<RealT>};
      RealT Se1_{ZERO<RealT>};
      RealT E2_{ZERO<RealT>};
      RealT Se2_{ZERO<RealT>};
      IdxT  UEL_{0};

      RealT spd_on_{ZERO<RealT>}; ///< One with the speed multiplier
      RealT uel_on_{ZERO<RealT>}; ///< One if the limiter signal is added to the voltage error
      RealT lim_on_{ZERO<RealT>}; ///< One with the lower limit of the field-voltage state
      RealT SA_{ZERO<RealT>};     ///< Field voltage at which saturation begins
      RealT SB_{ZERO<RealT>};     ///< Gain of the saturation
    };
  } // namespace EMT
} // namespace GridKit
