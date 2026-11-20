/**
 * @file Ieeet1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT IEEE type 1 excitation system.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeet1/Ieeet1Data.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of an `Ieeet1`, each of size 1 and in per unit
    enum class Ieeet1InternalVariables : size_t
    {
      VTS,  ///< \f$V_{\mathrm{ts}}\f$ Sensed terminal voltage
      VR,   ///< \f$V_R\f$ Voltage regulator output
      EFDP, ///< \f$E_{\mathrm{fd}}'\f$ Field voltage before the speed multiplier
      VFX,  ///< \f$V_{\mathrm{fx}}\f$ Stabilizing feedback state
      VTR,  ///< \f$V_{\mathrm{tr}}\f$ Terminal voltage error
      VF,   ///< \f$V_F\f$ Stabilizing feedback
      VE,   ///< \f$V_E\f$ Exciter control voltage
      EFD,  ///< \f$E_{\mathrm{fd}}\f$ Field voltage
      KSAT, ///< \f$E_{\mathrm{fd}}'S(E_{\mathrm{fd}}')\f$ Saturation contribution
    };

    /// IEEE type 1 excitation system
    template <typename scalar_type, typename index_type>
    class Ieeet1 final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Ieeet1Data<RealT, IdxT>;
      using InternalVariablesT = Ieeet1InternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Ieeet1Inputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Ieeet1";

      static constexpr std::array<IdxT, Utilities::enum_size<Ieeet1InternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>              NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<Ieeet1Inputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                NUM_OUTPUT_VARIABLES{};

      explicit Ieeet1(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in Ieeet1InternalVariables order
      enum : size_t
      {
        VTS,
        VR,
        EFDP,
        VFX,
        VTR,
        VF,
        VE,
        EFD,
        KSAT,
      };

      // Offsets of the external values, in Ieeet1Inputs order
      enum : size_t
      {
        VT,
        OMEGA,
        VREF,
        VS,
        VUEL,
        VOEL,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `Tr`, `Ta`, `Te` and `Tf` [s]

      RealT V_{ZERO<RealT>};
      RealT Tr_{ZERO<RealT>};
      RealT Ka_{ZERO<RealT>};
      RealT Ta_{ZERO<RealT>};
      RealT Ke_{ZERO<RealT>};
      RealT Te_{ZERO<RealT>};
      RealT Kf_{ZERO<RealT>};
      RealT Tf_{ZERO<RealT>};
      RealT Vrmin_{ZERO<RealT>};
      RealT Vrmax_{ZERO<RealT>};
      RealT E1_{ZERO<RealT>};
      RealT E2_{ZERO<RealT>};
      RealT Se1_{ZERO<RealT>};
      RealT Se2_{ZERO<RealT>};
      RealT Ispdlim_{ZERO<RealT>};

      RealT SA_{ZERO<RealT>}; ///< Field voltage at which saturation begins
      RealT SB_{ZERO<RealT>}; ///< Gain of the saturation
    };
  } // namespace EMT
} // namespace GridKit
