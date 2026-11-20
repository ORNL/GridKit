/**
 * @file Bridge.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT converter bridge.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Converter/Bridge/BridgeData.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/Pwm.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Bridge`
    enum class BridgeInternalVariables : size_t
    {
      i,    ///< \f$\mathbf{i}\f$ Converter-side currents, size 3 [A]
      vo,   ///< \f$\mathbf{v}_\mathrm{o}\f$ Capacitor voltages, size 3 [V]
      ig,   ///< \f$\mathbf{i}_g\f$ Grid-side currents, size 3 [A]
      idq,  ///< \f$\mathbf{i}_\mathrm{dq}\f$ Converter-side currents, d then q [A]
      vodq, ///< \f$\mathbf{v}_\mathrm{o,dq}\f$ Capacitor voltages, d then q [V]
    };

    /// Two-level three-leg converter with an LCL filter to its bus
    template <typename scalar_type, typename index_type>
    class Bridge final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using PwmT       = Pwm<RealT, IdxT>;

      static constexpr size_t LEGS{3};

    public:
      using ModelDataT         = BridgeData<RealT, IdxT>;
      using InternalVariablesT = BridgeInternalVariables;
      using TerminalsT         = BridgeTerminals;
      using InputsT            = BridgeInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Bridge";

      static constexpr std::array<IdxT, Utilities::enum_size<BridgeInternalVariables>()> NUM_INTERNAL_VARIABLES{LEGS, LEGS, LEGS, 2, 2};
      static constexpr std::array<IdxT, Utilities::enum_size<BridgeTerminals>()>         NUM_TERMINAL_VARIABLES{LEGS, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<BridgeInputs>()>            NUM_INPUT_VARIABLES{1, 2};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                NUM_OUTPUT_VARIABLES{};

      explicit Bridge(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

    private:
      template <PwmSeries series, typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <PwmSeries series, typename T>
      [[gnu::always_inline]] void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      /// Kernels of a bridge whose switching functions are evaluated by `series`
      template <PwmSeries series>
      struct Kernels final
      {
        using InternalVariablesT = BridgeInternalVariables;
        using TerminalsT         = BridgeTerminals;
        using InputsT            = BridgeInputs;
        using OutputsT           = Utilities::EmptyEnum;

        static constexpr auto NUM_INTERNAL_VARIABLES = Bridge::NUM_INTERNAL_VARIABLES;
        static constexpr auto NUM_TERMINAL_VARIABLES = Bridge::NUM_TERMINAL_VARIABLES;
        static constexpr auto NUM_INPUT_VARIABLES    = Bridge::NUM_INPUT_VARIABLES;
        static constexpr auto NUM_OUTPUT_VARIABLES   = Bridge::NUM_OUTPUT_VARIABLES;

        template <typename T>
        [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const
        {
          bridge.template internalResidual<series>(y, yp, w, wp, f);
        }

        template <typename T>
        [[gnu::always_inline]] void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const
        {
          bridge.template externalResidual<series>(y, yp, w, wp, f);
        }

        const Bridge& bridge;
      };

      template <typename EvaluateT>
      int evaluate(EvaluateT evaluate) const;

      // Local indices of the internal variables, in BridgeInternalVariables order
      enum : size_t
      {
        I  = 0,
        VO = LEGS,
        IG = 2 * LEGS,
        ID = 3 * LEGS,
        IQ,
        VOD,
        VOQ,
      };

      // Offsets of the external values: the terminals, then inputs `theta` and `mdq`
      enum : size_t
      {
        BUS  = 0,
        VPOS = LEGS,
        VNEG,
        THETA,
        MD,
        MQ,
      };

      /// Cosine and sine of the angle of each leg, and its switching function
      template <typename T>
      struct Legs
      {
        std::array<T, LEGS> cosine;
        std::array<T, LEGS> sine;
        std::array<T, LEGS> s;
      };

      template <PwmSeries series, typename T>
      [[gnu::always_inline]] Legs<T> legs(const T* w) const;

      IdxT  N_{0};
      RealT Rs_{ZERO<RealT>};
      RealT Ls_{ZERO<RealT>};
      RealT C_{ZERO<RealT>};
      RealT Rg_{ZERO<RealT>};
      RealT Lg_{ZERO<RealT>};
      PwmT  pwm_;

      std::array<RealT, LEGS> gamma_{}; ///< Phase offset of each leg [rad]
    };
  } // namespace EMT
} // namespace GridKit
