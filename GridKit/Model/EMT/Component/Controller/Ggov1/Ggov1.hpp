/**
 * @file Ggov1.hpp
 * @brief General turbine governor.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ggov1/Ggov1Data.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Ggov1InternalVariables : size_t
    {
      PE,
      XD,
      XI,
      VALVE,
      TURB,
      TEMP,
      LOAD,
      XL,
      XMW,
      SPEED,
      DEMAND,
      FSR,
      PM,
      ERROR,
      MEMORY,
      FSRN,
      FSRT,
      MEMORY_DEV, ///< Differential fuel-memory deviation from the initial MEMORY.
    };

    /// General turbine governor
    template <typename ScalarT, typename IdxT>
    class Ggov1 final : public Component<ScalarT, IdxT>
    {
      using ComponentT = Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Ggov1Data<RealT, IdxT>;
      using InternalVariablesT = Ggov1InternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Ggov1Inputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view                                             CLASS_NAME = "Ggov1";
      static constexpr std::array<IdxT, Utilities::enum_size<InternalVariablesT>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_OUTPUT_VARIABLES{};

      explicit Ggov1(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;
      int                        tagDifferentiable() override;
      int                        evaluateResidual() override;
      int                        evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      enum : size_t
      {
        PE,
        XD,
        XI,
        VALVE,
        TURB,
        TEMP,
        LOAD,
        XL,
        XMW,
        SPEED,
        DEMAND,
        FSR,
        PM,
        ERROR,
        MEMORY,
        FSRN,
        FSRT,
        MEMORY_DEV,
      };

      enum : size_t
      {
        OMEGA_IN  = 0,
        P_IN      = 1,
        PREF_IN   = 2,
        PAUX_IN   = 3,
        PMWSET_IN = 4,
      };

      ModelDataT data_;
      RealT      initial_memory_{};
      RealT      inv_rating_{}, to_machine_{};
      RealT      governor_tracking_, load_tracking_, temperature_gain_, turbine_gain_, acceleration_gain_, derivative_gain_;
    };
  } // namespace EMT
} // namespace GridKit
