/**
 * @file Ieeeg1.hpp
 * @brief IEEE type 1 steam turbine governor.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeeg1/Ieeeg1Data.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Ieeeg1InternalVariables : size_t
    {
      SPEED,
      VALVE,
      BOWL,
      REHEAT,
      CROSSOVER,
      REHEAT2,
      PM,
      PM_LP,
    };

    /// IEEE type 1 steam turbine governor
    template <typename ScalarT, typename IdxT>
    class Ieeeg1 final : public Component<ScalarT, IdxT>
    {
      using ComponentT = Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Ieeeg1Data<RealT, IdxT>;
      using InternalVariablesT = Ieeeg1InternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Ieeeg1Inputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view                                             CLASS_NAME = "Ieeeg1";
      static constexpr std::array<IdxT, Utilities::enum_size<InternalVariablesT>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>            NUM_INPUT_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_OUTPUT_VARIABLES{};

      explicit Ieeeg1(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;
      int                        evaluateResidual() override;
      int                        evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      enum : size_t
      {
        SPEED,
        VALVE,
        BOWL,
        REHEAT,
        CROSSOVER,
        REHEAT2,
        PM,
        PM_LP,
      };

      enum : size_t
      {
        OMEGA_IN = 0,
        PREF_IN  = 1,
        PAUX_IN  = 2,
      };

      ModelDataT           data_;
      RealT                to_machine_{}, lead_gain_{}, straight_valve_{};
      std::array<RealT, 5> slope_{};
    };
  } // namespace EMT
} // namespace GridKit
