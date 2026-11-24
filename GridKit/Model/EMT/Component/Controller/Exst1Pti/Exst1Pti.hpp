/**
 * @file Exst1Pti.hpp
 * @brief PTI static excitation system.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Exst1Pti/Exst1PtiData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Exst1PtiInternalVariables : size_t
    {
      VC,
      LL,
      VA,
      WASH,
      EFD,
    };

    /// PTI static excitation system
    template <typename ScalarT, typename IdxT>
    class Exst1Pti final : public Component<ScalarT, IdxT>
    {
      using ComponentT = Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Exst1PtiData<RealT, IdxT>;
      using InternalVariablesT = Exst1PtiInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Exst1PtiInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view                                             CLASS_NAME = "Exst1Pti";
      static constexpr std::array<IdxT, Utilities::enum_size<InternalVariablesT>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_OUTPUT_VARIABLES{};

      explicit Exst1Pti(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;
      int                        evaluateResidual() override;
      int                        evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      enum : size_t
      {
        VC,
        LL,
        VA,
        WASH,
        EFD,
      };

      enum : size_t
      {
        VT_IN   = 0,
        IFD_IN  = 1,
        VREF_IN = 2,
        VS_IN   = 3,
        VUEL_IN = 4,
        VOEL_IN = 5,
      };

      ModelDataT data_;
      RealT      washout_gain_, lead_gain_;
    };
  } // namespace EMT
} // namespace GridKit
