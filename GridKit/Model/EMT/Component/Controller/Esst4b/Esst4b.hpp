/**
 * @file Esst4b.hpp
 * @brief IEEE ST4B static excitation system.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esst4b/Esst4bData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Esst4bInternalVariables : size_t
    {
      VC,
      XVR,
      VA,
      XVM,
      EFD,
    };

    /// IEEE ST4B static excitation system
    template <typename ScalarT, typename IdxT>
    class Esst4b final : public Component<ScalarT, IdxT>
    {
      using ComponentT = Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Esst4bData<RealT, IdxT>;
      using InternalVariablesT = Esst4bInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Esst4bInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view                                             CLASS_NAME = "Esst4b";
      static constexpr std::array<IdxT, Utilities::enum_size<InternalVariablesT>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>            NUM_INPUT_VARIABLES{1, 3, 3, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_OUTPUT_VARIABLES{};

      explicit Esst4b(const ModelDataT& data);

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
        XVR,
        VA,
        XVM,
        EFD,
      };

      enum : size_t
      {
        VT_IN   = 0,
        V_IN    = 1,
        I_IN    = 4,
        IFD_IN  = 7,
        VREF_IN = 8,
        VS_IN   = 9,
        VUEL_IN = 10,
        VOEL_IN = 11,
      };

      ModelDataT data_;
      RealT      inv_v_peak_{}, inv_i_peak_{}, kp_real_{}, kp_imag_{};
    };
  } // namespace EMT
} // namespace GridKit
