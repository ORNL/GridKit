/**
 * @file SexsPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT simplified excitation system.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/SexsPti/SexsPtiData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `SexsPti`, each of size 1 and in per unit
    enum class SexsPtiInternalVariables : size_t
    {
      VR,  ///< \f$V_R\f$ Lead-lag block state
      EFD, ///< \f$E_{\mathrm{fd}}\f$ Field voltage
      VTR, ///< \f$V_{\mathrm{tr}}\f$ Terminal voltage error
    };

    /// Simplified excitation system
    template <typename scalar_type, typename index_type>
    class SexsPti final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = SexsPtiData<RealT, IdxT>;
      using InternalVariablesT = SexsPtiInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = SexsPtiInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "SexsPti";

      static constexpr std::array<IdxT, Utilities::enum_size<SexsPtiInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>               NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<SexsPtiInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                 NUM_OUTPUT_VARIABLES{};

      explicit SexsPti(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in SexsPtiInternalVariables order
      enum : size_t
      {
        VR,
        EFD,
        VTR,
      };

      // Offsets of the external values, in SexsPtiInputs order
      enum : size_t
      {
        VT,
        VREF,
        VS,
        VUEL,
        VOEL,
      };

      RealT Ta_{ZERO<RealT>};
      RealT Tb_{ZERO<RealT>};
      RealT Te_{ZERO<RealT>};
      RealT K_{ZERO<RealT>};
      RealT Efdmax_{ZERO<RealT>};
      RealT Efdmin_{ZERO<RealT>};
      RealT V_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
