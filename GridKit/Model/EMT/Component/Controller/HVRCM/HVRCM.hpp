/**
 * @file HVRCM.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of EMT high-voltage reactive current management.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/HVRCM/HVRCMData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of an `HVRCM`, each of size 1 and in per unit
    enum class HVRCMInternalVariables : size_t
    {
      iq, ///< \f$i_q\f$ Reactive-current order after high-voltage management
    };

    /// Reduction of a reactive-current order at high terminal voltage
    template <typename scalar_type, typename index_type>
    class HVRCM final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = HVRCMData<RealT, IdxT>;
      using InternalVariablesT = HVRCMInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = HVRCMInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "HVRCM";

      static constexpr std::array<IdxT, Utilities::enum_size<HVRCMInternalVariables>()> NUM_INTERNAL_VARIABLES{1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>             NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<HVRCMInputs>()>            NUM_INPUT_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit HVRCM(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      RealT Vhvmax_{ZERO<RealT>};
      RealT Khv_{ZERO<RealT>};
      RealT Qmin_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
