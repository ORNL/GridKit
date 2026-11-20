/**
 * @file LVRT.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of EMT low-voltage ride-through.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/LVRT/LVRTData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of an `LVRT`, each of size 1 and in per unit
    enum class LVRTInternalVariables : size_t
    {
      il, ///< \f$i_l\f$ Recovering active-current ceiling
      ip, ///< \f$i_p\f$ Limited active-current order
    };

    /// Limit of an active-current order at low terminal voltage
    template <typename scalar_type, typename index_type>
    class LVRT final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = LVRTData<RealT, IdxT>;
      using InternalVariablesT = LVRTInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = LVRTInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "LVRT";

      static constexpr std::array<IdxT, Utilities::enum_size<LVRTInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>            NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<LVRTInputs>()>            NUM_INPUT_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>              NUM_OUTPUT_VARIABLES{};

      explicit LVRT(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      RealT TM_{ZERO<RealT>};
      RealT VL0_{ZERO<RealT>};
      RealT VL1_{ZERO<RealT>};
      RealT IL1_{ZERO<RealT>};
      RealT Rpmax_{ZERO<RealT>};
      RealT VA0_{ZERO<RealT>};
      RealT VA1_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
