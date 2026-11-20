/**
 * @file GastPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT gas turbine governor.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/GastPti/GastPtiData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `GastPti`, each of size 1 and in per unit
    enum class GastPtiInternalVariables : size_t
    {
      XVALVE, ///< \f$x_V\f$ Fuel-valve state on the turbine rating
      XFLOW,  ///< \f$x_F\f$ Fuel-flow state on the turbine rating
      XTEMP,  ///< \f$x_T\f$ Exhaust-temperature feedback state on the turbine rating
      VLOAD,  ///< \f$V_D\f$ Speed and load fuel demand on the turbine rating
      VTEMP,  ///< \f$V_T\f$ Temperature-limit demand on the turbine rating
      VLV,    ///< \f$V\f$ Smooth low-value selector output on the turbine rating
      PMECH,  ///< \f$P_\mathrm{m}\f$ Mechanical power on the machine rating
    };

    /// Gas turbine governor
    template <typename scalar_type, typename index_type>
    class GastPti final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = GastPtiData<RealT, IdxT>;
      using InternalVariablesT = GastPtiInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = GastPtiInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "GastPti";

      static constexpr std::array<IdxT, Utilities::enum_size<GastPtiInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>               NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<GastPtiInputs>()>            NUM_INPUT_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                 NUM_OUTPUT_VARIABLES{};

      explicit GastPti(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in GastPtiInternalVariables order
      enum : size_t
      {
        XVALVE,
        XFLOW,
        XTEMP,
        VLOAD,
        VTEMP,
        VLV,
        PMECH,
      };

      // Offsets of the external values, in GastPtiInputs order
      enum : size_t
      {
        OMEGA,
        PREF,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `T1`, `T2` and `T3` [s]

      RealT R_{ZERO<RealT>};
      RealT T1_{ZERO<RealT>};
      RealT T2_{ZERO<RealT>};
      RealT T3_{ZERO<RealT>};
      RealT At_{ZERO<RealT>};
      RealT Kt_{ZERO<RealT>};
      RealT Vmax_{ZERO<RealT>};
      RealT Vmin_{ZERO<RealT>};
      RealT Dturb_{ZERO<RealT>};
      RealT Trate_{ZERO<RealT>};
      RealT S_{ZERO<RealT>};

      RealT machine_to_turbine_{ZERO<RealT>}; ///< Machine rating over turbine rating
      RealT valve_on_{ZERO<RealT>};           ///< One unless the valve limits coincide
    };
  } // namespace EMT
} // namespace GridKit
