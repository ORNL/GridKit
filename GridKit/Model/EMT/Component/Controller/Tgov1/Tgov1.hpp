/**
 * @file Tgov1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT steam turbine governor.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Tgov1/Tgov1Data.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Tgov1`, each of size 1 and in per unit
    enum class Tgov1InternalVariables : size_t
    {
      PTX, ///< \f$P_\mathrm{t}\f$ Turbine-block output on the turbine rating
      PV,  ///< \f$P_\mathrm{v}\f$ Valve position on the turbine rating
      PM,  ///< \f$P_\mathrm{m}\f$ Mechanical power on the machine rating
    };

    /// Steam turbine governor
    template <typename scalar_type, typename index_type>
    class Tgov1 final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = Tgov1Data<RealT, IdxT>;
      using InternalVariablesT = Tgov1InternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Tgov1Inputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Tgov1";

      static constexpr std::array<IdxT, Utilities::enum_size<Tgov1InternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>             NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<Tgov1Inputs>()>            NUM_INPUT_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit Tgov1(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in Tgov1InternalVariables order
      enum : size_t
      {
        PTX,
        PV,
        PM,
      };

      // Offsets of the external values, in Tgov1Inputs order
      enum : size_t
      {
        OMEGA,
        PREF,
      };

      static constexpr RealT TIME_CONSTANT_MINIMUM{1.0e-3}; ///< Floor of `T1` and `T3` [s]

      RealT Trate_{ZERO<RealT>};
      RealT R_{ZERO<RealT>};
      RealT Pvmin_{ZERO<RealT>};
      RealT Pvmax_{ZERO<RealT>};
      RealT T1_{ZERO<RealT>};
      RealT T2_{ZERO<RealT>};
      RealT T3_{ZERO<RealT>};
      RealT Dt_{ZERO<RealT>};
      RealT S_{ZERO<RealT>};

      RealT machine_to_turbine_{ZERO<RealT>}; ///< Machine rating over turbine rating
    };
  } // namespace EMT
} // namespace GridKit
