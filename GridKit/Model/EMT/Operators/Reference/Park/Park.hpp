/**
 * @file Park.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT Park transform.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/ParkData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Park`, each of size 2
    enum class ParkInternalVariables : size_t
    {
      dq,   ///< \f$\mathbf{u}_\mathrm{dq}\f$ Input in rotating coordinates, d then q, in the units of the input
      axis, ///< \f$(c, s)\f$ Unit direction of the input in rotating coordinates
    };

    /// Three-phase quantities in rotating coordinates, with their direction
    template <typename scalar_type, typename index_type>
    class Park final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = ParkData<RealT, IdxT>;
      using InternalVariablesT = ParkInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = ParkInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Park";

      static constexpr std::array<IdxT, Utilities::enum_size<ParkInternalVariables>()> NUM_INTERNAL_VARIABLES{2, 2};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>            NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<ParkInputs>()>            NUM_INPUT_VARIABLES{3, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>              NUM_OUTPUT_VARIABLES{};

      explicit Park(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      RealT base_{ONE<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
