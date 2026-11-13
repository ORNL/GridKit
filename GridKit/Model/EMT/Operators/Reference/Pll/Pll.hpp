/**
 * @file Pll.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT phase-locked loop.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/PllData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Pll`, each of size 1
    enum class PllInternalVariables : size_t
    {
      theta, ///< \f$\theta\f$ Angle [rad]
      x,     ///< \f$x\f$ Integrator state [rad/s]
      omega, ///< \f$\omega\f$ Frequency [rad/s]
    };

    /// Synchronous-reference-frame PLL on a three-phase bus
    template <typename scalar_type, typename index_type>
    class Pll final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = PllData<RealT, IdxT>;
      using InternalVariablesT = PllInternalVariables;
      using PortsT             = PllPorts;

      static constexpr std::string_view CLASS_NAME = "Pll";

      static constexpr std::array<IdxT, Utilities::enum_size<PllInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<PllPorts>()>             NUM_PORT_VARIABLES{3};

      explicit Pll(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void externalResidual(const T*, const T*, const T*, const T*, T*) const
      {
      }

    private:
      RealT Kp_{ZERO<RealT>};
      RealT Ki_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
