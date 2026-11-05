/**
 * @file Switch.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT switch.
 */

#pragma once

#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Switch`, each of size N
    enum class SwitchInternalVariables : size_t
    {
      i12, ///< \f$\mathbf{i}_{12}\f$ Series current from terminal 1 to terminal 2 [A]
    };

    /// Ideal N-phase switch between two buses, opened and closed by events
    template <typename scalar_type, typename index_type>
    class Switch final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = SwitchData<RealT, IdxT>;
      using InternalVariablesT = SwitchInternalVariables;
      using PortsT             = SwitchPorts;

      static constexpr std::string_view CLASS_NAME = "Switch";

      explicit Switch(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;

      /// Opens or closes every phase, with one Jacobian pattern for both
      void setOpen(bool open);

      template <typename T>
      void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      IdxT  N_{0};
      RealT open_{ZERO<RealT>};
      RealT closed_{ONE<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
