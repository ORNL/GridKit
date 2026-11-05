/**
 * @file VoltageSource.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT voltage source.
 */

#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `VoltageSource`, each of size N
    enum class VoltageSourceInternalVariables : size_t
    {
      u, ///< \f$\mathbf{u}\f$ Branch voltage e - v [V]
      i, ///< \f$\mathbf{i}\f$ Current injection into the bus [A]
    };

    /// Submodels of a `VoltageSource`
    enum class VoltageSourceSubmodels : size_t
    {
      y, ///< \f$\mathbf{y}\f$ Terminal admittance
    };

    /// N-phase sinusoidal source behind its terminal admittance
    template <typename scalar_type, typename index_type>
    class VoltageSource final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using VectorFitT = VectorFit<scalar_type, index_type>;

    public:
      using ModelDataT         = VoltageSourceData<RealT, IdxT>;
      using InternalVariablesT = VoltageSourceInternalVariables;
      using PortsT             = VoltageSourcePorts;
      using SubmodelsT         = VoltageSourceSubmodels;

      static constexpr std::string_view CLASS_NAME = "VoltageSource";

      explicit VoltageSource(const ModelDataT& data);

      VectorFitT& submodel(SubmodelsT symbol);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;

      template <typename T>
      void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      RealT source(IdxT n) const;

      IdxT               N_{0};
      std::vector<RealT> E_;
      std::vector<RealT> phi_;
      RealT              omega_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
