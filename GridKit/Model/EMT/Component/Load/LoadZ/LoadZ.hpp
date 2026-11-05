/**
 * @file LoadZ.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT impedance load.
 */

#pragma once

#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LoadZ`, each of size N
    enum class LoadZInternalVariables : size_t
    {
      i, ///< \f$\mathbf{i}\f$ Current injection into the bus [A]
    };

    /// Submodels of a `LoadZ`
    enum class LoadZSubmodels : size_t
    {
      z, ///< \f$\mathbf{z}\f$ Impedance
    };

    /// N-phase impedance load
    template <typename scalar_type, typename index_type>
    class LoadZ final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using VectorFitT = VectorFit<scalar_type, index_type>;

    public:
      using ModelDataT         = LoadZData<RealT, IdxT>;
      using InternalVariablesT = LoadZInternalVariables;
      using PortsT             = LoadZPorts;
      using SubmodelsT         = LoadZSubmodels;

      static constexpr std::string_view CLASS_NAME = "LoadZ";

      explicit LoadZ(const ModelDataT& data);

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
      IdxT N_{0};
    };
  } // namespace EMT
} // namespace GridKit
