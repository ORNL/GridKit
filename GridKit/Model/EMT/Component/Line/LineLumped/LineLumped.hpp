/**
 * @file LineLumped.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT lumped line.
 */

#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LineLumped`, each of size K
    enum class LineLumpedInternalVariables : size_t
    {
      i12,   ///< \f$\mathbf{i}_{12}\f$ Series current from terminal 1 to terminal 2 [A]
      i_sh1, ///< \f$\mathbf{i}_1^\mathrm{sh}\f$ Shunt current at terminal 1 [A]
      i_sh2, ///< \f$\mathbf{i}_2^\mathrm{sh}\f$ Shunt current at terminal 2 [A]
    };

    /// Submodels of a `LineLumped`
    enum class LineLumpedSubmodels : size_t
    {
      z,  ///< \f$\mathbf{z}\f$ Series impedance per unit length
      y1, ///< \f$\mathbf{y}_1\f$ Shunt admittance per unit length at terminal 1
      y2, ///< \f$\mathbf{y}_2\f$ Shunt admittance per unit length at terminal 2
    };

    /// N-phase, K-conductor lumped line, with conductor k on phase c_k
    template <typename scalar_type, typename index_type>
    class LineLumped final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using VectorFitT = VectorFit<scalar_type, index_type>;

    public:
      using ModelDataT         = LineLumpedData<RealT, IdxT>;
      using InternalVariablesT = LineLumpedInternalVariables;
      using PortsT             = LineLumpedPorts;
      using SubmodelsT         = LineLumpedSubmodels;

      static constexpr std::string_view CLASS_NAME = "LineLumped";

      explicit LineLumped(const ModelDataT& data);

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
      IdxT              N_{0};
      IdxT              K_{0};
      std::vector<IdxT> c_; ///< Phase of each conductor, from 0
      RealT             dx_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
