/**
 * @file LoadZImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT impedance load.
 */

#pragma once

#include <stdexcept>

#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(const ModelDataT& data)
      : N_(data.N)
    {
      this->setLayout({N_}, {N_});

      const IdxT I = this->local(LoadZInternalVariables::i);
      auto&      z = this->template addSubmodel<VectorFitT>(LoadZSubmodels::z, data.Z, N_, N_);
      if (z.connect(VectorFitPorts::input, *this, I) != 0 || z.connect(VectorFitPorts::out, *this, I) != 0)
      {
        throw std::logic_error("LoadZ: cannot wire the impedance");
      }
    }

    template <typename scalar_type, typename index_type>
    typename LoadZ<scalar_type, index_type>::VectorFitT& LoadZ<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<VectorFitT>(symbol);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LoadZ<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("LoadZ");
      checks.check(N_ > 0, "N must be positive");
      return checks;
    }

    /// Currents at the system power base over the bus base voltage
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i = this->absoluteTolerance(LoadZInternalVariables::i);
      for (IdxT n = 0; n < N_; ++n)
      {
        const RealT v_base = this->absoluteTolerance(LoadZPorts::bus, n) / rel_tol;
        i[n]               = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadZ<scalar_type, index_type>::internalResidual(const T*, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I   = this->local(LoadZInternalVariables::i);
      const IdxT BUS = this->external(LoadZPorts::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[I + n] = w[BUS + n];
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadZ<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I   = this->local(LoadZInternalVariables::i);
      const IdxT BUS = this->external(LoadZPorts::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS + n] = y[I + n];
      }
    }
  } // namespace EMT
} // namespace GridKit
