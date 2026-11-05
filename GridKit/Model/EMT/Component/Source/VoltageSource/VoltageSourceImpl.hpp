/**
 * @file VoltageSourceImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT voltage source.
 */

#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSource.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    VoltageSource<scalar_type, index_type>::VoltageSource(const ModelDataT& data)
      : N_(data.N),
        E_(data.E),
        phi_(data.phi),
        omega_(data.omega)
    {
      if (phi_.empty())
      {
        phi_.assign(N_, ZERO<RealT>);
      }
      this->setLayout({N_, N_}, {N_});

      auto& y = this->template addSubmodel<VectorFitT>(VoltageSourceSubmodels::y, data.Y, N_, N_);
      if (y.connect(VectorFitPorts::input, *this, this->local(VoltageSourceInternalVariables::u)) != 0
          || y.connect(VectorFitPorts::out, *this, this->local(VoltageSourceInternalVariables::i)) != 0)
      {
        throw std::logic_error("VoltageSource: cannot wire the admittance");
      }
    }

    template <typename scalar_type, typename index_type>
    typename VoltageSource<scalar_type, index_type>::VectorFitT& VoltageSource<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<VectorFitT>(symbol);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks VoltageSource<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("VoltageSource");
      checks.check(N_ > 0, "N must be positive");
      checks.check(E_.size() == N_ && phi_.size() == N_, "E and phi must have N values");
      checks.check(omega_ > ZERO<RealT>, "omega must be positive");
      for (const RealT magnitude : E_)
      {
        checks.check(magnitude >= ZERO<RealT>, "E must be nonnegative");
      }
      return checks;
    }

    /// Voltages at the bus voltage tolerance, and currents at the system power base over the bus base voltage
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto u = this->absoluteTolerance(VoltageSourceInternalVariables::u);
      const auto i = this->absoluteTolerance(VoltageSourceInternalVariables::i);
      for (IdxT n = 0; n < N_; ++n)
      {
        const RealT v_tol  = this->absoluteTolerance(VoltageSourcePorts::bus, n);
        const RealT v_base = v_tol / rel_tol;
        u[n]               = v_tol;
        i[n]               = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void VoltageSource<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const IdxT U   = this->local(VoltageSourceInternalVariables::u);
      const IdxT I   = this->local(VoltageSourceInternalVariables::i);
      const IdxT BUS = this->external(VoltageSourcePorts::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[U + n] = -y[U + n] + source(n) - w[BUS + n];
        f[I + n] = -y[I + n];
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void VoltageSource<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I   = this->local(VoltageSourceInternalVariables::i);
      const IdxT BUS = this->external(VoltageSourcePorts::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS + n] = y[I + n];
      }
    }

    template <typename scalar_type, typename index_type>
    typename VoltageSource<scalar_type, index_type>::RealT VoltageSource<scalar_type, index_type>::source(IdxT n) const
    {
      return std::numbers::sqrt2_v<RealT> * E_[n] * std::cos(omega_ * this->time_ + phi_[n]);
    }
  } // namespace EMT
} // namespace GridKit
