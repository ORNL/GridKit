/**
 * @file SwitchImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT switch.
 */

#pragma once

#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>::Switch(const ModelDataT& data)
      : N_(data.N)
    {
      setOpen(data.open);
      this->setLayout({N_}, {N_, N_});
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Switch<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Switch");
      checks.check(N_ > 0, "N must be positive");
      return checks;
    }

    /// Currents at the system power base over the base voltage of terminal 1
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i12 = this->absoluteTolerance(SwitchInternalVariables::i12);
      for (IdxT n = 0; n < N_; ++n)
      {
        const RealT v_base = this->absoluteTolerance(SwitchPorts::bus1, n) / rel_tol;
        i12[n]             = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Switch<scalar_type, index_type>::setOpen(bool open)
    {
      open_   = ZERO<RealT>;
      closed_ = ONE<RealT>;
      if (open)
      {
        open_   = ONE<RealT>;
        closed_ = ZERO<RealT>;
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void Switch<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I12  = this->local(SwitchInternalVariables::i12);
      const IdxT BUS1 = this->external(SwitchPorts::bus1);
      const IdxT BUS2 = this->external(SwitchPorts::bus2);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[I12 + n] = open_ * y[I12 + n] + closed_ * (w[BUS2 + n] - w[BUS1 + n]);
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void Switch<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I12  = this->local(SwitchInternalVariables::i12);
      const IdxT BUS1 = this->external(SwitchPorts::bus1);
      const IdxT BUS2 = this->external(SwitchPorts::bus2);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS1 + n] = -y[I12 + n];
        f[BUS2 + n] = y[I12 + n];
      }
    }
  } // namespace EMT
} // namespace GridKit
