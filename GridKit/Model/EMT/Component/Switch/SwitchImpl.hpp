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
    /**
     * @brief Construct an ideal switch with `N` phases.
     *
     * Layout: internal variable `i12` of size `N`, terminals `bus1` and
     * `bus2` of size `N`, and input `open` of size 1.
     */
    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>::Switch(const ModelDataT& data)
      : N_(data.N)
    {
      this->setLayout({N_}, {N_, N_}, {N_}, {});
    }

    /**
     * @brief Check that `N` is positive.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Switch<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Switch");
      checks.check(N_ > 0, "N must be positive");
      return checks;
    }

    /**
     * @brief Set the current tolerances to `rel_tol` times the system power
     *        base over the base voltage of terminal 1.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The bus at `bus1` has set its tolerances with the same `rel_tol`.
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i12  = this->absoluteTolerance(SwitchInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(SwitchTerminals::bus1);
      for (IdxT n = 0; n < N_; ++n)
      {
        const RealT v_base = this->externalTolerance(BUS1 + n) / rel_tol;
        i12[n]             = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    /**
     * @brief Write the contact equation of each phase.
     *
     * An open contact carries no current and a closed contact no voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void Switch<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I12  = this->local(SwitchInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(SwitchTerminals::bus1);
      const IdxT BUS2 = this->terminalOffset(SwitchTerminals::bus2);
      const IdxT OPEN = this->inputOffset(SwitchInputs::open);

      for (IdxT n = 0; n < N_; ++n)
      {
        const T open = w[OPEN + n];
        f[I12 + n]   = open * y[I12 + n] + (ONE<RealT> - open) * (w[BUS2 + n] - w[BUS1 + n]);
      }
    }

    /**
     * @brief Write the series current to the residuals of the bus voltages.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void Switch<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I12  = this->local(SwitchInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(SwitchTerminals::bus1);
      const IdxT BUS2 = this->terminalOffset(SwitchTerminals::bus2);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS1 + n] = -y[I12 + n];
        f[BUS2 + n] = y[I12 + n];
      }
    }
  } // namespace EMT
} // namespace GridKit
