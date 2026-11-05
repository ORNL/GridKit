/**
 * @file PllImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT phase-locked loop.
 */

#pragma once

#include <cmath>
#include <numbers>
#include <type_traits>

#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/Pll.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Pll<scalar_type, index_type>::Pll(const ModelDataT& data)
      : ComponentT(std::type_identity<Pll>{}),
        Kp_(data.Kp),
        Ki_(data.Ki)
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Pll<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Pll");
      checks.check(Kp_ > ZERO<RealT>, "Kp must be positive");
      checks.check(Ki_ > ZERO<RealT>, "Ki must be positive");
      return checks;
    }

    /// Angle in radians, integrator state and frequency at the base frequency
    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      this->absoluteTolerance(PllInternalVariables::theta)[0] = rel_tol;
      this->absoluteTolerance(PllInternalVariables::x)[0]     = rel_tol * omega_base;
      this->absoluteTolerance(PllInternalVariables::omega)[0] = rel_tol * omega_base;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Pll<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto THETA = static_cast<size_t>(PllInternalVariables::theta);
      const auto X     = static_cast<size_t>(PllInternalVariables::x);
      const auto OMEGA = static_cast<size_t>(PllInternalVariables::omega);

      constexpr RealT TWO_THIRDS = TWO<RealT> / THREE<RealT>;
      constexpr RealT SHIFT      = TWO<RealT> * std::numbers::pi_v<RealT> / THREE<RealT>;

      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      const T vq = -TWO_THIRDS
                   * (w[0] * std::sin(y[THETA])
                      + w[1] * std::sin(y[THETA] - SHIFT)
                      + w[2] * std::sin(y[THETA] + SHIFT));

      f[THETA] = -yp[THETA] + y[OMEGA];
      f[X]     = -yp[X] + Ki_ * vq;
      f[OMEGA] = -y[OMEGA] + omega_base + Kp_ * vq + y[X];
    }
  } // namespace EMT
} // namespace GridKit
