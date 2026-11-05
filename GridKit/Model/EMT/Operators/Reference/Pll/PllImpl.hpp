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
#include <GridKit/Model/EMT/Operators/Reference/ThreePhase.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a PLL on a three-phase voltage.
     *
     * Layout: internal variables `theta` and `x` of size 1, input `v` of
     * size 3, and output `omega` of size 1.
     */
    template <typename scalar_type, typename index_type>
    Pll<scalar_type, index_type>::Pll(const ModelDataT& data)
      : ComponentT(std::type_identity<Pll>{}),
        Kp_(data.Kp),
        Ki_(data.Ki)
    {
    }

    /**
     * @brief Check that the gains are positive.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Pll<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Pll");
      checks.check(Kp_ > ZERO<RealT>, "Kp must be positive");
      checks.check(Ki_ > ZERO<RealT>, "Ki must be positive");
      return checks;
    }

    /**
     * @brief Set the tolerance of the angle to `rel_tol` radians, and those
     *        of the integrator state and of the frequency to `rel_tol` times
     *        the system base frequency.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     */
    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      this->absoluteTolerance(PllInternalVariables::theta)[0] = rel_tol;
      this->absoluteTolerance(PllInternalVariables::x)[0]     = rel_tol * omega_base;
      this->outputTolerance(PllOutputs::omega)[0]             = rel_tol * omega_base;
      return 0;
    }

    /**
     * @brief Write the residuals of the angle and of the integrator state.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Pll<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto THETA = static_cast<size_t>(PllInternalVariables::theta);
      const auto X     = static_cast<size_t>(PllInternalVariables::x);

      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;
      const T     vq         = park(w, y[THETA]).q;

      f[THETA] = -yp[THETA] + omega_base + Kp_ * vq + y[X];
      f[X]     = -yp[X] + Ki_ * vq;
    }

    /**
     * @brief Write the frequency `omega`.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Pll<scalar_type, index_type>::output(OutputsT, const T* y, const T*, const T* w, const T*, T* f) const
    {
      const auto THETA = static_cast<size_t>(PllInternalVariables::theta);
      const auto X     = static_cast<size_t>(PllInternalVariables::x);

      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      f[0] = omega_base + Kp_ * park(w, y[THETA]).q + y[X];
    }
  } // namespace EMT
} // namespace GridKit
