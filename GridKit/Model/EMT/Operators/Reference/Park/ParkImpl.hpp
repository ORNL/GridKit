/**
 * @file ParkImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT Park transform.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/Park.hpp>
#include <GridKit/Model/EMT/Operators/Reference/ThreePhase.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a Park transform.
     *
     * Layout: internal variables `dq` and `axis` of size 2, input `abc` of
     * size 3 and input `theta` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Park<scalar_type, index_type>::Park(const ModelDataT& data)
      : ComponentT(std::type_identity<Park>{}),
        base_(data.base)
    {
    }

    /**
     * @brief Check that `base` is positive and finite.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Park<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Park");
      checks.check(base_ > ZERO<RealT> && std::isfinite(base_), "base must be positive and finite");
      return checks;
    }

    /**
     * @brief Set the tolerance of `dq` to `rel_tol` times the rated peak
     *        magnitude, and that of `axis` to `rel_tol`.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     */
    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::ranges::fill(this->absoluteTolerance(ParkInternalVariables::dq), rel_tol * base_);
      std::ranges::fill(this->absoluteTolerance(ParkInternalVariables::axis), rel_tol);
      return 0;
    }

    /**
     * @brief Write the residuals of the rotating coordinates and of their
     *        direction.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Park<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      // Local indices of `dq` and `axis`, each d then q
      constexpr size_t D      = 0;
      constexpr size_t Q      = 1;
      constexpr size_t AXIS_D = 2;
      constexpr size_t AXIS_Q = 3;

      // Offset of input `theta`, after the three phases of `abc`
      constexpr size_t THETA = 3;

      // Squared magnitude, in per unit of the rated peak, at or below which the d axis gives the direction
      constexpr RealT MIN_SQUARED = 1.0e-12;

      const auto dq = park(w, w[THETA]);
      f[D]          = -y[D] + dq.d;
      f[Q]          = -y[Q] + dq.q;

      const T d       = y[D] / base_;
      const T q       = y[Q] / base_;
      const T squared = d * d + q * q;

      // The mask keeps the dependencies on `dq` where the direction is not measured
      T measured = ZERO<RealT>;
      if (squared > MIN_SQUARED)
      {
        measured = ONE<RealT>;
      }
      const T inverse = measured / std::sqrt(std::max(squared, T{MIN_SQUARED}));

      f[AXIS_D] = -y[AXIS_D] + (ONE<RealT> - measured) + inverse * d;
      f[AXIS_Q] = -y[AXIS_Q] + inverse * q;
    }
  } // namespace EMT
} // namespace GridKit
