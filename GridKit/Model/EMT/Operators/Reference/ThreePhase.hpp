/**
 * @file ThreePhase.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Transforms of three-phase quantities to stationary and rotating coordinates.
 */

#pragma once

#include <cmath>
#include <numbers>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Projections of three-phase quantities on the d and q axes
    template <typename T>
    struct DQ
    {
      T d;
      T q;
    };

    /// Projections of three-phase quantities on the stationary alpha and beta axes
    template <typename T>
    struct AlphaBeta
    {
      T alpha;
      T beta;
    };

    /**
     * @brief alpha and beta projections of three-phase quantities.
     *
     * The transform keeps the amplitude of a balanced set, and the alpha axis
     * is the axis of phase a.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] abc Phases a, b and c.
     * @return The projections.
     */
    template <typename T>
    [[gnu::always_inline]] inline AlphaBeta<T> clarke(const T* abc)
    {
      using RealT = typename ScalarTraits<T>::RealT;

      return {(TWO<RealT> * abc[0] - abc[1] - abc[2]) / THREE<RealT>,
              (abc[1] - abc[2]) / std::numbers::sqrt3_v<RealT>};
    }

    /**
     * @brief d and q projections of three-phase quantities at angle `theta`.
     *
     * The transform keeps the amplitude of a balanced set, and q leads d by
     * pi/2.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] abc Phases a, b and c.
     * @param[in] theta Angle of the d axis from the axis of phase a [rad].
     * @return The projections.
     */
    template <typename T>
    [[gnu::always_inline]] inline DQ<T> park(const T* abc, T theta)
    {
      const AlphaBeta<T> u = clarke(abc);
      const T            c = std::cos(theta);
      const T            s = std::sin(theta);
      return {u.alpha * c + u.beta * s, u.beta * c - u.alpha * s};
    }
  } // namespace EMT
} // namespace GridKit
