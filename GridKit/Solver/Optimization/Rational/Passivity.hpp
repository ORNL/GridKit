/**
 * @file Passivity.hpp
 *
 * @brief Passivity assessment of fitted rational models.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <vector>

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /// Outcome of a passivity assessment.
    template <typename scalar_type, typename index_type>
    struct PassivityReport
    {
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      /// One contiguous frequency band where passivity is violated. A
      /// zero start is the DC edge itself; an infinite end marks an
      /// asymptotic violation carried by the constant term.
      struct Violation
      {
        RealT omega_start{0.0}; ///< Lower band edge.
        RealT omega_end{0.0};   ///< Upper band edge.
      };

      bool passive{false};

      /// Every pole strictly inside the left half plane. An unstable
      /// model is never passive and its conductance screen is skipped.
      bool stable{false};

      std::vector<Violation> violations;
    };

    /**
     * @brief Assess passivity of a fitted admittance model.
     *
     * Gates on stability, then screens the conductance eigenvalue floor
     * at DC, over a logarithmic grid spanning the union of the fitted
     * band and the pole magnitudes, and at the asymptotic limit through
     * the constant term. Violation band edges inside the grid are
     * refined by bisection.
     *
     * @param[in] omega_low  Lower edge of the fitted band, positive
     * @param[in] omega_high Upper edge of the fitted band, above the
     *                       lower edge
     *
     * @return 0 on success, negative on invalid inputs or an
     *         eigensolver failure
     */
    template <typename scalar_type, typename index_type>
    [[nodiscard("May fail. Check error code.")]]
    int assessPassivity(
        const RationalModel<scalar_type, index_type>&      model,
        PassivityReport<scalar_type, index_type>&          report,
        typename GridKit::ScalarTraits<scalar_type>::RealT omega_low,
        typename GridKit::ScalarTraits<scalar_type>::RealT omega_high);
  } // namespace Optimization
} // namespace GridKit
