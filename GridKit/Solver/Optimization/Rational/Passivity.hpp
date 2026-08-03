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

      /// One contiguous frequency band where passivity is violated.
      struct Violation
      {
        RealT omega_start{0.0}; ///< Lower band edge.
        RealT omega_end{0.0};   ///< Upper band edge.
      };

      bool passive{false};

      std::vector<Violation> violations;
    };

    /**
     * @brief Assess passivity of a fitted admittance model.
     *
     * Screens with a singular-value sweep and certifies violation band
     * edges with the Hamiltonian eigenvalue test.
     *
     * @return 0 on success, negative on an eigensolver failure
     */
    template <typename scalar_type, typename index_type>
    [[nodiscard("May fail. Check error code.")]]
    int assessPassivity(const RationalModel<scalar_type, index_type>& model,
                        PassivityReport<scalar_type, index_type>& report);
  } // namespace Optimization
} // namespace GridKit
