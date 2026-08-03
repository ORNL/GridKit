/**
 * @file StateSpaceRealization.hpp
 *
 * @brief Exact-rank factorized realization of a rational model.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <complex>
#include <vector>

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Factorized form H(s) = D + sE + C (sI - P)^-1 B.
     *
     * Each residue matrix is factored at its numerical rank, so the
     * realization is exact and reduces to the rank-one form for scalar and
     * vector responses. Matches the StateSpace model coefficient contract.
     */
    template <typename scalar_type, typename index_type>
    struct StateSpaceRealization
    {
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ComplexT = std::complex<RealT>;

      std::vector<ComplexT> poles; ///< One entry per state.
      std::vector<ComplexT> c;     ///< rows x states, row-major.
      std::vector<ComplexT> b;     ///< states x cols, row-major.
      std::vector<RealT>    d;     ///< rows * cols, empty when absent.
      std::vector<RealT>    e;     ///< rows * cols, empty when absent.
      IdxT                  states{0};
      IdxT                  rows{1};
      IdxT                  cols{1};
    };

    /**
     * @brief Factor a rational model into its state-space realization.
     *
     * @return 0 on success, negative when a residue factorization fails
     */
    template <typename scalar_type, typename index_type>
    [[nodiscard("May fail. Check error code.")]]
    int toStateSpace(const RationalModel<scalar_type, index_type>& model,
                     StateSpaceRealization<scalar_type, index_type>& realization);
  } // namespace Optimization
} // namespace GridKit
