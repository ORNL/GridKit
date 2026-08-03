```cpp
/**
 * @file RationalModel.hpp
 *
 * @brief Rational approximation of a sampled frequency response.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <complex>
#include <vector>

#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Rational model H(s) = D + sE + sum of R_q / (s - p_q).
     *
     * Coefficient conventions and conjugate-pair ordering follow the
     * VectorFit model contract: real poles carry real residues, and every
     * nonreal pole and residue is immediately followed by its conjugate.
     * Fitted instances satisfy the contract by construction; validate()
     * gates instances arriving from any other source.
     */
    template <typename scalar_type, typename index_type>
    struct RationalModel
    {
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ComplexT = std::complex<RealT>;

      std::vector<ComplexT> poles;    ///< Q poles, conjugate pairs adjacent.
      std::vector<ComplexT> residues; ///< Q * rows * cols, pole-major.
      std::vector<RealT>    d;        ///< rows * cols, empty when absent.
      std::vector<RealT>    e;        ///< rows * cols, empty when absent.
      IdxT                  rows{1};  ///< Output dimension N.
      IdxT                  cols{1};  ///< Input dimension K.

      /**
       * @brief One residue matrix element for pole q.
       */
      ComplexT residue(IdxT q, IdxT row, IdxT col) const;

      /**
       * @brief Evaluate one matrix element of the model at j omega.
       */
      ComplexT evaluate(RealT omega, IdxT row, IdxT col) const;

      /**
       * @brief Whether every pole satisfies Re(p) <= -margin.
       */
      bool isStable(RealT margin = 0.0) const;

      /**
       * @brief Check the coefficient contract.
       *
       * @param[in] tolerance Relative comparison tolerance; zero demands
       *                      exact conjugacy
       * @return 0 when valid, -1 on inconsistent sizes, -2 on a nonfinite
       *         value, -3 on a conjugate-pair violation
       *
       * @pre tolerance >= 0
       */
      [[nodiscard("May fail. Check error code.")]]
      int validate(RealT tolerance = 0.0) const;
    };
  } // namespace Optimization
} // namespace GridKit
```
