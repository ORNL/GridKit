/**
 * @file SampledResponse.hpp
 *
 * @brief Sampled frequency-response carrier for rational approximation.
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
     * @brief Sampled frequency response H(j omega_m).
     *
     * Angular frequencies are strictly increasing and positive. The
     * response stores omega.size() * rows * cols entries, sample-major,
     * row-major within each sample; use the element accessor rather than
     * indexing the flat storage.
     */
    template <typename scalar_type, typename index_type>
    struct SampledResponse
    {
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ComplexT = std::complex<RealT>;

      std::vector<RealT>    omega;    ///< Angular frequency samples.
      std::vector<ComplexT> response; ///< Sampled response values.
      IdxT                  rows{1};  ///< Output dimension N.
      IdxT                  cols{1};  ///< Input dimension K.

      IdxT sampleCount() const
      {
        return static_cast<IdxT>(omega.size());
      }

      ComplexT operator()(IdxT sample, IdxT row, IdxT col) const
      {
        return response[(sample * rows + row) * cols + col];
      }

      ComplexT& operator()(IdxT sample, IdxT row, IdxT col)
      {
        return response[(sample * rows + row) * cols + col];
      }
    };
  } // namespace Optimization
} // namespace GridKit
