/**
 * @file MinimumPhase.hpp
 *
 * @brief Delay extraction and minimum-phase shifting of sampled responses.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Optimization/Rational/SampledResponse.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Smallest delay carried by a sampled delay trace,
     *        tau = min over every channel and every sample.
     */
    template <typename scalar_type, typename index_type>
    typename GridKit::ScalarTraits<scalar_type>::RealT
    minimumDelay(const SampledResponse<scalar_type, index_type>& tau);

    /**
     * @brief Shift every channel of the sampled response to minimum
     *        phase by multiplying with exp(j omega tau).
     */
    template <typename scalar_type, typename index_type>
    void applyDelayShift(SampledResponse<scalar_type, index_type>&          samples,
                         typename GridKit::ScalarTraits<scalar_type>::RealT tau);
  } // namespace Optimization
} // namespace GridKit
