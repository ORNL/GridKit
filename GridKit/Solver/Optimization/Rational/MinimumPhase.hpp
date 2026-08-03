/**
 * @file MinimumPhase.hpp
 *
 * @brief Delay extraction and minimum-phase shifting of sampled responses.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <vector>

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Optimization/Rational/SampledResponse.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Smallest delay over the sweep, tau = min over omega of
     *        tau(omega), per column of the sampled delay trace.
     */
    template <typename scalar_type, typename index_type>
    std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>
    minimumDelay(const SampledResponse<scalar_type, index_type>& tau);

    /**
     * @brief Shift each column of the sampled response to minimum phase by
     *        multiplying with exp(j omega tau) for that column's delay.
     */
    template <typename scalar_type, typename index_type>
    void applyDelayShift(
        SampledResponse<scalar_type, index_type>& samples,
        const std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>&
            tau);
  } // namespace Optimization
} // namespace GridKit
