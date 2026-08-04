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
     * @brief Per-mode transport delays, each the smallest sampled
     *        delay over the frequencies where the mode still carries
     *        magnitude, |h_m| >= mag_floor * max |h_m|.
     *
     * Below its alive band a mode's phase delay keeps sliding toward
     * the lossless front, which carries no energy, so samples the fit
     * treats as zeros must not decide the delay.
     *
     * @param[in] tau       Delay traces, one row per mode, values in
     *                      the real part
     * @param[in] h         Modal propagation on the same grid
     * @param[in] mag_floor Per-mode peak fraction below which samples
     *                      are ignored; zero keeps every sample
     *
     * @pre tau and h share their grid and mode count, with at least
     *      one sample
     */
    template <typename scalar_type, typename index_type>
    std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>
    modalDelays(const SampledResponse<scalar_type, index_type>&    tau,
                const SampledResponse<scalar_type, index_type>&    h,
                typename GridKit::ScalarTraits<scalar_type>::RealT mag_floor);

    /**
     * @brief Shift every channel of the sampled response to minimum
     *        phase by multiplying with exp(j omega tau).
     */
    template <typename scalar_type, typename index_type>
    void applyDelayShift(SampledResponse<scalar_type, index_type>&          samples,
                         typename GridKit::ScalarTraits<scalar_type>::RealT tau);
  } // namespace Optimization
} // namespace GridKit
