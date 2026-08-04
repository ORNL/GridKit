/**
 * @file MinimumPhase.cpp
 *
 * @brief Delay extraction and minimum-phase shifting.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "MinimumPhase.hpp"

#include <algorithm>
#include <cmath>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Per-mode minimum over the alive samples of the delay
     *        trace. The per-mode peak sample is always alive, so every
     *        mode yields a delay.
     */
    template <typename scalar_type, typename index_type>
    std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>
    modalDelays(const SampledResponse<scalar_type, index_type>&    tau,
                const SampledResponse<scalar_type, index_type>&    h,
                typename GridKit::ScalarTraits<scalar_type>::RealT mag_floor)
    {
      using RealT = typename GridKit::ScalarTraits<scalar_type>::RealT;

      const auto modes        = tau.rows;
      const auto sample_count = tau.sampleCount();

      std::vector<RealT> delays(static_cast<size_t>(modes));
      for (index_type mode = 0; mode < modes; ++mode)
      {
        RealT peak = RealT{0};
        for (index_type m = 0; m < sample_count; ++m)
        {
          peak = std::max(peak, std::abs(h(m, mode, 0)));
        }
        const RealT floor_value = mag_floor * peak;

        bool  seen     = false;
        RealT smallest = RealT{0};
        for (index_type m = 0; m < sample_count; ++m)
        {
          if (std::abs(h(m, mode, 0)) < floor_value)
          {
            continue;
          }
          const RealT value = tau(m, mode, 0).real();
          smallest          = seen ? std::min(smallest, value) : value;
          seen              = true;
        }
        delays[static_cast<size_t>(mode)] = smallest;
      }
      return delays;
    }

    /**
     * @brief Shift every channel to minimum phase with one delay,
     *        removing the transport delay of the response before
     *        fitting.
     */
    template <typename scalar_type, typename index_type>
    void applyDelayShift(SampledResponse<scalar_type, index_type>&          samples,
                         typename GridKit::ScalarTraits<scalar_type>::RealT tau)
    {
      using RealT    = typename GridKit::ScalarTraits<scalar_type>::RealT;
      using ComplexT = std::complex<RealT>;

      const auto channels     = static_cast<size_t>(samples.rows) * static_cast<size_t>(samples.cols);
      const auto sample_count = samples.omega.size();

      for (size_t m = 0; m < sample_count; ++m)
      {
        const RealT    angle = samples.omega[m] * tau;
        const ComplexT shift{std::cos(angle), std::sin(angle)};
        for (size_t ch = 0; ch < channels; ++ch)
        {
          samples.response[m * channels + ch] *= shift;
        }
      }
    }

    template std::vector<double>
    modalDelays<double, long int>(const SampledResponse<double, long int>&,
                                  const SampledResponse<double, long int>&,
                                  double);
    template std::vector<double>
    modalDelays<double, size_t>(const SampledResponse<double, size_t>&,
                                const SampledResponse<double, size_t>&,
                                double);

    template void
    applyDelayShift<double, long int>(SampledResponse<double, long int>&,
                                      double);
    template void
    applyDelayShift<double, size_t>(SampledResponse<double, size_t>&, double);
  } // namespace Optimization
} // namespace GridKit
