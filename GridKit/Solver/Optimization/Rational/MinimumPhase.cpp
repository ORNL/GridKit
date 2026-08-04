/**
 * @file MinimumPhase.cpp
 *
 * @brief Delay extraction and minimum-phase shifting.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "MinimumPhase.hpp"

#include <cmath>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Smallest delay over the sweep per channel of the sampled
     *        delay trace, tau = min over omega of tau(omega). Delay
     *        samples carry their values in the real part.
     *
     * @pre The trace holds at least one sample
     */
    template <typename scalar_type, typename index_type>
    std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>
    minimumDelay(const SampledResponse<scalar_type, index_type>& tau)
    {
      using RealT = typename GridKit::ScalarTraits<scalar_type>::RealT;

      const auto channels     = static_cast<size_t>(tau.rows) * static_cast<size_t>(tau.cols);
      const auto sample_count = tau.omega.size();

      std::vector<RealT> delays(channels, RealT{0});
      for (size_t ch = 0; ch < channels; ++ch)
      {
        RealT smallest = tau.response[ch].real();
        for (size_t m = 1; m < sample_count; ++m)
        {
          smallest =
              std::min(smallest, tau.response[m * channels + ch].real());
        }
        delays[ch] = smallest;
      }
      return delays;
    }

    /**
     * @brief Shift each channel to minimum phase by multiplying with
     *        exp(j omega tau) for that channel's delay, removing the
     *        bulk transport delay before fitting.
     *
     * @pre tau holds rows * cols entries
     */
    template <typename scalar_type, typename index_type>
    void applyDelayShift(
        SampledResponse<scalar_type, index_type>& samples,
        const std::vector<typename GridKit::ScalarTraits<scalar_type>::RealT>&
            tau)
    {
      using RealT    = typename GridKit::ScalarTraits<scalar_type>::RealT;
      using ComplexT = std::complex<RealT>;

      const auto channels     = static_cast<size_t>(samples.rows) * static_cast<size_t>(samples.cols);
      const auto sample_count = samples.omega.size();

      for (size_t m = 0; m < sample_count; ++m)
      {
        for (size_t ch = 0; ch < channels; ++ch)
        {
          const RealT angle = samples.omega[m] * tau[ch];
          samples.response[m * channels + ch] *=
              ComplexT{std::cos(angle), std::sin(angle)};
        }
      }
    }

    /**
     * @brief Shift every channel to minimum phase with one shared
     *        delay, removing the bulk transport delay of the whole
     *        response before fitting.
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
    minimumDelay<double, long int>(const SampledResponse<double, long int>&);
    template std::vector<double>
    minimumDelay<double, size_t>(const SampledResponse<double, size_t>&);

    template void
    applyDelayShift<double, long int>(SampledResponse<double, long int>&,
                                      const std::vector<double>&);
    template void
    applyDelayShift<double, size_t>(SampledResponse<double, size_t>&,
                                    const std::vector<double>&);

    template void
    applyDelayShift<double, long int>(SampledResponse<double, long int>&,
                                      double);
    template void
    applyDelayShift<double, size_t>(SampledResponse<double, size_t>&, double);
  } // namespace Optimization
} // namespace GridKit
