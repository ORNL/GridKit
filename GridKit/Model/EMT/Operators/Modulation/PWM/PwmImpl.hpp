/**
 * @file PwmImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of EMT pulse-width modulation.
 *
 * The model equations are in the `README.md` of this directory.
 */

#pragma once

#include <cmath>
#include <limits>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/Pwm.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct the modulation and select the series for its carrier.
     *
     * The edges are resolved by the pulse sum when the smoothing scale times
     * the carrier period is at least sqrt(2) pi. Below that the carrier
     * harmonics are summed, and the duty ratios are used alone once the
     * ripple is below round-off. `verify` reports invalid data.
     */
    template <typename RealT, typename IdxT>
    Pwm<RealT, IdxT>::Pwm(const ModelDataT& data)
      : N_(data.N),
        fc_(data.fc),
        alignment_(data.alignment),
        injection_(data.injection),
        k0_(data.k0),
        mu_(data.mu)
    {
      if (!(mu_ > ZERO<RealT> && fc_ > ZERO<RealT>) )
      {
        return;
      }

      constexpr RealT PI = std::numbers::pi_v<RealT>;

      c_ = mu_ / fc_;
      if (c_ >= std::numbers::sqrt2_v<RealT> * PI)
      {
        series_ = PwmSeries::pulses;
        return;
      }

      // With z = 2 pi^2 / c, the ripple and its duty derivative are within 8 z exp(-z), kept below a quarter of round-off
      constexpr RealT BOUND = 32.0;

      const RealT z = TWO<RealT> * PI * PI / c_;
      if (std::exp(-z) <= std::numeric_limits<RealT>::epsilon() / (BOUND * z))
      {
        series_ = PwmSeries::duty;
        return;
      }

      series_ = PwmSeries::harmonics;
      for (int h = 1; h <= HARMONICS; ++h)
      {
        const RealT zh                    = z * static_cast<RealT>(h);
        gain_[static_cast<size_t>(h - 1)] = zh / std::sinh(zh) / (PI * static_cast<RealT>(h));
      }
    }

    /**
     * @brief Check the channel count, the carrier frequency, the smoothing
     *        scale, the alignment and the zero-state partition.
     *
     * @return The failed checks.
     */
    template <typename RealT, typename IdxT>
    Model::ConfigurationChecks Pwm<RealT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Pwm");
      checks.check(N_ == CHANNELS, "N must be 3");
      checks.check(std::isfinite(fc_) && fc_ > ZERO<RealT>, "fc must be positive");
      checks.check(std::isfinite(mu_) && mu_ > ZERO<RealT>, "mu must be finite and positive");
      checks.check(alignment_ >= ZERO<RealT> && alignment_ <= ONE<RealT>, "alignment must be between 0 and 1");
      checks.check(k0_ >= ZERO<RealT> && k0_ <= ONE<RealT>, "k0 must be between 0 and 1");
      return checks;
    }

    /**
     * @brief Switching functions of the modulation signals at time `t`,
     *        evaluated by `series`.
     *
     * @tparam series Series that evaluates the switching functions.
     * @tparam T Scalar type of the kernel.
     * @param[in] m Modulation signal of each channel.
     * @param[in] t Time [s].
     * @param[out] s Switching function of each channel.
     *
     * @pre The modulation signals give duty ratios between 0 and 1.
     */
    template <typename RealT, typename IdxT>
    template <PwmSeries series, typename T>
    [[gnu::always_inline]] void Pwm<RealT, IdxT>::switching(const T* m, RealT t, T* s) const
    {
      constexpr RealT TWO_PI = TWO<RealT> * std::numbers::pi_v<RealT>;

      T m0 = ZERO<RealT>;
      if (injection_)
      {
        const T m_max = Math::max(Math::max(m[0], m[1]), m[2]);
        const T m_min = Math::min(Math::min(m[0], m[1]), m[2]);
        m0            = -(ONE<RealT> - TWO<RealT> * k0_ + k0_ * m_max + (ONE<RealT> - k0_) * m_min);
      }

#pragma clang loop unroll(full)
      for (size_t n = 0; n < CHANNELS; ++n)
      {
        s[n] = HALF<RealT> * (ONE<RealT> + m[n] + m0);
      }
      if constexpr (series == PwmSeries::duty)
      {
        return;
      }

      const RealT phase = t * fc_ - std::floor(t * fc_);
#pragma clang loop unroll(full)
      for (size_t n = 0; n < CHANNELS; ++n)
      {
        const T d = s[n];
        const T a = alignment_ * (ONE<RealT> - d);

        if constexpr (series == PwmSeries::pulses)
        {
          // tanh(x) - tanh(y) = tanh(x - y) (1 - tanh(x) tanh(y))
          const RealT scale  = HALF<RealT> * c_;
          T           pulses = ZERO<RealT>;
#pragma clang loop unroll(full)
          for (int k = -RADIUS; k <= RADIUS; ++k)
          {
            const RealT start  = phase - static_cast<RealT>(k);
            pulses            += ONE<RealT> - std::tanh(scale * (start - a)) * std::tanh(scale * (start - a - d));
          }
          s[n] = HALF<RealT> * std::tanh(scale * d) * pulses;
        }
        else
        {
          s[n] = d + (harmonicSum(TWO_PI * (phase - a)) - harmonicSum(TWO_PI * (phase - a - d)));
        }
      }
    }

    /**
     * @brief Switching functions of the modulation signals at time `t`,
     *        evaluated by the series that the carrier selects.
     *
     * @tparam T Scalar type of the values.
     * @param[in] m Modulation signal of each channel.
     * @param[in] t Time [s].
     * @param[out] s Switching function of each channel.
     *
     * @pre The modulation signals give duty ratios between 0 and 1.
     */
    template <typename RealT, typename IdxT>
    template <typename T>
    void Pwm<RealT, IdxT>::switching(const T* m, RealT t, T* s) const
    {
      switch (series_)
      {
      case PwmSeries::pulses:
        return switching<PwmSeries::pulses>(m, t, s);
      case PwmSeries::harmonics:
        return switching<PwmSeries::harmonics>(m, t, s);
      case PwmSeries::duty:
        return switching<PwmSeries::duty>(m, t, s);
      }
    }

    /**
     * @brief Sum over the carrier harmonics h of gain[h - 1] sin(h angle),
     *        by the Clenshaw recurrence.
     *
     * @tparam T Scalar type of the kernel.
     */
    template <typename RealT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] T Pwm<RealT, IdxT>::harmonicSum(T angle) const
    {
      const T c    = TWO<RealT> * std::cos(angle);
      T       next = ZERO<RealT>;
      T       last = ZERO<RealT>;
#pragma clang loop unroll(full)
      for (int h = HARMONICS; h > 0; --h)
      {
        const T current = gain_[static_cast<size_t>(h - 1)] + c * next - last;
        last            = next;
        next            = current;
      }
      return std::sin(angle) * next;
    }
  } // namespace EMT
} // namespace GridKit
