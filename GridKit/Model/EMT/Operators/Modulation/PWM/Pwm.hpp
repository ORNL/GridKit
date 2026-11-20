/**
 * @file Pwm.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of EMT pulse-width modulation.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/PwmData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Evaluation of the switching functions, chosen by the smoothing scale times the carrier period
    enum class PwmSeries
    {
      pulses,    ///< Sum of the nearest pulses, for resolved edges
      harmonics, ///< Carrier harmonics about the duty ratios
      duty,      ///< Duty ratios, where the ripple is below round-off
    };

    /**
     * @brief Pulse-width modulation of three channels on one carrier
     *
     * The model has no variables, so it is no `Component`. Its parent
     * evaluates the switching functions inside its own residual kernels, with
     * the series as a constant: Enzyme does not differentiate a choice among
     * the series.
     */
    template <typename RealT, typename IdxT>
    class Pwm
    {
    public:
      using ModelDataT = PwmData<RealT, IdxT>;

      explicit Pwm(const ModelDataT& data);

      Model::ConfigurationChecks verify() const;

      template <PwmSeries series, typename T>
      [[gnu::always_inline]] void switching(const T* m, RealT t, T* s) const;

      template <typename T>
      void switching(const T* m, RealT t, T* s) const;

      /// Series that the smoothing scale and the carrier select
      PwmSeries series() const
      {
        return series_;
      }

    private:
      template <typename T>
      [[gnu::always_inline]] T harmonicSum(T angle) const;

      /// Number of channels
      static constexpr size_t CHANNELS{3};

      /// Carrier periods summed on each side of the current one
      static constexpr int RADIUS = 10;

      /// Carrier harmonics summed
      static constexpr int HARMONICS = 10;

      IdxT  N_{0};
      RealT fc_{ZERO<RealT>};
      RealT alignment_{HALF<RealT>};
      bool  injection_{true};
      RealT k0_{HALF<RealT>};
      RealT mu_{ZERO<RealT>};

      RealT                        c_{ZERO<RealT>}; ///< Smoothing scale times the carrier period
      PwmSeries                    series_{PwmSeries::harmonics};
      std::array<RealT, HARMONICS> gain_{}; ///< Amplitude of each carrier harmonic
    };
  } // namespace EMT
} // namespace GridKit
