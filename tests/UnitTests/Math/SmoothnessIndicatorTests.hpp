#pragma once

#include <type_traits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT>
    class SmoothnessIndicatorTests
    {
    private:
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

      static constexpr RealT kSmoothTolerance   = 1.0e-2;
      static constexpr RealT kNearOne           = 1.0 - kSmoothTolerance;
      static constexpr RealT kRoundoffTolerance = 100.0 * std::numeric_limits<RealT>::epsilon();

      static ScalarT scalar(const RealT value)
      {
        return static_cast<ScalarT>(value);
      }

      static bool within(const ScalarT value,
                         const ScalarT expected,
                         const ScalarT tolerance)
      {
        return std::abs(value - expected) < tolerance;
      }

    public:
      SmoothnessIndicatorTests()  = default;
      ~SmoothnessIndicatorTests() = default;

      TestOutcome clamp()
      {
        TestStatus success = true;

        const ScalarT lower = scalar(-0.25);
        const ScalarT upper = scalar(0.75);
        const ScalarT tol   = scalar(kSmoothTolerance);

        success *= (Math::clamp(scalar(-1.0), lower, upper) < lower + tol);
        success *= (Math::clamp(scalar(0.4), lower, upper) > lower);
        success *= (Math::clamp(scalar(0.4), lower, upper) < upper);
        success *= (Math::clamp(scalar(1.5), lower, upper) > upper - tol);

        return success.report(__func__);
      }

      TestOutcome deadband1()
      {
        TestStatus success = true;

        const ScalarT lower = scalar(-0.05);
        const ScalarT upper = scalar(0.10);
        const ScalarT tol   = scalar(kSmoothTolerance);

        // Type 1 (no-offset) deadband returns ~0 inside the band ...
        const ScalarT inside_input  = scalar(0.02);
        success                    *= (std::abs(Math::deadband1(inside_input, lower, upper)) < tol * tol);

        // ... and passes the input through unchanged (no offset) outside the band.
        const ScalarT far_above_input = scalar(4.0);
        const ScalarT far_below_input = scalar(-4.0);

        success *= std::isfinite(Math::deadband1(far_above_input, lower, upper));
        success *= within(Math::deadband1(far_above_input, lower, upper), far_above_input, tol);
        success *= std::isfinite(Math::deadband1(far_below_input, lower, upper));
        success *= within(Math::deadband1(far_below_input, lower, upper), far_below_input, tol);

        return success.report(__func__);
      }

      TestOutcome deadband2()
      {
        TestStatus success = true;

        const ScalarT lower = scalar(-0.05);
        const ScalarT upper = scalar(0.10);
        const ScalarT tol   = scalar(kSmoothTolerance);

        const ScalarT below_input    = scalar(-1.0);
        const ScalarT inside_input   = scalar(0.02);
        const ScalarT above_input    = scalar(1.0);
        const ScalarT expected_below = below_input - lower;
        const ScalarT expected_above = above_input - upper;

        success *= within(Math::deadband2(below_input, lower, upper), expected_below, tol);
        success *= (std::abs(Math::deadband2(inside_input, lower, upper)) < tol * tol);
        success *= within(Math::deadband2(above_input, lower, upper), expected_above, tol);

        const ScalarT lower_breakpoint = Math::deadband2(lower, lower, upper);
        const ScalarT upper_breakpoint = Math::deadband2(upper, lower, upper);

        success *= (lower_breakpoint < scalar(0.0));
        success *= (upper_breakpoint > scalar(0.0));
        success *= (std::abs(lower_breakpoint) < tol);
        success *= (std::abs(upper_breakpoint) < tol);

        const ScalarT x  = scalar(-0.4);
        success         *= (std::abs(Math::deadband2(x, lower, upper)
                             - (x - Math::clamp(x, lower, upper)))
                    < scalar(kRoundoffTolerance));

        const ScalarT far_above_input    = scalar(4.0);
        const ScalarT far_below_input    = scalar(-4.0);
        const ScalarT expected_far_above = far_above_input - upper;
        const ScalarT expected_far_below = far_below_input - lower;

        success *= std::isfinite(Math::deadband2(far_above_input, lower, upper));
        success *= within(Math::deadband2(far_above_input, lower, upper), expected_far_above, tol);
        success *= std::isfinite(Math::deadband2(far_below_input, lower, upper));
        success *= within(Math::deadband2(far_below_input, lower, upper), expected_far_below, tol);

        const ScalarT point        = scalar(0.25);
        const ScalarT above_point  = scalar(0.75);
        const ScalarT below_point  = scalar(-0.25);
        success                   *= (std::abs(Math::deadband2(above_point, point, point) - (above_point - point))
                    < scalar(kRoundoffTolerance));
        success                   *= (std::abs(Math::deadband2(below_point, point, point) - (below_point - point))
                    < scalar(kRoundoffTolerance));

        return success.report(__func__);
      }

      TestOutcome limitIndicators()
      {
        TestStatus success = true;

        const ScalarT limit_min = scalar(0.0);
        const ScalarT limit_max = scalar(3.0);
        const ScalarT near_zero = scalar(kSmoothTolerance);
        const ScalarT near_one  = scalar(kNearOne);

        success *= (Math::above(scalar(1.0), limit_min) > near_one);
        success *= (Math::above(scalar(-0.2), limit_min) < near_zero);
        success *= (Math::below(scalar(1.0), limit_max) > near_one);
        success *= (Math::below(scalar(3.2), limit_max) < near_zero);

        success *= (Math::inside(scalar(1.5), limit_min, limit_max) > near_one);
        success *= (Math::inside(scalar(-0.2), limit_min, limit_max) < near_zero);
        success *= (Math::inside(scalar(3.2), limit_min, limit_max) < near_zero);

        success *= (Math::outside(scalar(1.5), limit_min, limit_max) < near_zero);
        success *= (Math::outside(scalar(-0.2), limit_min, limit_max) > near_one);
        success *= (Math::outside(scalar(3.2), limit_min, limit_max) > near_one);

        const ScalarT x  = scalar(1.5);
        success         *= (std::abs(Math::inside(x, limit_min, limit_max)
                             + Math::outside(x, limit_min, limit_max)
                             - scalar(1.0))
                    < scalar(kRoundoffTolerance));

        return success.report(__func__);
      }

      TestOutcome slew()
      {
        TestStatus success = true;

        const ScalarT rate = scalar(0.5);
        const ScalarT tol  = scalar(kSmoothTolerance);

        success *= within(Math::slew(scalar(2.0), rate), rate, tol);
        success *= within(Math::slew(scalar(-2.0), rate), -rate, tol);
        success *= within(Math::slew(scalar(0.2), rate), scalar(0.2), tol);
        success *= within(Math::slew(scalar(-0.2), rate), scalar(-0.2), tol);

        return success.report(__func__);
      }

      TestOutcome linseg()
      {
        TestStatus success = true;

        const ScalarT lower  = scalar(0.0);
        const ScalarT upper  = scalar(2.0);
        const ScalarT height = scalar(4.0);
        const ScalarT tol    = scalar(kSmoothTolerance);

        const ScalarT below_input      = scalar(-1.0);
        const ScalarT midpoint_input   = scalar(1.0);
        const ScalarT above_input      = scalar(3.0);
        const ScalarT expected_mid     = height * (midpoint_input - lower) / (upper - lower);
        const ScalarT negative_height  = -height;
        const ScalarT expected_mid_neg = negative_height * (midpoint_input - lower) / (upper - lower);

        success *= (Math::linseg(below_input, lower, upper, height) < tol);
        success *= within(Math::linseg(midpoint_input, lower, upper, height), expected_mid, tol);
        success *= within(Math::linseg(above_input, lower, upper, height), height, tol);
        success *= within(Math::linseg(midpoint_input, lower, upper, negative_height), expected_mid_neg, tol);

        return success.report(__func__);
      }

      TestOutcome ramp()
      {
        TestStatus success = true;

        const ScalarT tol         = scalar(kSmoothTolerance);
        const ScalarT roundoff    = scalar(kRoundoffTolerance);
        const ScalarT tau         = scalar(1.0 / Math::MU<RealT>);
        const ScalarT at_zero     = tau * std::log(scalar(2.0));
        const ScalarT far_above   = scalar(4.0);
        const ScalarT far_below   = scalar(-4.0);
        const ScalarT lower       = scalar(-0.25);
        const ScalarT upper       = scalar(0.75);
        const ScalarT x           = scalar(0.4);
        const ScalarT smooth_clip = lower + Math::ramp(x - lower)
                                    - Math::ramp(x - upper);

        success *= (Math::ramp(scalar(1.0)) > scalar(kNearOne));
        success *= (Math::ramp(scalar(-1.0)) < tol);
        success *= (std::abs(Math::ramp(scalar(0.0)) - at_zero) < roundoff);
        success *= (Math::ramp(-tol) > scalar(0.0));
        success *= (Math::ramp(tol) > Math::ramp(scalar(0.0)));
        success *= std::isfinite(Math::ramp(far_above));
        success *= (Math::ramp(far_above) > far_above - tol);
        success *= std::isfinite(Math::ramp(far_below));
        success *= (Math::ramp(far_below) < roundoff);

        success *= (smooth_clip > lower);
        success *= (smooth_clip < upper);
        success *= std::isfinite(Math::clamp(far_above, lower, upper));
        success *= (Math::clamp(far_above, lower, upper) < upper + roundoff);
        success *= std::isfinite(Math::clamp(far_below, lower, upper));
        success *= (Math::clamp(far_below, lower, upper) > lower - roundoff);

        return success.report(__func__);
      }

      TestOutcome minMax()
      {
        TestStatus success = true;

        const ScalarT high = scalar(2.0);
        const ScalarT low  = scalar(-1.0);
        const ScalarT tol  = scalar(kSmoothTolerance);

        success *= within(Math::max(high, low), high, tol);
        success *= within(Math::max(low, high), high, tol);

        success *= within(Math::min(high, low), low, tol);
        success *= within(Math::min(low, high), low, tol);

        using Variable = GridKit::DependencyTracking::Variable;

        // This mirrors model equations where a differentiable state is
        // limited by a plain real parameter.
        const RealT    parameter_bound = kSmoothTolerance;
        const Variable state_below_bound{-1.0, 0};
        const Variable state_above_bound{1.0, 1};
        const auto     lower_bounded_state = Math::max(state_below_bound, parameter_bound);
        const auto     upper_bounded_state = Math::min(state_above_bound, parameter_bound);

        static_assert(std::is_same<typename std::decay<decltype(lower_bounded_state)>::type, Variable>::value,
                      "Mixed-type max should keep the differentiable scalar type.");
        static_assert(std::is_same<typename std::decay<decltype(upper_bounded_state)>::type, Variable>::value,
                      "Mixed-type min should keep the differentiable scalar type.");

        success *= (std::abs(lower_bounded_state.getValue() - parameter_bound) < kSmoothTolerance);
        success *= (std::abs(upper_bounded_state.getValue() - parameter_bound) < kSmoothTolerance);

        const ScalarT x  = scalar(0.4);
        const ScalarT y  = scalar(-0.7);
        success         *= (std::abs(Math::min(x, y) + Math::max(x, y) - (x + y))
                    < scalar(kRoundoffTolerance));

        const ScalarT point = scalar(0.25);
        const ScalarT bias  = std::log(scalar(2.0)) / scalar(Math::MU<RealT>);

        success *= (std::abs(Math::max(point, point) - (point + bias))
                    < scalar(kRoundoffTolerance));
        success *= (std::abs(Math::min(point, point) - (point - bias))
                    < scalar(kRoundoffTolerance));

        return success.report(__func__);
      }

      TestOutcome antiWindupIndicator()
      {
        TestStatus success = true;

        const ScalarT limit_min   = scalar(0.0);
        const ScalarT limit_max   = scalar(3.0);
        const ScalarT inside      = scalar(1.5);
        const ScalarT above_limit = scalar(3.2);
        const ScalarT below_limit = scalar(-0.2);
        const ScalarT f_positive  = scalar(0.03);
        const ScalarT f_negative  = -f_positive;
        const ScalarT near_zero   = scalar(kSmoothTolerance);
        const ScalarT near_one    = scalar(kNearOne);

        // Inside the limits, both increasing and decreasing dynamics pass.
        success *= (Math::indicator(inside, f_positive, limit_min, limit_max) > near_one);
        success *= (Math::indicator(inside, f_negative, limit_min, limit_max) > near_one);

        // Above the upper limit, increasing dynamics are blocked.
        success *= (Math::indicator(above_limit, f_positive, limit_min, limit_max) < near_zero);

        // Above the upper limit, decreasing dynamics restore the state.
        success *= (Math::indicator(above_limit, f_negative, limit_min, limit_max) > near_one);

        // Below the lower limit, decreasing dynamics are blocked.
        success *= (Math::indicator(below_limit, f_negative, limit_min, limit_max) < near_zero);

        // Below the lower limit, increasing dynamics restore the state.
        success *= (Math::indicator(below_limit, f_positive, limit_min, limit_max) > near_one);

        return success.report(__func__);
      }

      TestOutcome antiWindup()
      {
        TestStatus success = true;

        const ScalarT limit_min   = scalar(0.0);
        const ScalarT limit_max   = scalar(3.0);
        const ScalarT inside      = scalar(1.5);
        const ScalarT above_limit = scalar(3.2);
        const ScalarT below_limit = scalar(-0.2);

        const ScalarT f_positive        = scalar(0.03);
        const ScalarT f_negative        = -f_positive;
        const ScalarT blocked_tolerance = f_positive * scalar(kSmoothTolerance);

        success *= (std::abs(Math::antiwindup(inside, f_positive, limit_min, limit_max)
                             - Math::indicator(inside, f_positive, limit_min, limit_max) * f_positive)
                    < scalar(kRoundoffTolerance));

        success *= (std::abs(Math::antiwindup(above_limit, f_positive, limit_min, limit_max)) < blocked_tolerance);
        success *= within(Math::antiwindup(above_limit, f_negative, limit_min, limit_max),
                          f_negative,
                          blocked_tolerance);
        success *= (std::abs(Math::antiwindup(below_limit, f_negative, limit_min, limit_max)) < blocked_tolerance);
        success *= within(Math::antiwindup(below_limit, f_positive, limit_min, limit_max),
                          f_positive,
                          blocked_tolerance);

        return success.report(__func__);
      }

      TestOutcome dynamicAntiWindupBounds()
      {
        TestStatus success = true;

        using Variable = GridKit::DependencyTracking::Variable;

        const Variable state{0.0, 0};
        const Variable rate{0.03, 1};
        const Variable lower{-0.05, 2};
        const Variable upper{0.05, 3};

        const auto gate    = Math::indicator(state, rate, lower, upper);
        const auto limited = Math::antiwindup(state, rate, lower, upper);

        static_assert(std::is_same<typename std::decay<decltype(gate)>::type,
                                   Variable>::value,
                      "Dynamic-bound indicator should retain dependency tracking.");
        static_assert(std::is_same<typename std::decay<decltype(limited)>::type,
                                   Variable>::value,
                      "Dynamic-bound antiwindup should retain dependency tracking.");

        success                          *= (gate.getValue() > kNearOne);
        success                          *= within(limited.getValue(), rate.getValue(), kSmoothTolerance);
        const auto& gate_dependencies     = gate.getDependencies();
        const auto& limited_dependencies  = limited.getDependencies();
        for (size_t variable = 0; variable < 4; ++variable)
        {
          success *= gate_dependencies.contains(variable);
          success *= limited_dependencies.contains(variable);
        }
        for (const size_t bound : {size_t{2}, size_t{3}})
        {
          const auto gate_bound    = gate_dependencies.find(bound);
          const auto limited_bound = limited_dependencies.find(bound);
          if (gate_bound != gate_dependencies.end())
          {
            success *= std::abs(gate_bound->second) > 0.0;
          }
          if (limited_bound != limited_dependencies.end())
          {
            success *= std::abs(limited_bound->second) > 0.0;
          }
        }

        const RealT    real_lower{-0.05};
        const Variable mixed_upper{0.05, 4};
        const auto     mixed_gate = Math::indicator(state, rate, real_lower, mixed_upper);
        static_assert(std::is_same<typename std::decay<decltype(mixed_gate)>::type,
                                   Variable>::value,
                      "A real lower bound and dynamic upper bound should retain the scalar type.");
        success *= mixed_gate.getDependencies().contains(4);

        const Variable equal_lower{0.0, 5};
        const Variable equal_upper{0.0, 6};
        const auto     equal_gate = Math::indicator(state, rate, equal_lower, equal_upper);
        const auto     equal_limited =
            Math::antiwindup(state, rate, equal_lower, equal_upper);
        success *= within(equal_gate.getValue(), 0.75, kRoundoffTolerance);
        success *= within(equal_limited.getValue(), 0.75 * rate.getValue(), kRoundoffTolerance);
        success *= equal_gate.getDependencies().contains(5);
        success *= equal_gate.getDependencies().contains(6);
        success *= equal_limited.getDependencies().contains(5);
        success *= equal_limited.getDependencies().contains(6);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
