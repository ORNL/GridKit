#pragma once

#include <GridKit/CommonMath.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT>
    class SmoothnessIndicatorTests
    {
    public:
      SmoothnessIndicatorTests()  = default;
      ~SmoothnessIndicatorTests() = default;

      TestOutcome antiWindupIndicator()
      {
        TestStatus success = true;

        const ScalarT limit_min = 0.0;
        const ScalarT limit_max = 3.0;

        // Inside the limits the indicator passes dynamics through, regardless
        // of the sign of f: value is close to 1.
        success *= (Math::indicator(limit_min, limit_max, static_cast<ScalarT>(1.5), static_cast<ScalarT>(0.01)) > static_cast<ScalarT>(0.99));

        // Above the upper limit with f pushing further out: blocked (≈ 0).
        success *= (Math::indicator(limit_min, limit_max, static_cast<ScalarT>(3.2), static_cast<ScalarT>(0.01)) < static_cast<ScalarT>(0.1));

        // Above the upper limit but f pulling back in: passed (≈ 1).
        success *= (Math::indicator(limit_min, limit_max, static_cast<ScalarT>(3.2), static_cast<ScalarT>(-0.01)) > static_cast<ScalarT>(0.9));

        // Below the lower limit with f pushing further out: blocked (≈ 0).
        success *= (Math::indicator(limit_min, limit_max, static_cast<ScalarT>(-0.2), static_cast<ScalarT>(-0.01)) < static_cast<ScalarT>(0.1));

        // Below the lower limit but f pulling back in: passed (≈ 1).
        success *= (Math::indicator(limit_min, limit_max, static_cast<ScalarT>(-0.2), static_cast<ScalarT>(0.01)) > static_cast<ScalarT>(0.9));

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
