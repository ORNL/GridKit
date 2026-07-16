/**
 * @file EnzymeTests.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class EnzymeTests
    {
    public:
      EnzymeTests()  = default;
      ~EnzymeTests() = default;

      TestOutcome scalar_square()
      {
        TestStatus success = true;


        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
