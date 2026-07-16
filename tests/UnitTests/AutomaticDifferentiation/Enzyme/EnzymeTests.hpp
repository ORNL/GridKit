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
  
        ScalarT var        = 5.0;
        ScalarT sq         = square(var);
        ScalarT dsq_ref    = dsquare(var);
        ScalarT dsq_enzyme = GridKit::Enzyme::Scalar::__enzyme_autodiff<ScalarT>(square, var);

        success *= GridKit::Testing::isEqual(dsq_enzyme, dsq_ref);

        if (!success)
        {
          std::cout << "x = " << var << "\n"
                    << "x^2 = " << sq << "\n"
                    << "Reference d(x^2)/dx = " << dsq_ref << "\n"
                    << "Enzyme d(x^2)/dx = " << dsq_enzyme << "\n";
        }

        return success.report(__func__);
      }

    private:
      static ScalarT square(ScalarT x)
      {
        return x * x;
      }

      static ScalarT dsquare(ScalarT x)
      {
        return 2.0 * x;
      }
    };

  } // namespace Testing
} // namespace GridKit
