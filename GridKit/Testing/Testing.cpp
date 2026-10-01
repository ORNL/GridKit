
#include <iostream>

#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Colors.hpp>

namespace GridKit
{
  namespace Testing
  {

    TestStatus::TestStatus()
      : outcome_(TestOutcome::PASS)
    {
    }

    TestStatus::TestStatus(bool success)
      : outcome_(success ? TestOutcome::PASS : TestOutcome::FAIL)
    {
    }

    TestStatus::TestStatus(const char* funcname)
      : outcome_(TestOutcome::PASS),
        funcname_(funcname)
    {
    }

    TestStatus::~TestStatus()
    {
    }

    TestOutcome TestStatus::report(const char* funcname)
    {
      using namespace Utilities::Colors;

      if (expectFailure_)
      {
        if ((outcome_ == FAIL) || (outcome_ == EXPECTED_FAIL))
          outcome_ = EXPECTED_FAIL;
        else if ((outcome_ == PASS) || (outcome_ == UNEXPECTED_PASS))
          outcome_ = UNEXPECTED_PASS;
        else
          outcome_ = SKIP;
      }

      switch (outcome_)
      {
      case PASS:
        std::cout << "--- " << GREEN << "PASS" << CLEAR
                  << ": Test " << funcname << "\n";
        break;
      case FAIL:
        std::cout << "--- " << RED << "FAIL" << CLEAR
                  << ": Test " << funcname << "\n";
        break;
      case SKIP:
        std::cout << "--- " << YELLOW << "SKIP" << CLEAR
                  << ": Test " << funcname << CLEAR << "\n";
        break;
      case EXPECTED_FAIL:
        std::cout << "--- " << ORANGE << "FAIL" << CLEAR
                  << " (EXPECTED)"
                  << ": Test " << funcname << "\n";
        break;
      case UNEXPECTED_PASS:
        std::cout << "--- " << BLUE << "PASS" << CLEAR
                  << " (UNEXPECTED)"
                  << ": Test " << funcname << "\n";
        break;
      default:
        std::cout << "--- " << RED << "FAIL" << CLEAR
                  << "Unrecognized test result " << outcome_
                  << " for test " << funcname << "\n";
      }
      return outcome_;
    }

    TestingResults::TestingResults()
    {
    }

    TestingResults::~TestingResults()
    {
    }

    TestingResults::TestingResults(const TestingResults& r)
    {
      this->success            = r.success;
      this->failure            = r.failure;
      this->skip               = r.skip;
      this->expected_failure   = r.expected_failure;
      this->unexpected_success = r.unexpected_success;
    }

    void TestingResults::init()
    {
      this->success            = 0;
      this->failure            = 0;
      this->skip               = 0;
      this->expected_failure   = 0;
      this->unexpected_success = 0;
    }

    TestingResults& TestingResults::operator+=(const TestingResults& rhs)
    {
      this->success            += rhs.success;
      this->failure            += rhs.failure;
      this->skip               += rhs.skip;
      this->expected_failure   += rhs.expected_failure;
      this->unexpected_success += rhs.unexpected_success;

      return *this;
    }

    TestingResults& TestingResults::operator+=(const TestOutcome outcome)
    {
      switch (outcome)
      {
      case PASS:
        this->success++;
        break;
      case FAIL:
        this->failure++;
        break;
      case SKIP:
        this->skip++;
        break;
      case EXPECTED_FAIL:
        this->expected_failure++;
        break;
      case UNEXPECTED_PASS:
        this->unexpected_success++;
        break;
      default:
        std::cout << "Warning: Unrecognized test outcome code "
                  << outcome << ". Assuming failure ...\n";
        this->failure++;
      }
      return *this;
    }

    int TestingResults::summary()
    {
      std::cout << "\n\nTest Summary\n";
      std::cout << "----------------------------\n";
      std::cout << "\tSuccessful tests:     " << success << "\n";
      std::cout << "\tFailed test:          " << failure << "\n";
      std::cout << "\tSkipped tests:        " << skip << "\n";
      std::cout << "\tExpected failures:    " << expected_failure << "\n";
      std::cout << "\tUnexpected successes: " << unexpected_success << "\n";
      std::cout << "\n";

      return failure;
    }

  } // namespace Testing
} // namespace GridKit
