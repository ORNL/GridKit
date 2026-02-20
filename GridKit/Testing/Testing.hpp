#pragma once

/**
 * @file TestingCore.hpp
 *
 * @author Slaven Peles <peless@ornl.gov>, ORNL
 *
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <utility>

namespace GridKit
{
  namespace Testing
  {

    enum TestOutcome
    {
      PASS = 0,
      FAIL,
      SKIP,
      EXPECTED_FAIL,
      UNEXPECTED_PASS
    };

    class TestStatus
    {
    public:
      TestStatus();

      TestStatus(bool success);

      TestStatus(const char* funcname);

      ~TestStatus();

      TestStatus& operator=(const bool isPass)
      {
        if (isPass)
          outcome_ = TestOutcome::PASS;
        else
          outcome_ = TestOutcome::FAIL;
        return *this;
      }

      TestStatus& operator*=(const bool isPass)
      {
        if (!isPass)
          outcome_ = TestOutcome::FAIL;
        return *this;
      }

      void skipTest()
      {
        outcome_ = TestOutcome::SKIP;
      }

      void expectFailure()
      {
        expectFailure_ = true;
      }

      TestOutcome report()
      {
        return report(funcname_);
      }

      TestOutcome report(const char* funcname);

    private:
      TestOutcome outcome_;
      const char* funcname_;
      bool        expectFailure_ = false;
    };

    struct TestingResults
    {
      int success            = 0;
      int failure            = 0;
      int skip               = 0;
      int expected_failure   = 0;
      int unexpected_success = 0;

      TestingResults();

      ~TestingResults();

      TestingResults(const TestingResults& r);

      void init();

      TestingResults& operator+=(const TestingResults& rhs);

      TestingResults& operator+=(const TestOutcome outcome);

      int summary();
    };

    inline TestingResults operator+(const TestingResults& lhs, const TestingResults& rhs)
    {
      return TestingResults(lhs) += rhs;
    }

    inline TestingResults operator+(const TestingResults& lhs, const TestOutcome outcome)
    {
      return TestingResults(lhs) += outcome;
    }

    inline TestingResults operator+(const TestOutcome outcome, const TestingResults& rhs)
    {
      return TestingResults(rhs) += outcome;
    }

    template <typename T>
    bool isEqual(const T value,
                 const T ref,
                 const T tol = std::numeric_limits<T>::epsilon())
    {
      T error = std::abs(value - ref) / (1.0 + std::abs(ref));
      return (error < tol);
    }

    /**
     * @brief Equality comparison between maps with a tolerance for the scalar value
     *
     * @tparam IdxT - Integer data type
     * @tparam RealT - Real type
     *
     * @param[in] a - first map to compare
     * @param[in] b - second map to compare
     * @param[in] tol - tolerance
     * @return bool - true if the maps are equal; false otherwise
     */
    template <typename IdxT = size_t, typename RealT = double>
    inline bool isEqual(std::map<IdxT, RealT> a,
                        std::map<IdxT, RealT> b,
                        const RealT           tol = std::numeric_limits<RealT>::epsilon())
    {
      int fail = 0;

      if (a.size() != b.size())
      {
        fail++;
        std::cerr << "Containers do not have the same size! "
                  << "a.size() = " << a.size() << ", and "
                  << "b.size() = " << b.size() << "\n";
      }

      for (const auto& pair_a : a)
      {
        auto it = b.find(pair_a.first);
        if (it != b.end())
        {
          if (!isEqual(pair_a.second, it->second, tol))
          {
            fail++;
            std::cerr << "Mismatching map values! "
                      << "a.first = " << pair_a.first << ", "
                      << "a.second = " << std::setprecision(16) << pair_a.second << ", and "
                      << "b.second = " << std::setprecision(16) << it->second << "\n";
          }
        }
        else
        {
          fail++;
          std::cerr << "Entry not found in the second container! "
                    << "a.first = " << pair_a.first << "\n";
        }
      }

      return fail == 0;
    }

    namespace Detail
    {
      template <typename TErr>
      struct ThrowsHelper
      {
        template <typename F>
        bool operator()(F f)
        {
          try
          {
            f();
          }
          catch (const TErr&)
          {
            return true;
          }
          return false;
        }
      };

      template <>
      struct ThrowsHelper<void>
      {
        template <typename F>
        bool operator()(F f)
        {
          try
          {
            f();
          }
          catch (...)
          {
            return true;
          }
          return false;
        }
      };
    } // namespace Detail

    template <typename TErr = void, typename TFunc, typename... TArgs>
    inline bool
    throws(TFunc func, TArgs&&... args)
    {
      return Detail::ThrowsHelper<TErr>{}(
          [&]()
          { func(std::forward<TArgs>(args)...); });
    }

  } // namespace Testing
} // namespace GridKit
