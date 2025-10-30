/**
 * @file Testing.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains utilies for testing.
 *
 */
#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

namespace GridKit
{
  namespace Testing
  {

    template <typename T>
    bool isEqual(const T value,
                 const T ref,
                 const T tol = std::numeric_limits<T>::epsilon())
    {
      T error = std::abs(value - ref) / (1.0 + std::abs(ref));
      return (error < tol);
    }

    /**
     * @brief Equatlity comparison between maps with a tolerance for the scalar value
     *
     * @tparam IdxT
     * @tparam RealT
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

  } // namespace Testing

} // namespace GridKit
