/**
 * @file EnzymeDefinitions.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

/**
 * @brief Enzyme constants for activity analysis
 *
 */
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme fwddiff template
       *
       * @tparam T - return type
       * @tparam ModelT - model type
       */
      template <typename T, typename... ModelT>
      extern T __enzyme_fwddiff(void*, ModelT...) noexcept;
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
