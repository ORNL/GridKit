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
       * @details This is a core templated intrinsic that the Enzyme pass will
       * use to perform automatic differenciation. We define the template here
       * so it can later be used in different places.
       *
       * @tparam T - return type
       * @tparam ModelT - model type
       */
      template <typename T, typename... ModelT>
      extern T __enzyme_fwddiff(void*, ModelT...) noexcept;
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
