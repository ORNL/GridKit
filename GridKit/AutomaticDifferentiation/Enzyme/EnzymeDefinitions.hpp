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
    namespace Scalar
    {
      /**
       * @brief Enzyme autodiff template for scalar functions
       *
       * @tparam scalar_type - scalar type
       */
      template <typename scalar_type>
      extern scalar_type __enzyme_autodiff(scalar_type (*)(scalar_type), ...);
    } // namespace Scalar

    namespace Sparse
    {
      /**
       * @brief Enzyme fwddiff template for GridKit models
       *
       * @details This is a core templated intrinsic that the Enzyme pass will
       * use to perform automatic differenciation. We define the template here
       * so it can later be used in different places.
       *
       * @tparam T - return type
       * @tparam model_type - model type
       */
      template <typename T, typename... model_type>
      extern T __enzyme_fwddiff(void*, model_type...) noexcept;
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
