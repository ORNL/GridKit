/**
 * @file CarsonData.hpp
 *
 * @brief Static earth material data for the Carson model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Constants.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      struct CarsonData
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        ScalarT earth_sigma{0.0};
        ScalarT earth_eps{Constants::epsilon0<ScalarT>()};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
