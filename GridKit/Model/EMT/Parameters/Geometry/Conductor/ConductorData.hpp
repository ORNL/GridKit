/**
 * @file ConductorData.hpp
 *
 * @brief Static conductor geometry, material, weight, and phase data for the
 * Conductor model.
 *
 */

#pragma once

#include <string>
#include <vector>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      struct ConductorData
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        IdxT                     K{0};
        std::vector<ScalarT>     radius;
        std::vector<ScalarT>     inner_radius;
        std::vector<ScalarT>     sigma;
        std::vector<ScalarT>     mu;
        std::vector<ScalarT>     weight;
        std::vector<std::string> phase;
        /// Circuit index on multi-circuit structures; meaningful for phases a, b, and c
        std::vector<IdxT>        circuit;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
