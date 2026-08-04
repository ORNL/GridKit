/**
 * @file VectorFittingData.hpp
 *
 * @brief Input data for the VectorFitting application.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <filesystem>
#include <vector>

#include <GridKit/Solver/Optimization/VectorFitting/VectorFitting.hpp>

namespace GridKit
{
  namespace Optimization
  {
    namespace Application
    {
      /// Fitting request applied to the sampled response.
      struct FitSettings
      {
        size_t        poles{10};
        size_t        min_poles{2}; ///< Order search bounds.
        size_t        max_poles{0}; ///< Order search when > 0.
        double        target_rel_rms{1.0e-3};
        double        min_mag{0.0}; ///< Zero disables.
        RationalTerms terms{RationalTerms::CONSTANT};
        Weighting     weighting{Weighting::UNIFORM};
      };

      struct VectorFittingData
      {
        std::filesystem::path samples; ///< Sampled-response CSV.
        size_t                rows{1}; ///< Output dimension N.
        size_t                cols{1}; ///< Input dimension K.
        FitSettings           fit;
        std::filesystem::path output_model; ///< Rational model JSON.
      };
    } // namespace Application
  } // namespace Optimization
} // namespace GridKit
