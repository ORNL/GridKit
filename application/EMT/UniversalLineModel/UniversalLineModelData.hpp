/**
 * @file UniversalLineModelData.hpp
 *
 * @brief Input data for the EMT UniversalLineModel application.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <filesystem>

#include "FrequencyResponseData.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Application
    {
      /// Lowest-order fitting request for one target quantity.
      struct FitTargetSettings
      {
        size_t min_poles{2};
        size_t max_poles{30};
        double target_rel_rms{1.0e-3};
      };

      struct UniversalLineModelData
      {
        std::filesystem::path model;     ///< Line description JSON.
        FrequencyGrid         frequency; ///< Sweep grid.
        IdaSettings           ida;       ///< Sweep solver settings.

        FitTargetSettings yc; ///< Characteristic admittance fit.
        FitTargetSettings h;  ///< Minimum-phase propagation fit.

        std::filesystem::path output_directory; ///< Model JSON destination.
      };
    } // namespace Application
  } // namespace EMT
} // namespace GridKit
