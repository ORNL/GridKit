/**
 * @file VectorFittingJSONParser.hpp
 *
 * @brief JSON parsing for the VectorFitting application specification.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <filesystem>

#include "VectorFittingData.hpp"

namespace GridKit
{
  namespace Optimization
  {
    namespace Application
    {
      /**
       * @brief Parse a VectorFitting solver specification file.
       *
       * Follows the FrequencyResponse parser conventions, including
       * rejection of unsupported vocabulary.
       */
      VectorFittingData
      parseVectorFittingData(const std::filesystem::path& solver_file);
    } // namespace Application
  } // namespace Optimization
} // namespace GridKit
