/**
 * @file UniversalLineModelJSONParser.hpp
 *
 * @brief JSON parsing for the UniversalLineModel application specification.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <filesystem>

#include "UniversalLineModelData.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Application
    {
      /**
       * @brief Parse a UniversalLineModel solver specification file.
       *
       * Follows the FrequencyResponse parser conventions, including
       * rejection of unsupported vocabulary.
       */
      UniversalLineModelData
      parseUniversalLineModelData(const std::filesystem::path& solver_file);
    } // namespace Application
  } // namespace EMT
} // namespace GridKit
