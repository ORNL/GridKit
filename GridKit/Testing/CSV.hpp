#pragma once

#include <fstream>
#include <string>

#include <GridKit/Testing/AggregateErrors.hpp>

namespace GridKit
{
  namespace Testing
  {

    /**
     * @brief Throw error if file does not exist or if file fails to open (read)
     */
    std::ifstream checkOpenFile(const std::string& f);

    /**
     * @brief Generate ErrorSet for variables in two files
     *
     * @note This assumes a "time" variable is always in the first column.
     */
    std::unique_ptr<ErrorSet> compareCSV(
        const std::string& test_file,
        const std::string& reference_file,
        ErrorType          error_type    = ErrorType::RELATIVE,
        double             abs_threshold = DEFAULT_ABS_ERROR_THRESHOLD);

  } // namespace Testing
} // namespace GridKit
