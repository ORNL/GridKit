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
    ErrorSet compareCSV(const std::string& f_a, const std::string& f_b);

  } // namespace Testing
} // namespace GridKit
