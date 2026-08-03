/**
 * @file FrequencyResponse.cpp
 *
 * @brief Evaluate EMT line-parameter models over frequency.
 *
 */

#include <chrono>
#include <filesystem>
#include <iostream>

#include <GridKit/Model/EMT/Parameters/OverheadDataJSONParser.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "FrequencyResponseJSONParser.hpp"
#include "FrequencySweep.hpp"

using scalar_type = double;
using index_type  = size_t;
using Clock       = std::chrono::high_resolution_clock;

namespace
{
  namespace fs = std::filesystem;

  using Log = GridKit::Utilities::Logger;

  int usage()
  {
    std::cout << "\n"
              << "Usage:\n"
              << "       FrequencyResponse <solver-json-file>\n"
              << "\n"
              << "Please provide a FrequencyResponse solver JSON file.\n"
              << "\n";
    return 1;
  }

  int runFrequencyResponse(const fs::path& solver_file)
  {
    using namespace GridKit::EMT::Application;
    using namespace GridKit::EMT::Parameters;

    const auto spec = parseFrequencyResponseData(solver_file);
    auto       data = parseOverheadData<scalar_type, index_type>(spec.model);

    data.monitored_variables = spec.variables;
    data.monitor_sink        = {{GridKit::Model::VariableMonitorFormat::CSV,
                                 spec.output_file.string()}};

    return runFrequencySweep<scalar_type, index_type>(data,
                                                      spec.frequency,
                                                      spec.ida);
  }
} // namespace

int main(int argc, const char* argv[])
{
  if (argc != 2)
  {
    return usage();
  }

  const auto start = Clock::now();
  try
  {
    const int  retval = runFrequencyResponse(argv[1]);
    const auto stop   = Clock::now();
    const auto dur    = std::chrono::duration<double>(stop - start);
    std::cout << "\n\nComplete in " << dur << "\n";
    // Internal negative failure codes would wrap modulo 256 into large
    // shell statuses; 1 is the documented hard-failure exit.
    return retval < 0 ? 1 : retval;
  }
  catch (const std::exception& e)
  {
    Log::error() << e.what() << std::endl;
    return 1;
  }
}
