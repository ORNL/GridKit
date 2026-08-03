/**
 * @file VectorFitting.cpp
 *
 * @brief Fit a rational model to a sampled frequency response.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include <chrono>
#include <filesystem>
#include <iostream>

#include <GridKit/Solver/Optimization/Rational/StateSpaceRealization.hpp>
#include <GridKit/Solver/Optimization/VectorFitting/VectorFitting.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "VectorFittingJSONParser.hpp"

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
              << "       VectorFitting <solver-json-file>\n"
              << "\n"
              << "Please provide a VectorFitting solver JSON file.\n"
              << "\n";
    return 1;
  }

  int runVectorFitting(const fs::path& solver_file)
  {
    // Pipeline (implementation pending):
    //
    //   1. Parse the solver specification.
    //   2. Read the sampled response CSV into SampledResponse.
    //   3. Fit with GridKit::Optimization::VectorFitting.
    //   4. Write the rational model JSON; optionally the state-space
    //      realization JSON.
    //   5. Print the fit statistics report.

    (void)solver_file;
    return 0;
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
    const int  retval = runVectorFitting(argv[1]);
    const auto stop   = Clock::now();
    const auto dur    = std::chrono::duration<double>(stop - start);
    std::cout << "\n\nComplete in " << dur << "\n";
    return retval;
  }
  catch (const std::exception& e)
  {
    Log::error() << e.what() << std::endl;
    return 1;
  }
}
