/**
 * @file VectorFitting.cpp
 *
 * @brief Fit a rational model to a sampled frequency response.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>
#include <GridKit/Solver/Optimization/Rational/SampledResponse.hpp>
#include <GridKit/Solver/Optimization/VectorFitting/VectorFitting.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "RationalModelJSON.hpp"
#include "VectorFittingJSONParser.hpp"

using scalar_type = double;
using index_type  = size_t;
using Clock       = std::chrono::high_resolution_clock;

namespace
{
  namespace fs = std::filesystem;

  using Log       = GridKit::Utilities::Logger;
  using ResponseT = GridKit::Optimization::SampledResponse<scalar_type,
                                                           index_type>;
  using ModelT    = GridKit::Optimization::RationalModel<scalar_type,
                                                         index_type>;
  using FitterT   = GridKit::Optimization::VectorFitting<scalar_type,
                                                         index_type>;

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

  /**
   * @brief Read a sampled response CSV with header
   *        omega, re_0_0, im_0_0, re_0_1, ... (channels row-major).
   */
  ResponseT readSamplesCsv(const fs::path& file_path,
                           index_type      rows,
                           index_type      cols)
  {
    std::ifstream stream(file_path);
    if (!stream)
    {
      throw std::runtime_error("Failed to open samples CSV: " + file_path.string());
    }

    const auto channels =
        static_cast<size_t>(rows) * static_cast<size_t>(cols);

    ResponseT samples;
    samples.rows = rows;
    samples.cols = cols;

    std::string line;
    bool        header = true;
    while (std::getline(stream, line))
    {
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      if (line.empty())
      {
        continue;
      }
      if (header)
      {
        header = false;
        continue;
      }

      std::vector<double> values;
      std::stringstream   splitter(line);
      std::string         token;
      while (std::getline(splitter, token, ','))
      {
        values.push_back(std::stod(token));
      }
      if (values.size() != 1 + 2 * channels)
      {
        throw std::runtime_error("Samples CSV row width does not match "
                                 "rows * cols in "
                                 + file_path.string());
      }

      samples.omega.push_back(values[0]);
      for (size_t ch = 0; ch < channels; ++ch)
      {
        samples.response.push_back(
            {values[1 + 2 * ch], values[2 + 2 * ch]});
      }
    }
    if (samples.omega.size() < 2)
    {
      throw std::runtime_error("Samples CSV holds fewer than two rows: " + file_path.string());
    }
    return samples;
  }

  int runVectorFitting(const fs::path& solver_file)
  {
    using namespace GridKit::Optimization::Application;

    const auto spec = parseVectorFittingData(solver_file);
    const auto samples =
        readSamplesCsv(spec.samples, spec.rows, spec.cols);

    FitterT fitter(samples);

    typename FitterT::Parameters options;
    options.pole_count = spec.fit.poles;
    options.terms      = spec.fit.terms;
    options.weighting  = spec.fit.weighting;
    if (spec.fit.max_poles > 0)
    {
      options.order_search.enabled        = true;
      options.order_search.min_poles      = spec.fit.min_poles;
      options.order_search.max_poles      = spec.fit.max_poles;
      options.order_search.target_rel_rms = spec.fit.target_rel_rms;
    }

    ModelT    model;
    const int status = fitter.fit(model, options);
    if (status < 0)
    {
      Log::error() << "VectorFitting failed with code " << status
                   << std::endl;
      return status;
    }
    std::cout << fitter.getStats().report() << "\n";

    std::ofstream stream(spec.output_model);
    if (!stream)
    {
      throw std::runtime_error("Failed to write " + spec.output_model.string());
    }
    stream << modelToJson(model).dump(2) << "\n";

    if (!spec.output_state_space.empty())
    {
      std::cout << "Note: state-space export arrives with the "
                   "StateSpaceRealization implementation; skipped.\n";
    }

    // A relocation iteration that stops at its pass limit still yields a
    // usable model; only a missed order-search target is a failure here.
    return status == 2 ? 2 : 0;
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
