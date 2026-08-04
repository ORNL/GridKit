/**
 * @file UniversalLineModel.cpp
 *
 * @brief Build universal line model coefficients from a line description.
 *
 * Sweeps the line-parameter model over frequency and fits Yc(s) to
 * yc.model.json and the propagation function to propagation.model.json
 * as the modal sum H(s) = sum_m Hmin_m(s) exp(-s tau_m), keys K,
 * modes [{Hmin, delay {tau}}]; see README.md for the fitting targets.
 * All rational fits are stable with no term linear in s.
 *
 * Exit codes: 0 when every error target was met and the Yc fit is
 * passive, 1 on a hard failure, 2 when an error target was missed,
 * 3 when the targets were met but the Yc fit is nonpassive (the
 * passivity report also travels inside yc.model.json).
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Parameters/OverheadDataJSONParser.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Solver/Optimization/Rational/MinimumPhase.hpp>
#include <GridKit/Solver/Optimization/Rational/Passivity.hpp>
#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>
#include <GridKit/Solver/Optimization/Rational/SampledResponse.hpp>
#include <GridKit/Solver/Optimization/VectorFitting/VectorFitting.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "FrequencySweep.hpp"
#include "RationalModelJSON.hpp"

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
  using ComplexT  = ResponseT::ComplexT;

  using GridKit::Optimization::Application::modelToJson;

  constexpr double BIORTHOGONALITY_WARNING = 1.0e-5;
  constexpr double BIORTHOGONALITY_ERROR   = 1.0e-2;

  /// Application settings; every field has a CLI override.
  struct Settings
  {
    struct Target
    {
      index_type min_poles{2};
      index_type max_poles{30};
      double     target_rel_rms{1.0e-3};

      /// Fraction of the peak sample magnitude below which values are
      /// fit as exact zeros; zero disables the cleaning.
      double min_mag{0.0};
    };

    double fmin{10.0};
    double fmax{1.0e8};
    size_t points{401};

    GridKit::EMT::Application::IdaSettings ida;

    fs::path   output{"output"};
    Target     yc;
    Target     h;
    index_type restarts{3};
    bool       refine{false};
  };

  int usage()
  {
    std::cout << "\n"
              << "Usage:\n"
              << "       UniversalLineModel <line-json-file> [options]\n"
              << "\n"
              << "Please provide an overhead line description JSON file.\n"
              << "Pass --help after the line file to list the options.\n"
              << "\n";
    return 1;
  }

  /// CLI table; the line description is positional as argv[1] until
  /// CliArgs grows positional-argument support (Option.hpp todo).
  GridKit::Utilities::CliArgs makeArgs()
  {
    using GridKit::Utilities::ArgType;

    return GridKit::Utilities::CliArgs{
        {.name     = {"--fmin"},
         .help     = "Sweep start frequency in Hz",
         .type     = ArgType::Real,
         .defaults = 10.0},

        {.name     = {"--fmax"},
         .help     = "Sweep stop frequency in Hz",
         .type     = ArgType::Real,
         .defaults = 1.0e8},

        {.name     = {"--points"},
         .help     = "Sweep sample count on the log grid",
         .type     = ArgType::Integer,
         .defaults = 401},

        {.name     = {"--rtol"},
         .help     = "Sweep relative tolerance",
         .type     = ArgType::Real,
         .defaults = 1.0e-7},

        {.name     = {"--max-steps"},
         .help     = "Maximum internal integrator steps for the sweep",
         .type     = ArgType::Integer,
         .defaults = 200000},

        {.name     = {"--output", "-o"},
         .help     = "Output directory for the model files",
         .type     = ArgType::String,
         .defaults = "output"},

        {.name     = {"--yc-min-poles"},
         .help     = "Lowest pole count tried for the Yc fit",
         .type     = ArgType::Integer,
         .defaults = 2},

        {.name     = {"--yc-max-poles"},
         .help     = "Highest pole count tried for the Yc fit",
         .type     = ArgType::Integer,
         .defaults = 30},

        {.name     = {"--yc-target"},
         .help     = "Relative RMS error target for the Yc fit",
         .type     = ArgType::Real,
         .defaults = 1.0e-3},

        {.name     = {"--h-min-poles"},
         .help     = "Lowest pole count tried for the propagation fit",
         .type     = ArgType::Integer,
         .defaults = 2},

        {.name     = {"--h-max-poles"},
         .help     = "Highest pole count tried for the propagation fit",
         .type     = ArgType::Integer,
         .defaults = 30},

        {.name     = {"--h-target"},
         .help     = "Relative RMS error target for the propagation fit",
         .type     = ArgType::Real,
         .defaults = 1.0e-3},

        {.name     = {"--min-mag"},
         .help     = "Fraction of the peak magnitude below which propagation "
                     "values are fit as exact zeros",
         .type     = ArgType::Real,
         .defaults = 1.0e-4},

        {.name     = {"--restarts"},
         .help     = "Perturbation restarts per fit when the error stays "
                     "above the fitter's restart threshold",
         .type     = ArgType::Integer,
         .defaults = 3},

        {.name = {"--refine"},
         .help = "Polish each fit with the constrained refinement stage",
         .flag = true}};
  }

  /// One fit target from its option prefix, validated before the raw
  /// integers are widened to index_type.
  Settings::Target makeTarget(const GridKit::Utilities::CliArgs& args,
                              const std::string&                 prefix)
  {
    const int    min_poles = args.get<int>(prefix + "-min-poles");
    const int    max_poles = args.get<int>(prefix + "-max-poles");
    const double target    = args.get<double>(prefix + "-target");

    if (min_poles < 1 || max_poles < min_poles)
    {
      throw std::runtime_error("--" + prefix + "-min-poles and --" + prefix
                               + "-max-poles must satisfy 1 <= min <= max");
    }
    if (!(target > 0.0))
    {
      throw std::runtime_error("--" + prefix + "-target must be positive");
    }

    return {.min_poles      = static_cast<index_type>(min_poles),
            .max_poles      = static_cast<index_type>(max_poles),
            .target_rel_rms = target};
  }

  /// Settings from the parsed options; the JSON validation layer is
  /// bypassed here, so the sweep and fitter preconditions are enforced
  /// on the raw values.
  Settings makeSettings(const GridKit::Utilities::CliArgs& args)
  {
    Settings settings;
    settings.fmin = args.get<double>("fmin");
    settings.fmax = args.get<double>("fmax");
    if (!(settings.fmin > 0.0) || !(settings.fmax > settings.fmin))
    {
      throw std::runtime_error("Sweep range requires 0 < --fmin < --fmax");
    }

    const int points = args.get<int>("points");
    if (points < 2)
    {
      throw std::runtime_error("--points must be at least 2");
    }
    settings.points = static_cast<size_t>(points);

    settings.ida.tolerance = args.get<double>("rtol");
    if (!(settings.ida.tolerance > 0.0) || settings.ida.tolerance >= 1.0e-1)
    {
      throw std::runtime_error("--rtol must be in (0, 0.1)");
    }

    const int max_steps = args.get<int>("max-steps");
    if (max_steps < 1)
    {
      throw std::runtime_error("--max-steps must be positive");
    }
    settings.ida.max_steps = static_cast<size_t>(max_steps);

    settings.output = args.get("output");
    settings.yc     = makeTarget(args, "yc");
    settings.h      = makeTarget(args, "h");

    // The magnitude floor applies to the propagation fit only; Yc has
    // no decayed entries and keeps the fitter default.
    settings.h.min_mag = args.get<double>("min-mag");
    if (!(settings.h.min_mag >= 0.0
          && settings.h.min_mag < 1.0))
    {
      throw std::runtime_error("--min-mag must lie in [0, 1)");
    }

    const int restarts = args.get<int>("restarts");
    if (restarts < 0)
    {
      throw std::runtime_error("--restarts must be nonnegative");
    }
    settings.restarts = static_cast<index_type>(restarts);

    settings.refine = args.get<bool>("refine");
    return settings;
  }

  /// One monitor CSV loaded whole: header tokens and value columns.
  struct MonitorTable
  {
    std::vector<std::string>         header;
    std::vector<std::vector<double>> columns;

    size_t rowCount() const
    {
      return columns.empty() ? 0 : columns.front().size();
    }

    /// Column index of an exact header token; -1 when absent.
    long find(const std::string& name) const
    {
      for (size_t k = 0; k < header.size(); ++k)
      {
        if (header[k] == name)
        {
          return static_cast<long>(k);
        }
      }
      return -1;
    }
  };

  MonitorTable readMonitorCsv(const fs::path& file_path)
  {
    std::ifstream stream(file_path);
    if (!stream)
    {
      throw std::runtime_error("Failed to open monitor CSV: " + file_path.string());
    }

    MonitorTable table;
    std::string  line;
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

      std::vector<std::string> tokens;
      std::stringstream        splitter(line);
      std::string              token;
      while (std::getline(splitter, token, ','))
      {
        tokens.push_back(token);
      }

      if (table.header.empty())
      {
        table.header = tokens;
        table.columns.assign(tokens.size(), {});
        continue;
      }
      if (tokens.size() != table.header.size())
      {
        throw std::runtime_error("Ragged monitor CSV row in " + file_path.string());
      }
      for (size_t k = 0; k < tokens.size(); ++k)
      {
        table.columns[k].push_back(std::stod(tokens[k]));
      }
    }
    if (table.header.empty() || table.rowCount() < 2)
    {
      throw std::runtime_error("Monitor CSV holds no samples: " + file_path.string());
    }
    return table;
  }

  /// Conductor count from the largest contiguous Overhead_Yc index.
  index_type conductorCount(const MonitorTable& table)
  {
    index_type count = 0;
    while (table.find("Overhead_Yc_real_" + std::to_string(count) + "_" + std::to_string(count)) >= 0)
    {
      ++count;
    }
    if (count == 0)
    {
      throw std::runtime_error(
          "Monitor CSV holds no Overhead_Yc columns; "
          "the sweep must monitor Yc, H, Tau, Tv, and Ti");
    }
    return count;
  }

  std::vector<double> column(const MonitorTable& table,
                             const std::string&  name)
  {
    const long index = table.find(name);
    if (index < 0)
    {
      throw std::runtime_error("Monitor CSV misses column: " + name);
    }
    return table.columns[static_cast<size_t>(index)];
  }

  /// Assemble the complex matrix response prefix_real/imag_row_col; the
  /// monitor writes both indices at every conductor count.
  ResponseT gatherMatrix(const MonitorTable& table,
                         const std::string&  prefix,
                         index_type          k)
  {
    ResponseT samples;
    samples.rows  = k;
    samples.cols  = k;
    samples.omega = column(table, "omega");
    samples.response.assign(samples.omega.size() * static_cast<size_t>(k) * static_cast<size_t>(k),
                            {});

    for (index_type row = 0; row < k; ++row)
    {
      for (index_type col = 0; col < k; ++col)
      {
        const std::string suffix =
            "_" + std::to_string(row) + "_" + std::to_string(col);
        const auto real_part = column(table, prefix + "_real" + suffix);
        const auto imag_part = column(table, prefix + "_imag" + suffix);

        for (size_t m = 0; m < samples.omega.size(); ++m)
        {
          samples(m, row, col) = {real_part[m], imag_part[m]};
        }
      }
    }
    return samples;
  }

  /// Assemble the complex per-mode response prefix_real/imag_mode as a
  /// one-column sample set.
  ResponseT gatherModes(const MonitorTable& table,
                        const std::string&  prefix,
                        index_type          modes)
  {
    ResponseT samples;
    samples.rows  = modes;
    samples.cols  = 1;
    samples.omega = column(table, "omega");
    samples.response.assign(samples.omega.size() * static_cast<size_t>(modes),
                            {});

    for (index_type mode = 0; mode < modes; ++mode)
    {
      const std::string suffix    = "_" + std::to_string(mode);
      const auto        real_part = column(table, prefix + "_real" + suffix);
      const auto        imag_part = column(table, prefix + "_imag" + suffix);

      for (size_t m = 0; m < samples.omega.size(); ++m)
      {
        samples(m, mode, 0) = {real_part[m], imag_part[m]};
      }
    }
    return samples;
  }

  /// Real delay trace Overhead_Tau_m as a one-column response per mode.
  ResponseT gatherDelays(const MonitorTable& table, index_type modes)
  {
    ResponseT samples;
    samples.rows  = modes;
    samples.cols  = 1;
    samples.omega = column(table, "omega");
    samples.response.assign(samples.omega.size() * static_cast<size_t>(modes),
                            {});

    for (index_type mode = 0; mode < modes; ++mode)
    {
      const auto values =
          column(table, "Overhead_Tau_" + std::to_string(mode));
      for (size_t m = 0; m < samples.omega.size(); ++m)
      {
        samples(m, mode, 0) = {values[m], 0.0};
      }
    }
    return samples;
  }

  /**
   * @brief Stop on corrupt transform data: the modal decomposition
   *        guarantees Ti^H Tv = I, so a residual beyond the threshold
   *        means the monitored data is corrupt.
   */
  void validateTransforms(const ResponseT& tv, const ResponseT& ti)
  {
    const auto k            = tv.rows;
    const auto sample_count = tv.omega.size();

    double worst = 0.0;
    for (size_t m = 0; m < sample_count; ++m)
    {
      for (index_type row = 0; row < k; ++row)
      {
        for (index_type col = 0; col < k; ++col)
        {
          ComplexT entry{0.0, 0.0};
          for (index_type i = 0; i < k; ++i)
          {
            entry += std::conj(ti(m, i, row)) * tv(m, i, col);
          }
          if (row == col)
          {
            entry -= 1.0;
          }
          const double magnitude = std::abs(entry);
          if (!std::isfinite(magnitude))
          {
            throw std::runtime_error(
                "Non-finite transform sample in the monitored sweep data");
          }
          worst = std::max(worst, magnitude);
        }
      }
    }

    if (worst > BIORTHOGONALITY_ERROR)
    {
      std::ostringstream message;
      message << "Transform biorthogonality residual reaches " << worst
              << ", beyond " << BIORTHOGONALITY_ERROR
              << "; the modal decomposition guarantees this identity, so "
                 "the monitored data is corrupt.";
      throw std::runtime_error(message.str());
    }
    if (worst > BIORTHOGONALITY_WARNING)
    {
      std::cout << "Warning: transform biorthogonality residual reaches "
                << worst << "\n";
    }
  }

  /**
   * @brief One mode's phase-domain contribution, the rank-one dyad
   *        conj(ti_mode) modal_mode tv_mode^T per sample.
   *
   * The monitored transforms satisfy the adjoint pairing Ti^H Tv = I
   * while the physical current transformation follows the transpose
   * pairing Ti = Tv^-T, so the current-form assembly consumes
   * conj(Ti). The eigenvector normalization cancels inside each dyad,
   * and the dyads over all modes sum to the full propagation function.
   */
  ResponseT modeDyad(const ResponseT& modal,
                     const ResponseT& tv,
                     const ResponseT& ti,
                     index_type       mode)
  {
    const auto k            = tv.rows;
    const auto sample_count = tv.omega.size();

    ResponseT dyad;
    dyad.rows  = k;
    dyad.cols  = k;
    dyad.omega = tv.omega;
    dyad.response.assign(tv.response.size(), {});

    for (size_t m = 0; m < sample_count; ++m)
    {
      for (index_type row = 0; row < k; ++row)
      {
        for (index_type col = 0; col < k; ++col)
        {
          dyad(m, row, col) =
              std::conj(ti(m, row, mode)) * modal(m, mode, 0) * tv(m, col, mode);
        }
      }
    }
    return dyad;
  }

  void writeJson(const fs::path& file_path, const nlohmann::json& j)
  {
    std::ofstream stream(file_path);
    if (!stream)
    {
      throw std::runtime_error("Failed to write " + file_path.string());
    }
    stream << j.dump(2) << "\n";
  }

  /**
   * @brief Lowest-order fit of one target with a constant term and no
   *        term linear in s, per the Propagation submodel validation.
   *
   * A relocation iteration that stops at its pass limit still yields a
   * usable model, so only the error target decides success here.
   *
   * @return 0 when the target was met, 2 when the order search ended
   *         above it (best model returned), negative on a hard failure
   */
  int fitTarget(const ResponseT&        samples,
                const Settings::Target& target,
                index_type              restarts,
                bool                    refine,
                const std::string&      label,
                ModelT&                 model)
  {
    FitterT fitter(samples);

    typename FitterT::Parameters options;
    options.terms = GridKit::Optimization::RationalTerms::CONSTANT;
    options.weighting =
        GridKit::Optimization::Weighting::INVERSE_MAGNITUDE;
    options.min_mag                     = target.min_mag;
    options.order_search.enabled        = true;
    options.order_search.min_poles      = target.min_poles;
    options.order_search.max_poles      = target.max_poles;
    options.order_search.target_rel_rms = target.target_rel_rms;
    options.restarts.max_restarts       = restarts;
    options.refine                      = refine;

    const int status = fitter.fit(model, options);
    if (status < 0)
    {
      Log::error() << label << " fit failed with code " << status
                   << std::endl;
      return status;
    }
    std::cout << label << ": " << fitter.getStats().report() << "\n";
    return status == 2 ? 2 : 0;
  }

  /**
   * @brief Each mode's rank-one dyad unwound by its own delay and fit
   *        individually; the propagation function is the sum of the
   *        fitted modes behind their delays.
   */
  int fitPropagation(const ResponseT& tv,
                     const ResponseT& ti,
                     const ResponseT& h,
                     const ResponseT& tau,
                     const Settings&  settings,
                     nlohmann::json&  propagation)
  {
    // The fitter's magnitude floor also gates the delay extraction, so
    // delays are identified from exactly the data the fit sees.
    const auto delays =
        GridKit::Optimization::modalDelays<scalar_type, index_type>(
            tau, h, settings.h.min_mag);

    constexpr double to_microseconds = 1.0e6;
    auto             modes           = nlohmann::json::array();
    int              worst           = 0;
    for (index_type mode = 0; mode < tv.rows; ++mode)
    {
      auto target = modeDyad(h, tv, ti, mode);
      GridKit::Optimization::applyDelayShift<scalar_type, index_type>(
          target, delays[mode]);

      std::cout << "Propagation mode " << mode << ": delay "
                << to_microseconds * delays[mode] << " us\n";

      ModelT    model;
      const int status = fitTarget(target,
                                   settings.h,
                                   settings.restarts,
                                   settings.refine,
                                   "Hmin" + std::to_string(mode),
                                   model);
      if (status < 0)
      {
        return status;
      }
      worst = std::max(worst, status);

      // The delay coefficient follows the Delay operator's parameter set.
      modes.push_back({{"Hmin", modelToJson(model)},
                       {"delay", {{"tau", delays[mode]}}}});
    }

    propagation["K"]     = tv.rows;
    propagation["modes"] = std::move(modes);
    return worst;
  }

  int runUniversalLineModel(const fs::path& line_file,
                            const Settings& settings)
  {
    using namespace GridKit::EMT::Application;
    using namespace GridKit::EMT::Parameters;
    using Variable =
        OverheadData<scalar_type, index_type>::MonitorableVariables;

    auto data = parseOverheadData<scalar_type, index_type>(line_file);

    fs::create_directories(settings.output);
    const auto response_csv = settings.output / "response.csv";

    data.monitored_variables = {Variable::Yc, Variable::H, Variable::Tau, Variable::Tv, Variable::Ti};
    data.monitor_sink        = {{GridKit::Model::VariableMonitorFormat::CSV,
                                 response_csv.string()}};

    const FrequencyGrid frequency{settings.fmin, settings.fmax, settings.points, "log"};

    const int sweep = runFrequencySweep<scalar_type, index_type>(
        data, frequency, settings.ida);
    if (sweep != 0)
    {
      Log::error() << "Frequency sweep failed with code " << sweep
                   << std::endl;
      return sweep;
    }

    const auto table = readMonitorCsv(response_csv);
    const auto k     = conductorCount(table);

    const auto yc  = gatherMatrix(table, "Overhead_Yc", k);
    const auto h   = gatherModes(table, "Overhead_H", k);
    const auto tv  = gatherMatrix(table, "Overhead_Tv", k);
    const auto ti  = gatherMatrix(table, "Overhead_Ti", k);
    const auto tau = gatherDelays(table, k);

    validateTransforms(tv, ti);

    ModelT    yc_model;
    const int yc_status =
        fitTarget(yc, settings.yc, settings.restarts, settings.refine, "Yc", yc_model);
    if (yc_status < 0)
    {
      return yc_status;
    }

    nlohmann::json propagation;
    const int      h_status =
        fitPropagation(tv, ti, h, tau, settings, propagation);
    if (h_status < 0)
    {
      return h_status;
    }

    GridKit::Optimization::PassivityReport<scalar_type, index_type> report;
    const int                                                       passivity =
        GridKit::Optimization::assessPassivity(yc_model,
                                               report,
                                               yc.omega.front(),
                                               yc.omega.back());
    if (passivity != 0)
    {
      Log::error() << "Passivity assessment failed with code " << passivity
                   << std::endl;
      return passivity;
    }
    if (report.passive)
    {
      std::cout << "Yc fit is passive over the screened band\n";
    }
    else
    {
      constexpr double two_pi =
          2.0 * GridKit::EMT::Constants::pi<double>();
      std::cout << "Warning: Yc fit violates passivity in "
                << report.violations.size() << " band(s):\n";
      for (const auto& band : report.violations)
      {
        std::cout << "  " << band.omega_start / two_pi << " Hz to "
                  << band.omega_end / two_pi << " Hz\n";
      }
    }

    // The verdict travels with the artifact, so a consumer can reject a
    // nonpassive model without reparsing this run's output.
    nlohmann::json yc_json = modelToJson(yc_model);
    nlohmann::json passivity_json;
    passivity_json["passive"] = report.passive;
    passivity_json["stable"]  = report.stable;
    auto violations           = nlohmann::json::array();
    for (const auto& band : report.violations)
    {
      nlohmann::json entry;
      entry["omega_start"] = band.omega_start;
      if (std::isfinite(band.omega_end))
      {
        entry["omega_end"] = band.omega_end;
      }
      violations.push_back(entry);
    }
    passivity_json["violations"] = violations;
    yc_json["passivity"]         = passivity_json;
    writeJson(settings.output / "yc.model.json", yc_json);

    writeJson(settings.output / "propagation.model.json", propagation);

    // Zero only when every error target was met and the Yc fit is
    // passive: 2 flags a missed error target, 3 a nonpassive Yc fit
    // whose targets were met.
    const int fit_status = std::max(yc_status, h_status);
    if (fit_status == 0 && !report.passive)
    {
      return 3;
    }
    return fit_status;
  }
} // namespace

int main(int argc, const char* argv[])
{
  if (argc < 2)
  {
    return usage();
  }

  auto args = makeArgs();
  if (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")
  {
    args.printHelp();
    return 0;
  }

  const auto start = Clock::now();
  try
  {
    args.parseArgs(argc - 1, argv + 1);

    const int retval =
        runUniversalLineModel(argv[1], makeSettings(args));
    const auto stop = Clock::now();
    const auto dur  = std::chrono::duration<double>(stop - start);
    std::cout << "\n\nComplete in " << dur << "\n";
    // Internal negative failure codes would wrap modulo 256 into
    // large shell statuses; 1 is the documented hard-failure exit.
    return retval < 0 ? 1 : retval;
  }
  catch (const std::exception& e)
  {
    Log::error() << e.what() << std::endl;
    return 1;
  }
}
