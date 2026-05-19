#include "PDSim.hpp"

#include <cxxabi.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <optional>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/IdaDiagnostics.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

using Log = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;
using namespace AnalysisManager::Sundials;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

namespace
{
  /// Demangle a typeid name. Returns the mangled string if demangling fails.
  std::string demangleName(const char* mangled)
  {
    int         status = 0;
    char*       buf    = abi::__cxa_demangle(mangled, nullptr, nullptr, &status);
    std::string out    = (status == 0 && buf != nullptr) ? std::string(buf) : std::string(mangled);
    std::free(buf);
    return out;
  }

  /**
   * @brief Best-effort post-mortem: capture IDA's current internal state at the
   * point of failure, evaluate the residual once into the model's existing
   * buffers, identify the top-K largest |residual| entries, and write a JSON
   * dump to @p out_path.
   *
   * Zero allocation during normal operation: this function only runs on the
   * failure path and reuses the model's already-allocated y_/yp_/f_ vectors.
   * Mutates model state (acceptable since PDSim re-throws and exits right after).
   * Wrapped in a try/catch so any nested failure does not displace the
   * original SUNDIALS exception.
   */
  void dumpFailingResiduals(Ida<scalar_type, index_type>&         ida,
                            SystemModel<scalar_type, index_type>& sys,
                            const std::filesystem::path&          out_path,
                            const std::string&                    what,
                            std::size_t                           top_k = 10)
  {
    try
    {
      // 1. Borrow IDA's current internal yy/yp pointers (no copy on IDA's side).
      const real_type* yy = ida.yyData();
      const real_type* yp = ida.ypData();
      if (yy == nullptr || yp == nullptr)
        return;

      auto& my  = sys.y();
      auto& myp = sys.yp();
      if (my.empty() || myp.empty() || my.size() != myp.size())
        return;

      // 2. Copy into the model's existing buffers (one memcpy each — no allocation).
      std::memcpy(my.data(), yy, my.size() * sizeof(real_type));
      std::memcpy(myp.data(), yp, myp.size() * sizeof(real_type));

      // 3. Evaluate the residual at the failure state into model's existing f_.
      sys.evaluateResidual();
      const auto& f = sys.getResidual();
      if (f.empty())
        return;

      // 4. Build owner lookup: global residual index -> owning bus or component.
      struct OwnerInfo
      {
        std::string kind; // "Bus" or "Component"
        long        owner_id = -1;
        std::string owner_class;
        long        local_index = -1;
      };

      std::unordered_map<std::size_t, OwnerInfo> owners;
      owners.reserve(f.size());

      for (auto* bus : sys.getBuses())
      {
        if (bus == nullptr)
          continue;
        const auto class_name = demangleName(typeid(*bus).name());
        for (index_type j = 0; j < bus->size(); ++j)
        {
          owners.emplace(static_cast<std::size_t>(bus->getResidualIndex(j)),
                         OwnerInfo{"Bus", static_cast<long>(bus->busID()), class_name, static_cast<long>(j)});
        }
      }
      for (auto* comp : sys.getComponents())
      {
        if (comp == nullptr)
          continue;
        const auto class_name = demangleName(typeid(*comp).name());
        for (index_type j = 0; j < comp->size(); ++j)
        {
          owners.emplace(static_cast<std::size_t>(comp->getResidualIndex(j)),
                         OwnerInfo{"Component", static_cast<long>(comp->getGridKitComponentID()), class_name, static_cast<long>(j)});
        }
      }

      // 5. Find top-K by |f|. partial_sort is O(N log K).
      std::vector<std::size_t> idx(f.size());
      std::iota(idx.begin(), idx.end(), 0);
      const std::size_t take = std::min(top_k, f.size());
      std::partial_sort(idx.begin(),
                        idx.begin() + static_cast<std::ptrdiff_t>(take),
                        idx.end(),
                        [&f](std::size_t a, std::size_t b)
                        {
                          return std::abs(f[a]) > std::abs(f[b]);
                        });

      // 6. Compute residual norms for context.
      double rnorm_inf  = 0.0;
      double rnorm_2_sq = 0.0;
      for (std::size_t i = 0; i < f.size(); ++i)
      {
        const double v = std::abs(f[i]);
        if (v > rnorm_inf)
          rnorm_inf = v;
        rnorm_2_sq += v * v;
      }

      // 7. Serialize JSON.
      nlohmann::json j;
      j["what"]              = what;
      j["time"]              = ida.currentInternalTime();
      j["state_size"]        = f.size();
      j["residual_inf_norm"] = rnorm_inf;
      j["residual_2_norm"]   = std::sqrt(rnorm_2_sq);
      j["top_residuals"]     = nlohmann::json::array();
      for (std::size_t i = 0; i < take; ++i)
      {
        const std::size_t gi = idx[i];
        nlohmann::json    entry;
        entry["global_index"] = gi;
        entry["value"]        = f[gi];
        entry["abs_value"]    = std::abs(f[gi]);
        entry["y"]            = my[gi];
        entry["yp"]           = myp[gi];
        const auto it         = owners.find(gi);
        if (it != owners.end())
        {
          entry["owner_kind"]        = it->second.kind;
          entry["owner_id"]          = it->second.owner_id;
          entry["owner_class"]       = it->second.owner_class;
          entry["owner_local_index"] = it->second.local_index;
        }
        else
        {
          entry["owner_kind"] = "unknown";
        }
        j["top_residuals"].push_back(entry);
      }

      // 8. Write file (best-effort).
      std::ofstream stream(out_path);
      if (stream)
      {
        stream << j.dump(2) << '\n';
      }
    }
    catch (...)
    {
      // Defensive: dump must not displace the original SUNDIALS exception.
      // Silently abandon the post-mortem if anything throws (e.g. residual
      // evaluation itself fails at the failure point).
    }
  }
} // namespace

int main(int argc, const char* argv[])
{
  // Study file
  if (argc < 2)
  {
    Log::error() << "No input file provided" << std::endl;
    std::cout << "\n"
                 "Usage:\n"
                 "       pdsim <json-input-file>\n"
                 "\n"
                 "Please provide a json input file for the study to run.\n"
                 "\n";
    exit(1);
  }

  auto study = parseStudyData(argv[1]);

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  real_type dt = study.dt;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  IdaStatsRecorder             ida_stats(!study.ida_stats_file.empty());
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  ida.initializeSimulation(0.0, false);
  real_type                  curr_time    = 0.0;
  int                        solve_status = 0;
  std::exception_ptr         pending_exception;
  std::optional<std::string> pending_what;
  for (const auto& event : study.events)
  {
    // Run to event time
    int nout = static_cast<int>(std::round((event.time - curr_time) / dt));
    if (nout > 0)
    {
      ida_stats.beginSegment(ida);
      try
      {
        solve_status = ida.runSimulation(event.time, nout);
      }
      catch (const std::exception& ex)
      {
        solve_status      = -1;
        pending_exception = std::current_exception();
        pending_what      = std::string(ex.what());
      }
      // Capture stats even if the segment failed — IdaGetX counters remain valid.
      ida_stats.endSegment(ida, curr_time, event.time, nout);
    }
    if (solve_status != 0)
    {
      break;
    }

    // Set up run for event (to start at event time)
    if (event.type == EventType::FAULT_ON)
    {
      sys.getBusFault(event.element_id)->setStatus(true);
    }
    else if (event.type == EventType::FAULT_OFF)
    {
      sys.getBusFault(event.element_id)->setStatus(false);
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time, true);
    curr_time = event.time;
  }

  // Run to final time
  int nout = static_cast<int>(std::round((study.tmax - curr_time) / dt));
  if (solve_status == 0 && nout > 0)
  {
    ida_stats.beginSegment(ida);
    try
    {
      solve_status = ida.runSimulation(study.tmax, nout);
    }
    catch (const std::exception& ex)
    {
      solve_status      = -1;
      pending_exception = std::current_exception();
      pending_what      = std::string(ex.what());
    }
    ida_stats.endSegment(ida, curr_time, study.tmax, nout);
  }

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  if (!study.ida_stats_file.empty())
  {
    // Post-mortem residual dump (best-effort, only on failure).
    if (pending_exception)
    {
      const auto residual_path = study.ida_stats_file.parent_path() / "ida_residual_dump.json";
      dumpFailingResiduals(ida, sys, residual_path, pending_what.value_or(""));
    }
    writeIdaStatsJson(ida_stats.report(), {study.ida_stats_file, std::nullopt});
  }

  // Preserve original behaviour: surface the SUNDIALS exception (and its message)
  // by re-throwing now that diagnostics have been flushed to disk.
  if (pending_exception)
  {
    if (pending_what.has_value())
    {
      Log::error() << *pending_what << std::endl;
    }
    std::rethrow_exception(pending_exception);
  }

  // Generate aggregate errors comparing variable output to reference solution
  std::string func{"monitor file vs reference file"};
  TestStatus  status{func.c_str()};
  if (!study.output_file.empty() && !study.reference_file.empty())
  {
    auto errorSet = compareCSV(study.output_file, study.reference_file);

    // Print the errors
    errorSet.display();

    // Check against specified tolerance
    status *= errorSet.total.max_error < study.error_tol;

    status.report();
  }

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return solve_status == 0 ? status.get() : solve_status;
}
