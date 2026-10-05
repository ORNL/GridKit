#include "MultiStageCosimulation.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Interpolation.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

#include "GridKit/Solver/Dynamic/Ida.hpp"

namespace AnalysisManager
{

  template <class ScalarT, typename IdxT>
  MultiStageCosimulation<ScalarT, IdxT>::MultiStageCosimulation(
      size_t      max_iterations,
      size_t      num_stages,
      MemorySpace memspace)
    : num_stages_(num_stages),
      max_iterations_(max_iterations),
      memspace_(memspace)
  {
    if (num_stages_ == 0)
    {
      throw std::invalid_argument("MultiStageCosimulation requires at least one stage.");
    }

    if (max_iterations_ == 0)
    {
      throw std::invalid_argument("MultiStageCosimulation requires at least one iteration.");
    }
  }

  template <class ScalarT, typename IdxT>
  MultiStageCosimulation<ScalarT, IdxT>::~MultiStageCosimulation() = default;

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::addSubsystems(Subsystem* subsys)
  {
    partitions_.push_back(subsys);
    num_partitions_ = partitions_.size();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::allocate()
  {
    ida_solvers_.resize(num_partitions_);
    coupling_mat_y_.resize(num_partitions_);
    coupling_mat_yp_.resize(num_partitions_);
    interpolants_y_.resize(num_partitions_);
    interpolants_yp_.resize(num_partitions_);
    next_stages_y_.resize(num_stages_ + 1);
    curr_stages_y_.resize(num_stages_ + 1);
    next_stages_yp_.resize(num_stages_ + 1);
    curr_stages_yp_.resize(num_stages_ + 1);

    IdxT global_size = 0;

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto* partition = partitions_[p];

      partition->allocate();
      const auto coupling_size  = partition->getExternSize();
      const auto partition_size = partition->getInternalSize();

      // One coupling state per stage (index 0 = start of step).
      coupling_mat_y_[p].resize(num_stages_ + 1);
      coupling_mat_yp_[p].resize(num_stages_ + 1);

      for (size_t stage = 0; stage <= num_stages_; ++stage)
      {
        coupling_mat_y_[p][stage] = std::make_shared<State>(coupling_size);
        coupling_mat_y_[p][stage]->allocate(memspace_);
        coupling_mat_y_[p][stage]->setToZero(memspace_);

        coupling_mat_yp_[p][stage] = std::make_shared<State>(coupling_size);
        coupling_mat_yp_[p][stage]->allocate(memspace_);
        coupling_mat_yp_[p][stage]->setToZero(memspace_);
      }

      interpolants_y_[p].resize(num_stages_ + 1);
      interpolants_yp_[p].resize(num_stages_ + 1);

      for (size_t stage = 1; stage <= num_stages_; ++stage)
      {
        CouplingMatrix cols_y(coupling_mat_y_[p].begin(), coupling_mat_y_[p].begin() + static_cast<std::ptrdiff_t>(stage));
        CouplingMatrix cols_yp(coupling_mat_yp_[p].begin(), coupling_mat_yp_[p].begin() + static_cast<std::ptrdiff_t>(stage));

        // Placeholder interval [0, 1]; overwritten at each stage solve.
        interpolants_y_[p][stage]  = std::make_unique<Interpolant>(cols_y, static_cast<RealT>(0), static_cast<RealT>(1), GridKit::NodeType::Chebyshev, memspace_);
        interpolants_yp_[p][stage] = std::make_unique<Interpolant>(cols_yp, static_cast<RealT>(0), static_cast<RealT>(1), GridKit::NodeType::Chebyshev, memspace_);
      }

      global_size += static_cast<IdxT>(partition_size);

      ida_solvers_[p] = std::make_unique<IdaSolver>(partitions_[p]);
    }

    for (size_t stage = 0; stage <= num_stages_; ++stage)
    {
      curr_stages_y_[stage] = std::make_unique<State>(global_size);
      next_stages_y_[stage] = std::make_unique<State>(global_size);

      curr_stages_y_[stage]->allocate(memspace_);
      curr_stages_y_[stage]->setToZero(memspace_);
      next_stages_y_[stage]->allocate(memspace_);
      next_stages_y_[stage]->setToZero(memspace_);

      curr_stages_yp_[stage] = std::make_unique<State>(global_size);
      next_stages_yp_[stage] = std::make_unique<State>(global_size);

      curr_stages_yp_[stage]->allocate(memspace_);
      curr_stages_yp_[stage]->setToZero(memspace_);
      next_stages_yp_[stage]->allocate(memspace_);
      next_stages_yp_[stage]->setToZero(memspace_);
    }

    y_  = std::make_shared<State>(global_size);
    yp_ = std::make_shared<State>(global_size);

    y_->allocate(memspace_);
    yp_->allocate(memspace_);

    y_->setToZero(memspace_);
    yp_->setToZero(memspace_);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::initializeSimulation(std::vector<ScalarT>& time_points, State* y, State* yp)
  {
    for (size_t p = 0; p < num_partitions_; ++p)
    {
      partitions_[p]->initialize();
    }

    if (y != nullptr)
    {
      y_->copyFromExternal(*y, memspace_, memspace_);
      yp_->copyFromExternal(*yp, memspace_, memspace_);
    }
    else
    {
      for (size_t p = 0; p < num_partitions_; ++p)
      {
        gatherLocalVariables(*y_, *yp_, partitions_[p]);
      }
    }

    y_->setDataUpdated(memspace_);
    yp_->setDataUpdated(memspace_);

    time_points_ = time_points;

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      ida_solvers_[p]->setMaxSteps(100000);
      ida_solvers_[p]->setTolerance(1e-4);
      ida_solvers_[p]->configureSimulation();
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::runSimulation(std::optional<std::function<void(const RealT t, const State&)>> step_callback)
  {
    for (size_t t_idx = 0; t_idx + 1 < time_points_.size(); t_idx++)
    {
      const RealT t  = time_points_[t_idx];
      const RealT dt = time_points_[t_idx + 1] - time_points_[t_idx];

      curr_stages_y_[0] = y_;
      next_stages_y_[0] = y_;

      curr_stages_yp_[0] = yp_;
      next_stages_yp_[0] = yp_;

      for (size_t curr_stage = 1; curr_stage <= num_stages_; ++curr_stage)
      {
        const size_t next_stage = curr_stage + 1;

        partitionStageSolve(curr_stage, next_stage, t, dt);

        std::swap(curr_stages_y_, next_stages_y_);
        std::swap(curr_stages_yp_, next_stages_yp_);
      }

      std::cout << "Time point: " << time_points_[t_idx + 1] << std::endl;
      y_->copyFromExternal(*curr_stages_y_.back(), memspace_);
      yp_->copyFromExternal(*curr_stages_yp_.back(), memspace_);

      if (step_callback)
      {
        (*step_callback)(t + dt, *y_);
      }
    }
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeLocalVariables(const State& global_y, const State& global_yp, Subsystem* partition)
  {
    const auto* global_y_data  = global_y.getData(memspace_);
    const auto* global_yp_data = global_yp.getData(memspace_);

    auto* y_data  = partition->y().getData(memspace_);
    auto* yp_data = partition->yp().getData(memspace_);

    for (size_t i = 0; i < partition->getInternalSize(); ++i)
    {
      const auto g = partition->getNodeConnection(i);
      y_data[i]    = global_y_data[g];
      yp_data[i]   = global_yp_data[g];
    }

    partition->y().setDataUpdated(memspace_);
    partition->yp().setDataUpdated(memspace_);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::gatherLocalVariables(State& global_y, State& global_yp, Subsystem* partition)
  {
    auto* global_y_data  = global_y.getData(memspace_);
    auto* global_yp_data = global_yp.getData(memspace_);

    const auto* y_data  = partition->y().getData(memspace_);
    const auto* yp_data = partition->yp().getData(memspace_);

    for (size_t i = 0; i < partition->getInternalSize(); ++i)
    {
      const auto g      = partition->getNodeConnection(i);
      global_y_data[g]  = y_data[i];
      global_yp_data[g] = yp_data[i];
    }

    global_y.setDataUpdated(memspace_);
    global_yp.setDataUpdated(memspace_);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeCouplingVariables(size_t partition_idx, size_t num_stages)
  {
    for (size_t stage_idx = 0; stage_idx < num_stages; ++stage_idx)
    {
      distributeCouplingStage(partition_idx, stage_idx);
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeCouplingStage(
      size_t partition_idx,
      size_t stage_idx)
  {
    const auto& external_indices = partitions_[partition_idx]->getExternalDataIndices();
    const auto* source_data_y    = curr_stages_y_[stage_idx]->getData(memspace_);
    const auto* source_data_yp   = curr_stages_yp_[stage_idx]->getData(memspace_);

    auto& coupling_stage_y  = *(coupling_mat_y_[partition_idx][stage_idx]);
    auto& coupling_stage_yp = *(coupling_mat_yp_[partition_idx][stage_idx]);

    auto* coupling_y_data  = coupling_stage_y.getData(memspace_);
    auto* coupling_yp_data = coupling_stage_yp.getData(memspace_);

    for (size_t external_idx = 0; external_idx < external_indices.size(); ++external_idx)
    {
      const auto global_idx = external_indices[external_idx];

      coupling_y_data[external_idx]  = source_data_y[global_idx];
      coupling_yp_data[external_idx] = source_data_yp[global_idx];
    }
    coupling_stage_y.setDataUpdated(memspace_);
    coupling_stage_yp.setDataUpdated(memspace_);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::partitionStageSolve(size_t curr_stage, size_t next_stage, RealT t0, RealT dt)
  {
    assert(curr_stage >= 1 && curr_stage <= num_stages_);
    assert(next_stage == curr_stage + 1);

    const auto curr_stage_output_times = Interpolant::getPhysicalOutputTimes(std::max(curr_stage, static_cast<size_t>(2)), t0, dt);
    const auto next_stage_output_times = Interpolant::getPhysicalOutputTimes(next_stage, t0, dt);

    const size_t num_iterations = curr_stage == 1 ? 1 : max_iterations_;

// Per-stage setup (not per-iteration): rescale this stage's prebuilt
// interpolants to the current time interval. Node locations and weights
// are already computed; this only stores t0 and dt.
#pragma omp parallel for 
    for (size_t p = 0; p < num_partitions_; ++p)
    {
      interpolants_y_[p][curr_stage]->setTimeInterval(t0, dt);
      interpolants_yp_[p][curr_stage]->setTimeInterval(t0, dt);
    }

    for (size_t iteration = 0; iteration < num_iterations; iteration++)
    {
      const bool  is_last_iteration = (iteration + 1 == num_iterations);
      const auto& output_times      = is_last_iteration ? next_stage_output_times : curr_stage_output_times;

#pragma omp parallel for schedule(dynamic)
      for (size_t p = 0; p < num_partitions_; p++)
      {
        auto* partition = partitions_[p];

        distributeLocalVariables(*y_, *yp_, partition);
        distributeCouplingVariables(p, curr_stage);

        Interpolant* poly_y  = interpolants_y_[p][curr_stage].get();
        Interpolant* poly_yp = interpolants_yp_[p][curr_stage].get();

        partition->setForcingFunction(
            [poly_y, poly_yp, partition](ScalarT t) mutable
            {
              poly_y->evaluate(t, partition->getExternalDataY());
              poly_yp->evaluate(t, partition->getExternalDataYP());
            });
      }

      // Per-partition status so integrator failures are not silently dropped
      // inside the parallel region.
      std::vector<int> status(num_partitions_, 0);

#pragma omp parallel for schedule(dynamic)
      for (size_t p = 0; p < num_partitions_; p++)
      {
        status[p] = integrate(ida_solvers_[p].get(), partitions_[p], output_times, is_last_iteration);
      }

      for (size_t p = 0; p < num_partitions_; ++p)
      {
        if (status[p] != 0)
        {
          return status[p];
        }
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::integrate(
      IdaSolver*                solver,
      Subsystem*                model,
      const std::vector<RealT>& output_times,
      const bool                is_last_iteration)
  {
    if (output_times.empty())
    {
      return 0;
    }

    ScalarT t = output_times[0];

    solver->configureSimulation();
    solver->initializeSimulation(t, true);

    for (size_t time_idx = 1; time_idx < output_times.size(); time_idx++)
    {
      t = output_times[time_idx];

      const int retval = solver->runSimulation(t, 0.0);

      if (retval != 0)
      {
        return retval;
      }

      if (is_last_iteration)
      {
        gatherLocalVariables(*next_stages_y_[time_idx], *next_stages_yp_[time_idx], model);
      }
      else
      {
        gatherLocalVariables(*curr_stages_y_[time_idx], *curr_stages_yp_[time_idx], model);
      }
    }

    return 0;
  }

  template class MultiStageCosimulation<double, int>;

} // namespace AnalysisManager