#include "MultiStageCosimulation.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Interpolation.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

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

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::allocate()
  {

    linear_workspaces_.resize(num_partitions_);
    resolve_lin_solver_.resize(num_partitions_);
    gridkit_lin_solver_.resize(num_partitions_);
    solvers_.resize(num_partitions_);
    ida_solvers_.resize(num_partitions_);
    coupling_mat_y_.resize(num_partitions_);
    coupling_mat_yp_.resize(num_partitions_);
    // error_norms_.resize(num_partitions_);
    next_stages_y_.resize(num_stages_ + 1);
    curr_stages_y_.resize(num_stages_ + 1);
    next_stages_yp_.resize(num_stages_ + 1);
    curr_stages_yp_.resize(num_stages_ + 1);

    IdxT global_size = 0;
    // Allocate one coupling state for each stage of every partition.
    for (size_t p = 0; p < num_partitions_; ++p)
    {

      auto* partition = partitions_[p];

      partition->allocate();
      const auto coupling_size  = partition->getExternSize();
      const auto partition_size = partition->getInternalSize();

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

    // Allocate global solution vectors.
    y_.resize(global_size);
    yp_.resize(global_size);

    y_.allocate(memspace_);
    yp_.allocate(memspace_);

    y_.setToZero(memspace_);
    yp_.setToZero(memspace_);

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
      y_.copyFromExternal(*y, memspace_, memspace_);
      yp_.copyFromExternal(*yp, memspace_, memspace_);
    }
    else
    {
      for (size_t p = 0; p < num_partitions_; ++p)
      {
        gatherLocalVariables(y_, yp_, partitions_[p]);
      }
    }

    y_.setDataUpdated(memspace_);
    yp_.setDataUpdated(memspace_);

    time_points_ = time_points;

    return 0;
  }

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::runSimulation(std::optional<std::function<void(const State&)>> step_callback)
  {
    for (size_t t_idx = 0; t_idx + 1 < time_points_.size(); t_idx++)
    {
      const RealT t  = time_points_[t_idx];
      const RealT dt = time_points_[t_idx + 1] - time_points_[t_idx];

      curr_stages_y_[0]->copyFromExternal(y_, memspace_);
      next_stages_y_[0]->copyFromExternal(y_, memspace_);

      curr_stages_yp_[0]->copyFromExternal(yp_, memspace_);
      next_stages_yp_[0]->copyFromExternal(yp_, memspace_);

      for (size_t curr_stage = 1; curr_stage <= num_stages_; ++curr_stage)
      {
        const size_t next_stage = curr_stage + 1;

        partitionStageSolve(
            curr_stage,
            next_stage,
            t,
            dt);

        std::swap(curr_stages_y_, next_stages_y_);
        std::swap(curr_stages_yp_, next_stages_yp_);
      }

      // // std::cout << "End of step: " << t << std::endl;

      y_.copyFromExternal(*curr_stages_y_.back(), memspace_);
      yp_.copyFromExternal(*curr_stages_yp_.back(), memspace_);

      // ZOH(t,
      //     dt);
      // for (size_t partition_idx = 0; partition_idx < num_partitions_; partition_idx++)
      // {
      //   gatherLocalVariables(y_, yp_, partitions_[partition_idx]);
      // }
      // printState("Final y", y_);

      if (step_callback)
      {
        (*step_callback)(y_);
      }
    }

    // printState("Final solution: ", y_);

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
      y_data[i]  = global_y_data[partition->getNodeConnection(i)];
      yp_data[i] = global_yp_data[partition->getNodeConnection(i)];
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
      global_y_data[partition->getNodeConnection(i)]  = y_data[i];
      global_yp_data[partition->getNodeConnection(i)] = yp_data[i];
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
  typename MultiStageCosimulation<ScalarT, IdxT>::MergedTimePoints
  MultiStageCosimulation<ScalarT, IdxT>::mergeTimePoints(
      const std::vector<ScalarT>& first,
      const std::vector<ScalarT>& second)
  {
    constexpr RealT tolerance = 1e-14;

    MergedTimePoints result;

    result.times.reserve(first.size() + second.size());
    result.from_first.reserve(first.size() + second.size());
    result.from_second.reserve(first.size() + second.size());

    size_t first_idx  = 0;
    size_t second_idx = 0;

    while (first_idx < first.size() && second_idx < second.size())
    {
      const RealT first_time  = first[first_idx];
      const RealT second_time = second[second_idx];

      if (std::abs(first_time - second_time) <= tolerance)
      {
        result.times.push_back(first_time);
        result.from_first.push_back(true);
        result.from_second.push_back(true);

        ++first_idx;
        ++second_idx;
      }
      else if (first_time < second_time)
      {
        result.times.push_back(first_time);
        result.from_first.push_back(true);
        result.from_second.push_back(false);

        ++first_idx;
      }
      else
      {
        result.times.push_back(second_time);
        result.from_first.push_back(false);
        result.from_second.push_back(true);

        ++second_idx;
      }
    }

    while (first_idx < first.size())
    {
      result.times.push_back(first[first_idx]);
      result.from_first.push_back(true);
      result.from_second.push_back(false);

      ++first_idx;
    }

    while (second_idx < second.size())
    {
      result.times.push_back(second[second_idx]);
      result.from_first.push_back(false);
      result.from_second.push_back(true);

      ++second_idx;
    }

    return result;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::partitionStageSolve(size_t curr_stage, size_t next_stage, RealT t0, RealT dt)
  {
    assert(curr_stage >= 1);
    assert(next_stage >= 2);

    const auto curr_stage_output_times = Interpolant::getPhysicalOutputTimes(curr_stage, t0, dt);
    const auto next_stage_output_times = Interpolant::getPhysicalOutputTimes(next_stage, t0, dt);

    MergedTimePoints output_data = mergeTimePoints(curr_stage_output_times, next_stage_output_times);

    output_data.times.erase(output_data.times.begin());
    output_data.from_first.erase(output_data.from_first.begin());
    output_data.from_second.erase(output_data.from_second.begin());

    const size_t num_iterations = curr_stage == 1 ? 1 : max_iterations_;

    for (size_t iteration = 0; iteration < num_iterations; iteration++)
    {
      for (size_t partition_idx = 0; partition_idx < num_partitions_; partition_idx++)
      {
        auto* partition = partitions_[partition_idx];
        distributeLocalVariables(y_, yp_, partition);
        distributeCouplingVariables(partition_idx, curr_stage);
      }

      for (size_t partition_idx = 0; partition_idx < num_partitions_; partition_idx++)
      {
        auto* partition = partitions_[partition_idx];

        CouplingMatrix& coupling_y  = coupling_mat_y_[partition_idx];
        CouplingMatrix& coupling_yp = coupling_mat_yp_[partition_idx];

        CouplingMatrix stage_coupling_y(coupling_y.begin(), coupling_y.begin() + curr_stage);
        CouplingMatrix stage_coupling_yp(coupling_yp.begin(), coupling_yp.begin() + curr_stage);

        auto external_polynomial_y  = std::make_shared<Interpolant>(stage_coupling_y, t0, dt);
        auto external_polynomial_yp = std::make_shared<Interpolant>(stage_coupling_yp, t0, dt);

        partition->setForcingFunction([&](ScalarT t) mutable
                                      {
          external_polynomial_y->evaluate(t, partition->getExternalDataY());
          external_polynomial_yp->evaluate(t, partition->getExternalDataYP()); });

        size_t output_idx     = 0;
        size_t curr_stage_idx = 1;
        size_t next_stage_idx = 1;

        auto sol_out = [&]([[maybe_unused]] double tt) -> void
        {
          if (std::abs(tt - output_data.times[output_idx]) > 1e-14)
          {
            throw std::runtime_error("Output time mismatch in stage output.");
          }

          const bool is_last_iteration = (iteration == (num_iterations - 1));

          // Non-final iterations: the solve is still provisional. Park the freshly
          // computed *next-stage* values (from_second schedule) in next_stages_, where
          // they feed the coupling interpolant on the following pass.
          if (output_data.from_second[output_idx] && is_last_iteration)
          {
            gatherLocalVariables(*next_stages_y_[next_stage_idx],
                                 *next_stages_yp_[next_stage_idx],
                                 partition);
            ++next_stage_idx;
          }

          // Final iteration: converged values. Commit the *current-stage* values
          // (from_first schedule) to curr_stages_.
          if (output_data.from_first[output_idx] && !is_last_iteration)
          {
            gatherLocalVariables(*curr_stages_y_[curr_stage_idx],
                                 *curr_stages_yp_[curr_stage_idx],
                                 partition);
            ++curr_stage_idx;
          }

          ++output_idx;
        };

        ida_solvers_[partition_idx] = std::make_unique<IdaSolver>(partitions_[partition_idx]);
        integrate(
            ida_solvers_[partition_idx].get(),
            partition,
            t0,
            output_data.times,
            sol_out);
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::ZOH(
      RealT t0,
      RealT dt)
  {
    const RealT t1 = t0 + dt;

    const auto* global_y  = y_.getData(memspace_);
    const auto* global_yp = yp_.getData(memspace_);

    // Distribute the state and freeze the coupling variables at t0.
    for (size_t partition_idx = 0; partition_idx < num_partitions_; ++partition_idx)
    {
      auto* partition = partitions_[partition_idx];

      distributeLocalVariables(y_, yp_, partition);

      auto& external_y  = partition->getExternalDataY();
      auto& external_yp = partition->getExternalDataYP();

      auto* external_y_data  = external_y.getData(memspace_);
      auto* external_yp_data = external_yp.getData(memspace_);

      const auto& external_indices = partition->getExternalDataIndices();

      for (size_t i = 0; i < external_indices.size(); ++i)
      {
        const auto global_idx = external_indices[i];

        external_y_data[i]  = global_y[global_idx];
        external_yp_data[i] = global_yp[global_idx];
      }

      external_y.setDataUpdated(memspace_);
      external_yp.setDataUpdated(memspace_);

      std::cout << "partition: " << partition_idx << " Size: " << partition->size() << std::endl;
    }

    // Integrate every partition independently over [t0, t1].
    for (size_t partition_idx = 0; partition_idx < num_partitions_; ++partition_idx)
    {
      auto* partition = partitions_[partition_idx];

      const auto before = computeMaxResidual(partition, t0);

      std::cout << std::scientific << std::setprecision(6)
                << "t0 = " << t0 << "  partition " << partition_idx
                << "\n  before: max|F| = " << before.max_abs
                << "  (local eq " << before.local_idx
                << ", global " << before.global_idx << ")\n";

      auto output = [](RealT) {};

      ida_solvers_[partition_idx] = std::make_unique<IdaSolver>(partitions_[partition_idx]);
      const int retval            = integrate(
          ida_solvers_[partition_idx].get(),
          partition,
          t0,
          {t1},
          output);

      if (retval != 0)
      {
        return retval;
      }

      // Residual AFTER: IDA's solution at t1, same frozen external data
      const auto after = computeMaxResidual(partition, t1);

      std::cout << "  after:  max|F| = " << after.max_abs
                << "  (local eq " << after.local_idx
                << ", global " << after.global_idx << ")\n";
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::integrate(
      IdaSolver*                  solver,
      [[maybe_unused]] Subsystem* model,
      RealT                       t0,
      const std::vector<RealT>&   output_times,
      std::function<void(RealT)>  sol_out)
  {

    if (output_times.empty())
    {
      return 0;
    }

    // Reinitialize IDA from the current partition state at the
    // beginning of this stage solve.

    // solver->setSuppressAlgebraicErrors(true);
    solver->setMaxSteps(5000);
    solver->setTolerance(1e-12, 1e-12);
    solver->configureSimulation();
    solver->getDefaultInitialCondition();
    solver->initializeSimulation(t0, true);
    for (const RealT output_time : output_times)
    {
      const int retval = solver->runSimulation(
          output_time,
          0.0,
          sol_out);

      if (retval != 0)
      {
        return retval;
      }
    }

    return 0;
  }

  template class MultiStageCosimulation<double, int>;

} // namespace AnalysisManager
