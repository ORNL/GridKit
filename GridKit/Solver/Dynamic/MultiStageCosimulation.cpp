
#include "MultiStageCosimulation.hpp";

#include <alloca.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "GridKit/Solver/Dynamic/Interpolation.hpp"
#include "GridKit/Solver/Dynamic/Rosenbrock.hpp"

namespace Integrator
{

  template <class ScalarT, typename IdxT>
  MultiStageCosimulation<ScalarT, IdxT>::MultiStageCosimulation(
      std::vector<PartitionEvaluator*> partitions,
      std::size_t                      num_stages,
      MemorySpace                      memspace)
    : num_partitions_(partitions.size()),
      num_stages_(num_stages),
      memspace_(memspace),
      partitions_(std::move(partitions))
  {
    if (num_partitions_ == 0)
    {
      throw std::invalid_argument("MultiStageCosimulation requires at least one partition.");
    }

    if (num_stages_ == 0)
    {
      throw std::invalid_argument("MultiStageCosimulation requires at least one stage.");
    }

    for (auto* partition : partitions_)
    {
      if (partition == nullptr)
      {
        throw std::invalid_argument("Partition pointer cannot be null.");
      }
    }
  }

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::allocate()
  {
    num_partitions_ = partitions_.size();

    x0_local_.resize(num_partitions_);
    coupling_mat_.resize(num_partitions_);
    linear_workspaces_.resize(num_partitions_);
    lin_solver_.resize(num_partitions_);
    solvers_.resize(num_partitions_);

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto* part = partitions_[p];
      assert(part != nullptr);

      const IdxT num_states   = part->getStateSize();
      const IdxT num_coupling = part->getCouplingSize();

      x0_local_[p] = std::make_unique<State>(num_states);
      x0_local_[p]->allocate(memspace_);

      coupling_mat_[p] = std::make_unique<State>(num_coupling);
      coupling_mat_[p]->allocate(memspace_);

      linear_workspaces_[p] = std::make_unique<ReSolve::LinAlgWorkspaceCpu>();

      lin_solver_[p] = std::make_unique<ReSolve::SystemSolver>(
          linear_workspaces_[p].get(),
          "klu",
          "klu",
          "klu");

      lin_solver_[p]->initialize();

      solvers_[p] = std::make_unique<Rosenbrock>(
          Rosenbrock::Tableau::rodas5p(),
          part,
          lin_solver_[p].get(),
          &vector_handler_);

      solvers_[p]->allocate();
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeLocal(const State& global_y, PartitionEvaluator* partition)
  {
    auto&  internal = partition->getInternalIndecies();
    size_t n        = partition->size();

    matrix_handler_.matvec();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeCoupling(const State& global_y, int stage)
  {

    auto* global_data = global_y.getData(memspace_);

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto*  part     = partitions_[p];
      auto&  external = part->getExternalIndices(); // TODO: add this function to partition evaluate
      size_t n        = part->size();

      auto* local_data = coupling_mat_[p].getData(stage, memspace_);

      for (size_t i = 0; i < n; ++i)
      {
        local_data[i] = global_data[external[i]];
      }
    }

    return 0;
  }

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::timestep(const State& global_y)
  {
    return 0;
  }

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::partitionStageSolve(int stage, double t0, double dt)
  {

    std::vector<double> out_times;

    for (size_t iter = 0; iter < max_iterations_; iter++)
    {
      for (size_t part = 0; part < num_partitions_; part++)
      {

        auto* partitionEvaluator = partitions_[part];
        auto* partition_solver   = solvers_[part];

        if (stage == 0 && iter == 0)
        {
          distributeLocal(y_, partitionEvaluator);
        }

        auto& coupling_stage_values = coupling_mat_[part];

        GridKit::LagrangeInterpolant interp(coupling_stage_values, t0, dt);
        partitionEvaluator->addForcing(interp);

        for (size_t i = 0; i < num_partitions_; i++)
        {
          partition_solution_[i]->setToZero(memspace_);
        }

        int counter = 1;

        auto sol_out = [&](double tt, ReSolve::vector::Vector yp) -> void
        {
          auto* internal_csr = partition_internal_csr_[part].get();
          auto* internal_sol = partition_solution_[counter].get();

          double alpha = 1.0;
          double beta  = 1.0;

          matrix_handler_.matvec(internal_csr, &yp, internal_sol, &alpha, &beta, memspace_);

          counter++;
        };

        partition_solver->integrate(out_times, stepController_, params_, sol_out);
      }

      if (stage == 0)
      {
        break;
      }
    }

    return 0;
  }

} // namespace Integrator
