
#include "MultiStageCosimulation.hpp";

#include <alloca.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace Integrator
{

  // TODO: Complete implementation
  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::allocate()
  {

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto* part = partitions_[p];

      int num_coupling = part->getCouplingSize();
      int num_states   = part->getStateSize();

      x0_local_[p]     = State(num_states);
      coupling_mat_[p] = State(num_coupling, num_stages_);

      x0_local_[p].allocate(memspace_);
      coupling_mat_[p].allocate(memspace_);
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeLocal(const State& global_y)
  {

    auto* global_data = global_y.getData(memspace_);

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto*  part     = partitions_[p];
      auto&  internal = part->getInternalIndecies(); // TODO: add this function to partition evaluate
      size_t n        = part->size();

      auto* local_data = x0_local_[p].getData(memspace_);

      for (size_t i = 0; i < n; ++i)
      {
        local_data[i] = global_data[internal[i]];
      }
    }

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
  int MultiStageCosimulation<ScalarT, IdxT>::partitionSolve()
  {

    std::vector<double> output_time = output_times;

    for (size_t i = 0; i < max_iterations_; i++)
    {
      for (auto* partitionEvaluator : partitions_)
      {

        // TODO: add a coupling function to the partition
        partitionEvaluator->addForcing();

        if (i + 1 < max_iterations_)
        {
          solver.integrate({t0}, partition_solution_.getData(stage));
        }
        else
        {
          solver.integrate(output_time, partition_solution_);
        }
      }

      updateCoupling(stage + 1)

          if (stage == 0)
      {
        break;
      }
    }

    return 0;
  }

} // namespace Integrator