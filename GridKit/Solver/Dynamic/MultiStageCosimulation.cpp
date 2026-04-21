
#include "MultiStageCosimulation.hpp";

#include <alloca.h>

#include <cmath>
#include <cstddef>
#include <format>
#include <iterator>
#include <vector>

namespace Integrator
{

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::allocate()
  {
    x0_local_     = std::make_unique<std::unique_ptr<State>[]>(num_partitions_);
    coupling_Mat_ = std::make_unique<std::unique_ptr<State>[]>(num_partitions_);

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto* part         = partitions_[p];
      int   num_coupling = part->getCouplingSize();
      int   num_states   = part->getStateSize();

      x0_local_[p]     = std::make_unique<State>(num_states);
      coupling_Mat_[p] = std::make_unique<State>(num_coupling, num_stages_);

      x0_local_[p]->allocate(memspace_);
      coupling_Mat_[p]->allocate(memspace_);
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeLocal(const State* global_y)
  {

    for (size_t p = 0; p < num_partitions_; ++p)
    {
      auto&       part     = partitions_[p];
      const auto& external = part.getExternalIndices();
      size_t      n        = part.size();

      auto* local_data = x0_local_[p].get()->getData(p, memspace_);

      for (size_t i = 0; i < n; ++i)
      {
        local_data[i] = global_y[external[i]];
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::distributeCoupling(const State* global_y, int stage)
  {
    const IdxT num_parts = static_cast<IdxT>(partitions_.size());

    for (IdxT p = 0; p < num_parts; ++p)
    {
      auto*      part     = partitions_[p];
      auto       external = part->getCoupling();
      const IdxT n        = part->size();

      auto* local_data = x0_local_[p].get()->getData(p, memspace_);

      for (IdxT i = 0; i < n; ++i)
      {
        local_data[i] = global_y[external[i]];
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::timestep(const State& global_y)
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int MultiStageCosimulation<ScalarT, IdxT>::partitionSolve()
  {

    std::vector<double> output_time = output_times;

    for (size_t i = 0; i < max_iterations_; i++)
    {
      for (auto* partitionEvaluator : partitions_)
      {

        partitionEvaluator->addForcing(coupling);

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

    std::swap(stages_, partition_solutions_);

    return 0;
  }

} // namespace Integrator